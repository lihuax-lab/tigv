"""Genome coordinates and file-backed tracks. Coordinates are zero-based internally."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable
from collections import defaultdict
import heapq
import math
import re

import pysam

from .gtf import Gene, GeneIndex

try:
    import pyBigWig
except ImportError:  # A useful error is raised only when a bigWig is opened.
    pyBigWig = None


@dataclass
class View:
    chrom: str
    start: int
    span: int

    @property
    def end(self) -> int:
        return self.start + self.span

    def clamp(self, chrom_size: int | None = None) -> None:
        self.span = max(1, self.span)
        if chrom_size is not None:
            self.span = min(self.span, max(1, chrom_size))
            self.start = min(max(0, self.start), max(0, chrom_size - self.span))
        else:
            self.start = max(0, self.start)

    def pan(self, fraction: float, chrom_size: int | None = None) -> None:
        self.start += round(self.span * fraction)
        self.clamp(chrom_size)

    def zoom(self, factor: float, min_span: int, chrom_size: int | None = None) -> None:
        center = self.start + self.span / 2
        self.span = max(min_span, round(self.span * factor))
        self.start = round(center - self.span / 2)
        self.clamp(chrom_size)

    def interval(self, column: int, width: int) -> tuple[int, int]:
        left = self.start + column * self.span // width
        right = self.start + (column + 1) * self.span // width
        return left, max(left + 1, right)

    def column(self, pos: int, width: int) -> int:
        return (pos - self.start) * width // self.span


@dataclass(frozen=True)
class ReadHit:
    track_index: int
    chrom: str
    start: int
    end: int
    name: str
    flag: int


def parse_region(region: str) -> tuple[str, int | None, int | None]:
    match = re.fullmatch(r"([^:\s]+)(?::([\d,]+)(?:-([\d,]+))?)?", region.strip())
    if not match:
        raise ValueError("Use chr, chr:position, or chr:start-end")
    chrom, first, last = match.groups()
    if first is None:
        return chrom, None, None
    start = int(first.replace(",", ""))
    end = int(last.replace(",", "")) if last else start
    if start < 1 or end < start:
        raise ValueError("Coordinates must be positive and start <= end")
    return chrom, start - 1, end


class Tracks:
    def __init__(self) -> None:
        self.bams: list[tuple[str, pysam.AlignmentFile]] = []
        self.bigwigs: list[tuple[str, object]] = []
        self.fasta: pysam.FastaFile | None = None
        self.fasta_path: str | None = None
        self.genes: GeneIndex | None = None
        self.gtf_path: str | None = None
        self.chrom_sizes: dict[str, int] = {}
        self.read_name_cache: dict[str, list[ReadHit]] = {}

    def close(self) -> None:
        for _, bam in self.bams:
            bam.close()
        for _, bigwig in self.bigwigs:
            bigwig.close()
        if self.fasta:
            self.fasta.close()

    def add_bam(self, path: str) -> None:
        path = str(Path(path).expanduser().resolve())
        bam = pysam.AlignmentFile(path, "rb")
        try:
            if not bam.has_index():
                raise ValueError(f"BAM index missing: run samtools index {path}")
            self.bams.append((Path(path).name, bam))
            self.read_name_cache.clear()
            for chrom, length in zip(bam.references, bam.lengths):
                self.chrom_sizes[chrom] = max(self.chrom_sizes.get(chrom, 0), length)
        except Exception:
            bam.close()
            raise

    def add_bigwig(self, path: str) -> None:
        if pyBigWig is None:
            raise RuntimeError("pyBigWig is required for bigWig tracks; install environment.yml")
        path = str(Path(path).expanduser().resolve())
        bigwig = pyBigWig.open(path)
        if bigwig is None:
            raise ValueError(f"Could not open bigWig: {path}")
        self.bigwigs.append((Path(path).name, bigwig))
        for chrom, length in bigwig.chroms().items():
            self.chrom_sizes[chrom] = max(self.chrom_sizes.get(chrom, 0), length)

    def set_fasta(self, path: str) -> None:
        path = str(Path(path).expanduser().resolve())
        fasta = pysam.FastaFile(path)
        if self.fasta:
            self.fasta.close()
        self.fasta = fasta
        self.fasta_path = path
        for chrom, length in zip(fasta.references, fasta.lengths):
            self.chrom_sizes[chrom] = max(self.chrom_sizes.get(chrom, 0), length)

    def set_gtf(self, path: str) -> None:
        path = str(Path(path).expanduser().resolve())
        genes = GeneIndex.from_file(path)
        self.genes = genes
        self.gtf_path = path
        for chrom, length in genes.chrom_sizes.items():
            self.chrom_sizes[chrom] = max(self.chrom_sizes.get(chrom, 0), length)

    def resolve_chrom(self, chrom: str) -> str:
        if chrom in self.chrom_sizes:
            return chrom
        alias = chrom[3:] if chrom.startswith("chr") else "chr" + chrom
        if alias in self.chrom_sizes:
            return alias
        raise ValueError(f"Unknown chromosome: {chrom}")

    def chrom_size(self, chrom: str) -> int | None:
        return self.chrom_sizes.get(chrom)

    def reference(self, chrom: str, start: int, end: int) -> str | None:
        if not self.fasta or chrom not in self.fasta.references:
            return None
        return self.fasta.fetch(chrom, max(0, start), min(end, self.fasta.get_reference_length(chrom))).upper()

    def find_gene(self, name: str) -> list[Gene]:
        if not self.genes:
            raise ValueError("Load a GTF before searching genes")
        matches = self.genes.search(name)
        if not matches:
            suggestions = self.genes.suggest(name)
            suffix = f"; prefixes: {', '.join(suggestions)}" if suggestions else ""
            raise ValueError(f"Gene not found: {name}{suffix}")
        return matches

    def find_read(self, name: str, progress: Callable[[int, int], None] | None = None) -> list[ReadHit]:
        """Search exact query names; coordinate BAM indexes cannot answer this directly."""
        if not self.bams:
            raise ValueError("Load a BAM before searching reads")
        if name in self.read_name_cache:
            return self.read_name_cache[name]
        hits: list[ReadHit] = []
        for track_index, (_, bam) in enumerate(self.bams):
            scanned = 0
            # A separate handle leaves the coordinate iterator used by rendering untouched.
            with pysam.AlignmentFile(bam.filename, "rb") as reader:
                for read in reader.fetch(until_eof=True):
                    scanned += 1
                    if read.query_name == name and not read.is_unmapped:
                        hits.append(ReadHit(track_index, reader.get_reference_name(read.reference_id),
                                            read.reference_start, read.reference_end or read.reference_start + 1,
                                            name, read.flag))
                    if progress and scanned % 100_000 == 0:
                        progress(track_index, scanned)
        if not hits:
            raise ValueError(f"Mapped read not found: {name}")
        self.read_name_cache[name] = hits
        return hits


def bigwig_bins(bigwig: object, view: View, width: int) -> list[float | None]:
    chrom_size = bigwig.chroms(view.chrom)
    if chrom_size is None or view.start >= chrom_size:
        return [None] * width
    end = min(view.end, chrom_size)
    # pyBigWig stats returns one value per screen column without fetching huge arrays.
    count = min(width, max(1, end - view.start))
    values = bigwig.stats(view.chrom, view.start, end, nBins=count, type="mean", exact=False)
    return [values[min(count - 1, column * count // width)] for column in range(width)]


def coverage_and_reads(
    bam: pysam.AlignmentFile, view: View, width: int, max_reads: int = 12_000,
    highlighted_name: str | None = None,
) -> tuple[list[int], list[pysam.AlignedSegment], bool]:
    chrom_size = bam.get_reference_length(view.chrom) if view.chrom in bam.references else 0
    if not chrom_size or view.start >= chrom_size:
        return [0] * width, [], False
    diff = [0] * (width + 1)
    reads: list[pysam.AlignedSegment] = []
    truncated = False
    for read in bam.fetch(view.chrom, view.start, min(view.end, chrom_size)):
        if read.is_unmapped:
            continue
        display_only = read.is_secondary or read.is_supplementary
        if read.is_secondary and read.query_name != highlighted_name:
            continue
        if len(reads) >= max_reads and read.query_name != highlighted_name:
            truncated = True
            if highlighted_name is None:
                break
            continue
        if len(reads) >= max_reads:
            truncated = True
        reads.append(read)
        if display_only:
            continue
        for start, end in read.get_blocks():
            left = max(0, view.column(max(start, view.start), width))
            right = min(width, math.ceil((min(end, view.end) - view.start) * width / view.span))
            if right > left:
                diff[left] += 1
                diff[right] -= 1
    running = 0
    coverage = []
    for delta in diff[:-1]:
        running += delta
        coverage.append(running)
    return coverage, reads, truncated


def read_cells(
    read: pysam.AlignedSegment, view: View, width: int, reference: str | None
) -> dict[int, tuple[str, bool]]:
    """Return column -> (glyph, mismatch), hiding indels until fully resolved."""
    cells: dict[int, tuple[str, bool]] = {}
    insertion_columns: list[int] = []
    md_bases: dict[int, str] = {}
    if reference is None and read.has_tag("MD"):
        try:
            md_bases = {
                ref_pos: base.upper()
                for _, ref_pos, base in read.get_aligned_pairs(with_seq=True)
                if ref_pos is not None and base is not None
            }
        except ValueError:
            pass
    seq = read.query_sequence or ""
    query_pos = 0
    ref_pos = read.reference_start
    base_resolution = view.span <= width
    for op, length in read.cigartuples or ():
        if op in (0, 7, 8):  # M, =, X
            for offset in range(length):
                pos = ref_pos + offset
                if view.start <= pos < view.end:
                    col = view.column(pos, width)
                    if 0 <= col < width:
                        base = seq[query_pos + offset].upper() if query_pos + offset < len(seq) else "N"
                        ref_base = (
                            reference[pos - view.start]
                            if reference and pos - view.start < len(reference)
                            else md_bases.get(pos)
                        )
                        mismatch = (base != ref_base if ref_base and ref_base != "N" else op == 8)
                        if base_resolution:
                            cells[col] = (base, mismatch)
                        elif mismatch:
                            cells[col] = ("*", True)
                        else:
                            cells.setdefault(col, ("#", False))
            query_pos += length
            ref_pos += length
        elif op == 1:  # insertion before the next reference base
            if (base_resolution and length <= view.span and query_pos + length <= len(seq)
                    and view.start <= ref_pos - 1 and ref_pos < view.end):
                col = view.column(ref_pos, width)
                if 0 <= col < width:
                    insertion_columns.append(col)
            query_pos += length
        elif op in (2, 3):  # deletion or skipped intron
            if op == 2 and base_resolution and view.start <= ref_pos and ref_pos + length <= view.end:
                for pos in range(ref_pos, ref_pos + length):
                    col = view.column(pos, width)
                    if 0 <= col < width:
                        cells[col] = ("-", True)
            ref_pos += length
        elif op == 4:  # soft clipping
            query_pos += length
    for col in insertion_columns:
        cells[col] = ("+", True)
    return cells


def read_groups(reads: list[pysam.AlignedSegment]) -> list[list[pysam.AlignedSegment]]:
    """Keep alignments carrying the same query name together on one row."""
    groups: dict[object, list[pysam.AlignedSegment]] = defaultdict(list)
    for read in reads:
        key: object = read.query_name if read.query_name is not None else id(read)
        groups[key].append(read)
    return sorted(groups.values(), key=lambda group: min(read.reference_start for read in group))


def read_fastq(read: pysam.AlignedSegment) -> tuple[str | None, str | None]:
    """Reconstruct four FASTQ lines from the sequence and qualities stored in BAM."""
    sequence = read.get_forward_sequence()
    if not sequence:
        return None, "Sequence is absent from BAM"
    if read.query_qualities is None:
        return None, "Base qualities are absent or invalid in BAM"
    qualities = read.get_forward_qualities()
    if qualities is None or len(qualities) != len(sequence) or any(q > 93 for q in qualities):
        return None, "Base qualities are absent or invalid in BAM"
    quality = "".join(chr(q + 33) for q in qualities)
    return f"@{read.query_name or 'unnamed'}\n{sequence}\n+\n{quality}", None


def _gap_cells(left: int, right: int, crosses_genes: bool) -> dict[int, str]:
    """Draw the open interval [left, right) between two aligned blocks."""
    size = right - left
    if size <= 0:
        return {}
    if not crosses_genes:
        return {column: "·" for column in range(left, right)}
    if size == 1:
        return {left: "["}
    return {left: "[", right - 1: "]"}


def _remote_targets(read: pysam.AlignedSegment) -> list[tuple[int, int]]:
    targets: list[tuple[int, int]] = []
    if read.is_paired and not read.mate_is_unmapped and read.next_reference_id >= 0:
        targets.append((read.next_reference_id, read.next_reference_start))
    if read.has_tag("SA"):
        header = read.header
        for entry in read.get_tag("SA").rstrip(";").split(";"):
            fields = entry.split(",")
            if len(fields) < 2:
                continue
            try:
                targets.append((header.get_tid(fields[0]), int(fields[1]) - 1))
            except (ValueError, KeyError):
                continue
    return [(tid, pos) for tid, pos in targets if tid >= 0 and pos >= 0]


def _gene_ids(genes: GeneIndex | None, chrom: str, start: int, end: int) -> set[tuple[str, str]]:
    if genes is None:
        return set()
    candidates = genes.overlap(chrom, start, end)
    exonic = {
        (gene.chrom, gene.gene_id)
        for gene in candidates
        if gene.exons and any(exon_start < end and exon_end > start for exon_start, exon_end in gene.exons)
    }
    return exonic or {(gene.chrom, gene.gene_id) for gene in candidates}


def _crosses_genes(
    genes: GeneIndex | None, left_chrom: str, left_start: int, left_end: int,
    right_chrom: str, right_start: int, right_end: int,
) -> bool:
    left_ids = _gene_ids(genes, left_chrom, left_start, left_end)
    right_ids = _gene_ids(genes, right_chrom, right_start, right_end)
    return bool(left_ids and right_ids and left_ids.isdisjoint(right_ids))


def _deletion_spans(read: pysam.AlignedSegment) -> list[tuple[int, int]]:
    spans: list[tuple[int, int]] = []
    ref_pos = read.reference_start
    for op, length in read.cigartuples or ():
        if op == 2:
            spans.append((ref_pos, ref_pos + length))
            ref_pos += length
        elif op in (0, 3, 7, 8):
            ref_pos += length
    return spans


def read_connection_cells(
    reads: list[pysam.AlignedSegment], view: View, width: int,
    genes: GeneIndex | None = None,
) -> dict[int, str]:
    """Join one read's pieces; bracket only transitions between annotated genes."""
    if not reads:
        return {}
    blocks = sorted(block for read in reads for block in read.get_blocks() if block[1] > block[0])
    if not blocks:
        return {}
    merged: list[list[int]] = []
    for start, end in blocks:
        if merged and start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    cells: dict[int, str] = {}
    deletions = [span for read in reads for span in _deletion_spans(read)]
    for (left_start, left_end), (right_start, right_end) in zip(merged, merged[1:]):
        if any(start <= left_end and end >= right_start for start, end in deletions):
            continue
        first = max(0, view.column(left_end - 1, width) + 1)
        last = min(width, view.column(right_start, width))
        crossing = _crosses_genes(genes, view.chrom, left_start, left_end,
                                  view.chrom, right_start, right_end)
        cells.update(_gap_cells(first, last, crossing))

    # An offscreen mate or supplementary segment still gets a visible open end.
    for read in reads:
        for target_tid, target_pos in _remote_targets(read):
            if any(other is not read and other.reference_id == target_tid
                   and other.reference_start <= target_pos < (other.reference_end or other.reference_start + 1)
                   for other in reads):
                continue
            try:
                target_chrom = read.header.get_reference_name(target_tid)
            except (ValueError, IndexError):
                continue
            read_blocks = read.get_blocks()
            if not read_blocks:
                continue
            if target_tid == read.reference_id and target_pos < read.reference_start:
                edge_start, edge_end = read_blocks[0]
                crossing = _crosses_genes(genes, view.chrom, edge_start, edge_end,
                                          target_chrom, target_pos, target_pos + 1)
                anchor = view.column(read.reference_start, width)
                glyphs = "]" if crossing else "·····"
                for column, glyph in zip(range(anchor - len(glyphs), anchor), glyphs):
                    if 0 <= column < width:
                        cells[column] = glyph
            else:
                edge_start, edge_end = read_blocks[-1]
                crossing = _crosses_genes(genes, view.chrom, edge_start, edge_end,
                                          target_chrom, target_pos, target_pos + 1)
                anchor = view.column((read.reference_end or read.reference_start + 1) - 1, width) + 1
                glyphs = "[" if crossing else "·····"
                for column, glyph in zip(range(anchor, anchor + len(glyphs)), glyphs):
                    if 0 <= column < width:
                        cells[column] = glyph
    return cells


def layout_reads(
    reads: list[pysam.AlignedSegment], view: View, width: int,
    genes: GeneIndex | None = None,
) -> list[list[pysam.AlignedSegment]]:
    """Place read groups on rows without putting other reads in their gaps."""
    rows: list[list[pysam.AlignedSegment]] = []
    row_ends: list[tuple[int, int]] = []
    for group in read_groups(reads):
        left = max(0, min(view.column(read.reference_start, width) for read in group))
        right = min(width, max(view.column((read.reference_end or read.reference_start + 1) - 1, width) + 1
                               for read in group))
        connector = read_connection_cells(group, view, width, genes)
        if connector:
            left = min(left, min(connector))
            right = max(right, max(connector) + 1)
        if right <= left:
            continue
        if row_ends and left > row_ends[0][0]:
            _, row_no = heapq.heappop(row_ends)
            rows[row_no].extend(group)
        else:
            row_no = len(rows)
            rows.append(list(group))
        heapq.heappush(row_ends, (right, row_no))
    return rows


def pack_reads(reads: list[pysam.AlignedSegment], view: View, width: int, limit: int, offset: int = 0,
               genes: GeneIndex | None = None) -> list[list[pysam.AlignedSegment]]:
    return layout_reads(reads, view, width, genes)[offset : offset + limit]
