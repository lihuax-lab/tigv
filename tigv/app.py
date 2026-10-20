"""Curses interface for tigv."""

from __future__ import annotations

import argparse
import curses
from dataclasses import dataclass
import math
import shlex
import sys

import pysam

from .core import (ReadHit, Tracks, View, bigwig_bins,
                   coverage_and_reads, layout_reads, parse_region, read_cells,
                   read_connection_cells, read_fastq, read_groups)
from .gtf import Gene


LABEL = 13
HELP = [
    "Navigation: h/l move 10%, H/L move 50%, +/- zoom, z base view",
    "            gg start, G end, [/] previous/next chromosome",
    "            J/K select track (including Reference), j/k scroll its rows",
    "Search:     /gene-name, /chr:start-end, /@read-name",
    "            :read NAME, n/N next/previous hit",
    "            :gene NAME, :goto chr:start-end",
    "Genes:      e toggles layout; :expand separates genes; :collapse packs",
    "Import:     :bam PATH, :bw PATH, :fa PATH, :gtf PATH",
    "Other:      :help, :quit, q quit, ? help, r redraw",
    "Colors:     red = mismatch / insertion / deletion; + insertion; - deletion",
    "            # = match, * = mismatch while zoomed out",
    "            +/- need base resolution and the complete event in view",
    "            · joins exons; cross-gene chRNA shows [     ]",
    "            GTF annotation is needed to identify cross-gene breaks",
    "Mouse:      click a read to open details and reconstructed FASTQ",
    "            detail window: j/k scroll, h/l move sequence, Esc closes",
    "Files:      BAM needs .bai; FASTA needs .fai; chromosome names must match",
    "            read-name search scans each BAM once; exact, case-sensitive",
    "Press any key to return",
]


@dataclass
class ReadSelection:
    track_name: str
    track_index: int
    read: pysam.AlignedSegment
    group: list[pysam.AlignedSegment]


def gene_lanes(genes: list[Gene], view: View, width: int, expanded: bool) -> list[list[Gene]]:
    """Expanded genes each own a row; packed genes share rows only when separated."""
    if expanded:
        return [[gene] for gene in genes]
    lanes: list[list[Gene]] = []
    ends: list[int] = []
    for gene in genes:
        left = max(0, view.column(gene.start, width))
        right = min(width, view.column(gene.end - 1, width) + 1)
        lane = next((index for index, end in enumerate(ends) if left > end + 2), len(ends))
        if lane == len(lanes):
            lanes.append([])
            ends.append(right)
        else:
            ends[lane] = right
        lanes[lane].append(gene)
    return lanes


class Browser:
    def __init__(self, tracks: Tracks, view: View, gene: str | None = None):
        self.tracks = tracks
        self.view = view
        self.status = "Press ? for help  |  / gene or locus search  |  : command"
        self.read_offsets: dict[int, int] = {}
        self.gene_offset = 0
        self.gene_expanded = False
        self.focus_index = 0
        self.track_top = 1
        self.last_search = gene
        self.gene_hits = []
        self.gene_hit_index = 0
        self.read_hits: list[ReadHit] = []
        self.read_hit_index = 0
        self.selected_read_name: str | None = None
        self.pending_read_reveal = False
        self.search_mode: str | None = None
        self.show_help = False
        self.last_g = False
        self.has_colors = False
        self.read_click_map: dict[tuple[int, int], ReadSelection] = {}
        self.read_detail: ReadSelection | None = None
        self.detail_lines: list[tuple[str, bool]] = []
        self.detail_scroll = 0
        self.detail_x = 0
        self.detail_fastq_length = 0

    def color(self, number: int, extra: int = 0) -> int:
        return (curses.color_pair(number) if self.has_colors else 0) | extra

    def put(self, screen, row: int, col: int, value: str, attr: int = 0) -> None:
        height, width = screen.getmaxyx()
        if 0 <= row < height and 0 <= col < width - 1:
            try:
                screen.addnstr(row, col, value, width - col - 1, attr)
            except curses.error:
                pass

    def label(self, screen, row: int, name: str, focused: bool = False) -> None:
        label = (">" if focused else " ") + name[: LABEL - 2]
        self.put(screen, row, 0, label.ljust(LABEL - 1), self.color(1 if focused else 3, curses.A_BOLD))

    def track_items(self) -> list[tuple[str, int, str, object]]:
        items: list[tuple[str, int, str, object]] = []
        if self.tracks.genes:
            items.append(("genes", 0, "Genes", self.tracks.genes))
        items.extend(("bigwig", index, name, bw) for index, (name, bw) in enumerate(self.tracks.bigwigs))
        items.extend(("bam", index, name, bam) for index, (name, bam) in enumerate(self.tracks.bams))
        return items

    def focus_track(self, kind: str, track_index: int = 0) -> None:
        for index, (item_kind, item_index, _, _) in enumerate(self.track_items(), 1):
            if item_kind == kind and item_index == track_index:
                self.focus_index = index
                self.track_top = index
                return

    def draw(self, screen) -> None:
        screen.erase()
        self.read_click_map.clear()
        height, cols = screen.getmaxyx()
        if height < 12 or cols < 45:
            self.put(screen, 0, 0, "Terminal needs at least 45 columns and 12 rows")
            screen.refresh()
            return
        if self.show_help:
            self.draw_help(screen)
            return
        width = cols - LABEL - 1
        self.view.clamp(self.tracks.chrom_size(self.view.chrom))
        # One column per base is the closest meaningful zoom level.
        if self.view.span < width:
            self.view.zoom(width / self.view.span, width, self.tracks.chrom_size(self.view.chrom))
        row = 0
        title = f" tigv  {self.view.chrom}:{self.view.start + 1:,}-{self.view.end:,}  ({self.view.span:,} bp; {self.view.span / width:.1f} bp/col) "
        self.put(screen, row, 0, title.ljust(cols - 1), self.color(1, curses.A_BOLD))
        row += 1
        self.label(screen, row, "Position")
        self.draw_ruler(screen, row, width)
        row += 1
        self.label(screen, row, "Reference", self.focus_index == 0)
        reference = self.tracks.reference(self.view.chrom, self.view.start, self.view.end)
        self.draw_reference(screen, row, width, reference)
        row += 1
        content_end = height - 2
        items = self.track_items()
        for position in range(max(0, self.track_top - 1), len(items)):
            if row >= content_end:
                break
            kind, track_index, name, handle = items[position]
            focused = self.focus_index == position + 1
            remaining = len(items) - position - 1
            if kind == "genes":
                budget = max(1, content_end - row - 1 - remaining * 3)
                row = self.draw_genes(screen, row, width, content_end, budget, focused)
            elif kind == "bigwig":
                panel_end = max(row + 1, content_end - remaining * 3)
                row = self.draw_signal(screen, row, width, min(panel_end, content_end), name, handle, focused)
            else:
                panel_end = max(row + 4, content_end - remaining * 2)
                row = self.draw_bam(screen, row, width, min(panel_end, content_end), name, handle,
                                    reference, track_index, focused)
        if not self.tracks.bams and not self.tracks.bigwigs and not self.tracks.genes:
            self.put(screen, row, 0, "Load data: :bam file.bam  :bw file.bw  :fa ref.fa  :gtf genes.gtf")
        self.put(screen, height - 2, 0, " h/l pan  +/- zoom  J/K track  / search  e expand  click read  ? help ".ljust(cols - 1), self.color(1))
        self.put(screen, height - 1, 0, self.status[: cols - 1].ljust(cols - 1), self.color(4))
        screen.refresh()
        if self.read_detail is not None:
            self.draw_read_detail(screen)

    def draw_help(self, screen) -> None:
        height, width = screen.getmaxyx()
        self.put(screen, 0, 0, " tigv help ".ljust(width - 1), self.color(1, curses.A_BOLD))
        for row, line in enumerate(HELP, 2):
            if row >= height:
                break
            self.put(screen, row, 2, line)
        screen.refresh()

    def draw_ruler(self, screen, row: int, width: int) -> None:
        for col in range(width):
            self.put(screen, row, LABEL + col, "·", self.color(4))
        step = max(1, width // 6)
        for col in range(0, width, step):
            pos, _ = self.view.interval(col, width)
            self.put(screen, row, LABEL + col, f"|{pos + 1:,}", self.color(3))

    def draw_reference(self, screen, row: int, width: int, reference: str | None) -> None:
        if reference is None:
            self.put(screen, row, LABEL, "(load FASTA to compare read bases)", self.color(4))
            return
        if self.view.span > width:
            self.put(screen, row, LABEL, "zoom in with + or z to see bases", self.color(4))
            return
        for offset, base in enumerate(reference):
            col = self.view.column(self.view.start + offset, width)
            if 0 <= col < width:
                self.put(screen, row, LABEL + col, base, self.color(2, curses.A_BOLD))

    def draw_genes(self, screen, row: int, width: int, bottom: int, budget: int, focused: bool) -> int:
        genes = self.tracks.genes.overlap(self.view.chrom, self.view.start, self.view.end)
        lanes = gene_lanes(genes, self.view, width, self.gene_expanded)
        self.gene_offset = min(self.gene_offset, max(0, len(lanes) - 1))
        self.label(screen, row, "Genes", focused)
        mode = "expanded" if self.gene_expanded else "packed"
        first = self.gene_offset + 1 if lanes else 0
        last = min(len(lanes), self.gene_offset + budget)
        self.put(screen, row, LABEL, f"{len(genes)} genes  {mode}  rows {first}-{last}/{len(lanes)}", self.color(4))
        row += 1
        if not genes and row < bottom:
            self.put(screen, row, LABEL, "No annotated genes in view", self.color(4))
            return row + 1
        visible = lanes[self.gene_offset : self.gene_offset + min(budget, bottom - row)]
        for line, genes_in_lane in enumerate(visible):
            y = row + line
            if self.gene_expanded:
                self.label(screen, y, genes_in_lane[0].name)
            for gene in genes_in_lane:
                left = max(0, self.view.column(gene.start, width))
                right = min(width, self.view.column(gene.end - 1, width) + 1)
                for col in range(left, right):
                    self.put(screen, y, LABEL + col, "→" if gene.strand == "+" else "←" if gene.strand == "-" else "─", self.color(4))
                for start, end in gene.merged_exons():
                    a = max(0, self.view.column(max(start, self.view.start), width))
                    b = min(width, math.ceil((min(end, self.view.end) - self.view.start) * width / self.view.span))
                    for col in range(a, max(a, b)):
                        self.put(screen, y, LABEL + col, "█", self.color(2, curses.A_BOLD))
                if not self.gene_expanded and right - left >= len(gene.name) + 2:
                    self.put(screen, y, LABEL + left + 1, gene.name, self.color(3, curses.A_BOLD))
        return row + len(visible)

    def draw_signal(self, screen, row: int, width: int, bottom: int, name: str, bigwig, focused: bool) -> int:
        if row >= bottom:
            return row
        values = bigwig_bins(bigwig, self.view, width)
        finite = [value for value in values if value is not None and math.isfinite(value)]
        max_value = max(finite, default=0)
        min_value = min(finite, default=0)
        self.label(screen, row, "bigWig", focused)
        self.put(screen, row, LABEL, f"{name}  min {min_value:.3g}  max {max_value:.3g}", self.color(4))
        row += 1
        plot_height = min(3, bottom - row)
        if min_value < 0:
            upper = 0 if max_value <= 0 else max(1, plot_height // 2)
            lower = plot_height - upper
            for col, value in enumerate(values):
                if value is None or not math.isfinite(value):
                    continue
                if value >= 0 and upper:
                    filled = math.ceil(value / max(max_value, 1e-9) * upper)
                    for band in range(filled):
                        self.put(screen, row + upper - 1 - band, LABEL + col, "█", self.color(5))
                elif value < 0 and lower:
                    filled = math.ceil(value / min_value * lower)
                    for band in range(filled):
                        self.put(screen, row + upper + band, LABEL + col, "█", self.color(6))
        else:
            positive = max(max_value, 1e-9)
            for col, value in enumerate(values):
                if value is None or not math.isfinite(value):
                    continue
                level = max(0, min(plot_height * 8, round(value / positive * plot_height * 8)))
                for band in range(plot_height):
                    units = level - (plot_height - band - 1) * 8
                    if units > 0:
                        glyph = " ▁▂▃▄▅▆▇█"[min(8, units)]
                        self.put(screen, row + band, LABEL + col, glyph, self.color(5))
        return row + plot_height

    def draw_bam(self, screen, row: int, width: int, bottom: int, name: str, bam,
                 reference: str | None, track_index: int, focused: bool) -> int:
        coverage, reads, truncated = coverage_and_reads(
            bam, self.view, width, highlighted_name=self.selected_read_name)
        peak = max(coverage, default=0)
        self.label(screen, row, "BAM", focused)
        more = " (read cap reached)" if truncated else ""
        self.put(screen, row, LABEL, f"{name}  {len(reads)} reads  coverage max {peak}{more}", self.color(4))
        row += 1
        plot_height = min(2, bottom - row)
        for col, value in enumerate(coverage):
            level = math.ceil(value * plot_height / peak) if peak else 0
            for band in range(plot_height):
                if level > plot_height - band - 1:
                    self.put(screen, row + band, LABEL + col, "█", self.color(5))
        row += plot_height
        all_rows = layout_reads(reads, self.view, width, self.tracks.genes)
        if self.pending_read_reveal and self.read_hits:
            hit = self.read_hits[self.read_hit_index]
            if hit.track_index == track_index and hit.chrom == self.view.chrom:
                for index, packed in enumerate(all_rows):
                    if any(read.query_name == hit.name and read.reference_start == hit.start
                           and (read.reference_end or read.reference_start + 1) == hit.end
                           for read in packed):
                        self.read_offsets[track_index] = index
                        self.pending_read_reveal = False
                        break
        offset = self.read_offsets.get(track_index, 0)
        read_rows = all_rows[offset : offset + max(0, bottom - row)]
        for packed in read_rows:
            if row >= bottom:
                break
            groups = read_groups(packed)
            group_by_read = {id(read): group for group in groups for read in group}
            for group in groups:
                highlighted = self.selected_read_name is not None and group[0].query_name == self.selected_read_name
                for col, glyph in read_connection_cells(group, self.view, width, self.tracks.genes).items():
                    pair = 7 if highlighted else (3 if glyph in "[]" else 4)
                    self.put(screen, row, LABEL + col, glyph,
                             self.color(pair, curses.A_BOLD if highlighted or glyph in "[]" else 0))
                    self.read_click_map[(row, LABEL + col)] = ReadSelection(name, track_index, group[0], group)
            if len(groups) == 1:
                selection = ReadSelection(name, track_index, groups[0][0], groups[0])
                for col in range(LABEL):
                    self.read_click_map[(row, col)] = selection
            for read in packed:
                cells = read_cells(read, self.view, width, reference)
                highlighted = self.selected_read_name is not None and read.query_name == self.selected_read_name
                if highlighted:
                    self.put(screen, row, 0, ("*" + read.query_name)[: LABEL - 1].ljust(LABEL - 1),
                             self.color(7, curses.A_BOLD | (0 if self.has_colors else curses.A_REVERSE)))
                for col, (glyph, mismatch) in cells.items():
                    pair = (8 if mismatch else 7) if highlighted else (6 if mismatch else 4)
                    self.put(screen, row, LABEL + col, glyph,
                             self.color(pair, curses.A_BOLD if mismatch or highlighted else 0))
                    group = group_by_read[id(read)]
                    self.read_click_map[(row, LABEL + col)] = ReadSelection(name, track_index, read, group)
            row += 1
        return row

    def open_read_detail(self, selection: ReadSelection) -> None:
        self.read_detail = selection
        self.detail_scroll = 0
        self.detail_x = 0
        read = selection.read
        group = selection.group
        selected_index = next(index for index, item in enumerate(group, 1) if item is read)
        lines: list[tuple[str, bool]] = [
            (f"Read: {read.query_name or '(unnamed)'}", False),
            (f"BAM: {selection.track_name}   alignment {selected_index}/{len(group)} selected", False),
            ("", False),
        ]
        for index, item in enumerate(group, 1):
            chrom = item.reference_name or "*"
            start = item.reference_start + 1
            end = item.reference_end or start
            role = "supplementary" if item.is_supplementary else "secondary" if item.is_secondary else "primary"
            mate = " read1" if item.is_read1 else " read2" if item.is_read2 else ""
            marker = ">" if item is read else " "
            lines.extend([
                (f"{marker}{index}. {chrom}:{start:,}-{end:,}  {'-' if item.is_reverse else '+'}  {role}{mate}", False),
                (f"    CIGAR {item.cigarstring or '*'}   MAPQ {item.mapping_quality}   FLAG {item.flag}", False),
                (f"    query {item.query_alignment_start + 1}-{item.query_alignment_end} / "
                 f"{item.query_length or 0} bp", False),
            ])
            if item.is_paired and not item.mate_is_unmapped:
                mate_chrom = item.next_reference_name or "*"
                lines.append((f"    mate {mate_chrom}:{item.next_reference_start + 1:,}  TLEN {item.template_length}", False))
            for tag, value in item.get_tags():
                lines.append((f"    {tag}: {value}", False))
            lines.append(("", False))

        same_mate = [item for item in group if item.is_read1 == read.is_read1
                     and item.is_read2 == read.is_read2]
        candidates = sorted(same_mate, key=lambda item: (
            item.is_secondary or item.is_supplementary, item is not read))
        fastq = None
        source = read
        reason = None
        for item in candidates:
            fastq, reason = read_fastq(item)
            if fastq is not None:
                source = item
                break
        lines.append(("FASTQ (reconstructed from BAM; original header may differ)", False))
        if fastq is None:
            lines.append((reason or "FASTQ is unavailable", False))
            self.detail_fastq_length = 0
        else:
            if source is not read:
                source_index = next(index for index, item in enumerate(group, 1) if item is source)
                lines.append((f"Using alignment {source_index} for sequence and qualities", False))
            if any(op == 5 for op, _ in source.cigartuples or ()):
                lines.append(("Hard-clipped bases are absent from BAM; sequence may be partial", False))
            header, sequence, plus, qualities = fastq.split("\n")
            lines.extend([(header, False), (sequence, True), (plus, False), (qualities, True)])
            self.detail_fastq_length = len(sequence)
        self.detail_lines = lines

    def draw_read_detail(self, screen) -> None:
        height, cols = screen.getmaxyx()
        popup_height = min(height - 2, max(12, height * 4 // 5))
        popup_width = min(cols - 4, max(50, cols * 4 // 5))
        if popup_height < 7 or popup_width < 20:
            return
        top = (height - popup_height) // 2
        left = (cols - popup_width) // 2
        popup = curses.newwin(popup_height, popup_width, top, left)
        popup.erase()
        popup.border()
        body_height = popup_height - 4
        body_width = popup_width - 4
        self.detail_scroll = max(0, min(self.detail_scroll, len(self.detail_lines) - body_height))
        self.detail_x = max(0, min(self.detail_x, max(0, self.detail_fastq_length - body_width)))
        title = f" Read details  [{self.read_detail.track_name}] "
        popup.addnstr(0, 2, title, body_width, self.color(3, curses.A_BOLD))
        for line_number, (line, sequence_line) in enumerate(
                self.detail_lines[self.detail_scroll:self.detail_scroll + body_height], 2):
            shown = line[self.detail_x:] if sequence_line else line
            popup.addnstr(line_number, 2, shown, body_width,
                          self.color(2 if sequence_line else 4))
        footer = " Esc/q close  j/k scroll  h/l FASTQ"
        if self.detail_fastq_length:
            footer += (f"  {self.detail_x + 1}-{min(self.detail_fastq_length, self.detail_x + body_width)}"
                       f"/{self.detail_fastq_length} bp")
        popup.addnstr(popup_height - 2, 2, footer, body_width, self.color(3))
        popup.refresh()

    def mouse(self) -> bool:
        try:
            _, x, y, _, buttons = curses.getmouse()
        except curses.error:
            return True
        if self.read_detail is not None:
            return True
        if buttons & (curses.BUTTON1_CLICKED | curses.BUTTON1_DOUBLE_CLICKED):
            selection = self.read_click_map.get((y, x))
            if selection is not None:
                self.open_read_detail(selection)
        return True

    def prompt(self, screen, prefix: str) -> str | None:
        chars: list[str] = []
        height, cols = screen.getmaxyx()
        curses.curs_set(1)
        try:
            while True:
                line = prefix + "".join(chars)
                self.put(screen, height - 1, 0, line[: cols - 1].ljust(cols - 1), self.color(1))
                screen.move(height - 1, min(len(line), cols - 2))
                screen.refresh()
                key = screen.get_wch()
                if key in ("\n", "\r", curses.KEY_ENTER):
                    return "".join(chars).strip()
                if key in ("\x1b", "\x03"):
                    return None
                if key in (curses.KEY_BACKSPACE, "\x7f", "\b"):
                    if chars:
                        chars.pop()
                elif isinstance(key, str) and key.isprintable():
                    chars.append(key)
        finally:
            curses.curs_set(0)

    def goto(self, region: str, width: int) -> None:
        chrom, start, end = parse_region(region)
        chrom = self.tracks.resolve_chrom(chrom)
        if start is None:
            self.view = View(chrom, 0, min(self.tracks.chrom_size(chrom) or 1000, 1000))
        elif end == start + 1:
            self.view = View(chrom, max(0, start - width // 2), width)
        else:
            padding = max(2, (end - start) // 20)
            self.view = View(chrom, max(0, start - padding), max(width, end - start + padding * 2))
        self.view.clamp(self.tracks.chrom_size(chrom))
        self.read_offsets.clear()
        self.gene_offset = 0
        self.pending_read_reveal = False
        self.search_mode = None
        self.status = f"Moved to {region}"

    def search_gene(self, name: str, width: int, next_hit: bool = False, direction: int = 1) -> None:
        if not next_hit:
            self.gene_hits = self.tracks.find_gene(name)
            self.gene_hit_index = 0
            self.last_search = name
            self.search_mode = "gene"
        elif not self.gene_hits:
            raise ValueError("No previous gene search")
        else:
            self.gene_hit_index = (self.gene_hit_index + direction) % len(self.gene_hits)
        gene = self.gene_hits[self.gene_hit_index]
        chrom = self.tracks.resolve_chrom(gene.chrom)
        pad = max(20, (gene.end - gene.start) // 10)
        self.view = View(chrom, max(0, gene.start - pad), max(width, gene.end - gene.start + 2 * pad))
        self.view.clamp(self.tracks.chrom_size(chrom))
        self.read_offsets.clear()
        self.gene_offset = 0
        self.pending_read_reveal = False
        self.focus_track("genes")
        self.status = f"{gene.name} ({gene.gene_id})  {chrom}:{gene.start + 1:,}-{gene.end:,}  hit {self.gene_hit_index + 1}/{len(self.gene_hits)}"

    def search_read(self, name: str, width: int, progress=None,
                    next_hit: bool = False, direction: int = 1) -> None:
        if not next_hit:
            self.read_hits = self.tracks.find_read(name, progress)
            self.read_hit_index = 0
            self.search_mode = "read"
        elif not self.read_hits:
            raise ValueError("No previous read search")
        else:
            self.read_hit_index = (self.read_hit_index + direction) % len(self.read_hits)
        hit = self.read_hits[self.read_hit_index]
        span = max(width, hit.end - hit.start + 40)
        start = hit.start - (span - (hit.end - hit.start)) // 2
        self.view = View(hit.chrom, max(0, start), span)
        self.view.clamp(self.tracks.chrom_size(hit.chrom))
        self.selected_read_name = hit.name
        self.read_offsets[hit.track_index] = 0
        self.pending_read_reveal = True
        self.focus_track("bam", hit.track_index)
        self.status = (f"Read {hit.name}  {hit.chrom}:{hit.start + 1:,}-{hit.end:,}  "
                       f"BAM {hit.track_index + 1}  hit {self.read_hit_index + 1}/{len(self.read_hits)}")

    def command(self, raw: str, width: int, progress=None) -> bool:
        parts = shlex.split(raw)
        if not parts:
            return True
        verb, *args = parts
        if verb in ("q", "quit", "exit"):
            return False
        if verb in ("help", "h"):
            self.show_help = True
        elif verb == "clear":
            self.selected_read_name = None
            self.read_hits = []
            self.pending_read_reveal = False
            if self.search_mode == "read":
                self.search_mode = None
            self.status = "Read highlight cleared"
        elif verb in ("expand", "collapse"):
            self.gene_expanded = verb == "expand"
            self.gene_offset = 0
            self.status = "Genes expanded" if self.gene_expanded else "Genes packed"
        elif verb in ("bam", "bw", "bigwig", "fa", "fasta", "gtf", "gene", "read", "reads", "readname", "goto", "chr"):
            if not args:
                raise ValueError(f":{verb} needs an argument")
            value = " ".join(args)
            if verb == "bam":
                self.tracks.add_bam(value)
            elif verb in ("bw", "bigwig"):
                self.tracks.add_bigwig(value)
            elif verb in ("fa", "fasta"):
                self.tracks.set_fasta(value)
            elif verb == "gtf":
                self.tracks.set_gtf(value)
            elif verb == "gene":
                self.search_gene(value, width)
            elif verb in ("read", "reads", "readname"):
                self.search_read(value, width, progress)
            else:
                self.goto(value, width)
            if verb not in ("gene", "read", "reads", "readname", "goto", "chr"):
                self.status = f"Loaded {value}"
                if self.view.chrom == "" and self.tracks.chrom_sizes:
                    self.goto(next(iter(self.tracks.chrom_sizes)), width)
                if verb == "bam":
                    self.focus_track("bam", len(self.tracks.bams) - 1)
                elif verb in ("bw", "bigwig"):
                    self.focus_track("bigwig", len(self.tracks.bigwigs) - 1)
                elif verb == "gtf":
                    self.focus_track("genes")
                elif verb in ("fa", "fasta"):
                    self.focus_index = 0
                    self.track_top = 1
        else:
            raise ValueError(f"Unknown command: {verb}")
        return True

    def key(self, screen, key) -> bool:
        height, cols = screen.getmaxyx()
        width = max(1, cols - LABEL - 1)
        if key == curses.KEY_MOUSE:
            return self.mouse()
        if self.read_detail is not None:
            if key in ("\x1b", "q", "\n", "\r", curses.KEY_ENTER):
                self.read_detail = None
                self.detail_lines = []
            elif key in ("j", curses.KEY_DOWN):
                self.detail_scroll += 1
            elif key in ("k", curses.KEY_UP):
                self.detail_scroll = max(0, self.detail_scroll - 1)
            elif key in ("h", curses.KEY_LEFT):
                self.detail_x = max(0, self.detail_x - 10)
            elif key in ("l", curses.KEY_RIGHT):
                self.detail_x += 10
            elif key == "g":
                self.detail_scroll = 0
            elif key == "G":
                self.detail_scroll = len(self.detail_lines)
            return True
        if self.show_help:
            self.show_help = False
            return True
        if key in ("q", "\x03"):
            return False
        if key in (":", "/"):
            answer = self.prompt(screen, key)
            if answer:
                try:
                    def progress(track_index: int, scanned: int) -> None:
                        self.status = f"Searching BAM {track_index + 1}: {scanned:,} reads scanned..."
                        self.put(screen, height - 1, 0, self.status[: cols - 1].ljust(cols - 1), self.color(4))
                        screen.refresh()
                    return self.command(answer, width, progress) if key == ":" else self._search(answer, width, progress)
                except (OSError, ValueError, RuntimeError) as exc:
                    self.status = f"Error: {exc}"
            return True
        size = self.tracks.chrom_size(self.view.chrom)
        if key == "g":
            if self.last_g:
                self.view.start = 0
                self.read_offsets.clear()
                self.last_g = False
                return True
            self.last_g = True
            self.status = "g: press g again for chromosome start"
            return True
        self.last_g = False
        if key in ("h", curses.KEY_LEFT):
            self.view.pan(-0.1, size)
        elif key in ("l", curses.KEY_RIGHT):
            self.view.pan(0.1, size)
        elif key == "H":
            self.view.pan(-0.5, size)
        elif key == "L":
            self.view.pan(0.5, size)
        elif key in ("+", "=", curses.KEY_UP):
            self.view.zoom(0.5, width, size)
        elif key in ("-", "_", curses.KEY_DOWN):
            self.view.zoom(2, width, size)
        elif key == "z":
            center = self.view.start + self.view.span // 2
            self.view.span = width
            self.view.start = center - width // 2
            self.view.clamp(size)
        elif key in ("J", "K"):
            count = 1 + len(self.track_items())
            self.focus_index = max(0, min(count - 1, self.focus_index + (1 if key == "J" else -1)))
            self.track_top = max(1, self.focus_index)
            focus = "Reference" if self.focus_index == 0 else self.track_items()[self.focus_index - 1][2]
            self.status = f"Selected {focus}  ({self.focus_index + 1}/{count})"
        elif key in ("j", "k", "\x04", "\x15"):
            step = 5 if key in ("\x04", "\x15") else 1
            step *= 1 if key in ("j", "\x04") else -1
            if self.focus_index:
                kind, track_index, _, _ = self.track_items()[self.focus_index - 1]
                if kind == "genes":
                    self.gene_offset = max(0, self.gene_offset + step)
                elif kind == "bam":
                    self.read_offsets[track_index] = max(0, self.read_offsets.get(track_index, 0) + step)
        elif key == "e":
            self.gene_expanded = not self.gene_expanded
            self.gene_offset = 0
            self.status = "Genes expanded" if self.gene_expanded else "Genes packed"
        elif key == "G":
            self.view.start = max(0, (size or self.view.end) - self.view.span)
            self.read_offsets.clear()
        elif key in ("[", "]"):
            chroms = list(self.tracks.chrom_sizes)
            if chroms:
                index = chroms.index(self.view.chrom) if self.view.chrom in chroms else 0
                self.view.chrom = chroms[(index + (-1 if key == "[" else 1)) % len(chroms)]
                self.view.start = 0
                self.view.clamp(self.tracks.chrom_size(self.view.chrom))
                self.read_offsets.clear()
        elif key in ("n", "N"):
            try:
                direction = 1 if key == "n" else -1
                if self.search_mode == "read":
                    self.search_read(self.selected_read_name or "", width, next_hit=True, direction=direction)
                elif self.search_mode == "gene":
                    self.search_gene(self.last_search or "", width, next_hit=True, direction=direction)
                else:
                    raise ValueError("No previous gene or read search")
            except ValueError as exc:
                self.status = str(exc)
        elif key == "?":
            self.show_help = True
        elif key == "r":
            self.status = "Redrawn"
        return True

    def _search(self, answer: str, width: int, progress=None) -> bool:
        if answer.startswith("@"):
            if not answer[1:]:
                raise ValueError("Use /@READ_NAME to search reads")
            self.search_read(answer[1:], width, progress)
        elif ":" in answer or answer in self.tracks.chrom_sizes:
            self.goto(answer, width)
        else:
            self.search_gene(answer, width)
        return True

    def run(self, screen) -> None:
        self.has_colors = curses.has_colors()
        if self.has_colors:
            curses.use_default_colors()
            curses.init_pair(1, curses.COLOR_BLACK, curses.COLOR_CYAN)
            curses.init_pair(2, curses.COLOR_GREEN, -1)
            curses.init_pair(3, curses.COLOR_YELLOW, -1)
            curses.init_pair(4, curses.COLOR_WHITE, -1)
            curses.init_pair(5, curses.COLOR_BLUE, -1)
            curses.init_pair(6, curses.COLOR_RED, -1)
            curses.init_pair(7, curses.COLOR_BLACK, curses.COLOR_YELLOW)
            curses.init_pair(8, curses.COLOR_RED, curses.COLOR_YELLOW)
        try:
            curses.curs_set(0)
        except curses.error:
            pass
        screen.keypad(True)
        try:
            _, old_mousemask = curses.mousemask(curses.BUTTON1_CLICKED | curses.BUTTON1_DOUBLE_CLICKED)
        except curses.error:
            old_mousemask = None
        try:
            while True:
                try:
                    self.draw(screen)
                    key = screen.get_wch()
                    if not self.key(screen, key):
                        break
                except curses.error:
                    continue
        finally:
            if old_mousemask is not None:
                curses.mousemask(old_mousemask)


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Terminal genome browser with Vim navigation")
    parser.add_argument("--bam", action="append", default=[], help="Indexed BAM (repeatable)")
    parser.add_argument("--bw", "--bigwig", action="append", default=[], help="bigWig (repeatable)")
    parser.add_argument("--fa", "--fasta", help="Indexed reference FASTA")
    parser.add_argument("--gtf", help="GTF or GTF.gz gene annotation")
    parser.add_argument("--region", help="Initial region, e.g. chr1:1000-2000")
    parser.add_argument("--gene", help="Initial gene name; requires --gtf")
    parser.add_argument("--read", help="Initial read name; scans loaded BAM files")
    args = parser.parse_args(argv)
    tracks = Tracks()
    try:
        if args.fa:
            tracks.set_fasta(args.fa)
        if args.gtf:
            tracks.set_gtf(args.gtf)
        for path in args.bam:
            tracks.add_bam(path)
        for path in args.bw:
            tracks.add_bigwig(path)
        chrom = next(iter(tracks.chrom_sizes), "")
        view = View(chrom, 0, min(tracks.chrom_size(chrom) or 1000, 1000))
        browser = Browser(tracks, view)
        if args.region:
            browser.goto(args.region, 80)
        if args.gene:
            browser.search_gene(args.gene, 80)
        if args.read:
            def report(track_index: int, scanned: int) -> None:
                print(f"Searching BAM {track_index + 1}: {scanned:,} reads scanned...", end="\r", file=sys.stderr, flush=True)
            browser.search_read(args.read, 80, report)
        curses.wrapper(browser.run)
    except (OSError, ValueError, RuntimeError) as exc:
        parser.exit(2, f"tigv: {exc}\n")
    finally:
        tracks.close()


if __name__ == "__main__":
    main()
