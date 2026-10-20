"""Small in-memory gene index for GTF and GFF-like attribute fields."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from urllib.parse import unquote
import gzip
import re


@dataclass
class Gene:
    chrom: str
    start: int  # zero-based, inclusive
    end: int  # zero-based, exclusive
    gene_id: str
    name: str
    strand: str = "."
    exons: list[tuple[int, int]] = field(default_factory=list)

    def merged_exons(self) -> list[tuple[int, int]]:
        intervals = sorted(self.exons)
        merged: list[list[int]] = []
        for start, end in intervals:
            if merged and start <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], end)
            else:
                merged.append([start, end])
        return [(start, end) for start, end in merged]


def parse_attributes(raw: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for part in raw.strip().rstrip(";").split(";"):
        part = part.strip()
        if not part:
            continue
        match = re.match(r'^([^\s=]+)\s+(?:"([^"]*)"|([^\s]+))$', part)
        if match:
            result[match.group(1)] = match.group(2) if match.group(2) is not None else match.group(3)
        elif "=" in part:
            key, value = part.split("=", 1)
            result[key] = unquote(value)
    return result


class GeneIndex:
    BIN_SIZE = 100_000

    def __init__(self, genes: list[Gene]):
        self.genes = genes
        self.by_name: dict[str, list[Gene]] = {}
        self.bins: dict[tuple[str, int], list[Gene]] = {}
        self.chrom_sizes: dict[str, int] = {}
        for gene in genes:
            for key in {gene.name.casefold(), gene.gene_id.casefold()}:
                self.by_name.setdefault(key, []).append(gene)
            self.chrom_sizes[gene.chrom] = max(self.chrom_sizes.get(gene.chrom, 0), gene.end)
            for bin_no in range(gene.start // self.BIN_SIZE, (gene.end - 1) // self.BIN_SIZE + 1):
                self.bins.setdefault((gene.chrom, bin_no), []).append(gene)

    @classmethod
    def from_file(cls, path: str | Path) -> "GeneIndex":
        genes: dict[tuple[str, str], Gene] = {}
        opener = gzip.open if str(path).endswith(".gz") else open
        with opener(path, "rt", encoding="utf-8") as handle:
            for line_no, line in enumerate(handle, 1):
                if line.startswith("#") or not line.strip():
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) != 9 or fields[2].lower() not in {"gene", "exon"}:
                    continue
                chrom, _, feature, raw_start, raw_end, _, strand, _, raw_attrs = fields
                try:
                    start, end = int(raw_start) - 1, int(raw_end)
                except ValueError as exc:
                    raise ValueError(f"{path}:{line_no}: invalid coordinates") from exc
                if start < 0 or end <= start:
                    continue
                attrs = parse_attributes(raw_attrs)
                gene_id = attrs.get("gene_id") or attrs.get("gene") or attrs.get("ID")
                if not gene_id:
                    continue
                gene_name = attrs.get("gene_name") or attrs.get("Name") or gene_id
                key = (chrom, gene_id)
                gene = genes.get(key)
                if gene is None:
                    gene = Gene(chrom, start, end, gene_id, gene_name, strand)
                    genes[key] = gene
                else:
                    gene.start = min(gene.start, start)
                    gene.end = max(gene.end, end)
                    if gene.name == gene.gene_id and gene_name != gene_id:
                        gene.name = gene_name
                if feature.lower() == "exon":
                    gene.exons.append((start, end))
        return cls(list(genes.values()))

    def search(self, query: str) -> list[Gene]:
        return self.by_name.get(query.casefold(), [])

    def suggest(self, query: str, limit: int = 8) -> list[str]:
        needle = query.casefold()
        names = {gene.name for gene in self.genes if gene.name.casefold().startswith(needle)}
        return sorted(names)[:limit]

    def overlap(self, chrom: str, start: int, end: int) -> list[Gene]:
        seen: set[str] = set()
        found: list[Gene] = []
        for bin_no in range(start // self.BIN_SIZE, (end - 1) // self.BIN_SIZE + 1):
            for gene in self.bins.get((chrom, bin_no), ()):
                key = gene.gene_id
                if key not in seen and gene.end > start and gene.start < end:
                    seen.add(key)
                    found.append(gene)
        return sorted(found, key=lambda gene: (gene.start, gene.end))
