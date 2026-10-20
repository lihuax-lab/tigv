# tigv

`tigv` is a terminal genome browser with Vim-style keys. It shows indexed BAM
alignments and coverage, bigWig signal, a reference sequence, and gene models
from GTF. At single-base zoom, aligned bases are printed on their reference
positions; mismatches are red. Insertions (`+`) and deletions (`-`) appear in
red when the complete event fits in the visible region.

## Install

Create the separate `tigv` conda environment and install the command:

```bash
conda env create -f environment.yml
conda run -n tigv python -m pip install -e . --no-deps
```

Alternatively, after installing the dependencies in `environment.yml`, run
`python -m tigv` from the repository directory.

## Input files

BAM files need a `.bai` or `.csi` index. Reference FASTA files need a `.fai`
index. Create them when needed:

```bash
samtools index sample.bam
samtools faidx reference.fa
```

Chromosome names should agree across files (`chr1` versus `1`). The coordinate
prompt accepts either form if one exists in the loaded files; tracks with other
chromosome names stay empty. GTF and GTF.gz are supported. GTF is read into a
gene index at startup, so large annotations can take a moment to load.

```bash
conda run -n tigv tigv \
  --bam sample.bam --bw signal.bw --fa reference.fa --gtf genes.gtf \
  --gene TP53
```

You can repeat `--bam` and `--bw` for multiple tracks. Start with any subset of
files, then load more from the command prompt. Quote paths containing spaces:

```text
:bam "/data/sample 2.bam"
:bw /data/signal.bw
:fa /data/reference.fa
:gtf /data/genes.gtf.gz
```

## Navigation

| Key | Action |
| --- | --- |
| `h` / `l` | Pan left / right by 10% |
| `H` / `L` | Pan left / right by 50% |
| `+` / `-` | Zoom in / out around the center |
| `z` | Zoom to one base per terminal column |
| `gg` / `G` | Go to chromosome start / end |
| `[` / `]` | Previous / next chromosome |
| `j` / `k` | Scroll rows within the selected BAM or gene track |
| `J` / `K` | Select next / previous track, including the reference row |
| `/TP53` | Search for a gene by name or gene ID |
| `/chr12:123-234` | Jump to a one-based chromosome range (also accepts `/chr12:123` or `/chr12`) |
| `/@read_name` or `:read read_name` | Find mapped reads by exact name, jump to and highlight the hit |
| Mouse click on a read | Open a detail window with alignment fields, BAM tags, and reconstructed FASTQ |
| `n` / `N` | Next / previous result of the most recent gene or read search |
| `:goto chr1:1000-2000` | Jump to a region (one-based, inclusive) |
| `:gene TP53` | Search for a gene |
| `e`, `:expand`, `:collapse` | Toggle, expand, or pack gene rows |
| `:clear` | Clear the read highlight |
| `?` | Show help |
| `q` | Quit |

The reference and read bases appear at single-base zoom. In a wider view,
matches become `#` and mismatches become `*`; insertion and deletion markers
are hidden. The `+` and `-` markers appear only at base resolution when the
complete event is visible. Coverage and bigWig tracks aggregate into terminal
columns. Coverage uses up to 12,000
primary reads per visible region; the panel reports when that limit is reached.
The bigWig plot scales each track to its visible range; negative signal is red.

Alignments with the same read name stay on one display row. Nearby aligned
pieces, including segments separated by CIGAR `N`, are joined with `·`:

```text
####········####
```

Distance does not change this symbol. If the two pieces map to **different
annotated genes**, a cross-gene chRNA junction is shown as `[     ]`.
If only one end is visible, a mate coordinate or `SA` tag can mark its
outward-facing bracket. Load a GTF to identify cross-gene junctions;
without annotation, pieces are joined with `·`.

`J` and `K` move the highlighted track selector. The reference row is the first
selectable row. On a selected BAM or gene track, `j` and `k` scroll the reads or
genes in that track. `e` expands genes so each occupies its own row; if they do
not all fit, select the gene track and use `j`/`k` to browse the remaining rows.
Read-name search scans each BAM file on the first query because BAM indexes are
by genomic position; repeated queries for the same name use cached hits. Use
`--read read_name` at launch to search immediately.

Click a read base or its `·` / `[` / `]` connector to open a bordered detail
window. It shows the clicked alignment and other visible alignments with the
same read name, including CIGAR, mapping quality, flags, and BAM tags. The
FASTQ sequence and qualities are reconstructed in the original read direction
from BAM. Use `j` / `k` to scroll the window, `h` / `l` to move across long FASTQ
lines, and `Esc`, `Enter`, or `q` to close it. BAM cannot restore an original
FASTQ header or hard-clipped bases; if sequence or qualities are missing, the
window states why a complete FASTQ cannot be shown.
