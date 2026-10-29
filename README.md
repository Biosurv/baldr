# BALDR

BALDR assigns individual sequencing reads from single amplicon protocols to known lineages using lineage-defining single nucleotide variants (SNVs) from a Freyja-style barcode file. It operates on a single BAM or a directory of BAMs and produces per-lineage read sets, optional BAM/FASTQ subsets, and per-sample mixture summaries.

The current intended use is single amplicon sequencing with single-contig references.

**This tool is still in active development**

---

## Installation

From source:

```bash
git clone https://github.com/Biosurv/baldr.git
cd baldr
conda env create -f environment.yml
conda activate baldr
pip install .
```

## Quick start

```bash
baldr \
  --bam-dir alignments \
  --barcode-csv barcodes/barcodes_lins.csv \
  --outdir out_example \
  --write-bams \
  --write-mix-summary
```

This runs assignment on every `*.bam` in `alignments/`, writing per-lineage outputs and a mixture summary into `out_example/<sample>/`.

---

## Inputs

### Barcode table (`--barcode-csv`)

Wide-format CSV (`.csv`) or TSV. The first column is the lineage name (header may be `lineage`, `lin`, or `name`). Remaining columns are markers like `A123G`, `C456T`, `G7890A`. Cell value `1` means the lineage carries that SNV; anything else means it does not.

```
lineage,A23403G,C14408T,G28883C
A,1,1,0
B,0,1,1
C,0,1,0
```

### BAM (`--bam` or `--bam-dir`)

Indexed BAMs aligned to the same single-contig reference the barcode positions refer to. If a `.bai` is missing BALDR will attempt to index in place.

---

## Method

Each read's barcode-position observations are scored against every lineage as a log-likelihood (using base qualities for per-base error). An EM loop estimates mixture weights `π`, posteriors `γ` are computed per read, and a read is hard-assigned to a lineage when `γ ≥ --post-min` (default `0.90`), otherwise marked ambiguous.

Lineages with no defining mutations (rows that are all `0` in the barcode) are reference-equivalent and are scored like any other lineage.

---

## Key parameters

| Flag | Default | Notes |
|---|---|---|
| `--bam` / `--bam-dir` | — | Mutually exclusive; one required. |
| `--barcode-csv` | — | Required. CSV or TSV. |
| `--lineage-include` | all | Optional whitelist of lineage names. |
| `--mapq-min` | `0` | Minimum MAPQ to consider a read. |
| `--baseq-min` | `7` | Minimum base quality to count a barcode site. |
| `--allow-secondary` | off | Include secondary/supplementary alignments. |
| `--post-min` | `0.90` | EM posterior threshold for hard assignment. |
| `--em-max-iters` | `50` | EM iteration cap. |
| `--em-tol` | `1e-4` | EM convergence tolerance on `π`. |
| `--min-sites` | `0` | Per-read floor on barcode positions covered. |
| `--write-bams` | off | Per-lineage BAMs (indexed). |
| `--write-fastq` | off | Per-lineage FASTQs. |
| `--write-unassigned` | off | Also emit BAM/FASTQ for ambiguous reads. |
| `--write-mix-summary` | off | Per-sample mix summary TSV. |
| `--threads` | `1` | Threads for pysam BAM I/O. |

Run `baldr --help` for the full list.