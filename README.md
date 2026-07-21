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

This runs the default EM mode on every `*.bam` in `alignments/`, writing per-lineage outputs and a mixture summary into `out_example/<sample>/`.

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

## Modes

### EM (default)

Each read's barcode-position observations are scored against every lineage as a log-likelihood (using base qualities for per-base error). An EM loop estimates mixture weights `π`, posteriors `γ` are computed per read, and a read is hard-assigned to a lineage when `γ ≥ --post-min` (default `0.90`), otherwise marked ambiguous.

### Voting (`--voting`)

For each read, every covered barcode site contributes `+1` to lineages whose expected ALT matches the read base, and `-1` (unless `--ignore-ref-mismatch`) to lineages whose expected ALT does not. The lineage with the highest score wins if it beats the runner-up by at least `--min-margin`.

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
| `--voting` | off | Use voting instead of EM. |
| `--post-min` | `0.90` | EM posterior threshold for hard assignment. |
| `--em-max-iters` | `50` | EM iteration cap. |
| `--em-tol` | `1e-4` | EM convergence tolerance on `π`. |
| `--min-sites` | `0` | EM: per-read floor on barcode positions covered. Voting: per-lineage floor on the winning lineage's defining sites covered. |
| `--min-margin` | `1` | Voting only. Score margin over second-best lineage. |
| `--ignore-ref-mismatch` | off | Voting only. Treat reference-base reads as 0 instead of −1. |
| `--write-bams` | off | Per-lineage BAMs (indexed). |
| `--write-fastq` | off | Per-lineage FASTQs. |
| `--write-unassigned` | off | Also emit BAM/FASTQ for ambiguous reads. |
| `--write-mix-summary` | off | Per-sample mix summary TSV. |
| `--threads` | `1` | Threads for pysam BAM I/O. |

Run `baldr --help` for the full list.

---

## Outputs

In `--bam` mode, files are written to `--outdir`. In `--bam-dir` mode, each sample gets its own subdirectory `<outdir>/<sample>/`.

| File | Always | Contents |
|---|---|---|
| `<prefix>.<lineage>.names.txt` | yes | Read names assigned to that lineage. |
| `<prefix>.ambiguous.names.txt` | if any | Read names that failed assignment. |
| `<prefix>.summary.txt` | yes | Counts: total/processed/assigned/ambiguous reads, skip reasons, per-lineage read counts. |
| `<prefix>.<lineage>.bam` (+ `.bai`) | `--write-bams` | Per-lineage BAM, indexed. |
| `<prefix>.<lineage>.fastq` | `--write-fastq` | Per-lineage FASTQ. |
| `<prefix>.ambiguous.bam` / `.fastq` | `--write-unassigned` | Same, for ambiguous reads. |
| `<prefix>.mix_summary.tsv` | `--write-mix-summary` | Per-lineage counts and mixture estimate. Columns: `lineage`, `reads`, `mix_weight`, `frac_assigned`, `frac_total`. |

`reads`, `frac_assigned`, and `frac_total` are hard counts in both modes. `mix_weight` is the EM mixture weight `π[lin]` (sums to ~1 across lineages, includes evidence from sub-threshold reads).
