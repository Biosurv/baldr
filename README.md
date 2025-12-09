# BALDR

BALDR assigns individual sequencing reads from a single amplicon sequencing protocol to known lineages using lineage-defining single nucleotide variants (SNVs) from a Freyja-style barcode file.  
It can operate on a single BAM file or an entire directory of BAMs, generating per-lineage read sets, BAM/FASTQ subsets, and mixture summaries.

---

## Installation

You can install BALDR from source:

```bash
git clone https://github.com/Biosurv/baldr.git
cd baldr
conda env create -f environment.yml
conda activate baldr
pip install .
```

## Example usage

```bash
baldr --bam-dir alignments --barcode-csv barcodes/barcodes_lins.csv --min-sites 1 --min-margin 2 --outdir out_example --write-bams --write-mix-summary
