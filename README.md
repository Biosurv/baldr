# BALDR

**BALDR** — *Barcode-based Assignment of Lineages by Demixing Reads*

BALDR assigns individual sequencing reads to known lineages using lineage-defining single nucleotide variants (SNVs) from a Freyja-style barcode file.  
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
