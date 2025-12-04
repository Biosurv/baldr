#!/usr/bin/env bash
set -euo pipefail

abspath() {
    case "$1" in
        /*)  printf '%s\n' "$1" ;;
        *)   printf '%s\n' "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")" ;;
    esac
}

if [[ $# -lt 3 ]]; then
    echo "Usage: $0 /path/to/barcode_directory /path/to/reference.fasta MIN_MAPPED_READS"
    exit 1
fi

BARCODE_ROOT="$(abspath "$1")"
REF="$(abspath "$2")"
MIN_MAPPED_READS="$3"

if ! [[ "$MIN_MAPPED_READS" =~ ^[0-9]+$ ]]; then
    echo "Error: MIN_MAPPED_READS must be an integer, got: $MIN_MAPPED_READS"
    exit 1
fi

if [[ ! -d "$BARCODE_ROOT" ]]; then
    echo "Error: Barcode directory does not exist: $BARCODE_ROOT"
    exit 1
fi

if [[ ! -f "$REF" ]]; then
    echo "Error: Reference fasta not found: $REF"
    exit 1
fi

for cmd in bcftools bgzip samtools; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "Error: $cmd not found in PATH."
        exit 1
    fi
done

echo "Barcode directory: $BARCODE_ROOT"
echo "Reference: $REF"
echo "MIN_MAPPED_READS: $MIN_MAPPED_READS (primary mapped only, no secondary/supplementary)"
echo

for barcode_dir in "$BARCODE_ROOT"/barcode*/; do
    [[ -d "$barcode_dir" ]] || continue

    echo "Processing $barcode_dir"
    cd "$barcode_dir"

    mkdir -p consensus

    bam_files=( *.bam )

    if [[ ${#bam_files[@]} -eq 1 && "${bam_files[0]}" == "*.bam" ]]; then
        echo "  No BAM files found, skipping."
        cd - >/dev/null
        continue
    fi

    for bam in "${bam_files[@]}"; do
        echo "  BAM: $bam"

        mapped_reads=$(samtools view -c -F 2308 "$bam" || echo 0)

        if (( mapped_reads < MIN_MAPPED_READS )); then
            echo "  Skipping $bam: only $mapped_reads mapped reads (<$MIN_MAPPED_READS)."
            continue
        fi

        base="${bam%.bam}"
        vcf_raw="${base}.vcf"
        vcf_gz="${vcf_raw}.gz"
        consensus_fa="${base}.consensus.fasta"

        bcftools mpileup -Ou -f "$REF" "$bam" \
          | bcftools call --ploidy 1 -mv -Ou \
          | bcftools filter -e 'DP<20 || QUAL<20' -s LowQual -Ov -o "$vcf_raw"

        bgzip -f "$vcf_raw"
        bcftools index -f --csi "$vcf_gz"
        bcftools consensus -f "$REF" -s - -o "consensus/$consensus_fa" "$vcf_gz"

        echo "    → consensus/$consensus_fa"
    done

    cd - >/dev/null
done

echo "Finished processing all barcode folders."