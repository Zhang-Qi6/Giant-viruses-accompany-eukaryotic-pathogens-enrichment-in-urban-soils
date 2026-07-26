#!/usr/bin/env bash

# Batch DIAMOND blastx annotation for metagenomic reads.
# This script builds a DIAMOND database from a protein FASTA file and scans
# R1 metagenomic reads against the database in parallel.
#
# Requirements:
#   DIAMOND >= 2.0
#   GNU find, grep, xargs

set -euo pipefail

usage() {
    echo "Usage: $0 -f <protein_fasta> -d <db_prefix> -o <output_dir> -i <input_dir1,input_dir2,...> [options]"
    echo
    echo "Required arguments:"
    echo "  -f   Protein FASTA file used to build the DIAMOND database"
    echo "  -d   DIAMOND database prefix"
    echo "  -o   Output directory"
    echo "  -i   Comma-separated input directories containing FASTQ files"
    echo
    echo "Optional arguments:"
    echo "  -p   Number of samples processed in parallel [default: 4]"
    echo "  -t   Threads per DIAMOND job [default: 8]"
    echo "  -e   E-value cutoff [default: 1e-5]"
    echo "  -h   Show this help message"
}

CONCURRENT_JOBS=4
THREADS_PER_JOB=8
EVALUE="1e-5"

while getopts ":f:d:o:i:p:t:e:h" opt; do
    case "$opt" in
        f) INPUT_FAA="$OPTARG" ;;
        d) DB_PREFIX="$OPTARG" ;;
        o) OUT_DIR="$OPTARG" ;;
        i) INPUT_DIRS="$OPTARG" ;;
        p) CONCURRENT_JOBS="$OPTARG" ;;
        t) THREADS_PER_JOB="$OPTARG" ;;
        e) EVALUE="$OPTARG" ;;
        h) usage; exit 0 ;;
        *) usage; exit 1 ;;
    esac
done

if [[ -z "${INPUT_FAA:-}" || -z "${DB_PREFIX:-}" || -z "${OUT_DIR:-}" || -z "${INPUT_DIRS:-}" ]]; then
    usage
    exit 1
fi

if ! command -v diamond >/dev/null 2>&1; then
    echo "Error: DIAMOND was not found in PATH." >&2
    exit 1
fi

if [[ ! -f "$INPUT_FAA" ]]; then
    echo "Error: protein FASTA file not found: $INPUT_FAA" >&2
    exit 1
fi

mkdir -p "$OUT_DIR"

export DB_PREFIX
export OUT_DIR
export THREADS_PER_JOB
export EVALUE

echo "Starting batch DIAMOND annotation"
echo "Protein FASTA:      $INPUT_FAA"
echo "DIAMOND DB prefix:  $DB_PREFIX"
echo "Output directory:   $OUT_DIR"
echo "Parallel jobs:      $CONCURRENT_JOBS"
echo "Threads per job:    $THREADS_PER_JOB"
echo "E-value cutoff:     $EVALUE"
echo

# Build DIAMOND database if it does not already exist
if [[ ! -f "${DB_PREFIX}.dmnd" ]]; then
    echo "Building DIAMOND database..."
    diamond makedb \
        --in "$INPUT_FAA" \
        -d "$DB_PREFIX" \
        --threads "$THREADS_PER_JOB"
else
    echo "Existing DIAMOND database detected: ${DB_PREFIX}.dmnd"
fi

run_diamond_single() {
    local read_file="$1"
    local filename
    local sample_name
    local output_tsv

    filename=$(basename "$read_file")

    case "$filename" in
        *.1.fastq)
            sample_name="${filename%.1.fastq}"
            ;;
        *_1.fastq)
            sample_name="${filename%_1.fastq}"
            ;;
        *_R1.fq.gz)
            sample_name="${filename%_R1.fq.gz}"
            ;;
        *_1.fastq.gz)
            sample_name="${filename%_1.fastq.gz}"
            ;;
        *)
            return 0
            ;;
    esac

    output_tsv="${OUT_DIR}/${sample_name}_diamond.tsv"

    if [[ -s "$output_tsv" ]]; then
        echo "Skipped existing result: $sample_name"
        return 0
    fi

    echo "Running DIAMOND blastx: $sample_name"

    diamond blastx \
        -d "$DB_PREFIX" \
        -q "$read_file" \
        -o "$output_tsv" \
        --threads "$THREADS_PER_JOB" \
        --evalue "$EVALUE" \
        --max-target-seqs 1 \
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp \
        --quiet
}

export -f run_diamond_single

# Convert comma-separated directories to an array
IFS=',' read -r -a DIR_ARRAY <<< "$INPUT_DIRS"

VALID_DIRS=()
for dir in "${DIR_ARRAY[@]}"; do
    if [[ -d "$dir" ]]; then
        VALID_DIRS+=("$dir")
    else
        echo "Warning: input directory not found: $dir" >&2
    fi
done

if [[ "${#VALID_DIRS[@]}" -eq 0 ]]; then
    echo "Error: no valid input directories found." >&2
    exit 1
fi

echo
echo "Searching FASTQ files and running DIAMOND..."
find -L "${VALID_DIRS[@]}" -type f | \
    grep -E '(\.1\.fastq|_1\.fastq|_R1\.fq\.gz|_1\.fastq\.gz)$' | \
    sort | \
    xargs -r -n 1 -P "$CONCURRENT_JOBS" bash -c 'run_diamond_single "$1"' _

echo
echo "Batch DIAMOND annotation completed."
echo "Results saved to: $OUT_DIR"