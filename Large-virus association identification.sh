#!/usr/bin/env bash

# Batch hmmsearch for scanning protein FASTA files against an HMM database.
# Requirements: HMMER (hmmsearch), GNU find, GNU sort, xargs

set -euo pipefail

usage() {
    echo "Usage: $0 -i <protein_dir> -d <hmm_db> -o <output_dir> [-p concurrent_jobs] [-t threads_per_job] [-e evalue]"
    echo
    echo "Options:"
    echo "  -i   Directory containing protein FASTA files (*.faa)"
    echo "  -d   HMM database file"
    echo "  -o   Output directory"
    echo "  -p   Number of files processed in parallel [default: 8]"
    echo "  -t   CPU threads per hmmsearch job [default: 4]"
    echo "  -e   E-value cutoff for hmmsearch [default: 1e-5]"
    echo "  -h   Show this help message"
}

CONCURRENT_JOBS=8
THREADS_PER_JOB=4
EVALUE="1e-5"

while getopts ":i:d:o:p:t:e:h" opt; do
    case "$opt" in
        i) PROT_DIR="$OPTARG" ;;
        d) HMM_DB="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        p) CONCURRENT_JOBS="$OPTARG" ;;
        t) THREADS_PER_JOB="$OPTARG" ;;
        e) EVALUE="$OPTARG" ;;
        h) usage; exit 0 ;;
        *) usage; exit 1 ;;
    esac
done

if [[ -z "${PROT_DIR:-}" || -z "${HMM_DB:-}" || -z "${OUTPUT_DIR:-}" ]]; then
    usage
    exit 1
fi

if ! command -v hmmsearch >/dev/null 2>&1; then
    echo "Error: hmmsearch was not found. Please install HMMER first." >&2
    exit 1
fi

if [[ ! -d "$PROT_DIR" ]]; then
    echo "Error: protein directory does not exist: $PROT_DIR" >&2
    exit 1
fi

if [[ ! -f "$HMM_DB" ]]; then
    echo "Error: HMM database does not exist: $HMM_DB" >&2
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

export HMM_DB
export OUTPUT_DIR
export THREADS_PER_JOB
export EVALUE

run_hmmsearch_single() {
    local faa_file="$1"
    local filename
    local prefix
    local out_tbl
    local tmp_tbl
    local hit_count

    filename=$(basename "$faa_file")
    prefix="${filename%_protein.faa}"
    prefix="${prefix%.faa}"

    out_tbl="${OUTPUT_DIR}/${prefix}_hmm_hits.tbl"
    tmp_tbl="${out_tbl}.tmp"

    if [[ -s "$out_tbl" ]]; then
        echo "Skipped existing result: $prefix"
        return 0
    fi

    hmmsearch \
        --cpu "$THREADS_PER_JOB" \
        -E "$EVALUE" \
        --tblout "$tmp_tbl" \
        "$HMM_DB" \
        "$faa_file" \
        >/dev/null

    mv "$tmp_tbl" "$out_tbl"

    hit_count=$(awk '!/^#/ && NF {count++} END {print count+0}' "$out_tbl")

    if [[ "$hit_count" -gt 0 ]]; then
        echo "Hits found: $prefix ($hit_count hits)"
    else
        echo "No hits: $prefix"
    fi
}

export -f run_hmmsearch_single

echo "Starting batch hmmsearch"
echo "Protein directory: $PROT_DIR"
echo "HMM database:      $HMM_DB"
echo "Output directory:  $OUTPUT_DIR"
echo "Parallel jobs:     $CONCURRENT_JOBS"
echo "Threads per job:   $THREADS_PER_JOB"
echo "E-value cutoff:    $EVALUE"
echo

find "$PROT_DIR" -type f -name "*.faa" -printf "%s %p\n" | \
    sort -nr | \
    cut -d' ' -f2- | \
    xargs -I {} -P "$CONCURRENT_JOBS" bash -c 'run_hmmsearch_single "$1"' _ {}

echo
echo "Batch hmmsearch completed."