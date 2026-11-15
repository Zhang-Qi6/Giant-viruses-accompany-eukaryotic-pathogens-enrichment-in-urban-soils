
#!/usr/bin/env bash

set -Eeuo pipefail

# Batch hmmsearch for scanning protein FASTA files against an HMM database.
# Requirements: HMMER, GNU find, GNU sort, GNU xargs, Bash >= 4.4

usage() {
    echo "Usage: $0 -i <protein_dir> -d <hmm_db> -o <output_dir> [-p jobs] [-t threads] [-e evalue] [-f]"
    echo
    echo "Options:"
    echo "  -i   Directory containing protein FASTA files (*.faa)"
    echo "  -d   HMM database file"
    echo "  -o   Output directory"
    echo "  -p   Number of parallel jobs [default: 8]"
    echo "  -t   Threads per hmmsearch job [default: 4]"
    echo "  -e   Sequence E-value cutoff [default: 1e-10]"
    echo "  -f   Force overwrite existing results"
    echo "  -h   Show this help message"
}

CONCURRENT_JOBS=8
THREADS_PER_JOB=4
EVALUE="1e-10"
FORCE=0

while getopts ":i:d:o:p:t:e:fh" opt; do
    case "$opt" in
        i) PROT_DIR="$OPTARG" ;;
        d) HMM_DB="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        p) CONCURRENT_JOBS="$OPTARG" ;;
        t) THREADS_PER_JOB="$OPTARG" ;;
        e) EVALUE="$OPTARG" ;;
        f) FORCE=1 ;;
        h) usage; exit 0 ;;
        *) usage; exit 1 ;;
    esac
done

# Check required arguments
if [[ -z "${PROT_DIR:-}" || -z "${HMM_DB:-}" || -z "${OUTPUT_DIR:-}" ]]; then
    usage
    exit 1
fi

# Check numeric parameters
if ! [[ "$CONCURRENT_JOBS" =~ ^[1-9][0-9]*$ ]] ||
   ! [[ "$THREADS_PER_JOB" =~ ^[1-9][0-9]*$ ]]; then
    echo "Error: jobs and threads must be positive integers." >&2
    exit 1
fi

if ! [[ "$EVALUE" =~ ^([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$ ]]; then
    echo "Error: invalid E-value: $EVALUE" >&2
    exit 1
fi

# Check HMMER
if ! command -v hmmsearch >/dev/null 2>&1; then
    echo "Error: hmmsearch was not found." >&2
    exit 1
fi

# Check input paths
if [[ ! -d "$PROT_DIR" ]]; then
    echo "Error: protein directory does not exist: $PROT_DIR" >&2
    exit 1
fi

if [[ ! -f "$HMM_DB" || ! -r "$HMM_DB" ]]; then
    echo "Error: HMM database does not exist or is unreadable: $HMM_DB" >&2
    exit 1
fi

# Convert paths to absolute paths
PROT_DIR="$(cd "$PROT_DIR" && pwd -P)"
HMM_DB="$(cd "$(dirname "$HMM_DB")" && pwd -P)/$(basename "$HMM_DB")"

mkdir -p "$OUTPUT_DIR"
OUTPUT_DIR="$(cd "$OUTPUT_DIR" && pwd -P)"

export PROT_DIR HMM_DB OUTPUT_DIR
export THREADS_PER_JOB EVALUE FORCE

run_hmmsearch_single() {

    local faa_file="$1"
    local rel_path
    local out_tbl
    local tmp_tbl
    local hit_count

    # Preserve relative directory structure to prevent filename collisions
    rel_path="${faa_file#"$PROT_DIR"/}"
    rel_path="${rel_path%.faa}"

    out_tbl="${OUTPUT_DIR}/${rel_path}_hmm_hits.tbl"

    mkdir -p "$(dirname "$out_tbl")"

    # Skip completed results unless force mode is enabled
    if [[ "$FORCE" -eq 0 && -s "$out_tbl" ]]; then
        echo "Skipped existing result: $rel_path"
        return 0
    fi

    # Generate a unique temporary file
    tmp_tbl="$(mktemp "${out_tbl}.tmp.XXXXXX")"
    rm -f "$tmp_tbl"

    # Run HMMER
    if ! hmmsearch \
        --cpu "$THREADS_PER_JOB" \
        -E "$EVALUE" \
        --noali \
        --tblout "$tmp_tbl" \
        "$HMM_DB" \
        "$faa_file" \
        >/dev/null
    then
        rm -f "$tmp_tbl"
        echo "ERROR: hmmsearch failed for $faa_file" >&2
        return 1
    fi

    # Move completed result into final output directory
    mv -f "$tmp_tbl" "$out_tbl"

    # Count reported sequence hits
    hit_count=$(awk '!/^#/ && NF {count++} END {print count+0}' "$out_tbl")

    if [[ "$hit_count" -gt 0 ]]; then
        echo "Hits found: $rel_path ($hit_count records)"
    else
        echo "No hits: $rel_path"
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
echo "Force overwrite:   $FORCE"
echo

# Collect input files safely and sort by file size (largest first)
mapfile -d '' -t INPUT_RECORDS < <(
    find "$PROT_DIR" -type f -name "*.faa" \
        -printf '%s\t%p\0' | sort -z -nr
)

if [[ "${#INPUT_RECORDS[@]}" -eq 0 ]]; then
    echo "Error: no .faa files found in $PROT_DIR" >&2
    exit 1
fi

FAA_FILES=()

for record in "${INPUT_RECORDS[@]}"; do
    FAA_FILES+=("${record#*$'\t'}")
done

# Parallel execution with strict error handling in child shells
printf '%s\0' "${FAA_FILES[@]}" |
    xargs -0 -r -I {} -P "$CONCURRENT_JOBS" \
    bash -euo pipefail -c 'run_hmmsearch_single "$1"' _ {}

echo
echo "Batch hmmsearch completed successfully."
