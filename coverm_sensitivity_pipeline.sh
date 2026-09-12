#!/bin/bash
# ==============================================================================
# Script: coverm_sensitivity_pipeline.sh
# Description: Automated CoverM mapping pipeline for threshold-sensitivity 
#              assessment across metagenomic datasets. Evaluates genome coverage 
#              and abundance at Stringent, Moderate, and Relaxed stringency levels.
# ==============================================================================

set -e # Exit immediately if a command exits with a non-zero status

# ==============================================================================
# 1. USER CONFIGURATION (Edit these paths before running)
# ==============================================================================
# Resource Allocation
export CONCURRENT_JOBS=6    # Number of samples to process simultaneously
export COVERM_THREADS=30    # Number of threads per CoverM job

# Input/Output Paths
export LIST_FILE="./path/to/your_target_genomes_list.txt"
export SOURCE_GENOME_DIR="./path/to/source_genomes_database"
export WORK_DIR="./CoverM_Sensitivity_Analysis"

# Metagenomic Data Directories (Add all directories containing your FASTQ files)
VALID_DIRS=(
    "/path/to/metagenomes_dir_1"
    "/path/to/metagenomes_dir_2"
)

# ==============================================================================
# 2. ENVIRONMENT SETUP
# ==============================================================================
export TARGET_GENOME_DIR="${WORK_DIR}/Target_Genomes_Symlinks"
export OUT_STRINGENT="${WORK_DIR}/01_Stringent_Results"
export OUT_MODERATE="${WORK_DIR}/02_Moderate_Results"
export OUT_RELAXED="${WORK_DIR}/03_Relaxed_Results"

mkdir -p "$TARGET_GENOME_DIR"
mkdir -p "$OUT_STRINGENT" "$OUT_MODERATE" "$OUT_RELAXED"

# ==============================================================================
# 3. DATA PREPARATION: Construct Custom Reference Library via Symlinks
# ==============================================================================
echo "[INFO] ========== Preparing target genome reference library =========="
# Clear existing symlinks to avoid conflicts from previous runs
rm -f "${TARGET_GENOME_DIR}"/*

# Read genome list and create symbolic links
while IFS= read -r raw_genome_id; do
    # Remove carriage returns and trailing spaces
    genome_id=$(echo "$raw_genome_id" | tr -d '\r' | xargs)
    
    [ -z "$genome_id" ] && continue 
    
    # Flexible matching for various FASTA extensions and naming conventions
    find "$SOURCE_GENOME_DIR" -maxdepth 1 -type f \( \
        -name "${genome_id}" -o \
        -name "${genome_id}.fa" -o \
        -name "${genome_id}.fasta" -o \
        -name "${genome_id}.fna" -o \
        -name "${genome_id}_Genome.fa" -o \
        -name "${genome_id}_Genome.fasta" -o \
        -name "${genome_id}_Genome.fna" \
    \) | while read -r match_file; do
        ln -s "$match_file" "${TARGET_GENOME_DIR}/$(basename "$match_file")"
    done
done < "$LIST_FILE"

genome_count=$(ls -1q "$TARGET_GENOME_DIR" | wc -l)
echo "[INFO] Successfully linked $genome_count genomes to $TARGET_GENOME_DIR"

if [ "$genome_count" -eq 0 ]; then
    echo "[ERROR] No matching genomes found. Please verify the IDs in your list file and the source directory path."
    exit 1
fi

# ==============================================================================
# 4. CORE FUNCTION: Multi-gradient CoverM Mapping
# ==============================================================================
run_coverm() {
    local fq1="$1"
    local fq2=""
    local sample_name=""

    # Infer R2 file and extract sample name based on common sequencing naming conventions
    if [[ "$fq1" == *".1.fastq" ]]; then fq2="${fq1%.1.fastq}.2.fastq"; sample_name=$(basename "${fq1%.1.fastq}");
    elif [[ "$fq1" == *"_1.fastq" ]]; then fq2="${fq1%_1.fastq}_2.fastq"; sample_name=$(basename "${fq1%_1.fastq}");
    elif [[ "$fq1" == *"_R1.fq.gz" ]]; then fq2="${fq1%_R1.fq.gz}_R2.fq.gz"; sample_name=$(basename "${fq1%_R1.fq.gz}");
    elif [[ "$fq1" == *"_1.fastq.gz" ]]; then fq2="${fq1%_1.fastq.gz}_2.fastq.gz"; sample_name=$(basename "${fq1%_1.fastq.gz}");
    else return; fi

    # Skip if R2 does not exist
    [ ! -f "$fq2" ] && return

    echo "[INFO] Processing sample: $sample_name"

    # Define base CoverM command with requested metrics
    local base_cmd="coverm genome --coupled $fq1 $fq2 --genome-fasta-directory $TARGET_GENOME_DIR -x fasta --methods mean covered_fraction tpm relative_abundance count --min-covered-fraction 0 --threads $COVERM_THREADS"
    
    # --- Gradient 1: Stringent (95% Identity, 80% Aligned) ---
    local out_1="${OUT_STRINGENT}/${sample_name}_coverage.tsv"
    if [ ! -f "$out_1" ]; then
        $base_cmd --min-read-percent-identity 95 --min-read-aligned-percent 80 --output-file "$out_1"
    fi

    # --- Gradient 2: Moderate (85% Identity, 70% Aligned) ---
    local out_2="${OUT_MODERATE}/${sample_name}_coverage.tsv"
    if [ ! -f "$out_2" ]; then
        $base_cmd --min-read-percent-identity 85 --min-read-aligned-percent 70 --output-file "$out_2"
    fi

    # --- Gradient 3: Relaxed (75% Identity, 50% Aligned) ---
    local out_3="${OUT_RELAXED}/${sample_name}_coverage.tsv"
    if [ ! -f "$out_3" ]; then
        $base_cmd --min-read-percent-identity 75 --min-read-aligned-percent 50 --output-file "$out_3"
    fi
}

# Export function and variables for xargs parallelization
export -f run_coverm
export TARGET_GENOME_DIR OUT_STRINGENT OUT_MODERATE OUT_RELAXED COVERM_THREADS

# ==============================================================================
# 5. PARALLEL EXECUTION
# ==============================================================================
echo "[INFO] ========== Initiating Multi-Gradient CoverM Analysis =========="

# Find all R1 fastq files in valid directories and pipe to xargs for parallel processing
find -L "${VALID_DIRS[@]}" -type f \
| grep -E "(\.1\.fastq|_1\.fastq|_R1\.fq\.gz|_1\.fastq\.gz)$" \
| tr '\n' '\0' \
| xargs -0 -n 1 -P "$CONCURRENT_JOBS" -I {} bash -c 'run_coverm "{}"'

echo "[INFO] ========== All samples processed successfully! =========="