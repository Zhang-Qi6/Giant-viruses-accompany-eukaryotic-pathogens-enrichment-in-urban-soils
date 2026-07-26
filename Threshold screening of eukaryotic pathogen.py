import pandas as pd
import numpy as np
import glob
import os
import sys
import re

# ==========================================
# 1. Configuration & Parameters
# ==========================================
RESULTS_DIR = "/public2/home/zq1/大病毒NC返修/真核致病菌分析/DIAMOND_Annotation"  
OUTPUT_FILE = "High_Confidence_Pathogen_Profiles.csv"

# Rigorous Biological Quality Control Thresholds
MIN_BUSCO_COUNT = 5
MAX_CV = 1.1

print("🚀 Starting the ultimate metagenomic competitive parsing engine...")
print(f"🔒 Locked Thresholds: Distinct Marker Genes >= {MIN_BUSCO_COUNT}, Abundance CV <= {MAX_CV}\n")

all_sample_results = []
tsv_files = glob.glob(os.path.join(RESULTS_DIR, "*_diamond.tsv"))

if not tsv_files:
    print(f"❌ No alignment result files found in {RESULTS_DIR}.")
    sys.exit(1)

print(f"📂 Found {len(tsv_files)} DIAMOND result files. Processing...")

# Regex pattern to extract Genome name and BUSCO ID
# It looks for anything before the BUSCO ID (ignoring optional "_Genome_") and captures the digits+at+digits pattern.
REGEX_PATTERN = r'^(?P<Genome>.*?)_?(?:Genome_)?(?P<BUSCO_ID>\d+at\d+)'

# ==========================================
# 2. Batch Processing for Each Metagenomic Sample
# ==========================================
for file in tsv_files:
    sample_name = os.path.basename(file).replace("_diamond.tsv", "")
    print(f"🧬 Parsing sample: {sample_name}")
    
    try:
        # We only need the first two columns (qseqid and sseqid) for counting. 
        # Using usecols=[0,1] prevents errors if DIAMOND output has unexpected extra columns.
        df = pd.read_csv(file, sep='\t', header=None, usecols=[0, 1], names=['qseqid', 'sseqid'])
        
        # Ensure one vote per read
        df = df.drop_duplicates(subset=['qseqid'], keep='first')
    except pd.errors.EmptyDataError:
        print(f"  ⚠️ Results for {sample_name} are empty, skipping.")
        continue
    except Exception as e:
        print(f"  ⚠️ Error reading file for {sample_name}: {e}. Skipping.")
        continue

    if df.empty:
        continue

    # ----------------------------------------------------
    # ⚠️ CRITICAL STEP: Robust Regex Parsing of sseqid
    # ----------------------------------------------------
    try:
        # Extract Genome and BUSCO_ID using the compiled regex pattern
        extracted_data = df['sseqid'].str.extract(REGEX_PATTERN)
        
        # Merge the extracted columns back to the dataframe
        df['Genome'] = extracted_data['Genome']
        df['BUSCO_ID'] = extracted_data['BUSCO_ID']
        
        # Drop rows where the regex failed to find a valid BUSCO ID
        initial_count = len(df)
        df = df.dropna(subset=['Genome', 'BUSCO_ID'])
        dropped_count = initial_count - len(df)
        
        if dropped_count > 0:
            print(f"    ℹ️ Notice: Dropped {dropped_count} alignments with unrecognizable FASTA headers.")
            
    except Exception as e:
        print(f"❌ Critical error applying Regex parsing in sample {sample_name}: {e}")
        sys.exit(1)

    if df.empty:
        print(f"  ⚠️ No valid sequences remained in {sample_name} after header parsing. Skipping.")
        continue

    # Count the abundance (number of reads) for each distinct BUSCO within each genome
    busco_counts = df.groupby(['Genome', 'BUSCO_ID']).size().reset_index(name='Read_Count')

    # Calculate dual metrics: Number of distinct BUSCOs and the Coefficient of Variation (CV)
    def calc_metrics(x):
        counts = x['Read_Count']
        mean_val = counts.mean()
        std_val = counts.std(ddof=1) if len(counts) > 1 else 0
        cv_val = std_val / mean_val if mean_val > 0 else 0
        
        return pd.Series({
            'Detected_BUSCOs': len(counts),
            'Mean_Read_Abundance': mean_val,
            'Abundance_CV': cv_val
        })

    genome_stats = busco_counts.groupby('Genome').apply(calc_metrics).reset_index()

    # ----------------------------------------------------
    # Apply Strict Dual-Filter Defense
    # ----------------------------------------------------
    high_confidence_genomes = genome_stats[
        (genome_stats['Detected_BUSCOs'] >= MIN_BUSCO_COUNT) & 
        (genome_stats['Abundance_CV'] <= MAX_CV)
    ].copy()

    # Append sample name if valid taxa are found
    if not high_confidence_genomes.empty:
        high_confidence_genomes.insert(0, 'Sample', sample_name)
        all_sample_results.append(high_confidence_genomes)

# ==========================================
# 3. Aggregation & Output
# ==========================================
if all_sample_results:
    final_df = pd.concat(all_sample_results, ignore_index=True)
    
    # Sort by sample name and number of detected BUSCOs
    final_df = final_df.sort_values(by=['Sample', 'Detected_BUSCOs'], ascending=[True, False])
    
    # Export results to the current working directory
    final_df.to_csv(OUTPUT_FILE, index=False)
    
    print("\n" + "="*50)
    print(" ✅ Parsing complete! Environmental noise has been rigorously filtered out.")
    print(f" 📊 High-confidence pathogenic profiles saved to: {OUTPUT_FILE}")
    print(f" 🛡️ Total high-confidence taxa retained across samples: {len(final_df)}")
    print("="*50)
else:
    print("\n⚠️ Parsing complete, but no taxa survived the stringent thresholds.")