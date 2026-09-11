#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import argparse
import time
import matplotlib.pyplot as plt
import seaborn as sns

def count_colocalization(df, window, target_types):
    total_count = 0
    for scaffold, group in df.groupby('Molecule'):
        stress_pos = group[group['Gene_Type'] == target_types['stress']]['position'].values
        ncldv_pos = group[group['Gene_Type'] == target_types['ncldv']]['position'].values
        mobile_pos = group[group['Gene_Type'] == target_types['mobile']]['position'].values
        
        if len(ncldv_pos) == 0 or len(mobile_pos) == 0:
            continue
            
        for pos in stress_pos:
            near_ncldv = np.any(np.abs(ncldv_pos - pos) <= window)
            near_mobile = np.any(np.abs(mobile_pos - pos) <= window)
            
            if near_ncldv and near_mobile:
                total_count += 1
    return total_count

def run_null_model(df, window, permutations, target_types):
    print(f"\n[{window//1000} kb Window] Calculating...")
    observed_count = count_colocalization(df, window, target_types)
    print(f"-> Observed co-localization count: {observed_count}")
    
    null_counts = []
    start_time = time.time()
    
    for i in range(permutations):
        shuffled_df = df.copy()
        # Randomly permute gene types within each molecule (scaffold/contig) while preserving spatial architecture
        shuffled_df['Gene_Type'] = shuffled_df.groupby('Molecule')['Gene_Type'].transform(lambda x: np.random.permutation(x.values))
        null_counts.append(count_colocalization(shuffled_df, window, target_types))
            
    end_time = time.time()
    null_counts = np.array(null_counts)
    exceed_count = np.sum(null_counts >= observed_count)
    empirical_p = (1 + exceed_count) / (1 + permutations)
    
    print(f"-> P-value = {empirical_p:.4f} (Time elapsed: {end_time - start_time:.1f} s)")
    return observed_count, null_counts, empirical_p

def plot_combined_distributions(results_dict):
    print("\nGenerating multi-window null model distribution plots...")
    sns.set_theme(style="ticks", font_scale=1.1)
    
    num_windows = len(results_dict)
    # Create a dynamic 1-row subplot layout based on the number of windows
    fig, axes = plt.subplots(1, num_windows, figsize=(6 * num_windows, 5.5))
    
    if num_windows == 1:
        axes = [axes]
        
    for ax, (window, data) in zip(axes, results_dict.items()):
        obs = data['obs']
        nulls = data['nulls']
        p_val = data['p']
        
        max_val = max(np.max(nulls), obs)
        bins = np.arange(-0.5, max_val + 1.5, 1)
        
        sns.histplot(nulls, bins=bins, color='lightgray', edgecolor='black', 
                     stat='density', alpha=0.8, ax=ax)
        
        ax.axvline(x=obs, color='red', linestyle='--', linewidth=2.5, 
                    label=f'Observed = {obs}')
        
        # Annotate P-value
        p_text = f'Empirical P = {p_val:.4f}'
        ax.text(0.95, 0.85, p_text, 
                 verticalalignment='top', horizontalalignment='right',
                 transform=ax.transAxes,
                 fontsize=11, fontweight='bold',
                 bbox=dict(facecolor='white', alpha=0.9, edgecolor='gray', boxstyle='round,pad=0.5'))
        
        ax.set_title(f'Null Model ({window//1000} kb Window)', fontsize=14, fontweight='bold', pad=15)
        ax.set_xlabel('Number of Co-localized Loci', fontsize=12)
        if ax == axes[0]:
            ax.set_ylabel('Density', fontsize=12)
        else:
            ax.set_ylabel('')
        
        # Add legend with an opaque background to prevent visual overlap
        ax.legend(loc='upper right', frameon=True, facecolor='white', edgecolor='gray', framealpha=1.0)
        sns.despine(ax=ax)
        
    plt.tight_layout()
    
    pdf_filename = 'Null_Model_Comparison_MultiWindow.pdf'
    png_filename = 'Null_Model_Comparison_MultiWindow.png'
    plt.savefig(pdf_filename, format='pdf', bbox_inches='tight')
    plt.savefig(png_filename, dpi=300, bbox_inches='tight')
    print(f"Comparison plots successfully saved as: {pdf_filename} and {png_filename}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="GVMAG Null Model Analysis (Multi-Window)")
    parser.add_argument("-i", "--input", required=True, help="Path to the input dataset (TSV format)")
    parser.add_argument("-p", "--permutations", type=int, default=1000, help="Number of random permutations (default: 1000)")
    parser.add_argument("--mobile_name", type=str, default="Mobile_Element", help="Identifier for mobile elements in the Gene_Type column")
    
    args = parser.parse_args()
    print(f"Loading data from: {args.input} ...")
    df = pd.read_csv(args.input, sep='\t')
    df['position'] = (df['Start'] + df['End']) / 2.0
    
    target_types = {
        'stress': 'Urban_Function', 
        'ncldv': 'NCLDV_Marker', 
        'mobile': args.mobile_name
    }
    
    # Core modification: Define the three spatial window sizes to test (5kb, 10kb, 20kb)
    windows_to_test = [5000, 10000, 20000]
    results = {}
    
    print(f"Running sensitivity analysis across {len(windows_to_test)} spatial windows...")
    for w in windows_to_test:
        obs_cnt, null_cnts, emp_p = run_null_model(df, window=w, permutations=args.permutations, target_types=target_types)
        results[w] = {'obs': obs_cnt, 'nulls': null_cnts, 'p': emp_p}
        
    plot_combined_distributions(results)