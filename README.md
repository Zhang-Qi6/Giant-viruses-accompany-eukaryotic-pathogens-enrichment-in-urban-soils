# Giant Viruses Accompany Eukaryotic Pathogen Enrichment in Urban Soils

This repository contains analysis scripts used to identify urban-associated eukaryotic pathogen signatures, detect large-virus-associated markers, and explore ecological drivers of urban interface pathogens (UIPs) in soil and soil-fauna metagenomes.

## Overview

The workflow includes five major analytical modules:

1. **Targeted annotation of eukaryotic pathogen markers**
2. **Threshold screening of eukaryotic pathogen signals**
3. **Identification of large-virus-associated markers**
4. **Ecological driver analysis of UIP abundance**
5. **Gene Co-localization Null Model Analysis**
6. **Global visualization of UIP distribution against Anthromes background**

These scripts were developed for reproducible analysis of metagenomic datasets and associated environmental metadata.

## Repository structure

| Script | Description |
|---|---|
| `eukaryotic pathogen_annotation.sh` | Performs targeted DIAMOND `blastx` annotation of metagenomic reads against a custom eukaryotic pathogen BUSCO protein database. |
| `Threshold screening of eukaryotic pathogen.py` | Screens eukaryotic pathogen signals using predefined detection and consistency thresholds. |
| `Large-virus association identification.sh` | Performs batch HMMER searches to identify large-virus-associated hallmark markers in predicted protein sequences. |
| `UIP-GAM.py` | Analyzes ecological drivers of UIP abundance using variance partitioning analysis, PLS-VIP, and GAM models. |
| `Null-Model.py` | Automatically generates publication-ready density plots comparing observed data to the null distribution. |
| `UIP-global-map.py` | Visualizes global UIP abundance on an Anthromes raster background and exports source data for mapping. |

## Requirements

### Command-line tools

- [DIAMOND](https://github.com/bbuchfink/diamond)
- [HMMER](http://hmmer.org/)

### Python packages

The Python scripts require the following packages:

```bash
pip install pandas numpy matplotlib scikit-learn pygam matplotlib-venn rasterio openpyxl

## Database Download

The custom competitive protein database (`Combined_Competitive_DB.faa`) required for the DIAMOND annotation step is hosted in the GitHub Releases section of this repository due to its large size.

You can download and extract the database (v1.0) using the following commands before running `eukaryotic pathogen_annotation.sh`:

```bash
# Download the database archive
wget https://github.com/Zhang-Qi6/Giant-viruses-accompany-eukaryotic-pathogens-enrichment-in-urban-soils/commits/v1.0

# Extract the fasta file
unzip Combined_Competitive_DB.zip

## Gene Co-localization Null Model Analysis

To evaluate the significance of gene clustering within the metagenomic assemblies—specifically the spatial proximity of urban stressor genes, Nucleocytoviricota (NCLDV) markers, and mobile genetic elements (MGEs)—we utilize a permutation-based null model. This script (`gvmag_null_model.py`) tests whether the observed co-localization frequencies across different proximity windows indicate potential horizontal gene transfer (HGT) events rather than random gene ordering.

### Features
* **Spatial Permutation:** Shuffles gene identities within individual contigs/scaffolds while preserving the underlying spatial architecture of the genomic fragments.
* **Sensitivity Analysis:** Simultaneously evaluates multiple proximity windows (5 kb, 10 kb, and 20 kb) to test the robustness of the physical linkage.
* **Empirical P-value Calculation:** Compares the observed number of tripartite co-localizations (Stress + NCLDV + MGE) against a simulated background distribution.
* **Visualization:** Automatically generates publication-ready density plots comparing observed data to the null distributions.

### Usage

```bash
python gvmag_null_model.py -i <annotated_genes.tsv> -p 1000 --mobile_name "Mobile_Element"
