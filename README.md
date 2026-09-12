# Giant Viruses Accompany Eukaryotic Pathogen Enrichment in Urban Soils

This repository contains analysis scripts used to identify urban-associated eukaryotic pathogen signatures, detect large-virus-associated markers, and explore ecological drivers of urban interface pathogens (UIPs) in soil and soil-fauna metagenomes.

## Overview

The workflow includes seven major analytical modules:

1. **Targeted annotation of eukaryotic pathogen markers**
2. **Threshold screening of eukaryotic pathogen signals**
3. **Identification of large-virus-associated markers**
4. **Ecological driver analysis of UIP abundance**
5. **Gene Co-localization Null Model Analysis**
6. **Global visualization of UIP distribution against Anthromes background**
7. **Threshold-sensitivity assessment of metagenomic read mapping**

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
| `coverm_sensitivity_pipeline.sh` | Automated CoverM mapping pipeline for threshold-sensitivity assessment across metagenomic datasets. |

## Requirements

### Command-line tools

- [DIAMOND](https://github.com/bbuchfink/diamond)
- [HMMER](http://hmmer.org/)
- [CoverM](https://github.com/wwood/CoverM)

### Python packages

The Python scripts require the following packages:

```bash
pip install pandas numpy matplotlib scikit-learn pygam matplotlib-venn rasterio openpyxl
