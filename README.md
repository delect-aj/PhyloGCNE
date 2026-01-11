# PhyloGCNE: Phylogenetic Graph Convolutional Network for Microbial Data Analysis

PhyloGCNE is a deep learning framework that integrates phylogenetic information with Graph Convolutional Network with edge features (GCNE) for microbial data analysis.

## Overview

This project implements a phylogenetic-aware GCNE model for analyzing microbial abundance data from various microbiome phenotype including:
- Synthetic data for method validation
- Colorectal Cancer (CRC) - 16S and WGS data
- Inflammatory Bowel Disease (IBD) - 16S data
- Type 2 Diabetes (T2D) - WGS data
- Dietary Fiber intervention studies

## Project Structure

```
PhyloGCNE/
├── GCN/                    # GCN model implementation
│   └── GCN_model.py       # Main GCN model implementation
├── script/                 # Analysis scripts
│   ├── get_datasets.ipynb # Data preprocessing notebook
│   ├── RF.ipynb           # Random Forest baseline
│   ├── Model_explain.ipynb # Model interpretation
│   └── *.R                # R scripts for visualization
├── data/                   # Dataset configurations
│   ├── CRC_16S/           # CRC 16S data
│   ├── CRC_WGS/           # CRC WGS data
│   ├── IBD_16S/           # IBD 16S data
│   ├── T2D_WGS/           # T2D WGS data
│   ├── dietary_fiber/     # Dietary fiber study
│   ├── synthetic_data/    # Synthetic data for validation
|   └── Multi_classification # multi-status data
└── results/               # Analysis results and scripts
```

## Features

1. **Phylogenetic Integration**: Incorporates phylogenetic tree information into GCN architecture
2. **Multiple Datasets**: Support for various microbial datasets (16S, WGS, synthetic)
3. **Ablation Studies**: Analysis of phylogenetic information importance
4. **Benchmarking**: Comparison with traditional methods (Random Forest)
5. **Model Interpretation**: Feature importance and biomarker discovery

## Requirements

- Python 3.8+
- PyTorch 1.9+
- NumPy, SciPy, Pandas
- scikit-learn
- Jupyter Notebook
- R (for visualization)

## Usage

### 1. Data Preparation
```bash
# Run data preprocessing
jupyter notebook script/get_datasets.ipynb
```

### 2. Run GCN Model
```bash
# Example for CRC 16S data
cd data/CRC_16S
bash run_GCN.sh
```

### 3. Run Ablation Studies
```bash
# Run phylogenetic ablation study
bash run_GCN_ablation.sh
```

### 4. Run Baseline Methods
```bash
# Run Random Forest baseline
jupyter notebook script/RF.ipynb
```

### 5. Visualization
```bash
# Generate AUC plots
Rscript script/plot_auc.R
```

## Configuration Files

Each dataset directory contains:
- `run_GCN.sh`: Main GCN execution script
- `run_GCN_ablation.sh`: Ablation study script
- `run_DeepPhylo.sh`: DeepPhylo comparison
- `run_phylospec.sh`: PhyloSpec comparison
- `PCA_32.txt`: PCA coordinates for DeepPhylo model
- Phylogenetic tree files (`.nwk`)

## Results

Key results include:
- AUC performance comparisons
- Ablation study results
- Biomarker discovery
- Model interpretation visualizations
