# PhyloGCNE: Phylogenetic Graph Convolutional Network for Microbial Data Analysis

PhyloGCNE is a deep learning framework that integrates phylogenetic information with Graph Convolutional Network with edge features (GCNE) for microbial data analysis.

## Overview

This project implements a phylogenetic-aware GCNE model for analyzing microbial abundance data from various microbiome phenotypes including:
- Synthetic data for method validation
- Colorectal Cancer (CRC) - 16S and WGS data
- Inflammatory Bowel Disease (IBD) - 16S data
- Type 2 Diabetes (T2D) - WGS data
- Gastritis - 16S data
- Oral Squamous Cell Carcinoma (OSCC) - 16S data
- Dietary Fiber intervention studies
- Multi-classification tasks

## Project Structure

```
PhyloGCNE/
├── GCN/
│   └── GCN_model.py           # Main GCN model implementation
├── script/                     # Analysis scripts
│   ├── RF.ipynb               # Random Forest baseline
│   ├── Multclass.ipynb        # Multi-classification analysis
│   ├── plot_auc.R             # AUC curve visualization
│   ├── plot_ablation.R        # Ablation study visualization
│   ├── plot_biomarker.R       # Biomarker importance plots
│   ├── plot_biomarker_tax.R   # Taxonomy-level biomarker plots
│   └── network_analysis.R     # Network analysis
├── data/                       # Dataset directories
│   ├── CRC_16S/               # Colorectal Cancer 16S data
│   ├── CRC_WGS/               # Colorectal Cancer WGS data
│   ├── IBD_16S/               # Inflammatory Bowel Disease 16S data
│   ├── T2D_WGS/               # Type 2 Diabetes WGS data
│   ├── Gastritis_16S/         # Gastritis 16S data
│   ├── OSCC_16S/              # Oral Squamous Cell Carcinoma 16S 
│   ├── dietary_fiber/         # Dietary fiber intervention study
│   ├── Multi_classification/  # Multi-class classification data
│   └── synthetic_data/        # Synthetic data for validation
├── results/                   # Analysis results
├── requirements.txt
├── INSTALL.md                 # Installation guide
└── INSTALL_TORCHGEOMETRIC.md  # PyTorch Geometric setup guide
```

Each dataset directory typically contains:
```
data/[dataset]/
├── data/                              # Train/test split BIOM files (train_1.biom ... test_5.biom)
├── output/                            # PhyloSpec output
├── results/                           # PhyloGCNE results
├── ablation_phylo/                    # Phylogeny ablation results
(some datasets)
├── phylogeny.nwk                      # Full phylogenetic tree
├── phylogeny_ablation.nwk             # Tree for phylogeny ablation
DeepPhylo
├── PCA_32.txt                         # PCA coordinates for phylogenetic distance 
run model
├── run_GCN.sh                         # PhyloGCNE execution script
├── run_GCN_ablation.sh                # Ablation study script
├── run_DeepPhylo.sh                   # DeepPhylo comparison script
└── run_phylospec.sh                   # PhyloSpec comparison script
```

## Features

1. **Phylogenetic Integration**: Incorporates phylogenetic tree structure and branch lengths into GCN architecture via Gaussian kernel edge weighting. Edge attributes are constructed via a fixed Gaussian kernel transform of phylogenetic distances ( $e_{uv} = \exp(-2\rho d^2_{uv}),\ \rho=2$ ), which are subsequently projected into high-dimensional embeddings via a learnable linear layer ( $W_\text{edge} \in \mathbb{R}^{1 \times d}$ ) within the message-passing framework.
2. **Multiple Datasets**: Supports 16S and WGS microbial datasets across 7+ disease phenotypes
3. **Multi-classification**: Handles multi-class classification
4. **Ablation Studies**: Systematic analysis of phylogenetic information importance (topology and branch-length ablation)
5. **Benchmarking**: Comparison with DeepPhylo and PhyloSpec
6. **Model Interpretation**: Feature importance for biomarker discovery

## Requirements

- Python 3.8+
- PyTorch 1.9+
- PyTorch Geometric 2.0+
- NumPy, SciPy, Pandas
- scikit-learn
- biom-format
- scikit-bio
- Jupyter Notebook
- R (for visualization: ggplot2, pROC, dplyr, tidyr, reshape2)

Install Python dependencies:
```bash
pip install -r requirements.txt
```

See [INSTALL.md](INSTALL.md) and [INSTALL_TORCHGEOMETRIC.md](INSTALL_TORCHGEOMETRIC.md) for detailed setup instructions.

## Usage

### 1. Run PhyloGCNE Model

```bash
cd data/CRC_16S
bash run_GCN.sh
```

Key model arguments:

| Argument | Description | Default |
|---|---|---|
| `--train_table` | Path to training BIOM file | required |
| `--test_table` | Path to test BIOM file | required |
| `--metadata_filename` | Path to metadata TSV file | required |
| `--phylogeny_file_path` | Path to Newick phylogenetic tree | required |
| `--fold` | Fold/study name for output file naming | `PRJDB11845` |
| `--output_dir` | Root output directory | `./graph_model_crc_workdir/v19` |
| `--label_column` | Label column name in metadata | `group` |
| `--hidden_channels` | Hidden units per GCN layer | `100` |
| `--num_layers` | Number of GCN layers | `4` |
| `--act` | Activation function (`relu`, `gelu`, etc.) | `relu` |
| `--ffn` | Enable FFN block in each GCN layer | `True` |
| `--residual` | Enable residual connections | `True` |
| `--norm_type` | Normalization type (`batch_norm`/`layer_norm`) | `batch_norm` |
| `--readout` | Graph readout (`mean`/`max`/`mean_max`/`attention`) | `mean` |
| `--dropout_rate` | Dropout rate | `0.0` |
| `--epochs` | Training epochs | `100` |
| `--batch_size` | Training batch size | `32` |
| `--learning_rate` | Learning rate | `1e-3` |
| `--weight_decay` | Weight decay | `1e-5` |
| `--optimizer` | Optimizer (`adamw`/`adam`) | `adamw` |
| `--lr_scheduler_type` | LR scheduler (`linear`/`cosine`/`none`) | `none` |
| `--lr_warmup_ratio` | Warmup ratio for LR scheduler | `0.0` |
| `--max_grad_norm` | Max gradient norm for clipping | `1.0` |
| `--seed` | Random seed | `42` |
| `--device` | Device (`cuda`/`cpu`) | `cuda` |
| `--workers` | DataLoader worker processes | `4` |
| `--feature_importance` | Compute feature importance instead of training | `False` |
| `--stdev_spread` | Noise std scale for feature importance | `0.15` |
| `--beta_alpha` | Beta distribution α for feature importance noise | `2` |
| `--beta_beta` | Beta distribution β for feature importance noise | `5` |

Example:
```bash
python GCN/GCN_model.py \
    --train_table data/CRC_16S/data/train_1.biom \
    --test_table data/CRC_16S/data/test_1.biom \
    --metadata_filename data/CRC_16S/metadata.tsv \
    --phylogeny_file_path data/CRC_16S/phylogeny.nwk \
    --fold 1 --output_dir results \
    --batch_size 32 --learning_rate 0.001 \
    --num_layers 4 --ffn True \
    --norm_type batch_norm --residual True \
    --hidden_channels 100 --dropout_rate 0.1 \
    --weight_decay 0.00001
```

### 2. Run Ablation Studies

```bash
cd data/CRC_16S
bash run_GCN_ablation.sh
```

Ablation variants use modified tree files to assess the contribution of phylogenetic topology and branch lengths.

### 3. Run Baseline Methods

**DeepPhylo**:
```bash
bash run_DeepPhylo.sh  # requires DeepPhylo installation
```

**PhyloSpec**:
```bash
bash run_phylospec.sh  # requires PhyloSpec installation
```

**Random Forest**:
```bash
jupyter notebook script/RF.ipynb
```

### 4. Visualization

```bash
# AUC comparison plots
Rscript script/plot_auc.R

# Ablation study results
Rscript script/plot_ablation.R

# Biomarker importance
Rscript script/plot_biomarker.R
Rscript script/plot_biomarker_tax.R
```

## Citation

If you use PhyloGCNE in your research, please cite our work.

## Issues

Please report bugs and feature requests at [GitHub Issues](https://github.com/delect-aj/PhyloGCNE/issues).
