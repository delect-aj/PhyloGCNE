# Installation and Setup Guide

## Prerequisites

- Python 3.8 or higher
- Git
- (Optional) Conda for environment management

## Installation Methods

### Method 1: Using Conda (Recommended)

```bash
# Clone the repository
git clone https://github.com/delect-aj/PhyloGCNE.git
cd PhyloGCNE

# Create and activate conda environment
conda create -n phylogcne python=3.9
conda activate phylogcne

# Install dependencies
pip install -r requirements.txt
```

### Method 2: Using pip

```bash
# Clone the repository
git clone https://github.com/delect-aj/PhyloGCNE.git
cd PhyloGCNE

# Create virtual environment (optional but recommended)
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
```

## R Dependencies (for visualization)

If you want to use the R visualization scripts:

```r
# Install required R packages
install.packages(c("ggplot2", "pROC", "dplyr", "tidyr", "reshape2"))
```

## Data Preparation

1. Place your BIOM format data files in the appropriate `data/[dataset]/` directory
2. Ensure phylogenetic tree files (`.nwk`) are in the same directory
3. Update configuration scripts if needed

## Quick Start

1. **Prepare environment**:
   ```bash
   conda activate phylogcne
   ```

2. **Run data preprocessing**:
   ```bash
   jupyter notebook script/get_datasets.ipynb
   ```

3. **Run GCN model on example data**:
   ```bash
   cd data/synthetic_data
   bash run_GCN.sh
   ```

## Testing Installation

To verify your installation is working:

```python
python -c "import torch; import numpy; print('PyTorch version:', torch.__version__); print('NumPy version:', numpy.__version__)"
```

## Troubleshooting

### Common Issues

1. **PyTorch installation fails**:
   - Visit https://pytorch.org/get-started/locally/ for platform-specific instructions
   - Try: `pip install torch torchvision torchaudio`

2. **BIOM format issues**:
   - Ensure biom-format is installed: `pip install biom-format`
   - Check BIOM file version compatibility

3. **R script errors**:
   - Ensure R is installed and in PATH
   - Install required R packages as shown above

### Getting Help

If you encounter issues:
1. Check the [GitHub Issues](https://github.com/delect-aj/PhyloGCNE/issues)
2. Review the README.md for additional information
3. Ensure all dependencies are correctly installed
