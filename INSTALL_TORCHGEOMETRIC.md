# PyTorch Geometric Installation Guide

PyTorch Geometric (torch-geometric) has specific installation requirements depending on your PyTorch and CUDA versions.

## Prerequisites

1. **PyTorch**: Must be installed first
2. **CUDA**: If using GPU (optional but recommended for training)
3. **System dependencies**: Some backends require system packages

## Installation Methods

### Method 1: Using pip (Simplest)

For most users, this method works well:

```bash
# First install PyTorch (choose appropriate version for your system)
# Visit https://pytorch.org/get-started/locally/ for PyTorch installation

# Then install torch-geometric and its dependencies
pip install torch-geometric
```

### Method 2: Specific version matching

If you need specific versions:

```bash
# Check your PyTorch and CUDA versions
python -c "import torch; print(f'PyTorch: {torch.__version__}'); print(f'CUDA: {torch.version.cuda}')"

# Install matching torch-geometric
pip install torch-scatter torch-sparse torch-cluster torch-spline-conv -f https://data.pyg.org/whl/torch-${TORCH}+${CUDA}.html
pip install torch-geometric
```

Replace `${TORCH}` and `${CUDA}` with your versions (e.g., `torch-1.13.0+cu117`).

### Method 3: Using conda (Recommended for complex environments)

```bash
conda install pyg -c pyg
```

## Common Issues and Solutions

### Issue 1: "No matching distribution found"
- Ensure your PyTorch version is compatible
- Check Python version (requires Python 3.7+)

### Issue 2: CUDA version mismatch
- Install PyTorch with matching CUDA version first
- Use `torch.cuda.is_available()` to check GPU availability

### Issue 3: Missing system dependencies
On Ubuntu/Debian:
```bash
sudo apt-get update
sudo apt-get install -y libopenblas-dev libomp-dev
```

On macOS:
```bash
brew install libomp
```

## Verification

After installation, verify with:

```python
import torch
import torch_geometric

print(f"PyTorch version: {torch.__version__}")
print(f"PyTorch Geometric version: {torch_geometric.__version__}")
print(f"CUDA available: {torch.cuda.is_available()}")

# Test basic functionality
from torch_geometric.data import Data
edge_index = torch.tensor([[0, 1, 1, 2], [1, 0, 2, 1]], dtype=torch.long)
x = torch.tensor([[-1], [0], [1]], dtype=torch.float)
data = Data(x=x, edge_index=edge_index)
print(f"Test graph created: {data}")
```

## For This Project

This project requires:
- `torch-geometric>=2.0.0`
- Compatible PyTorch version (>=1.9.0 recommended)

If you encounter issues, try:
```bash
# Minimal installation
pip install torch==1.13.0 torchvision==0.14.0 torchaudio==0.13.0
pip install torch-geometric==2.0.0
```

## References

- [PyTorch Geometric Documentation](https://pytorch-geometric.readthedocs.io/)
- [PyTorch Geometric GitHub](https://github.com/pyg-team/pytorch_geometric)
- [PyTorch Installation Guide](https://pytorch.org/get-started/locally/)
