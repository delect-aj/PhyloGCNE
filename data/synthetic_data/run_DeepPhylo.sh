#!/bin/bash
#SBATCH --job-name=deep_phylo
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -p gpu
#SBATCH --mem=20G
#SBATCH -o deepphylo.out
#SBATCH -e deepphylo.err

module load miniconda/4.9.2
source activate jupyter_notebook

# Path to DeepPhylo installation (https://github.com/liuchuwei/DeepPhylo)
SOFTWARE_PATH=/path/to/DeepPhylo

# Data directory (relative to this script)
DATA_DIR=$(dirname "$0")

for i in $(seq 1 5); do
    python ${SOFTWARE_PATH}/deepphylo_classification.py \
        --train_table ${DATA_DIR}/data/train_${i}.biom \
        --test_table ${DATA_DIR}/data/test_${i}.biom \
        --metadata_filename ${DATA_DIR}/metadata.txt \
        --pca_file_path ${DATA_DIR}/PCA_32.txt \
        --output_dir ${DATA_DIR}/results \
        --fold ${i} \
        --epochs 100 -hs 80 -kec 3 -l 0.0001 -bs 64 -kep 7 -act relu --hidden_size 30
done
