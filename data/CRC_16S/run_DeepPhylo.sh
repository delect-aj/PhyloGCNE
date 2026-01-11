#!/bin/bash
#SBATCH --job-name=DeepPhylo
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -p gpu
#SBATCH --mem=20G
#SBATCH -o deepphylo.out
#SBATCH -e deepphylo.err

module load miniconda/4.9.2 
source activate jupyter_notebook

software_path=/home/dongbiao/software/DeepPhylo

for i in $(seq 1 5); do
python ${software_path}/deepphylo_classification.py \
    --train_table data/train_${i}.biom \
    --test_table data/test_${i}.biom \
    --metadata_filename metadata.tsv \
    --pca_file_path PCA_32.txt --output_dir results --fold ${i} \
    --epochs 100 -hs 80 -kec 13 -l 0.0001 -bs 64 -kep 4 -act relu --hidden_size 32
done
