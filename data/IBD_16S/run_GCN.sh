#!/bin/bash
#SBATCH --job-name=GCN_model
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -p gpu
#SBATCH --mem=20G
#SBATCH -o GCN.out
#SBATCH -e GCN.err

module load miniconda/4.9.2 
source activate jupyter_notebook

model_file=../../GCN

for i in $(seq 1 5); do
CUDA_VISIBLE_DEVICES=7 python ${model_file}/GCN_model.py \
    --train_table data/train_${i}.biom \
    --test_table data/test_${i}.biom \
    --metadata_filename metadata.tsv \
    --phylogeny_file_path phylogeny_normalize.nwk \
    --fold ${i} --output_dir results --batch_size 128 \
    --learning_rate 0.001 --num_layers 4 --ffn True \
    --norm_type batch_norm --residual True \
    --hidden_channels 100 --dropout_rate 0.1 \
    --weight_decay 0.00001 # --feature_importance True --beta_alpha 2 --beta_beta 2 --stdev_spread 0.15
done
