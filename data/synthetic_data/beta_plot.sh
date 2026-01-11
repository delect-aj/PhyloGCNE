#!/bin/bash
#SBATCH --job-name=beta_plot
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -p cu
#SBATCH --mem=20G
#SBATCH -o beta.out
#SBATCH -e beta.err
#SBATCH --exclude=cu01

module load miniconda/4.9.2
source activate qiime2-amplicon-2024.2

qiime tools import \
  --input-path table.biom \
  --type 'FeatureTable[RelativeFrequency]' \
  --input-format BIOMV210Format \
  --output-path beta/table.qza

qiime tools import \
  --input-path Synthetic_Dataset_1.nwk \
  --output-path beta/rooted-tree.qza \
  --type 'Phylogeny[Rooted]'
  
qiime diversity beta \
  --i-table beta/table.qza \
  --p-metric braycurtis \
  --o-distance-matrix beta/bray_curtis_distance_matrix.qza

qiime diversity beta-phylogenetic \
  --i-table beta/table.qza \
  --i-phylogeny beta/rooted-tree.qza \
  --p-metric weighted_unifrac \
  --o-distance-matrix beta/weighted_unifrac_distance_matrix.qza

