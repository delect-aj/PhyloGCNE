#!/bin/bash
#SBATCH --job-name=GCN_model
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -p gpu
#SBATCH --mem=20G
#SBATCH -o phylo_spec.out
#SBATCH -e phylo_spec.err

module load miniconda/4.9.2 
source activate PhyloSpec

software_path=/path/to/Phylo-Spec

for i in $(seq 1 5); do
python ${software_path}/src/model/PhyloSpec_train_test.py \
    -t phylogeny.nwk -c data/train_${i}.csv --PhyloSpec train -bs 64 \
    -taxo greengene2/exported-taxonomy/taxonomy.tsv -fold ${i} \
    -o output/ -lr 0.0001

python ${software_path}/src/model/PhyloSpec_train_test.py \
    -t phylogeny.nwk -c data/test_${i}.csv --PhyloSpec test -bs 64 \
    -taxo greengene2/exported-taxonomy/taxonomy.tsv -fold ${i} \
    -o output/ -lr 0.0001
    
done
