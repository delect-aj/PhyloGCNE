#!/bin/bash
#SBATCH --job-name=phylo_spec
#SBATCH -N 1
#SBATCH -n 28
#SBATCH -p cu
#SBATCH --mem=250G
#SBATCH -o phylo_spec.out
#SBATCH -e phylo_spec.err

module load miniconda/4.9.2 
source activate PhyloSpec

software_path=/home/dongbiao/software/Phylo-Spec
file_path=/home/dongbiao/GCN/data/OSCC_16S

for i in $(seq 1 5); do
python ${software_path}/src/model/PhyloSpec_train_test.py \
    -t ${file_path}/phylogeny.nwk -c ${file_path}/data/train_${i}.csv --PhyloSpec train -bs 64 \
    -taxo /beegfs/dongbiao/greengene2/exported-taxonomy/taxonomy.tsv -ep 50 -fold ${i} \
    -o output/

python ${software_path}/src/model/PhyloSpec_train_test.py \
    -t ${file_path}/phylogeny.nwk -c ${file_path}/data/test_${i}.csv --PhyloSpec test -bs 64 \
    -taxo /beegfs/dongbiao/greengene2/exported-taxonomy/taxonomy.tsv -fold ${i} \
    -o output/
    
done
