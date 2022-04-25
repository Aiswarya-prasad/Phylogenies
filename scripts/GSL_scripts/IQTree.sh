#!/bin/bash



eval "$(conda shell.bash hook)"
conda activate PhyloTools

WorkingDir=/raw_data/garance1/Phylogenies/
pathToScript=${WorkingDir}/scripts
genus=("Gilliamella" "Firm5"


iqtree -s ${WorkingDir}/${genus}/05_Whole_Alignments/CoreGeneAlignment.fasta\
        -st AA -nt 16 -bb 1000 -seed 12345 -m TEST\
        -pre ${genus}_Phylogeny
