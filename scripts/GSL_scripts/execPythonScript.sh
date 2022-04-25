#!/bin/bash



#eval "$(conda shell.bash hook)"
#conda activate garancepip

WorkingDir=/Users/garancesarton-loheac/Desktop/single_copy_aligned
pathToScript=/Users/garancesarton-loheac/Documents/PhD/Phylogenies/ForPublication/scripts
#genus=Bifidobacteria


#python3 ${pathToScript}/SAGE_FetchAndConcat2.py ${WorkingDir}/${genus}/04_Aligned_OG ${WorkingDir}/${genus}/05_Whole_Alignments

python3 ${pathToScript}/SAGE_FetchAndConcat2.py ${WorkingDir}/ ${WorkingDir}/05_Whole_Alignments False
