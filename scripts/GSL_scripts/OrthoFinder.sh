#!/bin/bash

####--------------------------------------
##SLURM options
####--------------------------------------
#SBATCH --job-name OrthoFinder
#SBATCH --partition wally
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 16
#SBATCH --mem 12G
#SBATCH --time 04:30:00
#SBATCH --output /users/gsartonl/scratch/gsartonl/Phylogenies/logs/%x_%A_%a.out
#SBATCH --error /users/gsartonl/scratch/gsartonl/Phylogenies/logs/%x_%A_%a.err
#SBATCH --array 0-5
####--------------------------------------${sample}
##preparation
##set you bash variables in order to quickly call them in the script
####

module load Phylogeny/OrthoFinder/2.3.8 

username=gsartonl
Genus=(Acinetobacteria Bifidobacteria Enterobacteriaceae Firm5 Floricoccus Leuconostoc Snodgrassella)

working_directory=/users/${username}/scratch/${username}/Phylogenies/${Genus[$SLURM_ARRAY_TASK_ID]}
AnnotationsDir=${working_directory}/02_PROKKA
FilesLocation=${working_directory}/03_GenesFaa

mkdir -p ${FilesLocation}

for FILE in $(find ${AnnotationsDir} | grep -v 'Data' | grep -v 'intergenic' | grep -v 'genes' | grep -v 'Ga' | grep 'faa') ; do cp ${FILE} ${FilesLocation} ; done


orthofinder -f ${FilesLocation} -t 16 -n ${Genus[$SLURM_ARRAY_TASK_ID]}