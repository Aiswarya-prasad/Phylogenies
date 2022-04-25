#!/bin/bash

####--------------------------------------
##SLURM options
####--------------------------------------
#SBATCH --job-name Acineto.Prokka
#SBATCH --partition wally
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 6G
#SBATCH --time 00:10:00
#SBATCH --output /users/gsartonl/scratch/gsartonl/Phylogenies/logs/Acineto/%x_%A_%a.out
#SBATCH --error /users/gsartonl/scratch/gsartonl/Phylogenies/logs/Acineto/%x_%A_%a.err
#SBATCH --array 0-18
####--------------------------------------${sample}
##preparation
##set you bash variables in order to quickly call them in the script
####--------------------------------------

username=gsartonl
SAMPLE=(2576861830.fna 2579779043.fna 2582580904.fna 2660238279.fna 2671180229.fna 2690315601.fna 2716885010.fna 2718218092.fna 2767802241.fna 2775507289.fna 2791355508.fna 2847069759.fna 2851902142.fna 2860756204.fna 2860764190.fna 2864751016.fna 2870361953.fna Asp_C16S1.fna Acinetobacter_sp_ESL0695.fna)
LocusTag=(2pS NIPH_512 CIP_110549 KCTC_2357 ANC_4422 ANC_5114 BBS3 AntiMn-1 BRTC-1 408225 PBJ7 MC5 2010C01-170 ATCC_17978 LXL_C1 S00005 QZS01 C16S1 ESL0695
)
threads=8

Genus=Acinetobacteria

working_directory=/users/${username}/scratch/${username}/Phylogenies/${Genus}
FilesLocation=${working_directory}/GenomesFNA
Overall_output_directory=${working_directory}/02_PROKKA

mkdir -p ${working_directory}

####---------------------
##modules
####--------------------------------------

module load HPC/Software 
module load UHTS/Analysis/prokka/1.13

####--------------------------------------
##start of script
####--------------------------------------

start=$SECONDS
echo "Step: PROKKA annotation"

#done
echo "=========================================================="
date +"START : %a %b %e %Y %H:%M:%S "
echo -e "Sample Name: "${SAMPLE[$SLURM_ARRAY_TASK_ID]}


###===========================
##PROKKA
echo -e "-------1. PROKKA annotation"
###===========================   


prokka --compliant  --outdir ${Overall_output_directory}/${SAMPLE[$SLURM_ARRAY_TASK_ID]} --locustag ${LocusTag[$SLURM_ARRAY_TASK_ID]} --prefix ${LocusTag[$SLURM_ARRAY_TASK_ID]} --proteins proteins.faa --evalue 0.001 ${FilesLocation}/${SAMPLE[$SLURM_ARRAY_TASK_ID]}

####--------------------------------------
##End of script
####--------------------------------------
date +"END : %a %b %e %H:%M:%S %Z %Y"
echo "=========================================================="

duration=$(( SECONDS - start ))
echo -e "The script ran for "${duration} "seconds"