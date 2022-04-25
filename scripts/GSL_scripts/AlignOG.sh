#!/bin/bash

WorkingDir=/raw_data/garance1/Phylogenies/
Genus=$1
pathToOGs=${WorkingDir}/${Genus}/03_GenesFaa/OrthoFinder/Results_${Genus}/Single_Copy_Orthologue_Sequences
outDir=${WorkingDir}/${Genus}/04_Aligned_OG

mkdir -p ${outDir}

N=8

for OG in $(ls ${pathToOGs}); do
    (
        mafft --amino --inputorder --localpair --maxiterate 1000 ${pathToOGs}/${OG} > ${outDir}/${OG/.fa/_aligned.fa};
        echo "starting task $OG.."
        sleep $(( (RANDOM % 3) + 1))
    ) &

    # allow to execute up to $N jobs in parallel
    if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
        # now there are $N jobs already running, so wait here for any job
        # to be finished so there is a place to start next one.
        wait -n
    fi

done

# no more jobs to be started but wait for pending jobs
# (all need to be finished)
wait

echo "all done"