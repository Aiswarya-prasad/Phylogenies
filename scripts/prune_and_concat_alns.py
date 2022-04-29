#!/usr/bin/env python3

"""
Edited version of GSL's SAGE_FetchAndConcat2.py to work with snakemake
# parse table and extract genes ID per genomes
# dict to store genomes names as key and  dict key= all values sequences as value
#%%
"""

# modules

import os
from os import path
import sys
import shutil
import numpy as np
import pandas as pd

#%%
#functions

def fasta_to_dict(fasta, Path, pipeNames) :
    dict_OG = {}
    seq = ''
    pathToFaa = Path + '/' + fasta
    with open(pathToFaa) as f:
        OG_ID = LongName(f.readline(), pipeNames)
        for line in f :
            if line.startswith(">") :
                dict_OG[OG_ID] = list(seq)
                OG_ID = LongName(line, pipeNames)
                seq = ''
            else :
                seq += line.rstrip()
            dict_OG[OG_ID] = list(seq)

    return(dict_OG)

def LongName(string, pipeNames='False') :
    string = "GrY2_F2_MAG_00000011_01168"
    if len(string.rstrip().split('_'))>2 and pipeNames=='False' :
        OG_ID = "_".join(string.rstrip().split('_')[0:-1])
        # OG_ID = "_".join(string.rstrip().split('_')[0:2])
        return(OG_ID)
    elif pipeNames=='True' :
        OG_ID=string.rstrip().split('|')[0].split('_')[0]
        print(OG_ID)
        return(OG_ID)
    else :
        OG_ID = string.rstrip().split('_')[0]
        return(OG_ID)



def insert_newlines(string, every=64) :
    return ('\n'.join(string[i:i + every] for i in range(0, len(string), every)))

PathAlignedOG=snakemake.input["aligned_dir"]
OutDir=snakemake.output["pruned_dir"]
WriteAll=True
pipeNames=snakemake.params["pipe_names"]

if not os.path.exists(OutDir):
    os.makedirs(OutDir)

OGList = [item for item in os.listdir(PathAlignedOG) if not item.startswith(".")]
AlignedFiles=[OG for OG in OGList if os.path.isfile(PathAlignedOG+'/'+OG)]
MainDictionary={}
for OG in AlignedFiles :
    OGDictionary = fasta_to_dict(OG, PathAlignedOG, pipeNames)
    MainDictionary[OG]=OGDictionary

#%%
# 1 : MainDictionary key(OG number) : value(OG_ID) (dictionary)
# 2 : dict_OG = key(sample name) : value(sequence)

allPrunedSequences = {}
FIRST='True'
gap = {'-'}
ProteinTable = pd.DataFrame(AlignedFiles)
ProteinTable.insert(0, 'Genome',AlignedFiles , True)
i=0

for alnOG in sorted(MainDictionary.keys()) :
    i+=1
    seq_table = pd.DataFrame.from_dict(MainDictionary[alnOG], orient='columns')
    gap_count = seq_table.isin(gap).sum(1)/seq_table.shape[1]
    #print(gap_count)
    prunedSeq_table = seq_table[gap_count<0.5]
    pruned_seq = prunedSeq_table.to_dict('list')
    pruned_seq = {key : ''.join(pruned_seq[key]) for key in pruned_seq.keys()}

    if WriteAll==True:
        pruned_aln = OutDir + '/' + alnOG.split(".")[0] + '_pruned.fasta'
        with open(pruned_aln, 'w') as pruned_alignment :
            for key, value in pruned_seq.items() :
                pruned_alignment.write(key + '\n')
                pruned_alignment.write(insert_newlines(value) + '\n')


    for key, value in pruned_seq.items() :
        print(key)
        if FIRST=='True' :
            allPrunedSequences[key] = (insert_newlines(value) + '\n')

        else :
            print('else')
            allPrunedSequences[key] += insert_newlines(value) + '\n'
    FIRST='False'
#%%
full_aln =  OutDir + '/' + 'CoreGeneAlignment.fasta'
with open(full_aln, 'w') as pruned_alignment :
    for key, value in allPrunedSequences.items() :
        GenomeName = key
        pruned_alignment.write(GenomeName + '\n')
        pruned_alignment.write(value)
   # list_allseq.append(list_seq)
