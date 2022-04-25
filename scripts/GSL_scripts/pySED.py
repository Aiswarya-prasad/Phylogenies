#!/usr/bin/env python3
####################################################
				# @gsartonl
				# 08.06.2021
####################################################

############# python sed function #############
# input :
## argv[1] : If 'True' : change strings in file using correspondance file. If 'False' change a single
###			 string/pattern in file. Default is False
## argv[2] : file where names should be changed (include path ! )
## argv[3] : output name (only file name)
## argv[4] : table with names correspondance OR new string
## argv[5] : old pattern if False


import sys
import os
from os import path
import pandas as pd
import numpy as np
FP = sys.argv[2]
ChangeIn = sys.argv[1]
print(ChangeIn)

FilePath = os.path.dirname(FP)
OutName = FilePath +'/' +  sys.argv[3]

if ChangeIn == 'False' :
    oldStrings = sys.argv[5]
    if oldStrings=='empty' :
        oldStrings==''
    newStrings = sys.argv[4]
    with open (sys.argv[2], 'r') as fileToModify :
        FileContent = fileToModify.read()
        FileToWrite=FileContent.replace(oldStrings, newStrings)
        with open(FilePath +'/' +  sys.argv[3], 'w') as FW :
            FW.write(FileToWrite)


elif ChangeIn == 'True' :

    i=0
    data = pd.read_excel(sys.argv[4], dtype=str)
    #data = pd.read_excel('/Users/garancesarton-loheac/Documents/PhD/Phylogenies/ForPublication/Gilliamella/GilliamellaGenomesDB.xlsx' ,dtype=str)
    data['CompleteName'] = data['Species'] + ' ' + data['Strain_name']
    locusTags = data['Locus_tag']

   # with open (sys.argv[2], 'r') as FM :
    with open (FP, 'r') as FM :
        fileToModify = FM.read()
        first='True'
        for tag in locusTags :
            
            if first=='True' :
                FileToWrite = fileToModify.replace(tag,  data[data['Locus_tag']==tag] ['CompleteName'].to_string().split('    ')[1])
                first='False'
            else :

                FileToWrite = FileToWrite.replace(tag, data[data['Locus_tag']==tag] ['CompleteName'].to_string().split('    ')[1])





else :
 	print('Incorrect argument. Argument is False for single pattern replacementis True for multiple pattern replacement')

with open(OutName, 'w') as FW :
    FW.write(FileToWrite)
