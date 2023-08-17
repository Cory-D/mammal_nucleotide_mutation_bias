#!/usr/bin/env python

# coding: utf-8 

# Funding received from the European Research Council and the Sigrid Jusélius Foundation contributed to the development of this software.
# Author: Cory Dunn
# Institution: University of Helsinki
# Author Email: cory.david.dunn@gmail.com
# License: GPLv3

import os
import fileinput

# Set paths

inputpath = '/Users/corydunn/Library/CloudStorage/Dropbox/Python_software/orthomam_AUG_13_2023/OrthoMAM10c_CDS/ANCESTRAL/116_input/'
outputpath1 = '/Users/corydunn/Library/CloudStorage/Dropbox/Python_software/orthomam_AUG_13_2023/OrthoMAM10c_CDS/ANCESTRAL/116_prepared_for_fluct/'

# Find files in ancestral prediction directory

files = os.listdir(inputpath)

filestem_set = []

for file in files:
    if file[-22:] == '.raxml.ancestralStates':
        filestem = file.split("_")[1] + '_' + file.split("_")[2] + '_NT'
        filestem_set.append(filestem)

# Remove bad alignments

bad_alignments = []

for file in files:
    if file[-12:] == '.reduced.phy':
        bad = file.split("_")[1] + '_' + file.split("_")[2] + '_NT'
        bad_alignments.append(bad)
        
selected_alignments_list = list(set(filestem_set) - set(bad_alignments))

print(len(filestem_set))
print(len(bad_alignments))
print(len(selected_alignments_list))

# Add '>' to ancestral prediction nodes

for alignment in selected_alignments_list:
    for line in fileinput.input(outputpath1 + 'ANC_' + alignment + '.raxml.ancestralStates', inplace=True): #used outpath1, since I don't want to perturb the output from the last script
        print (line.replace("Node", ">Node")), # replace 'Node' with '>Node'

# Generate full FASTA files with ancestral and terminal nodes

    os.system("seqkit fx2tab " + outputpath1 + 'ANC_' + alignment + '.raxml.ancestralStates > ' + outputpath1 + 'ANC_' + alignment + '.tabular')
    os.system("seqkit fx2tab " + outputpath1 + alignment + '_116.fasta > ' + outputpath1 + alignment + '_116.tabular')
    os.system("cat " + outputpath1 + 'ANC_' + alignment + '.tabular ' + outputpath1 + alignment + '_116.tabular > ' + outputpath1 + alignment + '_full.tabular')
    os.system("seqkit tab2fx " +  outputpath1 + alignment + '_full.tabular > ' +  outputpath1 + alignment + '_full.fasta')