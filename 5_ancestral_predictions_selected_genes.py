#!/usr/bin/env python

# coding: utf-8 

# Funding received from the European Research Council and the Sigrid Jusélius Foundation contributed to the development of this software.
# Author: Cory Dunn
# Institution: University of Helsinki
# Author Email: cory.david.dunn@gmail.com
# License: GPLv3

# Load libraries

import pandas as pd
import numpy as np
import os
import random

# Load file with analysis of CDS alignments

gap_and_alignment_length_DICT = pd.read_csv('/Users/corydunn/Library/CloudStorage/Dropbox/Python_software/orthomam_AUG_13_2023/OrthoMam_CDS_gap_pc_alignment_length_sequence_length.csv')

# Select alignments with at least 111 OrthoMAM samples represented by a sequence

selected_alignments = gap_and_alignment_length_DICT[gap_and_alignment_length_DICT['ALIGN_DEPTH_W_DUPS'] == 112]

# Load existing list of accessions

orthomam_accessions = []
infile = open('/Users/corydunn/Library/CloudStorage/Dropbox/Python_software/orthomam_AUG_13_2023/OrthoMAM_accessions.txt','r')
for line in infile:
    orthomam_accessions.append(line.strip('\n').strip('>'))
infile.close()

selected_alignment_list = selected_alignments['FILE'].to_list()

print('Number of selected alignments: ', len(selected_alignment_list))

# Generate gene by taxon sequence dataframe

inputpath = '/Users/corydunn/Library/CloudStorage/Dropbox/Python_software/orthomam_AUG_13_2023/OrthoMAM10c_CDS/TABULAR/'
outputpath1 = '/Users/corydunn/Library/CloudStorage/Dropbox/Python_software/orthomam_AUG_13_2023/OrthoMAM10c_CDS/ANCESTRAL/116_input/'

# Load input FASTA file for each selected alignment

for index, row in selected_alignments.iterrows():
    
    accessions_versus_selected_sequences = pd.DataFrame(index=orthomam_accessions)
    file = row['FILE']
    
    file_tabular_load = pd.read_csv(inputpath+file,sep="\t",header=None)#,index_col='0')
    file_tabular_load.columns=['Accession',file,'X']
    del file_tabular_load['X']
    file_tabular_load = file_tabular_load.set_index('Accession')
    accessions_versus_selected_sequences = pd.concat([accessions_versus_selected_sequences, file_tabular_load], axis=1)
    
    # For those accessions without a sequence, enter gaps
    
    sequence_length = int(accessions_versus_selected_sequences[file].str.len().max())
    gap_filler = "-" * sequence_length
    accessions_versus_selected_sequences[file] = accessions_versus_selected_sequences[file].replace(np.nan, gap_filler)
    
    # Enter random characters into last six columns, so that duplicate calls are unlikely (even if the FASTA are actually identical)

    accessions_versus_selected_sequences['rand1'] = 'X'
    accessions_versus_selected_sequences['rand2'] = 'X'
    accessions_versus_selected_sequences['rand3'] = 'X'
    accessions_versus_selected_sequences['rand4'] = 'X'
    accessions_versus_selected_sequences['rand5'] = 'X'
    accessions_versus_selected_sequences['rand6'] = 'X'

    for i in orthomam_accessions:
        accessions_versus_selected_sequences.at[i,'rand1'] = random.choice('ACGT')
        accessions_versus_selected_sequences.at[i,'rand2'] = random.choice('ACGT')
        accessions_versus_selected_sequences.at[i,'rand3'] = random.choice('ACGT')
        accessions_versus_selected_sequences.at[i,'rand4'] = random.choice('ACGT')
        accessions_versus_selected_sequences.at[i,'rand5'] = random.choice('ACGT')
        accessions_versus_selected_sequences.at[i,'rand6'] = random.choice('ACGT')
    
    print(accessions_versus_selected_sequences)
    columns_to_combine = [file,'rand1','rand2','rand3','rand4','rand5','rand6']

    accessions_versus_selected_sequences['sequence_to_send'] = accessions_versus_selected_sequences[columns_to_combine].agg(''.join, axis=1)

    sequence_toward_FASTA = accessions_versus_selected_sequences['sequence_to_send']

    filestem = file.split(".")[0]

    # Save resulting FASTA

    ofile = open(outputpath1+filestem+'_116.fasta', "w")
    for index, value in sequence_toward_FASTA.items():
        ofile.write(">" + str(index) + "\n" + value + "\n")
    ofile.close()

    # Perform ancestral predictions on FASTA files

    os.system('raxml-ng --ancestral --msa ' + outputpath1 + filestem + '_116.fasta --tree /Users/corydunn/Library/CloudStorage/Dropbox/Python_software/orthomam_AUG_13_2023/T3_ORTHOMAM_bestTree_rooted_Ornithorhynchus_anatinus.nwk --model GTR+I+G4 --prefix ANC_' + filestem + ' --thread 7')

    # gzip ancestral probabilities to save disk space

    os.system('gzip -v9 ' + 'ANC_' + filestem + '.raxml.ancestralProbs')

    #os.system('rm ' + 'ANC_' + filestem + '.raxml.ancestralProbs')