#!/usr/bin/env python

# coding: utf-8 

# Funding received from the European Research Council and the Sigrid Jusélius Foundation contributed to the development of this software.
# Author: Cory Dunn
# Institution: University of Helsinki
# Author Email: cory.david.dunn@gmail.com
# License: GPLv3

# Load libraries

import pandas as pd
import time

# Load file with analysis of CDS alignments

gap_and_alignment_length_DICT = pd.read_csv('OrthoMam_CDS_gap_pc_alignment_length_sequence_length.csv')

# Select alignments for treebuilding based upon relatively low gap % in FASTA, relatively compact gene structure, and presence in all of the 116 OrthoMAM samples

selected_alignments = gap_and_alignment_length_DICT[(gap_and_alignment_length_DICT['ALIGN_DEPTH_WO_DUPS'] == 116) & \
    (gap_and_alignment_length_DICT['GAPFRACTION'] <= 0.10) & \
    (gap_and_alignment_length_DICT['SEQUENCE_LENGTH'] <= 5000)]
                                                    
# Load existing list of accessions

orthomam_accessions = []
infile = open('OrthoMAM_accessions.txt','r')
for line in infile:
    orthomam_accessions.append(line.strip('\n').strip('>'))
infile.close()

selected_alignment_list = selected_alignments['FILE'].to_list()

# Prepare dataframe to hold selected sequences

accessions_versus_selected_sequences = pd.DataFrame(index=orthomam_accessions)

# Generate gene by taxon sequence dataframe

for index, row in selected_alignments.iterrows():
    
    file = row['FILE']
    path = '/Users/corydunn/Library/CloudStorage/Dropbox/Python_software/orthomam_AUG_13_2023/OrthoMAM10c_CDS/TABULAR/'
    file_tabular_load = pd.read_csv(path+file,sep="\t",header=None)#,index_col='0')
    file_tabular_load.columns=['Accession',file,'X']
    del file_tabular_load['X']
    file_tabular_load = file_tabular_load.set_index('Accession')
    accessions_versus_selected_sequences = pd.concat([accessions_versus_selected_sequences, file_tabular_load], axis=1)

# Concatenate selected sequences and save as FASTA

accessions_versus_selected_sequences_columns = accessions_versus_selected_sequences.columns.to_list()
accessions_versus_selected_sequences['Combined_CDS_selected'] = accessions_versus_selected_sequences[accessions_versus_selected_sequences_columns].agg(''.join, axis=1)

sequence_toward_FASTA = accessions_versus_selected_sequences['Combined_CDS_selected']

ofile = open('Combined_CDS_selected.fasta', "w")
for index, value in sequence_toward_FASTA.items():
    print(index,value)
    time.sleep(0.05)
    ofile.write(">" + str(index) + "\n" + value + "\n")
ofile.close()