#!/usr/bin/env python

# coding: utf-8 

# Funding received from the European Research Council and the Sigrid Jusélius Foundation contributed to the development of this software.
# Author: Cory Dunn
# Institution: University of Helsinki
# Author Email: cory.david.dunn@gmail.com
# License: GPLv3

# Load libraries

import pandas as pd
import os

# Load existing list of accessions

orthomam_accessions = []
infile = open('OrthoMAM_accessions.txt','r')
for line in infile:
    orthomam_accessions.append(line.strip('\n').strip('>'))
infile.close()

# Determine the depth, length, gap fraction of each alignment file (tabular)

gap_and_alignment_length_DICT = {"FILE":[],"GAPFRACTION":[],"ALIGN_DEPTH_W_DUPS":[],"ALIGN_DEPTH_WO_DUPS":[],"SEQUENCE_LENGTH":[]}

path = '/Users/corydunn/Library/CloudStorage/Dropbox/Python_software/orthomam_AUG_13_2023/OrthoMAM10c_CDS/TABULAR/'
files = os.listdir(path)

for file in files:
    active_tabular = pd.read_csv(path+file,sep="\t",header=None)
    active_tabular.columns=['Accession','Sequence','X']
    del active_tabular['X']
    active_tabular = active_tabular.set_index('Accession')
    fraction_gapchar = active_tabular['Sequence'].str.count('-').sum() / active_tabular['Sequence'].str.count('').sum()
    
    gap_and_alignment_length_DICT["FILE"].append(file)
    gap_and_alignment_length_DICT["GAPFRACTION"].append(fraction_gapchar)
    
    alignment_depth_w_dups = len(active_tabular)
    gap_and_alignment_length_DICT["ALIGN_DEPTH_W_DUPS"].append(alignment_depth_w_dups)
    
    active_tabular_no_dups = active_tabular.drop_duplicates(subset='Sequence', keep='first')
    alignment_depth_WO_dups = len(active_tabular_no_dups)
    sequence_length_in_alignment = active_tabular['Sequence'].str.count('').sum() / len(active_tabular) - 1

    gap_and_alignment_length_DICT["ALIGN_DEPTH_WO_DUPS"].append(alignment_depth_WO_dups)
    
    gap_and_alignment_length_DICT["SEQUENCE_LENGTH"].append(sequence_length_in_alignment)
    
    print('File: ', file, 'Fraction gaps in alignment: ', fraction_gapchar, 'Depth of alignment with duplications: ',len(active_tabular), 'Depth of alignment without duplications: ',alignment_depth_WO_dups,'Sequence length in alignment: ', sequence_length_in_alignment)
    
gap_and_alignment_length_DICT = pd.DataFrame(gap_and_alignment_length_DICT,columns = ["FILE","GAPFRACTION","ALIGN_DEPTH_W_DUPS","ALIGN_DEPTH_WO_DUPS","SEQUENCE_LENGTH"])

gap_and_alignment_length_DICT.to_csv('OrthoMam_CDS_gap_pc_alignment_length_sequence_length.csv')