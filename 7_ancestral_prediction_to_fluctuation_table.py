#!/usr/bin/env python

# coding: utf-8 

# Funding received from the European Research Council and the Sigrid Jusélius Foundation contributed to the development of this software.
# Author: Cory Dunn
# Institution: University of Helsinki
# Author Email: cory.david.dunn@gmail.com
# License: GPLv3

import pandas as pd
import os
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import AlignIO

# Set paths

inputpath = '/Users/corydunn/Library/CloudStorage/Dropbox/Python_software/orthomam_AUG_13_2023/OrthoMAM10c_CDS/ANCESTRAL/116_prepared_for_fluct/'
outputpath1 = '/Users/corydunn/Library/CloudStorage/Dropbox/Python_software/orthomam_AUG_13_2023/OrthoMAM10c_CDS/ANCESTRAL/fluctuation_results/'

codon_table = 1 ## universal (standard) codon table
verbose_flag = 'y'

selected_accession = 'Homo_sapiens'

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

print('Alignments for potential study: ', len(filestem_set))
print('Alignments removed due to aberrations (eg. gap only columns): ', len(bad_alignments))
print('Alignments for final analysis of substitutions: ', len(selected_alignments_list))

# To load the FASTA alignment file

def read_fasta(alignment):  # To read the FASTA alignment file
    aa_dict = {}
    with open(alignment, mode='r') as handle:
        # Using Biopython's parse function to reduce memory footprint
        for record in SeqIO.parse(handle, 'fasta'):
            # Extract individual parts of the FASTA record
            identifier = record.id
            sequence = record.seq
            aa_dict[identifier] = sequence
    return aa_dict

# To retrieve clades by name

def lookup_by_names(tree):
    names = {}
    for clade in tree.find_clades():
        if clade.name:
            if clade.name in names:
                raise ValueError("Duplicate key: %s" % clade.name)
            names[clade.name] = clade
    return names

# To retrieve all edges in the tree

def all_edges(tree):

    alledges = []
    for parent in tree.find_clades(terminal=False, order='level'):
        for child in parent.clades:
            alledges.append(parent.name + '*' + child.name)

    return alledges

analyzed_four_fold_positions = []

report_everything_about_selected_positions_results = pd.DataFrame(columns = ['Protein','Degenerate_index','Edge_name','Edge_remark','Ancestral_node','Descendant_node_or_species','Ancestral_character','Descendant_character','Character_remark','Branch_length','Analyzed_strand_ancestral','Analyzed_strand_descendent'])

for gene in selected_alignments_list:
    print('Analyzing: ', gene)    
    # Initialize lists for FASTA information

    record_x_toward_seq_dataframe = []
    sequence_records = []
    alignment_record_name_list = []
    
    # Load input FASTA into dataframe
    
    alignfile = inputpath + gene + '_full.fasta'
    for record in SeqIO.parse(alignfile,"fasta"):
        alignment_record_name_list.append(record.name)
        record_x_toward_seq_dataframe = list(record.seq)
        record_x_toward_seq_dataframe_UPPER = [x.upper() for x in record_x_toward_seq_dataframe] 
        sequence_records.append(record_x_toward_seq_dataframe_UPPER)

    sequence_dataframe = pd.DataFrame(sequence_records,index=alignment_record_name_list)
    sequence_dataframe_joined = sequence_dataframe.apply(''.join, axis=1)
    reference_sequence_concatenate = sequence_dataframe_joined.at[selected_accession]
    reference_sequence_concatenate_Seq = Seq(reference_sequence_concatenate)
    try:
        translated_reference_sequence_concatenate_Seq = reference_sequence_concatenate_Seq.translate(table=codon_table)
    except:
        print('Fail to translate reference sequence for file:', gene)
        continue
    translated_reference_sequence_concatenate_STR = str(translated_reference_sequence_concatenate_Seq)

   # Determine degenerate sites among mammals (ignore Anolis punctatus when testing mammals)

    sequence_dataframe_WO_A_punctatus = sequence_dataframe.copy(deep=True)
    #sequence_dataframe_WO_A_punctatus.drop("NC_044125_1_Anolis_punctatus",axis="index")

    indices_by_amino_acid = []
    indices_by_nucleotide = []
    degenerate_count = 0
    
    for i in range(len(translated_reference_sequence_concatenate_STR)):
        amino_acid = translated_reference_sequence_concatenate_STR[i]

        # Filter out non-standard characters

        identity_test_codon_1 = sequence_dataframe_WO_A_punctatus[i*3]
        identity_test_codon_1_FILTER = identity_test_codon_1[(identity_test_codon_1 == 'A') | \
            (identity_test_codon_1 == 'C') | \
            (identity_test_codon_1 == 'G') | \
            (identity_test_codon_1 == 'T')]
        identity_test_codon_1_VC = identity_test_codon_1_FILTER.value_counts(normalize = True)

        identity_test_codon_2 = sequence_dataframe_WO_A_punctatus[i*3+1]
        identity_test_codon_2_FILTER = identity_test_codon_2[(identity_test_codon_2 == 'A') | \
            (identity_test_codon_2 == 'C') | \
            (identity_test_codon_2 == 'G') | \
            (identity_test_codon_2 == 'T')]
        identity_test_codon_2_VC = identity_test_codon_2_FILTER.value_counts(normalize = True)

        # Identify degenerate sites

        if ((identity_test_codon_1_VC == 1).any()) & ((identity_test_codon_2_VC == 1).any()) & (amino_acid in 'AGPTV'):
            
            #indices_by_amino_acid.append(i)
            j = i*3+2
            #indices_by_nucleotide.append(j)
            nucleotide_possibilities_at_degenerate = sequence_dataframe[j].value_counts(normalize=True)
            nucleotide_possibilities_at_degenerate_list = nucleotide_possibilities_at_degenerate.index.tolist()
            nucleotide_possibilities_at_degenerate_string = ''.join(nucleotide_possibilities_at_degenerate_list)
            nucleotide_possibilities_at_degenerate_string = sorted(re.sub('[^ACGT]', '', nucleotide_possibilities_at_degenerate_string))
            
            try:
                gap_pc = nucleotide_possibilities_at_degenerate['-']
            except:
                gap_pc = 0

            # Ensure site is truly degenerate and gap% < 5%
            
            if (len(nucleotide_possibilities_at_degenerate_string) ==  4) & (gap_pc < 0.05):
                degenerate_count += 1
                degenerate_index = gene + '_' + str(degenerate_count)
                label_of_site = gene, amino_acid, i+1, (i+1)*3,degenerate_index,nucleotide_possibilities_at_degenerate_string
                analyzed_four_fold_positions.append(label_of_site)
                indices_by_amino_acid.append(i)
                indices_by_nucleotide.append(j)
                print (label_of_site)
        
        elif ((identity_test_codon_1_VC == 1).any()) & ((identity_test_codon_2_VC == 1).any()) & (amino_acid in 'S'):
            if (identity_test_codon_1_VC[identity_test_codon_1_VC == 1].index[0] == 'T') & (identity_test_codon_2_VC[identity_test_codon_2_VC == 1].index[0] == 'C'):
                
                #indices_by_amino_acid.append(i)
                j = i*3+2
                #indices_by_nucleotide.append(j)
                nucleotide_possibilities_at_degenerate = sequence_dataframe[j].value_counts(normalize=True)
                nucleotide_possibilities_at_degenerate_list = nucleotide_possibilities_at_degenerate.index.tolist()
                nucleotide_possibilities_at_degenerate_string = ''.join(nucleotide_possibilities_at_degenerate_list)
                nucleotide_possibilities_at_degenerate_string = sorted(re.sub('[^ACGT]', '', nucleotide_possibilities_at_degenerate_string))
                
                try:
                    gap_pc = nucleotide_possibilities_at_degenerate['-']
                except:
                    gap_pc = 0

                # Ensure site is truly degenerate and gap% < 5%
            
                if (len(nucleotide_possibilities_at_degenerate_string) ==  4) & (gap_pc < 0.05):
                    degenerate_count += 1
                    degenerate_index = gene + '_' + str(degenerate_count)
                    label_of_site = gene, amino_acid, i+1, (i+1)*3,degenerate_index,nucleotide_possibilities_at_degenerate_string
                    indices_by_amino_acid.append(i)
                    indices_by_nucleotide.append(j)
                    analyzed_four_fold_positions.append(label_of_site)
                    print (label_of_site)
        
        elif ((identity_test_codon_1_VC == 1).any()) & ((identity_test_codon_2_VC == 1).any()) & (amino_acid in 'R'):
            if (identity_test_codon_1_VC[identity_test_codon_1_VC == 1].index[0] == 'C') & (identity_test_codon_2_VC[identity_test_codon_2_VC == 1].index[0] == 'G'):
                
                j = i*3+2

                nucleotide_possibilities_at_degenerate = sequence_dataframe[j].value_counts(normalize=True)
                nucleotide_possibilities_at_degenerate_list = nucleotide_possibilities_at_degenerate.index.tolist()
                nucleotide_possibilities_at_degenerate_string = ''.join(nucleotide_possibilities_at_degenerate_list)
                nucleotide_possibilities_at_degenerate_string = sorted(re.sub('[^ACGT]', '', nucleotide_possibilities_at_degenerate_string))
    
                try:
                    gap_pc = nucleotide_possibilities_at_degenerate['-']
                except:
                    gap_pc = 0
                
                # Ensure site is truly degenerate and gap% < 5%
            
                if (len(nucleotide_possibilities_at_degenerate_string) ==  4) & (gap_pc < 0.05):
                    degenerate_count += 1
                    degenerate_index = gene + '_' + str(degenerate_count)
                    label_of_site = gene, amino_acid, i+1, (i+1)*3,degenerate_index,nucleotide_possibilities_at_degenerate_string
                    indices_by_amino_acid.append(i)
                    indices_by_nucleotide.append(j)
                    analyzed_four_fold_positions.append(label_of_site)
                    print (label_of_site)

    four_fold_dataframe = sequence_dataframe[indices_by_nucleotide]
    four_fold_dataframe_joined = four_fold_dataframe.apply(''.join, axis=1)

    # Write degenerate sites to FASTA

    ofile = open(outputpath1 + 'degenerate_nucleotide_sites_' + gene + '.fasta', "w")
    for seqi in range(len(four_fold_dataframe_joined)):
        ofile.write(">" + alignment_record_name_list[seqi] + "\n" + four_fold_dataframe_joined[seqi] + "\n")
    ofile.close()
    
    alignment_file = outputpath1 + 'degenerate_nucleotide_sites_' + gene + '.fasta'
    
    # Read the ancestral tree and record all edges
    ancestral_tree = 'ANC_' + gene + '.raxml.ancestralTree'
    from Bio import Phylo
    my_tree = Phylo.read(inputpath+ancestral_tree, 'newick')
    edges = all_edges(my_tree)
    clades_dict = lookup_by_names(my_tree)

    # Read the four-fold degenerate FASTA file and record character values for each position

    all_values = read_fasta(alignment_file)

    report_everything_about_selected_positions_results_GENE = pd.DataFrame(columns = ['Protein','Alignment_position','Edge_name','Edge_remark','Ancestral_node','Descendant_node_or_species', 'Ancestral_character','Descendant_character','Character_remark','Branch_length','Analyzed_strand_ancestral','Analyzed_strand_descendent','Degenerate_index'])
    #length_of_report_array_count = 0

    for i in range(len(indices_by_nucleotide)):
        query_label = str(i+1)
        query = i
        for edge in edges:
            parent_edge = edge.split('*')[0]
            child_edge = edge.split('*')[1]
            edge_len = str(clades_dict[child_edge].branch_length)
            parent_value = all_values[parent_edge][query]
            child_value = all_values[child_edge][query]
            if edge.count('_') == 0:
                edge_remark = 'Internal'
            elif edge.count('_') > 0:
                edge_remark = 'All'
            if parent_value == child_value:
                value_remark = 'Conserved'
            else:
                value_remark = 'Fluctuating'
            
            Lanc = parent_value
            Ldes = child_value
            
            other_degenerate_index = gene + '_' + query_label
            to_append_report = ([gene,query_label, edge, edge_remark, parent_edge, child_edge, parent_value, child_value,
                             value_remark, edge_len, Lanc, Ldes, other_degenerate_index])
            report_everything_about_selected_positions_results_GENE.loc[len(report_everything_about_selected_positions_results_GENE)] = to_append_report
            #length_of_report_array_count += 1
    
    report_everything_about_selected_positions_results_GENE.to_csv(outputpath1+'fluctuation_analysis_degenerate_nucleotide_sites_' + gene + '.csv')
    report_everything_about_selected_positions_results
    report_everything_about_selected_positions_results = pd.concat([report_everything_about_selected_positions_results,report_everything_about_selected_positions_results_GENE], sort=False)

analyzed_four_fold_positions_DF = pd.DataFrame(analyzed_four_fold_positions, columns = ['Protein','Amino_acid','Amino_acid_position','CDS_nucleotide_position','Degenerate_index','Nucleotide_possibilities_at_degenerate (ACGT)'])

# Generate report with each site along each edge (substitution or conservation)

report_everything_after_merge = pd.merge(report_everything_about_selected_positions_results, analyzed_four_fold_positions_DF, on = ['Protein','Degenerate_index'], how = "left")
report_everything_after_merge.to_csv(outputpath1+'fluctuation_analysis_after_merge_conserved_and_substitution_all_chars.csv')

# Generate reports just referring to substitutions at each position

substitutions = report_everything_after_merge[report_everything_after_merge['Analyzed_strand_ancestral'] != report_everything_after_merge['Analyzed_strand_descendent']]
substitutions = substitutions[(substitutions['Analyzed_strand_ancestral'].isin(['A','C','G','T'])) & (substitutions['Analyzed_strand_descendent'].isin(['A','C','G','T']))] 
substitutions['Alignment_position'] = substitutions['Alignment_position'].astype(str)
substitutions['Degenerate_index'] = substitutions['Protein'] + '_' + substitutions['Alignment_position']

TSS_DF = substitutions['Degenerate_index'].value_counts()
TSS_DF = TSS_DF.reset_index()
TSS_DF.rename(columns = {'Degenerate_index':'TSS'}, inplace = True)
TSS_DF.rename(columns = {'index':'Degenerate_index'}, inplace = True)

FINAL_A = substitutions.set_index('Degenerate_index').join(TSS_DF.set_index('Degenerate_index'))
FINAL_A = FINAL_A.sort_values(by='TSS', ascending=True)

FINAL_A.to_csv(outputpath1+'fluctuation_analysis_with_TSS_no_odd_characters_or_gaps_only_substitutions.csv')

FINAL_B = analyzed_four_fold_positions_DF.set_index('Degenerate_index').join(TSS_DF.set_index('Degenerate_index'))
FINAL_B.to_csv(outputpath1+'analyzed_degenerate_positions_all_proteins_w_TSS.csv')