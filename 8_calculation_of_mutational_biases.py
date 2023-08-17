#!/usr/bin/env python

# coding: utf-8 

# Funding received from the European Research Council and the Sigrid Jusélius Foundation contributed to the development of this software.
# Author: Cory Dunn
# Institution: University of Helsinki
# Author Email: cory.david.dunn@gmail.com
# License: GPLv3

# Load dependencies

import pandas as pd
import numpy as np
import numpy.polynomial.polynomial as poly
import logging

# Input files

inputpath = '/Users/corydunn/Library/CloudStorage/Dropbox/Python_software/orthomam_AUG_13_2023/OrthoMAM10c_CDS/ANCESTRAL/fluctuation_results/'
outputpath = '/Users/corydunn/Library/CloudStorage/Dropbox/Python_software/orthomam_AUG_13_2023/OrthoMAM10c_CDS/ANCESTRAL/final_output/'

verbose_flag = 'y'

# Logging and streaming to console

mylogs = logging.getLogger(__name__)
mylogs.setLevel(logging.INFO)
stream = logging.StreamHandler()
stream.setLevel(logging.INFO)
streamformat = logging.Formatter("%(message)s")
stream.setFormatter(streamformat)
mylogs.addHandler(stream)

file = logging.FileHandler('orthomam_madprops_output.log')
mylogs.addHandler(file)

mylogs.info('___\n')

# Depending upon the nucleotide baseline determined above, calculate the relative fractions of total instantaneous nucleotide mutations. 

def total_inst_nuc_mut(nucleotide_baseline_input):

    sum_slopes_dN_F = slopes_compare_deplete_DF.loc['f_' + nucleotide_baseline_input + '_mut_X','f_A_mut_y':'f_T_mut_y'].sum()
    instant_fA_F = float(slopes_compare_deplete_DF.loc['f_' + nucleotide_baseline_input + '_mut_X','f_A_mut_y'] / sum_slopes_dN_F)
    instant_fC_F = float(slopes_compare_deplete_DF.loc['f_' + nucleotide_baseline_input + '_mut_X','f_C_mut_y'] / sum_slopes_dN_F)
    instant_fG_F = float(slopes_compare_deplete_DF.loc['f_' + nucleotide_baseline_input + '_mut_X','f_G_mut_y'] / sum_slopes_dN_F)
    instant_fT_F = float(slopes_compare_deplete_DF.loc['f_' + nucleotide_baseline_input + '_mut_X','f_T_mut_y'] / sum_slopes_dN_F)
    
    if verbose_flag == 'y':

        mylogs.info('\nInstantaneous fraction of total A mutated (' + nucleotide_baseline_input + ' baseline):' + str(round(instant_fA_F,2))+'\n')
        mylogs.info('\nInstantaneous fraction of total C mutated (' + nucleotide_baseline_input + ' baseline):' + str(round(instant_fC_F,2))+'\n')
        mylogs.info('\nInstantaneous fraction of total G mutated (' + nucleotide_baseline_input + ' baseline):' + str(round(instant_fG_F,2))+'\n')
        mylogs.info('\nInstantaneous fraction of total T mutated (' + nucleotide_baseline_input + ' baseline):' + str(round(instant_fT_F,2))+'\n')

    return(sum_slopes_dN_F,instant_fA_F,instant_fC_F,instant_fG_F,instant_fT_F)

def max_slope_out(i,j):
        
        edge_versus_mutation_counts_DF_column_X_label = edge_versus_mutation_counts_DF.columns[edge_versus_mutation_counts_DF_column_X]
        edge_versus_mutation_counts_DF_column_Y_label = edge_versus_mutation_counts_DF.columns[edge_versus_mutation_counts_DF_column_Y]
        
        fmut_X = edge_versus_mutation_counts_DF.iloc[:,edge_versus_mutation_counts_DF_column_X].to_numpy()
        fmut_y = edge_versus_mutation_counts_DF.iloc[:,edge_versus_mutation_counts_DF_column_Y].to_numpy()
        
        idx = np.isfinite(fmut_X) & np.isfinite(fmut_y)
        coeffs = np.polyfit(fmut_X[idx], fmut_y[idx], 3)

        coeff_x_power_3 = coeffs[0]
        coeff_x_power_2 = coeffs[1]
        coeff_x_power_1 = coeffs[2]
        coeff_0 = coeffs[3]

        deriv_power_2 = 3 * coeff_x_power_3
        deriv_power_1 = 2 * coeff_x_power_2
        deriv_0 = coeff_x_power_1

        test_slope_input = np.linspace (-0.025,0.1,100)
        calculate_slopes = lambda t: deriv_power_2*t*t + deriv_power_1*t + deriv_0
        slope_output = np.array([calculate_slopes(xi) for xi in test_slope_input])
        max_slope_output = np.amax(slope_output)
        
        slopes_compare_deplete_DF.iloc[i,j] = max_slope_output
        if verbose_flag == 'y':
        
            mylogs.info(edge_versus_mutation_counts_DF_column_X_label + ' versus ' + edge_versus_mutation_counts_DF_column_Y_label + \
                ' polynomial: ' + str(round(coeff_x_power_3,2)) + ' x^3 + '  + str(round(coeff_x_power_2,2)) + ' x^2 + ' + str(round(coeff_x_power_1,2)) + ' x + ' + str(round(coeff_0,2))+'\n')
            mylogs.info(edge_versus_mutation_counts_DF_column_X_label + ' versus ' + edge_versus_mutation_counts_DF_column_Y_label + ' derivative: ' + str(round(deriv_power_2,2))+' x^2 + ' + str(round(deriv_power_1,2)) + ' x + ' + str(round(deriv_0,2))+'\n')
            mylogs.info(edge_versus_mutation_counts_DF_column_X_label + ' versus ' + edge_versus_mutation_counts_DF_column_Y_label + ' maximum slope between X of -0.025 and 0.1 (used for calculations): ' + str(round(max_slope_output,2))+'\n')
            mylogs.info('\n')
        
        return

# Load fluctuation analyses and degenerate site TSS values)

report_everything_about_selected_positions_results_LOADDATAFRAME = pd.read_csv(inputpath + 'fluctuation_analysis_after_merge_conserved_and_substitution_all_chars.csv')
TSS_data = pd.read_csv(inputpath + 'analyzed_degenerate_positions_all_proteins_w_TSS.csv')

report_everything_about_selected_positions_results_MERGETSS = report_everything_about_selected_positions_results_LOADDATAFRAME.merge(TSS_data, how = 'outer',on=['Degenerate_index','Protein','Amino_acid',\
    'Amino_acid_position','CDS_nucleotide_position','CDS_nucleotide_position','Nucleotide_possibilities_at_degenerate (ACGT)'])

# Remove positions with TSS values ranking less than 50% to avoid positions under selection

report_everything_about_selected_positions_results_MERGETSS['Percentile_rank'] = report_everything_about_selected_positions_results_MERGETSS.TSS.rank(method = 'dense', pct = True)

report_everything_about_selected_positions_results = report_everything_about_selected_positions_results_MERGETSS[(report_everything_about_selected_positions_results_MERGETSS['Percentile_rank'] <= 1) & \
    (report_everything_about_selected_positions_results_MERGETSS['Percentile_rank'] >= 0.50)] # remove lowest <50%> of TSS values (may be under selection)

# Report upon TSS cut-offs in initial dataframe load versus positions selected for upper 50% of TSS values

minimum_TSS_of_input_positions = report_everything_about_selected_positions_results_MERGETSS['TSS'].min()
maximum_TSS_of_input_positions = report_everything_about_selected_positions_results_MERGETSS['TSS'].max()
minimum_TSS_of_analyzed_positions = report_everything_about_selected_positions_results['TSS'].min()
maximum_TSS_of_analyzed_positions = report_everything_about_selected_positions_results['TSS'].max()
if verbose_flag == 'y': 
    mylogs.info('Minimum TSS of input alignment positions: ' + str(minimum_TSS_of_input_positions))
    mylogs.info('\n')
    mylogs.info('Maximum TSS of input alignment  positions: ' + str(maximum_TSS_of_input_positions))
    mylogs.info('\n')
    mylogs.info('Minimum TSS of _analyzed_ positions: ' + str(minimum_TSS_of_analyzed_positions))
    mylogs.info('\n')
    mylogs.info('Maximum TSS of _analyzed_ positions: ' + str(maximum_TSS_of_analyzed_positions))
    mylogs.info('\n')
    mylogs.info('\n')

# Check number and status of mutated or static across each edge - instances where non-ACGT are found at ancestral or descendant node are ignored - build tables with mutational events

report_everything_about_selected_positions_results['Mutation'] = report_everything_about_selected_positions_results['Analyzed_strand_ancestral'] + report_everything_about_selected_positions_results['Analyzed_strand_descendent']
processed_DF = report_everything_about_selected_positions_results.groupby(['Edge_name','Mutation'])[['Analyzed_strand_ancestral']].count()
processed_DF = processed_DF.reset_index()

# Keep only those instances where ACGT are found at ancestral and descendant node

processed_DF_real_chars = processed_DF[(processed_DF['Mutation'] == 'AA')  |  \
    (processed_DF['Mutation'] == 'AC')  |   \
    (processed_DF['Mutation'] == 'AG')  |   \
    (processed_DF['Mutation'] == 'AT')  |   \
    (processed_DF['Mutation'] == 'CC')  |   \
    (processed_DF['Mutation'] == 'CA')  |   \
    (processed_DF['Mutation'] == 'CG')  |   \
    (processed_DF['Mutation'] == 'CT')  |   \
    (processed_DF['Mutation'] == 'GG')  |   \
    (processed_DF['Mutation'] == 'GA')  |   \
    (processed_DF['Mutation'] == 'GC')  |   \
    (processed_DF['Mutation'] == 'GT')  |   \
    (processed_DF['Mutation'] == 'TT')  |   \
    (processed_DF['Mutation'] == 'TA')  |   \
    (processed_DF['Mutation'] == 'TC')  |   \
    (processed_DF['Mutation'] == 'TG')]

# Organize counts based upon edge and mutation

processed_DF_real_chars = processed_DF_real_chars.sort_values(['Edge_name','Mutation'])
processed_DF_real_chars.columns = ['Edge_name','Mutation','Counts']
processed_DF_real_chars.reset_index(inplace = True)

print(processed_DF_real_chars)

# Generate dataframe with edge and mutation type and counts

edge_list = processed_DF_real_chars['Edge_name'].drop_duplicates().tolist()
mutation_list = ['AA','AC','AG','AT','CC','CA','CG','CT','GG','GA','GC','GT','TT','TA','TC','TG']
edge_versus_mutation_counts_DF = pd.DataFrame(index=[edge_list],columns=mutation_list)

for idx in range(len(processed_DF_real_chars)):
    mutation = processed_DF_real_chars.loc[idx,'Mutation']
    edge = processed_DF_real_chars.loc[idx,'Edge_name']
    count = processed_DF_real_chars.loc[idx,'Counts']
    edge_versus_mutation_counts_DF.loc[edge,mutation] = count

edge_versus_mutation_counts_DF = edge_versus_mutation_counts_DF.fillna(0)

# Generate columns in the correct locations to carry sums and fractions of nucleotides and changes

edge_versus_mutation_counts_DF['A_n'] = edge_versus_mutation_counts_DF['AA'] + edge_versus_mutation_counts_DF['AC'] + edge_versus_mutation_counts_DF['AG'] + edge_versus_mutation_counts_DF['AT']
edge_versus_mutation_counts_DF['C_n'] = edge_versus_mutation_counts_DF['CC'] + edge_versus_mutation_counts_DF['CA'] + edge_versus_mutation_counts_DF['CG'] + edge_versus_mutation_counts_DF['CT']
edge_versus_mutation_counts_DF['G_n'] = edge_versus_mutation_counts_DF['GG'] + edge_versus_mutation_counts_DF['GA'] + edge_versus_mutation_counts_DF['GC'] + edge_versus_mutation_counts_DF['GT']
edge_versus_mutation_counts_DF['T_n'] = edge_versus_mutation_counts_DF['TT'] + edge_versus_mutation_counts_DF['TA'] + edge_versus_mutation_counts_DF['TC'] + edge_versus_mutation_counts_DF['TG']
  
edge_versus_mutation_counts_DF['f_AA'] = edge_versus_mutation_counts_DF['AA'] / edge_versus_mutation_counts_DF['A_n']
edge_versus_mutation_counts_DF['f_AC'] = edge_versus_mutation_counts_DF['AC'] / edge_versus_mutation_counts_DF['A_n']
edge_versus_mutation_counts_DF['f_AG'] = edge_versus_mutation_counts_DF['AG'] / edge_versus_mutation_counts_DF['A_n']
edge_versus_mutation_counts_DF['f_AT'] = edge_versus_mutation_counts_DF['AT'] / edge_versus_mutation_counts_DF['A_n']

edge_versus_mutation_counts_DF['f_CC'] = edge_versus_mutation_counts_DF['CC'] / edge_versus_mutation_counts_DF['C_n']
edge_versus_mutation_counts_DF['f_CA'] = edge_versus_mutation_counts_DF['CA'] / edge_versus_mutation_counts_DF['C_n']
edge_versus_mutation_counts_DF['f_CG'] = edge_versus_mutation_counts_DF['CG'] / edge_versus_mutation_counts_DF['C_n']
edge_versus_mutation_counts_DF['f_CT'] = edge_versus_mutation_counts_DF['CT'] / edge_versus_mutation_counts_DF['C_n']

edge_versus_mutation_counts_DF['f_GG'] = edge_versus_mutation_counts_DF['GG'] / edge_versus_mutation_counts_DF['G_n']
edge_versus_mutation_counts_DF['f_GA'] = edge_versus_mutation_counts_DF['GA'] / edge_versus_mutation_counts_DF['G_n']
edge_versus_mutation_counts_DF['f_GC'] = edge_versus_mutation_counts_DF['GC'] / edge_versus_mutation_counts_DF['G_n']
edge_versus_mutation_counts_DF['f_GT'] = edge_versus_mutation_counts_DF['GT'] / edge_versus_mutation_counts_DF['G_n']

edge_versus_mutation_counts_DF['f_TT'] = edge_versus_mutation_counts_DF['TT'] / edge_versus_mutation_counts_DF['T_n']
edge_versus_mutation_counts_DF['f_TA'] = edge_versus_mutation_counts_DF['TA'] / edge_versus_mutation_counts_DF['T_n']
edge_versus_mutation_counts_DF['f_TC'] = edge_versus_mutation_counts_DF['TC'] / edge_versus_mutation_counts_DF['T_n']
edge_versus_mutation_counts_DF['f_TG'] = edge_versus_mutation_counts_DF['TG'] / edge_versus_mutation_counts_DF['T_n']


edge_versus_mutation_counts_DF['f_A_mut'] = edge_versus_mutation_counts_DF['f_AC'] + edge_versus_mutation_counts_DF['f_AG'] + edge_versus_mutation_counts_DF['f_AT']
edge_versus_mutation_counts_DF['f_C_mut'] = edge_versus_mutation_counts_DF['f_CA'] + edge_versus_mutation_counts_DF['f_CG'] + edge_versus_mutation_counts_DF['f_CT']
edge_versus_mutation_counts_DF['f_G_mut'] = edge_versus_mutation_counts_DF['f_GA'] + edge_versus_mutation_counts_DF['f_GC'] + edge_versus_mutation_counts_DF['f_GT']
edge_versus_mutation_counts_DF['f_T_mut'] = edge_versus_mutation_counts_DF['f_TA'] + edge_versus_mutation_counts_DF['f_TC'] + edge_versus_mutation_counts_DF['f_TG']

if verbose_flag == 'y':
    edge_versus_mutation_counts_DF.to_csv(outputpath + 'analysis_edge_versus_substitution.csv')

# Compare decay rates of nucleotides. Slopes are determined by fitting polynomials, then taking the maximum slope (values of first derivative formula ) close to the origin (X between -0.025 and 0.1).
# The initial slope, reflecting instantaneous mutational propensities, is the desired quantity.
# Some values behind 0 on the X axis are allowed in the test window, since there is some lag between when mutations occur and when they are detected.
# First-order linear regression is not used, because at extended times (nucleotide decays), multiple substitutions can cause plots to plateau, depending up on the regression input.

# Generate dataframe to store slopes

nucleotide_remain_list_X = ['f_A_mut_X','f_C_mut_X','f_G_mut_X','f_T_mut_X']
nucleotide_remain_list_y = ['f_A_mut_y','f_C_mut_y','f_G_mut_y','f_T_mut_y']
slopes_compare_deplete_DF = pd.DataFrame(index=nucleotide_remain_list_X,columns=nucleotide_remain_list_y)    

# Generate slopes from appropriate columns and 'max_slope_out' function

for i in range(4):
    for j in range(4):
        edge_versus_mutation_counts_DF_column_X = 36 + i
        edge_versus_mutation_counts_DF_column_Y = 36 + j
        max_slope_out(i,j)

# Determine which nucleotide is the most common at analyzed sites within the input sequences.
# The nucleotide with the greatest representation makes for the most appropriate independent variable across further comparisons (greater granularity over time).

med_A = edge_versus_mutation_counts_DF['A_n'].median()
med_C = edge_versus_mutation_counts_DF['C_n'].median()
med_G = edge_versus_mutation_counts_DF['G_n'].median()
med_T = edge_versus_mutation_counts_DF['T_n'].median()

if verbose_flag == 'y':
    mylogs.info("\nMedian number of 'A' found in ancestor at each edge : " + str(int(med_A))+'\n')
    mylogs.info("\nMedian number of 'C' found in ancestor at each edge : " + str(int(med_C))+'\n')
    mylogs.info("\nMedian number of 'G' found in ancestor at each edge : " + str(int(med_G))+'\n')
    mylogs.info("\nMedian number of 'T' found in ancestor at each edge : " + str(int(med_T))+'\n')

medians_N = []
medians_N = [med_A,med_C,med_G,med_T]

if max(medians_N) == medians_N[0]:
    nucleotide_baseline = 'A'
elif max(medians_N) == medians_N[1]:
    nucleotide_baseline = 'C'
elif max(medians_N) == medians_N[2]:
    nucleotide_baseline = 'G'
else:
    nucleotide_baseline = 'T'

if verbose_flag == 'y':
    mylogs.info('\nNucleotide baseline (used as independent variable in slope calculations) will be: ' + nucleotide_baseline+'\n')

slopes_compare_deplete_DF.to_csv(outputpath + 'slopes_compare_deplete_DF.csv')

sum_slopes_dN,instant_fA,instant_fC,instant_fG,instant_fT = total_inst_nuc_mut(nucleotide_baseline)

mylogs.info('\n')

# Following the above analysis of how often each nucleotide mutates compared to the others, one must determine what each nucleotide then tends to become.
# Again, the X-axis is chosen based upon the median number of each nucleotides across the analyzed positions of different sequences.

if nucleotide_baseline == 'A':
    m = 36
    what_A_becomes_list_X = ['f_A_mut_X_A']
    what_C_becomes_list_X = ['f_C_mut_X_A']
    what_G_becomes_list_X = ['f_G_mut_X_A']
    what_T_becomes_list_X = ['f_T_mut_X_A']

elif nucleotide_baseline == 'C':
    m = 37
    what_A_becomes_list_X = ['f_A_mut_X_C']
    what_C_becomes_list_X = ['f_C_mut_X_C']
    what_G_becomes_list_X = ['f_G_mut_X_C']
    what_T_becomes_list_X = ['f_T_mut_X_C']

elif nucleotide_baseline == 'G':
    m = 38
    what_A_becomes_list_X = ['f_A_mut_X_G']
    what_C_becomes_list_X = ['f_C_mut_X_G']
    what_G_becomes_list_X = ['f_G_mut_X_G']
    what_T_becomes_list_X = ['f_T_mut_X_G']

elif nucleotide_baseline == 'T':
    m = 39
    what_A_becomes_list_X = ['f_A_mut_X_T']
    what_C_becomes_list_X = ['f_C_mut_X_T']
    what_G_becomes_list_X = ['f_G_mut_X_T']
    what_T_becomes_list_X = ['f_T_mut_X_T']

what_A_becomes_list_y = ['to_A','to_C','to_G','to_T']
what_A_becomes_DF = pd.DataFrame(index=what_A_becomes_list_X,columns=what_A_becomes_list_y)

what_C_becomes_list_y = ['to_C','to_A','to_G','to_T']
what_C_becomes_DF = pd.DataFrame(index=what_C_becomes_list_X,columns=what_C_becomes_list_y)

what_G_becomes_list_y = ['to_G','to_A','to_C','to_T']
what_G_becomes_DF = pd.DataFrame(index=what_G_becomes_list_X,columns=what_G_becomes_list_y)

what_T_becomes_list_y = ['to_T','to_A','to_C','to_G']
what_T_becomes_DF = pd.DataFrame(index=what_T_becomes_list_X,columns=what_T_becomes_list_y)

def max_slope_out_nuc(k, anc_nuc, edge_versus_mutation_counts_DF_column_X, edge_versus_mutation_counts_DF_column_Y):
        
        edge_versus_mutation_counts_DF_column_X_label = edge_versus_mutation_counts_DF.columns[edge_versus_mutation_counts_DF_column_X]
        edge_versus_mutation_counts_DF_column_Y_label = edge_versus_mutation_counts_DF.columns[edge_versus_mutation_counts_DF_column_Y]
        
        fmut_X = edge_versus_mutation_counts_DF.iloc[:,edge_versus_mutation_counts_DF_column_X].to_numpy()
        fmut_y = edge_versus_mutation_counts_DF.iloc[:,edge_versus_mutation_counts_DF_column_Y].to_numpy()
        
        idx = np.isfinite(fmut_X) & np.isfinite(fmut_y)
        coeffs = np.polyfit(fmut_X[idx], fmut_y[idx], 3)

        coeff_x_power_3 = coeffs[0]
        coeff_x_power_2 = coeffs[1]
        coeff_x_power_1 = coeffs[2]
        coeff_0 = coeffs[3]

        deriv_power_2 = 3 * coeff_x_power_3
        deriv_power_1 = 2 * coeff_x_power_2
        deriv_0 = coeff_x_power_1

        test_slope_input = np.linspace (-0.025,0.1,100)
        calculate_slopes = lambda t: deriv_power_2*t*t + deriv_power_1*t + deriv_0
        slope_output = np.array([calculate_slopes(xi) for xi in test_slope_input])
        max_slope_output = np.amax(slope_output)

        if anc_nuc == 'A':
            what_A_becomes_DF.iloc[0,k] = max_slope_output
        elif anc_nuc == 'C':
            what_C_becomes_DF.iloc[0,k] = max_slope_output
        elif anc_nuc == 'G':
            what_G_becomes_DF.iloc[0,k] = max_slope_output
        elif anc_nuc == 'T':
            what_T_becomes_DF.iloc[0,k] = max_slope_output
        
        if verbose_flag == 'y':
            
            mylogs.info('')
            mylogs.info(edge_versus_mutation_counts_DF_column_X_label + ' versus ' + edge_versus_mutation_counts_DF_column_Y_label + ' polynomial: ' + str(round(coeff_x_power_3,2))+' x^3 + ' + str(round(coeff_x_power_2,2)) + ' x^2 + ' + str(round(coeff_x_power_1,2)) + ' x + ' + str(round(coeff_0,2))+'\n')
            mylogs.info(edge_versus_mutation_counts_DF_column_X_label + ' versus ' + edge_versus_mutation_counts_DF_column_Y_label + ' derivative polynomial: ' + str(round(deriv_power_2,2))+' x^2 + ' + str(round(deriv_power_1,2)) + ' x + ' + str(round(deriv_0,2))+'\n')
            mylogs.info(edge_versus_mutation_counts_DF_column_X_label + ' versus ' + edge_versus_mutation_counts_DF_column_Y_label + ' maximum slope between X of -0.025 and 0.1 (used for calculations): ' + str(round(max_slope_output,2))+'\n')
            mylogs.info('\n')

        return

list_of_p = [20, 24, 28, 32]
list_of_anc_nuc = ['A', 'C', 'G', 'T']

for n in range(0,4):
    p = list_of_p[n]
    anc_nuc = list_of_anc_nuc[n]

    for k in range(4):
        edge_versus_mutation_counts_DF_column_X = m
        edge_versus_mutation_counts_DF_column_Y = p + k
        max_slope_out_nuc(k, anc_nuc, edge_versus_mutation_counts_DF_column_X, edge_versus_mutation_counts_DF_column_Y)


# Combine output of N > X analyses into a single dataframe

what_N_becomes_slopes_DF = pd.concat([what_A_becomes_DF,what_C_becomes_DF,what_G_becomes_DF,what_T_becomes_DF],sort=True)

# Final calculations: multiply the fractions resulting from analysis of nucleotide decay by the fractions
# resulting from what those nucleotides become. Then, normalize so that the least common mutation is 
# represented by '1'.

if nucleotide_baseline == 'A':
    what_N_becomes_slopes_DF_X_label = ['f_A_mut_X_A','f_C_mut_X_A','f_G_mut_X_A','f_T_mut_X_A']
if nucleotide_baseline == 'C':
    what_N_becomes_slopes_DF_X_label = ['f_A_mut_X_C','f_C_mut_X_C','f_G_mut_X_C','f_T_mut_X_C']
if nucleotide_baseline == 'G':
    what_N_becomes_slopes_DF_X_label = ['f_A_mut_X_G','f_C_mut_X_G','f_G_mut_X_G','f_T_mut_X_G']
if nucleotide_baseline == 'T':
    what_N_becomes_slopes_DF_X_label = ['f_A_mut_X_T','f_C_mut_X_T','f_G_mut_X_T','f_T_mut_X_T']

sumAmut = what_N_becomes_slopes_DF.loc[what_N_becomes_slopes_DF_X_label[0],'to_C'] + what_N_becomes_slopes_DF.loc[what_N_becomes_slopes_DF_X_label[0],'to_G'] + what_N_becomes_slopes_DF.loc[what_N_becomes_slopes_DF_X_label[0],'to_T']

sumCmut = what_N_becomes_slopes_DF.loc[what_N_becomes_slopes_DF_X_label[1],'to_A'] + what_N_becomes_slopes_DF.loc[what_N_becomes_slopes_DF_X_label[1],'to_G'] + what_N_becomes_slopes_DF.loc[what_N_becomes_slopes_DF_X_label[1],'to_T']

sumGmut = what_N_becomes_slopes_DF.loc[what_N_becomes_slopes_DF_X_label[2],'to_A'] + what_N_becomes_slopes_DF.loc[what_N_becomes_slopes_DF_X_label[2],'to_C'] + what_N_becomes_slopes_DF.loc[what_N_becomes_slopes_DF_X_label[2],'to_T']

sumTmut = what_N_becomes_slopes_DF.loc[what_N_becomes_slopes_DF_X_label[3],'to_A'] + what_N_becomes_slopes_DF.loc[what_N_becomes_slopes_DF_X_label[3],'to_C'] + what_N_becomes_slopes_DF.loc[what_N_becomes_slopes_DF_X_label[3],'to_G']

fracAmut_S = what_N_becomes_slopes_DF.iloc[0] / sumAmut
fracCmut_S = what_N_becomes_slopes_DF.iloc[1] / sumCmut
fracGmut_S = what_N_becomes_slopes_DF.iloc[2] / sumGmut
fracTmut_S = what_N_becomes_slopes_DF.iloc[3] / sumTmut

final_A_to_C_not_normalized = fracAmut_S['to_C'] * instant_fA
final_A_to_G_not_normalized = fracAmut_S['to_G'] * instant_fA
final_A_to_T_not_normalized = fracAmut_S['to_T'] * instant_fA

final_C_to_A_not_normalized = fracCmut_S['to_A'] * instant_fC
final_C_to_G_not_normalized = fracCmut_S['to_G'] * instant_fC
final_C_to_T_not_normalized = fracCmut_S['to_T'] * instant_fC

final_G_to_A_not_normalized = fracGmut_S['to_A'] * instant_fG
final_G_to_C_not_normalized = fracGmut_S['to_C'] * instant_fG
final_G_to_T_not_normalized = fracGmut_S['to_T'] * instant_fG

final_T_to_A_not_normalized = fracTmut_S['to_A'] * instant_fT
final_T_to_C_not_normalized = fracTmut_S['to_C'] * instant_fT
final_T_to_G_not_normalized = fracTmut_S['to_G'] * instant_fT

final_N_to_N_not_normalized_list = [final_A_to_C_not_normalized,\
    final_A_to_G_not_normalized,\
    final_A_to_T_not_normalized,\
    final_C_to_A_not_normalized,\
    final_C_to_G_not_normalized,\
    final_C_to_T_not_normalized,\
    final_G_to_A_not_normalized,\
    final_G_to_C_not_normalized,\
    final_G_to_T_not_normalized,\
    final_T_to_A_not_normalized,\
    final_T_to_C_not_normalized,\
    final_T_to_G_not_normalized]

norm_factor = 1/min(final_N_to_N_not_normalized_list)


final_A_to_C_YES_normalized  = norm_factor * final_A_to_C_not_normalized
final_A_to_G_YES_normalized  = norm_factor * final_A_to_G_not_normalized
final_A_to_T_YES_normalized  = norm_factor * final_A_to_T_not_normalized

final_C_to_A_YES_normalized  = norm_factor * final_C_to_A_not_normalized
final_C_to_G_YES_normalized  = norm_factor * final_C_to_G_not_normalized
final_C_to_T_YES_normalized  = norm_factor * final_C_to_T_not_normalized

final_G_to_A_YES_normalized  = norm_factor * final_G_to_A_not_normalized
final_G_to_C_YES_normalized  = norm_factor * final_G_to_C_not_normalized
final_G_to_T_YES_normalized  = norm_factor * final_G_to_T_not_normalized

final_T_to_A_YES_normalized  = norm_factor * final_T_to_A_not_normalized
final_T_to_C_YES_normalized  = norm_factor * final_T_to_C_not_normalized
final_T_to_G_YES_normalized  = norm_factor * final_T_to_G_not_normalized

if verbose_flag == 'y':

    mylogs.info('\n________________________________________________\n')
    mylogs.info('Final mutational propensity approximations (strand analyzed):'+'\n')
    mylogs.info('\n')
    mylogs.info('Normalized A to C frequency: ' + '{:.{}f}'.format(final_A_to_C_YES_normalized,1)+'\n')
    mylogs.info('Normalized A to G frequency: ' + '{:.{}f}'.format(final_A_to_G_YES_normalized,1)+'\n')
    mylogs.info('Normalized A to T frequency: ' + '{:.{}f}'.format(final_A_to_T_YES_normalized,1)+'\n')

    mylogs.info('Normalized C to A frequency: ' + '{:.{}f}'.format(final_C_to_A_YES_normalized,1)+'\n')
    mylogs.info('Normalized C to G frequency: ' + '{:.{}f}'.format(final_C_to_G_YES_normalized,1)+'\n')
    mylogs.info('Normalized C to T frequency: ' + '{:.{}f}'.format(final_C_to_T_YES_normalized,1)+'\n')

    mylogs.info('Normalized G to A frequency: ' + '{:.{}f}'.format(final_G_to_A_YES_normalized,1)+'\n')
    mylogs.info('Normalized G to C frequency: ' + '{:.{}f}'.format(final_G_to_C_YES_normalized,1)+'\n')
    mylogs.info('Normalized G to T frequency: ' + '{:.{}f}'.format(final_G_to_T_YES_normalized,1)+'\n')

    mylogs.info('Normalized T to A frequency: ' + '{:.{}f}'.format(final_T_to_A_YES_normalized,1)+'\n')
    mylogs.info('Normalized T to C frequency: ' + '{:.{}f}'.format(final_T_to_C_YES_normalized,1)+'\n')
    mylogs.info('Normalized T to G frequency: ' + '{:.{}f}'.format(final_T_to_G_YES_normalized,1)+'\n')
    mylogs.info('\n')

    transition_over_transversion_ratio = (final_A_to_G_YES_normalized + \
        final_C_to_T_YES_normalized + \
        final_G_to_A_YES_normalized + \
        final_T_to_C_YES_normalized) / \
        \
        (final_A_to_C_YES_normalized + \
        final_A_to_T_YES_normalized + \
        final_C_to_A_YES_normalized + \
        final_C_to_G_YES_normalized + \
        final_G_to_C_YES_normalized + \
        final_G_to_T_YES_normalized + \
        final_T_to_A_YES_normalized + \
        final_T_to_G_YES_normalized)

    mylogs.info('\nTransition versus transversion ratio: ' + '{:.{}f}'.format(transition_over_transversion_ratio,1)+'\n')

output_table_list = []
output_table_list.append(final_A_to_C_YES_normalized)
output_table_list.append(final_A_to_G_YES_normalized)
output_table_list.append(final_A_to_T_YES_normalized)
output_table_list.append(final_C_to_A_YES_normalized)
output_table_list.append(final_C_to_G_YES_normalized)
output_table_list.append(final_C_to_T_YES_normalized)
output_table_list.append(final_G_to_A_YES_normalized)
output_table_list.append(final_G_to_C_YES_normalized)
output_table_list.append(final_G_to_T_YES_normalized)
output_table_list.append(final_T_to_A_YES_normalized)
output_table_list.append(final_T_to_C_YES_normalized)
output_table_list.append(final_T_to_G_YES_normalized)

output_table_SERIES = pd.Series(output_table_list,index = ['Normalized A to C frequency (Analyzed strand): ',\
    'Normalized A to G frequency (Analyzed strand): ',\
    'Normalized A to T frequency (Analyzed strand): ',\
    'Normalized C to A frequency (Analyzed strand): ',\
    'Normalized C to G frequency (Analyzed strand): ',\
    'Normalized C to T frequency (Analyzed strand): ',\
    'Normalized G to A frequency (Analyzed strand): ',\
    'Normalized G to C frequency (Analyzed strand): ',\
    'Normalized G to T frequency (Analyzed strand): ',\
    'Normalized T to A frequency (Analyzed strand): ',\
    'Normalized T to C frequency (Analyzed strand): ',\
    'Normalized T to G frequency (Analyzed strand): '])

output_table_SERIES.to_csv(outputpath + 'MADPROPS_final_output.csv', header=False)
