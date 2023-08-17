#!/usr/bin/env python

# coding: utf-8 

# Funding received from the European Research Council and the Sigrid Jusélius Foundation contributed to the development of this software.
# Author: Cory Dunn
# Institution: University of Helsinki
# Author Email: cory.david.dunn@gmail.com
# License: GPLv3

import os

inputpath = '/Users/corydunn/Library/CloudStorage/Dropbox/Python_software/orthomam_AUG_13_2023/OrthoMAM10c_CDS/FASTA/'
outputpath = '/Users/corydunn/Library/CloudStorage/Dropbox/Python_software/orthomam_AUG_13_2023/OrthoMAM10c_CDS/TABULAR/'
files = os.listdir(inputpath)

for FASTA in files:
    filestem = FASTA.rsplit( ".", 1 )[ 0 ]
    print('./seqkit fx2tab ' + inputpath + FASTA + ' > ' + outputpath + filestem + '.tabular')
    os.system('./seqkit fx2tab ' + inputpath + FASTA + ' > ' + outputpath + filestem + '.tabular')

