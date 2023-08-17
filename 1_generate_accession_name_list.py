#!/usr/bin/env python

# coding: utf-8 

# Funding received from the European Research Council and the Sigrid Jusélius Foundation contributed to the development of this software.
# Author: Cory Dunn
# Institution: University of Helsinki
# Author Email: cory.david.dunn@gmail.com
# License: GPLv3

# Load libraries

import os
import re

# Make accession list from OrthoMAM FASTAs

path = '/Users/corydunn/Library/CloudStorage/Dropbox/Python_software/orthomam_AUG_13_2023/OrthoMAM10c_CDS/FASTA/'
files = os.listdir(path)
accessions = set()

for FASTA in files:
    file = open(path + FASTA, "r")
    for line in file:
        if re.search('>', line):
            accessions.add(line)
            print(line)
            
with open('/Users/corydunn/Library/CloudStorage/Dropbox/Python_software/orthomam_AUG_13_2023/OrthoMAM_accessions.txt', 'w') as f:
    for line in accessions:
        f.write(line)
f.close()