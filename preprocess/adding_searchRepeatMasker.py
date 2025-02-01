#!/usr/bin/env python3

import pandas as pd
import glob
import numpy as np
from sys import argv
import os
import datetime
from pytz import timezone

print('\n**********************')
print('START TIME:', datetime.datetime.now(timezone('EST')))

# In this script, we take in a vcf/csv batch file and annotate send the SV to RepeatMasker
# This annotation requires the use of RepeatMasker

svs = pd.read_csv(str(argv[1]), comment='#', sep='\t')
project = str(argv[2])
loc = str(argv[3])
filename = str(argv[4])
chr = str(argv[5])

# We iterate through every row. We obtain the REF or ALT sequence (deletion vs insertion) to create a fasta file. We append each sequence from this vcf/csv into this one fasta.
# This fasta is saved within the chromosome and split's directory. This fasta file is then sent to RepeatMasker. We use the human database for annotation.

for idx, row in svs.iterrows():
	if (row.SV_Type == 'insertion') | (row.SV_Type == 'INS'):
		search = '>'+row['CHROM']+ '_' +str(row['POS'])+ '_' +str(row['SVlen'])+ '_' +row['SV_Type']+ '\n' +row['ALT'] + '\n'
		with open('/home/nboev/projects/def-sushant/nboev/preprocess/'+project+'/'+loc+'/searchRepeatMasker/'+chr+'/'+filename+'/'+filename+'.fa', 'a') as fh:
			fh.write(search)

	elif (row.SV_Type == 'deletion') | (row.SV_Type == 'DEL'):
		search = '>'+row['CHROM']+ '_' +str(row['POS'])+ '_' +str(row['SVlen'])+ '_' +row['SV_Type']+ '\n' +row['REF'] + '\n'
		with open('/home/nboev/projects/def-sushant/nboev/preprocess/'+project+'/'+loc+'/searchRepeatMasker/'+chr+'/'+filename+'/'+filename+'.fa', 'a') as fh:
			fh.write(search)

os.system('/bin/bash -c "RepeatMasker -species human -s /home/nboev/projects/def-sushant/nboev/preprocess/"'+project+'"/"'+loc+'"/searchRepeatMasker/"'+chr+'"/"'+filename+'"/"'+filename+'".fa"')

print('END TIME:', datetime.datetime.now(timezone('EST')))
