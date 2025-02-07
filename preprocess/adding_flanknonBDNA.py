#!/usr/bin/env python3

import pandas as pd
import glob
import numpy as np
from sys import argv
import datetime
from pytz import timezone

print('\n **********************')
print('START TIME:', datetime.datetime.now(timezone('EST')))

# In this script, we take the bed file counts from nonBDNA annotations, and merge with the vcf. Our goal is to determine the number of motifs obtained in the flanking region. Very similar to adding_flankRepeatMasker.py 

# Data inputs sent in from the execution files
svs = pd.read_csv(str(argv[1]), comment='#', sep='\t')
project = str(argv[2])
loc = str(argv[3])
filename = str(argv[4])
chr = str(argv[5])

ogrow = len(svs)
print('og rows:', ogrow)

# Adding an unique identifier that we can use to match across the flanking results and the bed file results
for idx, row in svs.iterrows():
	svs.loc[idx, 'unique_id'] = row['CHROM']+ '_' +str(row['POS'])+ '_' +str(row['SVlen'])+ '_' +row['SV_Type']

# The bed file itself along with the name of the nonBDNA class
annot = pd.read_csv(str(argv[6]), comment='#', sep='\t', header=None)
nonb = str(argv[7])

# Standardizing the naming convention and determining the presence/absence of the motif and the number of counts per unique SV (across 2 flanks)
annot.rename(columns={0: 'CHROM', 1: 'Start', 2:'End', 3:'unique_id', 4:'count'}, inplace=True)
naming = {'g-quadruplex_forming_repeats':'G4', 'z-dna_motifs':'ZDNA', 'a-phased_repeats' : 'Arep', 'direct_repeats':'DR', 'inverted_repeats':'IR', 'mirror_repeats':'MR', 'short_tandem_repeats': 'STR'}
counted = annot.groupby(['CHROM', 'unique_id'])['count'].agg(['sum']).add_prefix(naming[nonb]+'_').reset_index()

for idx, row in counted.iterrows():
	if row[naming[nonb] + '_sum'] != 0:
		counted.loc[idx, naming[nonb]] = True
	else:
		counted.loc[idx, naming[nonb]] = False

# Merging with the original vcf/csv file
merged = pd.merge(svs, counted, on='unique_id', how='left', indicator = True)
merged = merged.drop(['CHROM_y'], axis=1)
merged.rename(columns={'CHROM_x': 'CHROM'}, inplace=True)

newrow = len(merged)
print('new rows:', newrow)

# Checking if we lose any of the rows in this process. If the new and old number of rows match, we save the annotated file
if ogrow == newrow:
	merged.to_csv('/home/nboev/projects/def-sushant/nboev/preprocess/'+project+'/'+loc+'/nonBDNA/'+chr +'/postmerge/'+nonb+'.'+filename+'.flanks.csv', sep='\t', index=False)
else:
	print('poor merged')

print('end of script')
print('END TIME:', datetime.datetime.now(timezone('EST')))

