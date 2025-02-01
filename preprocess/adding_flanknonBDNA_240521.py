#!/usr/bin/env python3

import pandas as pd
import glob
import numpy as np
from sys import argv
import datetime
from pytz import timezone

print('\n **********************')
print('START TIME:', datetime.datetime.now(timezone('EST')))

# Getting the vcf and opening with pandas
print(argv[1])
svs = pd.read_csv(str(argv[1]), comment='#', sep='\t')
project = str(argv[2])
loc = str(argv[3])
filename = str(argv[4])
chr = str(argv[5])

annot = pd.read_csv(str(argv[6]), comment='#', sep='\t', header=None)
nonb = str(argv[7])
print(svs.head())

ogrow = len(svs)
print('og rows:', ogrow)

# we need to add the unique identifier
# we're going to save a version I can use with bedtools
for idx, row in svs.iterrows():
	svs.loc[idx, 'unique_id'] = row['CHROM']+ '_' +str(row['POS'])+ '_' +str(row['SVlen'])+ '_' +row['SV_Type']

annot.rename(columns={0: 'CHROM', 1: 'Start', 2:'End', 3:'unique_id', 4:'count'}, inplace=True)

naming = {'g-quadruplex_forming_repeats':'G4', 'z-dna_motifs':'ZDNA', 'a-phased_repeats' : 'Arep', 'direct_repeats':'DR', 'inverted_repeats':'IR', 'mirror_repeats':'MR', 'short_tandem_repeats': 'STR'}

counted = annot.groupby(['CHROM', 'unique_id'])['count'].agg(['sum']).add_prefix(naming[nonb]+'_').reset_index()

for idx, row in counted.iterrows():
	if row[naming[nonb] + '_sum'] != 0:
		counted.loc[idx, naming[nonb]] = True
	else:
		counted.loc[idx, naming[nonb]] = False


merged = pd.merge(svs, counted, on='unique_id', how='left', indicator = True)
merged = merged.drop(['CHROM_y'], axis=1)
merged.rename(columns={'CHROM_x': 'CHROM'}, inplace=True)

print(merged.head())

newrow = len(merged)
print('new rows:', newrow)

if ogrow == newrow:
	print('good merge')
	merged.to_csv('/home/nboev/projects/def-sushant/nboev/preprocess/'+project+'/'+loc+'/nonBDNA/'+chr +'/postmerge/'+nonb+'.'+filename+'.flanks.csv', sep='\t', index=False)
else:
	print('poor merged')

print('end of script')
print('END TIME:', datetime.datetime.now(timezone('EST')))

