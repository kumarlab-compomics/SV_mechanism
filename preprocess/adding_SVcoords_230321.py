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
svs = pd.read_csv(str(argv[1]), comment='#', sep='\t')
print('PRE ANNOTATION')
print(svs.head())
project = str(argv[2])
loc = str(argv[3])
filename = str(argv[4])
chr = str(argv[5])
size = str(argv[6])

# updated so we can alter the SV flank vals
def pulling(df, flank) :

	precoord = 'pre_'+str(flank)+'_coord'
	postcoord = 'post_'+str(flank)+'_coord'

# we're going to add coordinates we can add for our upstream and downstream flanks
# modifications after adding unresolved SVs --> ex. deletion and deletion_ME

	for idx, row in df.iterrows():
		if (row.SV_Type == 'deletion')|(row.SV_Type == 'DEL'):
			df.loc[idx, precoord] = int(row.POS - flank)
			df.loc[idx, postcoord] = int(row.POS + row.SVlen + flank)

		else :
			df.loc[idx, precoord] = int(row.POS - flank)
			df.loc[idx, postcoord] = int(row.POS + flank)

	return df


svs = pulling(svs, int(size))

# we're going to save a version I can use with bedtools
for idx, row in svs.iterrows():
	svs.loc[idx, 'unique_id'] = row['CHROM']+ '_' +str(row['POS'])+ '_' +str(row['SVlen'])+ '_' +row['SV_Type']

	if (row.SV_Type == 'deletion') | (row.SV_Type == 'DEL'):
		svs.loc[idx, 'preflank_start'] = row['pre_'+size+'_coord']
		svs.loc[idx, 'preflank_end'] = row.POS
		svs.loc[idx, 'postflank_start'] = row.POS + row.SVlen
		svs.loc[idx, 'postflank_end'] = row['post_'+size+'_coord']
	else:
		svs.loc[idx, 'preflank_start'] = row['pre_'+size+'_coord']
		svs.loc[idx, 'preflank_end'] = row.POS
		svs.loc[idx, 'postflank_start'] = row.POS
		svs.loc[idx, 'postflank_end'] = row['post_'+size+'_coord']


# we need this values to be integers
svs['preflank_start'] = svs['preflank_start'].apply(np.int64)
svs['preflank_end'] = svs['preflank_end'].apply(np.int64)
svs['postflank_start'] = svs['postflank_start'].apply(np.int64)
svs['postflank_end'] = svs['postflank_end'].apply(np.int64)


short = svs [['CHROM', 'preflank_start', 'preflank_end', 'postflank_start', 'postflank_end', 'unique_id']]
print(short.head())

starts = short[['CHROM','preflank_start', 'postflank_start', 'unique_id']]
ends = short[['CHROM','preflank_end', 'postflank_end', 'unique_id']]

starts = pd.melt(starts, id_vars=["CHROM", "unique_id"], value_name="Start")
ends = pd.melt(ends, id_vars=["CHROM", "unique_id"], value_name="End")
merged = pd.merge(starts, ends, on=['CHROM', 'unique_id'])


for idx, row in merged.iterrows():
	merged.loc[idx, 'difference'] = int(row.End - row.Start)

merged = merged[merged.difference == int(size)]
merged = merged[['CHROM', 'Start', 'End', 'unique_id']]

# getting rid of negatives
#for idc, row in merged.iterrows():
#	if int(row.Start) < 0:
#		merged.loc[idx, 'Start'] = 0
#	if int(row.End) < 0:
#		merged.loc[idx, 'End'] = 0

merged.Start = np.where(merged.Start < 0, 0, merged.Start)
merged.End = np.where(merged.End < 0, 0, merged.End)

# Saving the file itself --> this is the only one we'll save
merged.to_csv('/home/nboev/projects/def-sushant/nboev/preprocess/'+project+'/'+loc+'/SVcoords/' +chr+'/'+ filename+'.SVcoords' +size+ 'flank_bedtools.bed', sep='\t', index=False, header=False)


print('end of script')
print('END TIME:', datetime.datetime.now(timezone('EST')))
