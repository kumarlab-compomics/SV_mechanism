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

print(argv[1])
svs = pd.read_csv(str(argv[1]), comment='#', sep='\t')
project = str(argv[2])
loc = str(argv[3])
filename = str(argv[4])
chr = str(argv[5])
print(svs.head())
svs = svs[svs.CHROM == chr]

ogrow = len(svs)
print('og number of rows:', ogrow)

# we need to make columns in order to match the two csv regardless of alignment type
for idx, row in svs.iterrows():
        svs.loc[idx, 'ID'] = row['CHROM']+ '_' +str(row['POS'])+ '_' + str(row['SVlen'])+ '_' + row['SV_Type']

# Making a function for the loading in of a file:
# Having an issue where the Blast results return nothing back?
def loadin (file, cond):
	cols = ["qseqid" ,"sseqid", "pident", "length", "mismatch", "gapopen" ,"qstart", "qend", "sstart", "send", "evalue", "bitscore"]
	if os.stat(file).st_size == 0:
		print(file)
		df_filt = pd.DataFrame(columns=cols)
		df_filt['ID'] = 0

	else:
		df = pd.read_csv(str(file), sep='\t', header=None)
		print(df.head())
		trues = []
		for i in cols:
			trues.append(cond+i)

		print(trues)
		df.columns = trues
		df[['kind','ID_tmp']] = df[trues[0]].str.split('_chr',expand=True)

# taking the first entry / the "best" entry
		df_filt = df.drop_duplicates(subset=['ID_tmp'], keep='first')
		print(df_filt.head())

		for idx, row in df_filt.iterrows():
			df_filt.loc[idx, 'ID'] = 'chr'+str(row.ID_tmp)

		df_filt = df_filt.drop(['kind', 'ID_tmp'], axis=1)
	return df_filt

prepost = loadin(str(argv[6]), 'prepost_')
presv = loadin(str(argv[7]), 'presv_')
postsv = loadin(str(argv[8]), 'postsv_')


#we're going make a merge function
def merging (blastdf, cond):
	blastdf = blastdf.loc[:,~blastdf.columns.duplicated()].copy()

	merged = pd.merge(svs, blastdf, how = 'left', on = 'ID', indicator = True)

	for idx, row in merged.iterrows():
		if row['_merge'] == 'left_only':
			merged.loc[idx, cond+'BlastHomo'] = False
		elif row['_merge'] == 'both':
			merged.loc[idx, cond+'BlastHomo'] = True

	return merged

merged_prepost = merging(prepost, 'prepost_')
merged_presv = merging(presv, 'presv_')
merged_postsv = merging(postsv, 'postsv_')

if len(merged_prepost) > 0:
	merge1 = pd.merge(merged_prepost,merged_presv,how = 'outer', on='ID')
	merged = pd.merge(merge1, merged_postsv,how = 'outer', on='ID')
else:
	merge1 = pd.merge(merged_postsv,merged_presv,how = 'outer', on='ID')
	merged = pd.merge(merge1, merged_prepost,how = 'outer', on='ID')

#print(merged.head())
merged = merged.loc[:,~merged.columns.str.endswith('_x')]
merged = merged.loc[:,~merged.columns.str.endswith('_y')]

newrow = len(merged)
print('new number of rows:', newrow)

checkingblastcols = ['prepost_length', 'presv_qend', 'presv_qstart', 'presv_gapopen', 'postsv_qend', 'postsv_qstart', 'postsv_gapopen', 'prepost_pident', 'presv_pident', 'postsv_pident']

for i in checkingblastcols:
	if i not in merged.columns:
		merged[i] = 0

# checking if the rows match (they should)
#if ogrow == newrow:
#	print('good merge')
merged.to_csv('/home/nboev/projects/def-sushant/nboev/preprocess/'+project+'/'+loc+'/Blast/'+chr+'/postmerge/' +filename+ '.Blastmergesungap.csv', sep='\t', index=False)
#else:
#	print('poor merged')

print('end of script')
print('END TIME:', datetime.datetime.now(timezone('EST')))
