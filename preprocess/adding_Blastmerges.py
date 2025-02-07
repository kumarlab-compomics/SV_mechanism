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

# In this script, we take the blast results produced from adding_BlastDNAShape.py, identify the greatest homology results, and merge them with the original vcf/csv file

# Data inputs sent in from the execution files
svs = pd.read_csv(str(argv[1]), comment='#', sep='\t')
project = str(argv[2])
loc = str(argv[3])
filename = str(argv[4])
chr = str(argv[5])
svs = svs[svs.CHROM == chr]

# Standardizing the ID column for merging later
for idx, row in svs.iterrows():
        svs.loc[idx, 'ID'] = row['CHROM']+ '_' +str(row['POS'])+ '_' + str(row['SVlen'])+ '_' + row['SV_Type']

'''
The function, loadin, takes one of the .txt files produced by Blast and identifies the best entries 
The input: 
	file: The blast result file
 	cond: The "type" of alignment conducted. Recall, this could include comparisons between: Upstream-SV, Downstream-SV, Upstream-Downstream 

Procedure: 
1. We create standardized column names we will use later for merging
2. We check if a file is empty (more relevant for small files where no alignment results are produced, more common in somatic cases)
3. If the file is not empty, we load and standardize column names
4. We then drop all duplicated column using the standardized ID, whereby the top entry is the best entry, as per Blast
'''

def loadin (file, cond):
	cols = ["qseqid" ,"sseqid", "pident", "length", "mismatch", "gapopen" ,"qstart", "qend", "sstart", "send", "evalue", "bitscore"]
	if os.stat(file).st_size == 0:
		print(file)
		df_filt = pd.DataFrame(columns=cols)
		df_filt['ID'] = 0

	else:
		df = pd.read_csv(str(file), sep='\t', header=None)
		trues = []
		for i in cols:
			trues.append(cond+i)

		df.columns = trues
		df[['kind','ID_tmp']] = df[trues[0]].str.split('_chr',expand=True)
		df_filt = df.drop_duplicates(subset=['ID_tmp'], keep='first')

		for idx, row in df_filt.iterrows():
			df_filt.loc[idx, 'ID'] = 'chr'+str(row.ID_tmp)

		df_filt = df_filt.drop(['kind', 'ID_tmp'], axis=1)
	return df_filt

# Calling the loadin function, via the Upstream-Downstream, Upstream-SV, and Downstream-SV Blast results
prepost = loadin(str(argv[6]), 'prepost_')
presv = loadin(str(argv[7]), 'presv_')
postsv = loadin(str(argv[8]), 'postsv_')

'''
The function, merging, helps do a series of merge iterations
The input: 
	blastdf: Dataframe to merge with the original vcf/csv file
	cond: The blast condition (ie. prepost_, presv_, postsv_)
'''

def merging (blastdf, cond):
	blastdf = blastdf.loc[:,~blastdf.columns.duplicated()].copy()

	# We merge the files, Recall, not ALL SVs produce , therefore we use a "left" merge to allow for this. 
	merged = pd.merge(svs, blastdf, how = 'left', on = 'ID', indicator = True)

	# We create a binary column to describe the absense/ presence of an alignment result. 
	for idx, row in merged.iterrows():
		if row['_merge'] == 'left_only':
			merged.loc[idx, cond+'BlastHomo'] = False
		elif row['_merge'] == 'both':
			merged.loc[idx, cond+'BlastHomo'] = True

	return merged

# Calling the merging function for the 3 types of alignments. We then create a "merged" dataframe to represent the merged file where all 3 alignment types are represented
merged_prepost = merging(prepost, 'prepost_')
merged_presv = merging(presv, 'presv_')
merged_postsv = merging(postsv, 'postsv_')

if len(merged_prepost) > 0:
	merge1 = pd.merge(merged_prepost,merged_presv,how = 'outer', on='ID')
	merged = pd.merge(merge1, merged_postsv,how = 'outer', on='ID')
else:
	merge1 = pd.merge(merged_postsv,merged_presv,how = 'outer', on='ID')
	merged = pd.merge(merge1, merged_prepost,how = 'outer', on='ID')

merged = merged.loc[:,~merged.columns.str.endswith('_x')]
merged = merged.loc[:,~merged.columns.str.endswith('_y')]

# A quality control check for the columns
checkingblastcols = ['prepost_length', 'presv_qend', 'presv_qstart', 'presv_gapopen', 'postsv_qend', 'postsv_qstart', 'postsv_gapopen', 'prepost_pident', 'presv_pident', 'postsv_pident']
for i in checkingblastcols:
	if i not in merged.columns:
		merged[i] = 0

# Saving the annotated dataframe
merged.to_csv('/home/nboev/projects/def-sushant/nboev/preprocess/'+project+'/'+loc+'/Blast/'+chr+'/postmerge/' +filename+ '.Blastmergesungap.csv', sep='\t', index=False)

print('end of script')
print('END TIME:', datetime.datetime.now(timezone('EST')))
