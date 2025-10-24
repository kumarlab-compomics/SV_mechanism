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
5. We check the orientation of the alignments (they need to be in opposite dirs), and we check if the breakpoint junctions have been included (ie. +-2bp)
'''

def loadin (file, cond):
		cols = ["qseqid" ,"sseqid", "pident", "length", "mismatch", "gapopen" ,"qstart", "qend", "sstart", "send", "evalue", "bitscore"]
		if os.stat(file).st_size == 0:
				df_filt = pd.DataFrame(columns=cols)
				df_filt['ID'] = 0

		else:
				df = pd.read_csv(str(file), sep='\t', header=None)
				trues = []
				for i in cols:
					trues.append(cond+i)

				df.columns = trues
				df[['kind','ID_tmp']] = df[trues[0]].str.split('_chr',expand=True)

				if cond == 'iter1_':
						preval = [1973, 1974, 1975, 1976, 1977]
						postval = [23, 24,25, 26, 27]

				df[cond + 'orient'] = (((df[cond + 'qend'] - df[cond + 'qstart']) * (df[cond + 'send'] - df[cond + 'sstart'])) > 0)
				df[cond+'prebrkpt_loose'] = df.apply(lambda row: any(row[cond+'sstart'] <= val <= row[cond+'send'] for val in preval),axis=1)
				df[cond+'postbrkpt_loose'] = df.apply(lambda row: any(row[cond+'qstart'] <= val <= row[cond+'qend'] for val in postval),axis=1)
				df_filt_loose = df[(df[cond+'orient']==True)&(df[cond+'prebrkpt_loose']==True)&(df[cond+'postbrkpt_loose']==True)].drop_duplicates(subset=['ID_tmp'], keep='first')

				for idx, row in df_filt_loose.iterrows():
					df_filt_loose.loc[idx, 'ID'] = 'chr'+str(row.ID_tmp)
				df_filt_loose = df_filt_loose.drop(['kind', 'ID_tmp'], axis=1)

	return df_filt_loose

# Calling the loadin function, via the Blast results
iter1 = loadin(str(argv[6]), 'iter1_')

'''
The function, merging, helps do a series of merge iterations
The input: 
	blastdf: Dataframe to merge with the original vcf/csv file
	cond: The breakpoint junction condition. In this case, we've only kept iteration 1
'''

def merging (blastdf, cond):
		if blastdf.empty:
			blastdf['ID'] = pd.Series(dtype='object')

		blastdf = blastdf.loc[:,~blastdf.columns.duplicated()].copy()
		merged = pd.merge(svs, blastdf, how = 'left', on = 'ID', indicator = True)
		return merged

# We then create a "merged" dataframe to represent the merged file 
merged = merging(iter1, 'iter1_loose_')

merged = merged.loc[:,~merged.columns.str.endswith('_x')]
merged = merged.loc[:,~merged.columns.str.endswith('_y')]

merged.columns = merged.columns.str.replace('_x$', '', regex=True)
merged.columns = merged.columns.str.replace('_y$', '', regex=True)

# Saving the annotated dataframe
merged.to_csv('/home/nboev/projects/def-sushant/nboev/preprocess/'+project+'/'+loc+'/Blast/'+chr+'/postmerge/' +filename+ '.Blastmergesungap_iter1.csv', sep='\t', index=False)

print('end of script')
print('END TIME:', datetime.datetime.now(timezone('EST')))
