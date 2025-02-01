#!/usr/bin/env python3

import pandas as pd
import glob
import numpy as np
from sys import argv
import datetime
from pytz import timezone
import itertools

print('\n **********************')
print('START TIME:', datetime.datetime.now(timezone('EST')))

# In this script, we take the output from RepeatMasker (SV ref/alt sequence annotations) and append the results to the original vcf/csv

svs = pd.read_csv(str(argv[1]), comment='#', sep='\t')
project = str(argv[2])
loc = str(argv[3])
filename = str(argv[4])
chr = str(argv[5])

svs = svs[svs.CHROM == chr]
ogrow = len(svs)
print('og number of rows:', ogrow)

# Addition of a standardized ID we can use later to merge with
for idx, row in svs.iterrows():
	svs.loc[idx, 'ID_match'] = row['CHROM']+ '_' +str(row['POS'])+ '_' + str(row['SVlen'])+ '_' + row['SV_Type']

# Parsing the RepeatMasker file
rep2 = pd.read_csv(str(argv[6]), header = None, skiprows = [0,1,2])
rep2 = rep2.replace(r"^ +| +$", r"", regex=True)
rep2 = rep2.replace({' +':' '},regex=True)
rep2 = rep2[0].str.split(' ',  expand=True)
rep2 = rep2.iloc[:,0:15]
rep2.columns = ['swScore', 'milliDiv', 'milliDel', 'milliIns', 'ID_match', 'genoStart', 'genoEnd', 'genoLeft', 'strand', 'repName', 'repClass', 'repStart', 'repEnd', 'repLeft', 'id']

# Grab unique IDs from the RepeatMasker file.
# We will loop through these IDs to append a list of unique classes and names to each SV. This therefore is a record of which classes and names were present in the SV
unique_svs = list(rep2.ID_match.unique())

for i in unique_svs:
	grab_class = []
	grab_class.append(list(rep2[rep2.ID_match ==i ].repClass.values))
	merged_class = list(set(list(itertools.chain(*grab_class))))
	svs.loc[svs.ID_match ==i, 'RepMasker_repClass'] = str(merged_class)

	grab_names = []
	grab_names.append(list(rep2[rep2.ID_match ==i ].repName.values))
	merged_names = list(set(list(itertools.chain(*grab_names))))
	svs.loc[svs.ID_match ==i, 'RepMasker_repName'] = str(merged_names)

newrow = len(svs)
print('new number of rows:', newrow)

# If the new number of rows matches the original number of rows, save the annotated dataframe
if ogrow == newrow:
	svs.to_csv('/home/nboev/projects/def-sushant/nboev/preprocess/'+project+'/'+loc+'/searchRepeatMasker/'+chr+ '/' +filename+'/'+filename+'.searchRepMasker.csv', sep='\t', index=False)
else:
	print('poor merged')

print('END TIME:', datetime.datetime.now(timezone('EST')))


