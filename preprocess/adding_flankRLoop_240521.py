#!/usr/bin/env python3

import pandas as pd
import glob
import numpy as np
from sys import argv
import datetime
from pytz import timezone
import matplotlib.pyplot as plt
import matplotlib as mpl

print('\n**********************')
print('START TIME:', datetime.datetime.now(timezone('EST')))

# Getting the vcf and opening with pandas
print(argv[1])
svs = pd.read_csv(str(argv[1]), comment='#', sep='\t')
project = str(argv[2])
loc = str(argv[3])
filename = str(argv[4])
chr = str(argv[5])
print(svs.head())

ogrow = len(svs)
print('og rows:', ogrow)

# we need to add the unique identifier
# we're going to save a version I can use with bedtools
for idx, row in svs.iterrows():
	svs.loc[idx, 'unique_id'] = row['CHROM']+ '_' +str(row['POS'])+ '_' +str(row['SVlen'])+ '_' +row['SV_Type']

print(svs.head())

# Getting the file with flanked annotations and opening with pandas
level9 = pd.read_csv(str(argv[6]), comment='#', sep='\t', header=None)
level8 = pd.read_csv(str(argv[7]), comment='#', sep='\t', header=None)
level7 = pd.read_csv(str(argv[8]), comment='#', sep='\t', header=None)
level6 = pd.read_csv(str(argv[9]), comment='#', sep='\t', header=None)
level5 = pd.read_csv(str(argv[10]), comment='#', sep='\t', header=None)
level4 = pd.read_csv(str(argv[11]), comment='#', sep='\t', header=None)
level3 = pd.read_csv(str(argv[12]), comment='#', sep='\t', header=None)
level2 = pd.read_csv(str(argv[13]), comment='#', sep='\t', header=None)
level1 = pd.read_csv(str(argv[14]), comment='#', sep='\t', header=None)


def cleaning (rloop, lvl):
	rloop.rename(columns={0: 'CHROM', 1: 'Start', 2:'End', 3:'unique_id', 4:'count'}, inplace=True)
	counted = rloop.groupby(['CHROM', 'unique_id'])['count'].agg(['sum']).add_prefix(lvl+'_').reset_index()

	for idx, row in counted.iterrows():
		if row[lvl + '_sum'] != 0:
			counted.loc[idx, lvl] = True
		else:
			counted.loc[idx, lvl] = False

	return counted

count9 = cleaning(level9, 'lvl9')
count8 = cleaning(level8, 'lvl8')
count7 = cleaning(level7, 'lvl7')
count6 = cleaning(level6, 'lvl6')
count5 = cleaning(level5, 'lvl5')
count4 = cleaning(level4, 'lvl4')
count3 = cleaning(level3, 'lvl3')
count2 = cleaning(level2, 'lvl2')
count1 = cleaning(level1, 'lvl1')

def merging(df1, df2):
	df = pd.merge(df1, df2, on='unique_id', how='inner', suffixes=('', '_delme'))
	df = df[[c for c in df.columns if not c.endswith('_delme')]]

	return df

m1 = merging(count9, count8)
m2 = merging(m1, count7)
m3 = merging(m2, count6)
m4 = merging(m3, count5)
m5 = merging(m4, count4)
m6 = merging(m5, count3)
m7 = merging(m6, count2)
merged_df = merging(m7, count1)

for idx, row in merged_df.iterrows():
	if row.lvl9 == True:
		merged_df.loc[idx, 'RLoop'] = 9
	elif row.lvl8 == True:
		merged_df.loc[idx, 'RLoop'] = 8
	elif row.lvl7 == True:
		merged_df.loc[idx, 'RLoop'] = 7
	elif row.lvl6 == True:
		merged_df.loc[idx, 'RLoop'] = 6
	elif row.lvl5 == True:
		merged_df.loc[idx, 'RLoop'] = 5
	elif row.lvl4 == True:
		merged_df.loc[idx, 'RLoop'] = 4
	elif row.lvl3 == True:
		merged_df.loc[idx, 'RLoop'] = 3
	elif row.lvl2 == True:
		merged_df.loc[idx, 'RLoop'] = 2
	elif row.lvl1 == True:
		merged_df.loc[idx, 'RLoop'] = 1
	else:
		merged_df.loc[idx, 'RLoop'] = 0

print(merged_df.head())

# then we need to merge with the og svs then we need to check the lengths and save
newrow = len(merged_df)
print('new number of rows:', newrow)

# checking if the rows match (they should)
#if ogrow == newrow:
#	print('good merge')
merged_df.to_csv('/home/nboev/projects/def-sushant/nboev/preprocess/'+project+'/'+loc+'/RLoopForming/'+chr +'/postmerge/'+filename+'.RLoopflanks_merged.csv', sep='\t', index=False)
#else:
#	print('poor merged')

print('end of script')
print('END TIME:', datetime.datetime.now(timezone('EST')))

