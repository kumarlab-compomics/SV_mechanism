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

# now we're going to combine this with the merging steps


print('END TIME:', datetime.datetime.now(timezone('EST')))
