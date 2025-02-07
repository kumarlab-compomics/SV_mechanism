#!/usr/bin/env python3
import pandas as pd
import glob
import os
from sys import argv
import datetime
from pytz import timezone
import subprocess
import numpy as np

print('\n ********************************')
print('START TIME:', datetime.datetime.now(timezone('EST')))

svs = pd.read_csv(str(argv[1]), comment='#', sep='\t')
print(len(svs))
chromo = str(argv[2])
print('chr:', chromo)
project = str(argv[4])
loc = str(argv[5])

print(list(svs.SV_Type.unique()))
svs = svs[(svs.CHROM == chromo)]
print(len(svs))
print(list(svs.SV_Type.unique()))
print(svs.head())
print(svs.columns)

# this time, we need to make a parameter file for each SV
for idx, row in svs.iterrows():
	print(row.CHROM, row.POS+1, row.SVlen, row.SV_Type)

	if (row.SV_Type == 'deletion') |(row.SV_Type == 'DEL'):
		ID = row.CHROM +'-'+ str(row.POS+1) +'-'+ 'DEL' +'-'+ str(int(float(row.SVlen)))

		startlens = str(int(float(row.SVlen)))
		stoplens = str(int(float(row.SVlen))+1)

	elif (row.SV_Type == 'insertion') |(row.SV_Type == 'INS'):
		ID = row.CHROM +'-'+ str(row.POS+1) +'-'+ 'INS' +'-'+ str(int(float(row.SVlen)))

		startlens = str(int(float(row.SVlen))+1)
		stoplens = str(int(float(row.SVlen))+2)

#       output_file = f'/home/nboev/projects/def-sushant/nboev/data/SimulatedSVs/{project}/{loc}/parameterIND/{chromo}/{ID}.param.txt'
	output_file = f'/home/nboev/scratch/data/SimulatedSVs/{project}/{loc}/parameterIND/{chromo}/{ID}.param.txt'

# need to make this correct!!
	first = subprocess.run(['sed', f's/INDEL_minimum_length: 50/INDEL_minimum_length: {startlens}/', str(argv[3])], stdout=subprocess.PIPE)
	subprocess.run(['sed', f's/INDEL_maximum_length: 51/INDEL_maximum_length: {stoplens}/'], input=first.stdout, stdout=open(output_file, 'w'))
