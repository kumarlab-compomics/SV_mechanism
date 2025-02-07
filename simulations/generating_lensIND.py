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

# In this script, we take in a vcf/csv, and create matching parameter files, as per SURVIVOR requirements for SV simulation generation

svs = pd.read_csv(str(argv[1]), comment='#', sep='\t')
chromo = str(argv[2])
project = str(argv[4])
loc = str(argv[5])
svs = svs[(svs.CHROM == chromo)]

'''
The following loops through all the rows in the vcf/csv file to generate custome configuration files (parameter_file)
	Note: parameter_file was downloaded from: https://github.com/fritzsedlazeck/SURVIVOR/wiki
 1. We generate standardized IDs for each SV
 2. We generate start/stop SV lengths, as per the config file.
 3. We use a new location for the config file, to hold the correct INDEL_maximum_length and INDEL_maximum_length lengths, as per the SV. This is done by interacting with the shell
'''

for idx, row in svs.iterrows():
	if (row.SV_Type == 'deletion') |(row.SV_Type == 'DEL'):
		ID = row.CHROM +'-'+ str(row.POS+1) +'-'+ 'DEL' +'-'+ str(int(float(row.SVlen)))

		startlens = str(int(float(row.SVlen)))
		stoplens = str(int(float(row.SVlen))+1)

	elif (row.SV_Type == 'insertion') |(row.SV_Type == 'INS'):
		ID = row.CHROM +'-'+ str(row.POS+1) +'-'+ 'INS' +'-'+ str(int(float(row.SVlen)))

		startlens = str(int(float(row.SVlen))+1)
		stoplens = str(int(float(row.SVlen))+2)

	output_file = f'/home/nboev/scratch/data/SimulatedSVs/{project}/{loc}/parameterIND/{chromo}/{ID}.param.txt'
	first = subprocess.run(['sed', f's/INDEL_minimum_length: 50/INDEL_minimum_length: {startlens}/', str(argv[3])], stdout=subprocess.PIPE)
	subprocess.run(['sed', f's/INDEL_maximum_length: 51/INDEL_maximum_length: {stoplens}/'], input=first.stdout, stdout=open(output_file, 'w'))
