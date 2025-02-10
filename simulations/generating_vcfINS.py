#!/usr/bin/env python3

import pandas as pd
import glob
import os
from sys import argv
import datetime
from pytz import timezone
from Bio import SeqIO

print('\n ********************************')
print('START TIME:', datetime.datetime.now(timezone('EST')))

# In this script, we will take the simulated fasta and bed files from SURVIVOR to generate a vcf-like file, which we will use to annotate and compare to the real SV
# This script is only used when requesting simulations for insertions

# Data inputs sent in from the execution files
fastafile = str(argv[1])
bedfile = str(argv[2])
parent = os.path.dirname(bedfile)
chromo = str(argv[3])
naming = str(argv[4])
naming = naming.replace('.param.txt','')
maxnum = str(argv[5])

# We need to clean up the bed file provided and add some columns to make it obvious it was a simulation
bed = pd.read_csv(bedfile, header= None, index_col=None, sep='\t')
bed = bed.rename(columns={0: 'CHROM', 1: 'POS', 4:'SV_Type' })
bed = bed[bed.SV_Type=='INS']
bed['SVlen'] = bed[3]-bed['POS']-1
bed['SV_logic'] = True
bed['Sim'] = True
bed['Sim_ID'] = naming

# Grabbing the simulated inserted sequence from the fasta
identifiers = []
seq = []
for seq_record in SeqIO.parse(fastafile, 'fasta'):
	identifiers.append((seq_record.id).replace('_','-'))
	seq.append(str(seq_record.seq))

# Adding identifiers and standardizing the columns so it looks like a vcf file
for idx, row in bed.iterrows():
	bed.loc[idx, 'ID'] = row['CHROM']+'-'+str(row['POS'])+'-'+row['SV_Type']+'-'+str(row['SVlen'])
	bed.loc[idx, 'identifier'] = row['CHROM']+'-'+str(row['POS'])

for idx, row in bed.iterrows():
	for i in range(len(identifiers)):
		if row.identifier == identifiers[i]:
			bed.loc[idx, 'REF'] = seq[i][0]
			bed.loc[idx, 'ALT'] = seq[i]

# We are only requesting the maxnum (I chose 100), simulations for each real SV
bed = bed.head(int(maxnum))
bed = bed[['CHROM', 'POS', 'ID', 'REF', 'ALT', 'SV_Type', 'SVlen', 'SV_logic', 'Sim', 'Sim_ID']]

# Saving this set of simulations
bed.to_csv(parent+'/'+naming+'.csv', sep='\t', index=False)
