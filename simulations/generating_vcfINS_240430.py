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

# this is for the processing of insertions only
# well need the fasta and bed files

fastafile = str(argv[1])
bedfile = str(argv[2])
parent = os.path.dirname(bedfile)
print(parent)

chromo = str(argv[3])
naming = str(argv[4])
naming = naming.replace('.param.txt','')

# how many rows should we take at the end = we want a uniform number of simulations per real variant
maxnum = str(argv[5])

# cleaning up the bed file
bed = pd.read_csv(bedfile, header= None, index_col=None, sep='\t')
bed = bed.rename(columns={0: 'CHROM', 1: 'POS', 4:'SV_Type' })
bed = bed[bed.SV_Type=='INS']
bed['SVlen'] = bed[3]-bed['POS']-1
bed['SV_logic'] = True
bed['Sim'] = True
bed['Sim_ID'] = naming

print(bed.head())

identifiers = []
seq = []
for seq_record in SeqIO.parse(fastafile, 'fasta'):
	identifiers.append((seq_record.id).replace('_','-'))
	seq.append(str(seq_record.seq))


for idx, row in bed.iterrows():
	bed.loc[idx, 'ID'] = row['CHROM']+'-'+str(row['POS'])+'-'+row['SV_Type']+'-'+str(row['SVlen'])
	bed.loc[idx, 'identifier'] = row['CHROM']+'-'+str(row['POS'])

# adding insertion sequences
for idx, row in bed.iterrows():
	for i in range(len(identifiers)):
		if row.identifier == identifiers[i]:
			bed.loc[idx, 'REF'] = seq[i][0]
			bed.loc[idx, 'ALT'] = seq[i]

# we only want the the max number of entries
bed = bed.head(int(maxnum))
print(bed)

bed = bed[['CHROM', 'POS', 'ID', 'REF', 'ALT', 'SV_Type', 'SVlen', 'SV_logic', 'Sim', 'Sim_ID']]
bed.to_csv(parent+'/'+naming+'.csv', sep='\t', index=False)

