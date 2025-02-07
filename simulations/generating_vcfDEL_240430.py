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

# this is for the processing of deletions only
bedfile = str(argv[1])
parent = os.path.dirname(bedfile)
print(parent)

chromo = str(argv[2])
naming = str(argv[3])
naming = naming.replace('.param.txt','')

# how many rows should we take at the end = we want a uniform number of simulations per real variant
maxnum = str(argv[4])

# cleaning up the bed file
bed = pd.read_csv(bedfile, header= None, index_col=None, sep='\t')
bed = bed.rename(columns={0: 'CHROM', 1: 'POS', 4:'SV_Type' })
bed = bed[bed.SV_Type=='DEL']
bed['SVlen'] = bed[3]-bed['POS']
bed['SV_logic'] = True
bed['Sim'] = True
bed['Sim_ID'] = naming

print(bed.head())

for idx, row in bed.iterrows():
	bed.loc[idx, 'ID'] = row['CHROM']+'-'+str(row['POS'])+'-'+row['SV_Type']+'-'+str(row['SVlen'])
	bed.loc[idx, 'identifier'] = row['CHROM']+'-'+str(row['POS'])

# this is not going to work! we need to add the dir name etc..

# we need to pull sequences from faidx
for idx, row in bed.iterrows():
	os.system('/bin/bash -c "samtools faidx /home/nboev/projects/def-sushant/nboev/data/Genome/hg38.fa "' +row.CHROM+'":"' +str(row.POS)+ '"-"' +str(row.POS+row.SVlen)+'" -o "'+parent+'"/"'+str(row.CHROM)+str(row.POS)+str(row.SVlen)+'".fasta"')
	record = SeqIO.read(parent+"/"+str(row.CHROM)+str(row.POS)+str(row.SVlen)+".fasta", "fasta")
	bed.loc[idx, 'REF'] = str(record.seq).upper()
	bed.loc[idx, 'ALT'] = str(record.seq).upper()[0]
	os.system('/bin/bash -c "rm "'+parent+'"/"'+str(row.CHROM)+str(row.POS)+str(row.SVlen)+'".fasta"')

print(bed.head())

# we only want the the max number of entries
bed = bed.head(int(maxnum))
print(bed)

bed = bed[['CHROM', 'POS', 'ID', 'REF', 'ALT', 'SV_Type', 'SVlen', 'SV_logic', 'Sim', 'Sim_ID']]
bed.to_csv(parent+'/'+naming+'.csv', sep='\t', index=False)


