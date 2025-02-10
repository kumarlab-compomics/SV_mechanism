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

# In this script, we will take the simulated bed file from SURVIVOR to generate a vcf-like file, whereby we will need to refrieve the REF sequence, which we will use to annotate and compare to the real SV
# This script is only used when requesting simulations for deletions

# Data inputs sent in from the execution files
bedfile = str(argv[1])
parent = os.path.dirname(bedfile)
chromo = str(argv[2])
naming = str(argv[3])
naming = naming.replace('.param.txt','')
maxnum = str(argv[4])

# We need to clean up the bed file provided and add some columns to make it obvious it was a simulation
bed = pd.read_csv(bedfile, header= None, index_col=None, sep='\t')
bed = bed.rename(columns={0: 'CHROM', 1: 'POS', 4:'SV_Type' })
bed = bed[bed.SV_Type=='DEL']
bed['SVlen'] = bed[3]-bed['POS']
bed['SV_logic'] = True
bed['Sim'] = True
bed['Sim_ID'] = naming

# Grabbing the simulated inserted sequence from the fasta
for idx, row in bed.iterrows():
	bed.loc[idx, 'ID'] = row['CHROM']+'-'+str(row['POS'])+'-'+row['SV_Type']+'-'+str(row['SVlen'])
	bed.loc[idx, 'identifier'] = row['CHROM']+'-'+str(row['POS'])


# In this case, we need to pull the REF sequence using the reference genome and the locations from the bed file
# We interact with the shell and use samtools' faidx to pull sequences. We then add this to the REF column, as per a vcf-like format
for idx, row in bed.iterrows():
	os.system('/bin/bash -c "samtools faidx /home/nboev/projects/def-sushant/nboev/data/Genome/hg38.fa "' +row.CHROM+'":"' +str(row.POS)+ '"-"' +str(row.POS+row.SVlen)+'" -o "'+parent+'"/"'+str(row.CHROM)+str(row.POS)+str(row.SVlen)+'".fasta"')
	record = SeqIO.read(parent+"/"+str(row.CHROM)+str(row.POS)+str(row.SVlen)+".fasta", "fasta")
	bed.loc[idx, 'REF'] = str(record.seq).upper()
	bed.loc[idx, 'ALT'] = str(record.seq).upper()[0]
	os.system('/bin/bash -c "rm "'+parent+'"/"'+str(row.CHROM)+str(row.POS)+str(row.SVlen)+'".fasta"')

# We are only requesting the maxnum (I chose 100), simulations for each real SV
bed = bed.head(int(maxnum))
bed = bed[['CHROM', 'POS', 'ID', 'REF', 'ALT', 'SV_Type', 'SVlen', 'SV_logic', 'Sim', 'Sim_ID']]

# Saving this set of simulations
bed.to_csv(parent+'/'+naming+'.csv', sep='\t', index=False)
