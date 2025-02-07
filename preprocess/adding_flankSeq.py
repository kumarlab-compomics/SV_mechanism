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

# In this script, we will determine the coordinates of the flanks, retrieve the matching sequences from the reference genome, and append these sequences to the original vcf/csv

# Data inputs sent in from the execution files
svs = pd.read_csv(str(argv[1]), comment='#', sep='\t')
project = str(argv[2])
loc = str(argv[3])
filename = str(argv[4])
chr = str(argv[5])
flank = str(argv[6])

# Designating the appropriate column names based on the requested flank length
pre_str = "pre_flank_seq_"+str(flank)
post_str = "post_flank_seq_"+str(flank)

df = svs[svs.CHROM == chr]

'''
The following loops through all the rows in the vcf/csv file to pull the reference sequence
1. We identify the pre-flanking position using the requested flank length
2. We interact with the shell system to call the samtools faidx function. Whereby, we provide the pre flanking coordinates and the HG38 reference genome and save the sequence in a fasta file
3. We read this fasta file using Bio's SeqIO function
4. We save this sequence in the pre_str column in a standardized format (upper case)
5. We delete the fasta file to mitigate space + overwriting issues

We do these steps for both the up and downstream flanks, along for insertions and deletions. 
	Note : The difference in the DEL vs INS procedure lies in the downstream flank coordinates, since deletions require the consideration of the SV length
'''

for idx, row in df[df.SV_logic == True].iterrows():
	if ((row.SV_Type == 'deletion') | (row.SV_Type == 'DEL') ) :
		prestart = int(row.POS) - int(flank)
		os.system('/bin/bash -c "samtools faidx /home/nboev/projects/def-sushant/nboev/data/Genome/hg38.fa "' +chr+ '":"' +str(prestart)+ '"-"' +str(row.POS)+ '" -o /home/nboev/projects/def-sushant/nboev/preprocess/"'+project+'"/"'+loc+'"/flankSeq/"'+chr+'"/"'+filename+'"testing.fasta"')
		record = SeqIO.read("/home/nboev/projects/def-sushant/nboev/preprocess/"+project+"/"+loc+"/flankSeq/"+chr+'/'+filename+"testing.fasta", "fasta")
		df.loc[idx, pre_str] = str(record.seq).upper()
		os.system('/bin/bash -c "rm /home/nboev/projects/def-sushant/nboev/preprocess/"'+project+'"/"'+loc+'"/flankSeq/"'+chr+'"/"'+filename+'"testing.fasta"')

		poststart = int(row.POS) + int(row.SVlen)
		postend = int(row.POS) + int(row.SVlen) + int(flank)
		os.system('/bin/bash -c "samtools faidx /home/nboev/projects/def-sushant/nboev/data/Genome/hg38.fa "' +chr+ '":"' +str(poststart)+ '"-"' +str(postend)+ '" -o /home/nboev/projects/def-sushant/nboev/preprocess/"'+project+'"/"'+loc+'"/flankSeq/"'+chr+'"/"'+filename+'"testing.fasta"')
		record = SeqIO.read("/home/nboev/projects/def-sushant/nboev/preprocess/"+project+"/"+loc+"/flankSeq/"+chr+'/'+filename+"testing.fasta", "fasta")
		df.loc[idx, post_str] = str(record.seq).upper()
		os.system('/bin/bash -c "rm /home/nboev/projects/def-sushant/nboev/preprocess/"'+project+'"/"'+loc+'"/flankSeq/"'+chr+'"/"'+filename+'"testing.fasta"')

	elif ((row.SV_Type == 'insertion') | (row.SV_Type == 'INS') ):
		prestart = int(row.POS) - int(flank)
		os.system('/bin/bash -c "samtools faidx /home/nboev/projects/def-sushant/nboev/data/Genome/hg38.fa "' +chr+ '":"' +str(prestart)+ '"-"' +str(row.POS)+ '" -o /home/nboev/projects/def-sushant/nboev/preprocess/"'+project+'"/"'+loc+'"/flankSeq/"'+chr+'"/"'+filename+'"testing.fasta"')
		record = SeqIO.read("/home/nboev/projects/def-sushant/nboev/preprocess/"+project+"/"+loc+"/flankSeq/"+chr+'/'+filename+"testing.fasta", "fasta")
		df.loc[idx, pre_str] = str(record.seq).upper()
		os.system('/bin/bash -c "rm /home/nboev/projects/def-sushant/nboev/preprocess/"'+project+'"/"'+loc+'"/flankSeq/"'+chr+'"/"'+filename+'"testing.fasta"')
		
		postend = int(row.POS) + int(flank)
		os.system('/bin/bash -c "samtools faidx /home/nboev/projects/def-sushant/nboev/data/Genome/hg38.fa "' +chr+ '":"' +str(row.POS)+ '"-"' +str(postend)+ '" -o /home/nboev/projects/def-sushant/nboev/preprocess/"'+project+'"/"'+loc+'"/flankSeq/"'+chr+'"/"'+filename+'"testing.fasta"')
		record = SeqIO.read("/home/nboev/projects/def-sushant/nboev/preprocess/"+project+"/"+loc+"/flankSeq/"+chr+'/'+filename+"testing.fasta", "fasta")
		df.loc[idx, post_str] = str(record.seq).upper()
		os.system('/bin/bash -c "rm /home/nboev/projects/def-sushant/nboev/preprocess/"'+project+'"/"'+loc+'"/flankSeq/"'+chr+'"/"'+filename+'"testing.fasta"')


# Saving the annotated dataframe
df.to_csv('/home/nboev/projects/def-sushant/nboev/preprocess/'+project+'/'+loc+'/flankSeq/'+chr+'/' +filename+ '.with' +flank+'.flankseq_samtools.csv', sep='\t', index=False)

print('END TIME:', datetime.datetime.now(timezone('EST')))
