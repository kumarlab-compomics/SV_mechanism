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

# loading in SVs
print(argv[1])
svs = pd.read_csv(str(argv[1]), comment='#', sep='\t')
project = str(argv[2])
loc = str(argv[3])
filename = str(argv[4])
chr = str(argv[5])
flank = str(argv[6])

ogrow = len(svs)
print('og rows:', ogrow)

pre_str = "pre_flank_seq_"+str(flank)
post_str = "post_flank_seq_"+str(flank)

df = svs[svs.CHROM == chr]

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



# calling this function
print(df.head())
newrow = len(df)
print('new rows:', newrow)

# saving the csv
df.to_csv('/home/nboev/projects/def-sushant/nboev/preprocess/'+project+'/'+loc+'/flankSeq/'+chr+'/' +filename+ '.with' +flank+'.flankseq_samtools.csv', sep='\t', index=False)

print('END TIME:', datetime.datetime.now(timezone('EST')))
