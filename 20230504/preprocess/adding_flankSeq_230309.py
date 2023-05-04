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
print(svs.head())
print(svs.tail())

# ADDING JUST THE HEAD --> MUST BE REMOVED AFTER!!
svs = svs.head()

filename = str(argv[1]).split('/')[9]

ogrow = len(svs)
print('og rows:', ogrow)

# start with a subset of chrmo --> don't forget the Chromo Y!!
chromos = ['chr1', 'chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX', 'chrY']
print(chromos)

def addseq (dataframe, flank) :
	pre_str = "pre_flank_seq_"+str(flank)
	post_str = "post_flank_seq_"+str(flank)

# RECALL, WE'LL NEED TO CONCAT THESE TOGETHER FOR FINAL DF BEFORE RETURNING
	grab = []

	for k in chromos:
		print('chromosome:', k)
		df = dataframe[dataframe.CHROM == k]

		for idx, row in df.iterrows():
			if (row.SV_logic == True) & (row.SV_Type == 'deletion') :
				prestart = int(row.POS) - int(flank)
				os.system('/bin/bash -c "samtools faidx /home/nboev/projects/def-sushant/nboev/data/Genome/hg38.fa "' +k+ '":"' +str(prestart)+ '"-"' +str(row.POS)+ '" -o /home/nboev/projects/def-sushant/nboev/preprocess/HGSVC2_v2_integrated_callset/"'+str(argv[3])+'"/flankSeq/testing.fasta"')
				record = SeqIO.read("/home/nboev/projects/def-sushant/nboev/preprocess/HGSVC2_v2_integrated_callset/"+str(argv[3])+"/flankSeq/testing.fasta", "fasta")
#				print(record.id, record.seq)
				df.loc[idx, pre_str] = str(record.seq).upper()
				os.system('/bin/bash -c "rm /home/nboev/projects/def-sushant/nboev/preprocess/HGSVC2_v2_integrated_callset/"'+str(argv[3])+'"/flankSeq/testing.fasta"')

				poststart = int(row.POS) + int(row.SVlen)
				postend = int(row.POS) + int(row.SVlen) + int(flank)
				os.system('/bin/bash -c "samtools faidx /home/nboev/projects/def-sushant/nboev/data/Genome/hg38.fa "' +k+ '":"' +str(poststart)+ '"-"' +str(postend)+ '" -o /home/nboev/projects/def-sushant/nboev/preprocess/HGSVC2_v2_integrated_callset/"'+str(argv[3])+'"/flankSeq/testing.fasta"')
				record = SeqIO.read("/home/nboev/projects/def-sushant/nboev/preprocess/HGSVC2_v2_integrated_callset/"+str(argv[3])+"/flankSeq/testing.fasta", "fasta")
				df.loc[idx, post_str] = (record.seq).upper()
				os.system('/bin/bash -c "rm /home/nboev/projects/def-sushant/nboev/preprocess/HGSVC2_v2_integrated_callset/"'+str(argv[3])+'"/flankSeq/testing.fasta"')

			elif (row.SV_Type == 'insertion'):

				prestart = int(row.POS) - int(flank)
				os.system('/bin/bash -c "samtools faidx /home/nboev/projects/def-sushant/nboev/data/Genome/hg38.fa "' +k+ '":"' +str(prestart)+ '"-"' +str(row.POS)+ '" -o /home/nboev/projects/def-sushant/nboev/preprocess/HGSVC2_v2_integrated_callset/"'+str(argv[3])+'"/flankSeq/testing.fasta"')
				record = SeqIO.read("/home/nboev/projects/def-sushant/nboev/preprocess/HGSVC2_v2_integrated_callset/"+str(argv[3])+"/flankSeq/testing.fasta", "fasta")
				df.loc[idx, pre_str] = (record.seq).upper()
				os.system('/bin/bash -c "rm /home/nboev/projects/def-sushant/nboev/preprocess/HGSVC2_v2_integrated_callset/"'+str(argv[3])+'"/flankSeq/testing.fasta"')

				postend = int(row.POS) + int(flank)
				os.system('/bin/bash -c "samtools faidx /home/nboev/projects/def-sushant/nboev/data/Genome/hg38.fa "' +k+ '":"' +str(row.POS)+ '"-"' +str(postend)+ '" -o /home/nboev/projects/def-sushant/nboev/preprocess/HGSVC2_v2_integrated_callset/"'+str(argv[3])+'"/flankSeq/testing.fasta"')
				record = SeqIO.read("/home/nboev/projects/def-sushant/nboev/preprocess/HGSVC2_v2_integrated_callset/"+str(argv[3])+"/flankSeq/testing.fasta", "fasta")
				df.loc[idx, post_str] = (record.seq).upper()
				os.system('/bin/bash -c "rm /home/nboev/projects/def-sushant/nboev/preprocess/HGSVC2_v2_integrated_callset/"'+str(argv[3])+'"/flankSeq/testing.fasta"')

		grab.append(df)

	print('number of dfs', len(grab))
	graball = pd.concat(grab)
	return graball

# calling this function
newsvs = addseq(svs, str(argv[2]))
print(newsvs.head())
print(newsvs.tail())

newrow = len(newsvs)
print('new rows:', newrow)

# saving the csv
newsvs.to_csv('/home/nboev/projects/def-sushant/nboev/preprocess/HGSVC2_v2_integrated_callset/'+str(argv[3])+'/flankSeq/' +filename+ '.with' +str(argv[2])+ 'flankseq_samtools.csv', sep='\t', index=False)

print('END TIME:', datetime.datetime.now(timezone('EST')))
