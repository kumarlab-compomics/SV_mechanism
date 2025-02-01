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
size = str(argv[6])

def pullBlast (df, flank):
	prename = "pre_flank_seq_"+flank
	postname = "post_flank_seq_"+flank

	grabbing = []

	for idx, row in df[df.SV_logic == True].iterrows():
		ID = row['CHROM']+ '_' +str(row['POS'])+ '_' + str(row['SVlen'])+ '_' + row['SV_Type']

		pre = '>PRE_' +row['CHROM']+ '_' +str(row['POS'])+ '_' +str(row['SVlen'])+ '_' +row['SV_Type']+ '\n' +row[prename] + '\n'
		post = '>POST_' +row['CHROM']+ '_' +str(row['POS'])+ '_' +str(row['SVlen'])+ '_' +row['SV_Type']+ '\n' +row[postname]+ '\n'

		with open('/home/nboev/projects/def-sushant/nboev/preprocess/'+project+'/'+loc+'/Blast/'+chr+'/premerge/'+filename+'PRE_testing.fa', 'w+') as fh:
			fh.write(pre)
		with open('/home/nboev/projects/def-sushant/nboev/preprocess/'+project+'/'+loc+'/Blast/'+chr+'/premerge/'+filename+'POST_testing.fa', 'w+') as fh:
			fh.write(post)

		if (row.SV_Type == 'insertion') | (row.SV_Type == 'INS'):
			sv = '>SV_' +row['CHROM']+ '_' +str(row['POS'])+ '_' +str(row['SVlen'])+ '_' +row['SV_Type']+ '\n' +row['ALT'] + '\n'
		elif (row.SV_Type == 'deletion') | (row.SV_Type == 'DEL'):
			sv = '>SV_' +row['CHROM']+ '_' +str(row['POS'])+ '_' +str(row['SVlen'])+ '_' +row['SV_Type']+ '\n' +row['REF'] + '\n'

		with open('/home/nboev/projects/def-sushant/nboev/preprocess/'+project+'/'+loc+'/Blast/'+chr+'/premerge/'+filename+'SV_testing.fa', 'w+') as fh:
			fh.write(sv)

# NOTE ADDED ON 2023/10/24 = ADDED -ungapped in order to remove gapped alignments which may alter the way we calculate coverage
		os.system('/bin/bash -c "blastn -subject /home/nboev/projects/def-sushant/nboev/preprocess/"'+project+'"/"'+loc+ '"/Blast/"'+chr+'"/premerge/"'+filename+'"PRE_testing.fa -query /home/nboev/projects/def-sushant/nboev/preprocess/"'+project+'"/"'+loc+ '"/Blast/"'+chr+'"/premerge/"'+filename+'"SV_testing.fa -ungapped -outfmt 6 >> /home/nboev/projects/def-sushant/nboev/preprocess/"'+project+'"/"'+loc+ '"/Blast/"'+chr+'"/premerge/"'+filename+'"pre_sv_resultsungap.txt"')
		os.system('/bin/bash -c "blastn -subject /home/nboev/projects/def-sushant/nboev/preprocess/"'+project+'"/"'+loc+ '"/Blast/"'+chr+'"/premerge/"'+filename+'"POST_testing.fa -query /home/nboev/projects/def-sushant/nboev/preprocess/"'+project+'"/"'+loc+ '"/Blast/"'+chr+'"/premerge/"'+filename+'"SV_testing.fa -ungapped -outfmt 6 >> /home/nboev/projects/def-sushant/nboev/preprocess/"'+project+'"/"'+loc+ '"/Blast/"'+chr+'"/premerge/"'+filename+'"post_sv_resultsungap.txt"')
		os.system('/bin/bash -c "blastn -subject /home/nboev/projects/def-sushant/nboev/preprocess/"'+project+'"/"'+loc+ '"/Blast/"'+chr+'"/premerge/"'+filename+'"PRE_testing.fa -query /home/nboev/projects/def-sushant/nboev/preprocess/"'+project+'"/"'+loc+ '"/Blast/"'+chr+'"/premerge/"'+filename+'"POST_testing.fa -ungapped -outfmt 6 >> /home/nboev/projects/def-sushant/nboev/preprocess/"'+project+'"/"'+loc+ '"/Blast/"'+chr+'"/premerge/"'+filename+'"pre_post_resultsungap.txt"')

# in this version, were trying to run DNA shape at the same time using the fa we already have
		os.system('/bin/bash -c "Rscript adding_SVDNAshapeR_230630.R /home/nboev/projects/def-sushant/nboev/preprocess/"'+project+'"/"'+loc+'"/Blast/"'+chr+'"/premerge/"'+filename+'"SV_testing.fa /home/nboev/projects/def-sushant/nboev/preprocess/"'+project+'"/"'+loc+'"/Blast/"'+chr+'"/premerge/"'+filename+'"PRE_testing.fa /home/nboev/projects/def-sushant/nboev/preprocess/"'+project+'"/"'+loc+'"/Blast/"'+chr+'"/premerge/"'+filename+'"POST_testing.fa "'+ID+'" > /home/nboev/projects/def-sushant/nboev/preprocess/"'+project+'"/"'+loc+ '"/adding_SVDNAshapeR_230630.R.txt"')
		result = pd.read_csv('/home/nboev/projects/def-sushant/nboev/preprocess/'+project+'/'+loc+'/Blast/'+chr+'/premerge/'+filename+'SV_testing.csv', comment='#', sep=',')
		grabbing.append(result)

	for i in ['PRE', 'POST', 'SV']:
		os.system('/bin/bash -c "rm /home/nboev/projects/def-sushant/nboev/preprocess/"'+project+'"/"'+loc+'"/Blast/"'+chr+'"/premerge/"'+filename+i+'"_testing.fa"')
		os.system('/bin/bash -c "rm /home/nboev/projects/def-sushant/nboev/preprocess/"'+project+'"/"'+loc+'"/Blast/"'+chr+'"/premerge/"'+filename+i+'"_testing.fa.EP"')
		os.system('/bin/bash -c "rm /home/nboev/projects/def-sushant/nboev/preprocess/"'+project+'"/"'+loc+'"/Blast/"'+chr+'"/premerge/"'+filename+i+'"_testing.fa.HelT"')
		os.system('/bin/bash -c "rm /home/nboev/projects/def-sushant/nboev/preprocess/"'+project+'"/"'+loc+'"/Blast/"'+chr+'"/premerge/"'+filename+i+'"_testing.fa.MGW"')
		os.system('/bin/bash -c "rm /home/nboev/projects/def-sushant/nboev/preprocess/"'+project+'"/"'+loc+'"/Blast/"'+chr+'"/premerge/"'+filename+i+'"_testing.fa.ProT"')
		os.system('/bin/bash -c "rm /home/nboev/projects/def-sushant/nboev/preprocess/"'+project+'"/"'+loc+'"/Blast/"'+chr+'"/premerge/"'+filename+i+'"_testing.fa.Roll"')


	allresults = pd.concat(grabbing, axis=0, ignore_index=True)

	return allresults


# calling the Blast function
lol = pullBlast(svs, size)
print(lol.head())

lol.to_csv('/home/nboev/projects/def-sushant/nboev/preprocess/'+project+'/'+loc+'/DNAshape/'+chr+'/postmerge/' +filename+ '.DNAshape.csv', sep='\t', index=False)

print('END TIME:', datetime.datetime.now(timezone('EST')))

