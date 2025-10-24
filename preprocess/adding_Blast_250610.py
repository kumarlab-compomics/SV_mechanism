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

# iteration 1 = 200bp = 175flank/25SV
# iteration 2 = percentage of the SV itself, 200bp flank + 10% of the SV
# iteration 3 = 2010bp = 2000flank/10SV

		if (row.SV_Type == 'insertion') | (row.SV_Type == 'INS'):
#			preseq1 = row[prename][-1975:]+row['ALT'][:25]
#			postseq1 = row['ALT'][-25:]+row[postname][:1975]
			preseq1 = row[prename][-1975:]+row['ALT'][:25]
			postseq1 = row['ALT'][-25:]+row[postname][:1975]

			per10 = int(0.1*row['SVlen'])
			preseq2 = row[prename][-200:]+row['ALT'][:per10]
			postseq2 = row['ALT'][-per10:]+row[postname][:200]

			preseq3 = row[prename][-2000:]+row['ALT'][:10]
			postseq3 = row['ALT'][-10:]+row[postname][:2000]

		elif (row.SV_Type == 'deletion') | (row.SV_Type == 'DEL')| (row.SV_Type == 'inversion') | (row.SV_Type == 'INV'):
#			preseq1 = row[prename][-1975:]+row['REF'][:25]
#			postseq1 = row['REF'][-25:]+row[postname][:1975]
			preseq1 = row[prename][-1975:]+row['REF'][:25]
			postseq1 = row['REF'][-25:]+row[postname][:1975]

			per10 = int(0.1*row['SVlen'])
			preseq2 = row[prename][-200:]+row['REF'][:per10]
			postseq2 = row['REF'][-per10:]+row[postname][:200]

			preseq3 = row[prename][-2000:]+row['REF'][:10]
			postseq3 = row['REF'][-10:]+row[postname][:2000]


# need to do all the iterations!
		pre1 = '>PRE_' +row['CHROM']+ '_' +str(row['POS'])+ '_' +str(row['SVlen'])+ '_' +row['SV_Type']+ '\n' +preseq1+ '\n'
		post1 = '>POST_' +row['CHROM']+ '_' +str(row['POS'])+ '_' +str(row['SVlen'])+ '_' +row['SV_Type']+ '\n' +postseq1+ '\n'

		pre2 = '>PRE_' +row['CHROM']+ '_' +str(row['POS'])+ '_' +str(row['SVlen'])+ '_' +row['SV_Type']+ '\n' +preseq2+ '\n'
		post2 = '>POST_' +row['CHROM']+ '_' +str(row['POS'])+ '_' +str(row['SVlen'])+ '_' +row['SV_Type']+ '\n' +postseq2+ '\n'

		pre3 = '>PRE_' +row['CHROM']+ '_' +str(row['POS'])+ '_' +str(row['SVlen'])+ '_' +row['SV_Type']+ '\n' +preseq3+ '\n'
		post3 = '>POST_' +row['CHROM']+ '_' +str(row['POS'])+ '_' +str(row['SVlen'])+ '_' +row['SV_Type']+ '\n' +postseq3+ '\n'

		with open('/home/nboev/projects/def-sushant/nboev/preprocess/'+project+'/'+loc+'/Blast/'+chr+'/premerge/'+filename+'PREiter1_testing.fa', 'w+') as fh:
			fh.write(pre1)
		with open('/home/nboev/projects/def-sushant/nboev/preprocess/'+project+'/'+loc+'/Blast/'+chr+'/premerge/'+filename+'POSTiter1_testing.fa', 'w+') as fh:
			fh.write(post1)

		with open('/home/nboev/projects/def-sushant/nboev/preprocess/'+project+'/'+loc+'/Blast/'+chr+'/premerge/'+filename+'PREiter2_testing.fa', 'w+') as fh:
			fh.write(pre2)
		with open('/home/nboev/projects/def-sushant/nboev/preprocess/'+project+'/'+loc+'/Blast/'+chr+'/premerge/'+filename+'POSTiter2_testing.fa', 'w+') as fh:
			fh.write(post2)

		with open('/home/nboev/projects/def-sushant/nboev/preprocess/'+project+'/'+loc+'/Blast/'+chr+'/premerge/'+filename+'PREiter3_testing.fa', 'w+') as fh:
			fh.write(pre3)
		with open('/home/nboev/projects/def-sushant/nboev/preprocess/'+project+'/'+loc+'/Blast/'+chr+'/premerge/'+filename+'POSTiter3_testing.fa', 'w+') as fh:
			fh.write(post3)

# NOTE ADDED ON 2023/10/24 = ADDED -ungapped in order to remove gapped alignments which may alter the way we calculate coverage
		os.system('/bin/bash -c "blastn -subject /home/nboev/projects/def-sushant/nboev/preprocess/"'+project+'"/"'+loc+ '"/Blast/"'+chr+'"/premerge/"'+filename+'"PREiter1_testing.fa -query /home/nboev/projects/def-sushant/nboev/preprocess/"'+project+'"/"'+loc+ '"/Blast/"'+chr+'"/premerge/"'+filename+'"POSTiter1_testing.fa -word_size 5 -perc_identity 70 -ungapped -outfmt 6 >> /home/nboev/projects/def-sushant/nboev/preprocess/"'+project+'"/"'+loc+ '"/Blast/"'+chr+'"/premerge/"'+filename+'"iter1_resultsungap.txt"')
#		os.system('/bin/bash -c "blastn -subject /home/nboev/projects/def-sushant/nboev/preprocess/"'+project+'"/"'+loc+ '"/Blast/"'+chr+'"/premerge/"'+filename+'"PREiter2_testing.fa -query /home/nboev/projects/def-sushant/nboev/preprocess/"'+project+'"/"'+loc+ '"/Blast/"'+chr+'"/premerge/"'+filename+'"POSTiter2_testing.fa -ungapped -word_size 4 -perc_identity 70 -outfmt 6 >> /home/nboev/projects/def-sushant/nboev/preprocess/"'+project+'"/"'+loc+ '"/Blast/"'+chr+'"/premerge/"'+filename+'"iter2_resultsungap.txt"')
#		os.system('/bin/bash -c "blastn -subject /home/nboev/projects/def-sushant/nboev/preprocess/"'+project+'"/"'+loc+ '"/Blast/"'+chr+'"/premerge/"'+filename+'"PREiter3_testing.fa -query /home/nboev/projects/def-sushant/nboev/preprocess/"'+project+'"/"'+loc+ '"/Blast/"'+chr+'"/premerge/"'+filename+'"POSTiter3_testing.fa -ungapped -word_size 4 -perc_identity 70 -outfmt 6 >> /home/nboev/projects/def-sushant/nboev/preprocess/"'+project+'"/"'+loc+ '"/Blast/"'+chr+'"/premerge/"'+filename+'"iter3_resultsungap.txt"')

	return


# calling the Blast function
pullBlast(svs, size)

print('END TIME:', datetime.datetime.now(timezone('EST')))

