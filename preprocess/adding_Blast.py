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

# In this script, we're going to employ Blast to identify alignments between the 1975bp from the flanks + 25bp from the SVs. We do this using the up and downstream flanks 

svs = pd.read_csv(str(argv[1]), comment='#', sep='\t')
project = str(argv[2])
loc = str(argv[3])
filename = str(argv[4])
chr = str(argv[5])
size = str(argv[6])

'''
The function, pullBlast, pulls the sequences from the flanks and the SVs, and performs the Blast alignment through the command line. 
The Blast alignment is done using -word_size 5 -perc_identity 70 -ungapped -outfmt 6. Although in postprocessing we force a min of 80% identity
The input: 
	df: The file to annotate
	flank: the length of the total flank+SV (ie. 2000bp)
'''

def pullBlast (df, flank):
	prename = "pre_flank_seq_"+flank
	postname = "post_flank_seq_"+flank
	grabbing = []

	for idx, row in df[df.SV_logic == True].iterrows():
		ID = row['CHROM']+ '_' +str(row['POS'])+ '_' + str(row['SVlen'])+ '_' + row['SV_Type']

		if (row.SV_Type == 'insertion') | (row.SV_Type == 'INS'):
			preseq1 = row[prename][-1975:]+row['ALT'][:25]
			postseq1 = row['ALT'][-25:]+row[postname][:1975]

		elif (row.SV_Type == 'deletion') | (row.SV_Type == 'DEL')| (row.SV_Type == 'inversion') | (row.SV_Type == 'INV'):
			preseq1 = row[prename][-1975:]+row['REF'][:25]
			postseq1 = row['REF'][-25:]+row[postname][:1975]

		pre1 = '>PRE_' +row['CHROM']+ '_' +str(row['POS'])+ '_' +str(row['SVlen'])+ '_' +row['SV_Type']+ '\n' +preseq1+ '\n'
		post1 = '>POST_' +row['CHROM']+ '_' +str(row['POS'])+ '_' +str(row['SVlen'])+ '_' +row['SV_Type']+ '\n' +postseq1+ '\n'

		with open('/home/nboev/projects/def-sushant/nboev/preprocess/'+project+'/'+loc+'/Blast/'+chr+'/premerge/'+filename+'PREiter1_testing.fa', 'w+') as fh:
			fh.write(pre1)
		with open('/home/nboev/projects/def-sushant/nboev/preprocess/'+project+'/'+loc+'/Blast/'+chr+'/premerge/'+filename+'POSTiter1_testing.fa', 'w+') as fh:
			fh.write(post1)

		os.system('/bin/bash -c "blastn -subject /home/nboev/projects/def-sushant/nboev/preprocess/"'+project+'"/"'+loc+ '"/Blast/"'+chr+'"/premerge/"'+filename+'"PREiter1_testing.fa -query /home/nboev/projects/def-sushant/nboev/preprocess/"'+project+'"/"'+loc+ '"/Blast/"'+chr+'"/premerge/"'+filename+'"POSTiter1_testing.fa -word_size 5 -perc_identity 70 -ungapped -outfmt 6 >> /home/nboev/projects/def-sushant/nboev/preprocess/"'+project+'"/"'+loc+ '"/Blast/"'+chr+'"/premerge/"'+filename+'"iter1_resultsungap.txt"')
	return

# Calling the pullBlast function
pullBlast(svs, size)

print('END TIME:', datetime.datetime.now(timezone('EST')))

