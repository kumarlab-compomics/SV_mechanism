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

# In this script, we use the +-2000bp flanks and SVs to generate DNA Shape profiles using DNA Shape R

svs = pd.read_csv(str(argv[1]), comment='#', sep='\t')
project = str(argv[2])
loc = str(argv[3])
filename = str(argv[4])
chr = str(argv[5])
size = str(argv[6])

'''
The function, pullShape, generates fasta files from the flanking and SV sequences, and sends them to the script, adding_SVDNAshapeR.R, to run DNAShapeR. Finally, the function deletes any intermediate files
The input: 
	df: The file to annotate
	flank: the length of the total flank+SV (ie. 2000bp)
'''

def pullShape (df, flank):
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
		elif (row.SV_Type == 'deletion') | (row.SV_Type == 'DEL')| (row.SV_Type == 'inversion') | (row.SV_Type == 'INV'):
			sv = '>SV_' +row['CHROM']+ '_' +str(row['POS'])+ '_' +str(row['SVlen'])+ '_' +row['SV_Type']+ '\n' +row['REF'] + '\n'

		with open('/home/nboev/projects/def-sushant/nboev/preprocess/'+project+'/'+loc+'/Blast/'+chr+'/premerge/'+filename+'SV_testing.fa', 'w+') as fh:
			fh.write(sv)

		os.system('/bin/bash -c "Rscript adding_SVDNAshapeR.R /home/nboev/projects/def-sushant/nboev/preprocess/"'+project+'"/"'+loc+'"/Blast/"'+chr+'"/premerge/"'+filename+'"SV_testing.fa /home/nboev/projects/def-sushant/nboev/preprocess/"'+project+'"/"'+loc+'"/Blast/"'+chr+'"/premerge/"'+filename+'"PRE_testing.fa /home/nboev/projects/def-sushant/nboev/preprocess/"'+project+'"/"'+loc+'"/Blast/"'+chr+'"/premerge/"'+filename+'"POST_testing.fa "'+ID+'" > /home/nboev/projects/def-sushant/nboev/preprocess/"'+project+'"/"'+loc+ '"/adding_SVDNAshapeR_230630.R.txt"')
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


# Calling the pullShape function
df = pullShape(svs, size)
df.to_csv('/home/nboev/projects/def-sushant/nboev/preprocess/'+project+'/'+loc+'/DNAshape/'+chr+'/postmerge/' +filename+ '.DNAshape.csv', sep='\t', index=False)

print('END TIME:', datetime.datetime.now(timezone('EST')))

