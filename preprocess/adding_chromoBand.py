#!/usr/bin/env python3

import pandas as pd
import glob
from sys import argv
import datetime
from pytz import timezone

print('START TIME:', datetime.datetime.now(timezone('EST')))

# In this script, we will determine the chromoband state of a given breakpoint as per ChromosomeBand_hg38.txt
# This file was downloaded here: http://genome.ucsc.edu/cgi-bin/hgTables

# Data inputs sent in from the execution files
svs = pd.read_csv(str(argv[1]), comment='#', sep='\t')
project = str(argv[2])
loc = str(argv[3])
filename = str(argv[4])
chr = str(argv[5])

# The file we are using to annotate. Essentially has the start/stop coordinates and the chromoband state/ Giemsa stain.
chromos = pd.read_csv('/home/nboev/projects/def-sushant/nboev/data/chromoBand/ChromosomeBand_hg38.txt', delimiter = "\t")
chromos = chromos[chromos['#chrom'] == chr]

ogrow = len(svs)
print('og rows:', ogrow)

df_svs = svs[svs['CHROM'] == chr]
df_chromos = chromos[chromos['#chrom'] == chr]

# We use the IntervalIndex function to create intervals in the index. This is then compared with the SV position to produce the gieStain column in the vcf/csv file
idx = pd.IntervalIndex.from_arrays( df_chromos['chromStart'], df_chromos['chromEnd'], closed='left')
df_chromos.index=idx
df_svs['gieStain']=df_chromos.loc[df_svs.POS,'gieStain'].values

newrow = len(df_svs)
print('new rows:', newrow)

# Checking if we lose any of the rows in this process. If the new and old number of rows match, we save the annotated file
if ogrow == newrow:
	df_svs.to_csv('/home/nboev/projects/def-sushant/nboev/preprocess/'+project+'/'+loc+'/chromoBand/'+chr+'/' +filename+ '.SVlength.csv', sep='\t', index=False)
else:
	print('poor merged')

print('END TIME:', datetime.datetime.now(timezone('EST')))
