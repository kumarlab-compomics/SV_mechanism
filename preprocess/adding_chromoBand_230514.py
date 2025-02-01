#!/usr/bin/env python3

import pandas as pd
import glob
from sys import argv
import datetime
from pytz import timezone

print('START TIME:', datetime.datetime.now(timezone('EST')))

#loading in the chromosome band data, we don't need this yet... our goal is just to see if we can load this in...
# loading in the csv (which holds all the chromo)
print(argv[1])
svs = pd.read_csv(str(argv[1]), comment='#', sep='\t')
project = str(argv[2])
loc = str(argv[3])
filename = str(argv[4])
chr = str(argv[5])
print(svs.head())

# we're only going to do one chromo at a time
chromos = pd.read_csv('/home/nboev/projects/def-sushant/nboev/data/chromoBand/ChromosomeBand_hg38.txt', delimiter = "\t")
chromos = chromos[chromos['#chrom'] == chr]

ogrow = len(svs)
print('og rows:', ogrow)

df_svs = svs[svs['CHROM'] == chr]
df_chromos = chromos[chromos['#chrom'] == chr]

idx = pd.IntervalIndex.from_arrays( df_chromos['chromStart'], df_chromos['chromEnd'], closed='left')
df_chromos.index=idx

df_svs['gieStain']=df_chromos.loc[df_svs.POS,'gieStain'].values
newrow = len(df_svs)
print('new rows:', newrow)

# checking if the rows match (they should)
# we need to assure we don't have over writing
if ogrow == newrow:
	print('good merge')
	df_svs.to_csv('/home/nboev/projects/def-sushant/nboev/preprocess/'+project+'/'+loc+'/chromoBand/'+chr+'/' +filename+ '.SVlength.csv', sep='\t', index=False)
else:
	print('poor merged')


print('END TIME:', datetime.datetime.now(timezone('EST')))



