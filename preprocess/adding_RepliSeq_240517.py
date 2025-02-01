#!/usr/bin/env python3

import pandas as pd
import glob
from sys import argv
import datetime
from pytz import timezone

print('START TIME:', datetime.datetime.now(timezone('EST')))

# loading in the Repliseq data
repli = pd.read_csv(str(argv[1]), sep=',')
svs = pd.read_csv(str(argv[2]), comment='#', sep='\t')
project = str(argv[3])
loc = str(argv[4])
filename = str(argv[5])
chr = str(argv[6])

repli = repli[repli['chr'] == chr]
svs = svs[svs['CHROM'] == chr]

ogrow = len(svs)
print('og rows:', ogrow)

idx = pd.IntervalIndex.from_arrays(repli['start'], repli['stop'], closed='left')
repli.index=idx

svs['S50']=repli.loc[svs.POS,'S50'].values

newrow = len(svs)
print('new rows:', newrow)

# checking if the rows match (they should)
if ogrow == newrow:
	print('good merge')
	svs.to_csv('/home/nboev/projects/def-sushant/nboev/preprocess/'+project+'/'+loc+'/RepliSeq/'+chr+'/' +filename+ '.withRepliSeq.csv', sep='\t', index=False)
else:
	print('poor merged')


print('END TIME:', datetime.datetime.now(timezone('EST')))



