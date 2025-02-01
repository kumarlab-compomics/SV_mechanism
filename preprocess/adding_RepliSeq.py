#!/usr/bin/env python3

import pandas as pd
import glob
from sys import argv
import datetime
from pytz import timezone

print('START TIME:', datetime.datetime.now(timezone('EST')))

# In this script, we take in a vcf/csv batch file and annotate the SV's position with it's replication timing (as measured by S50 from a RepliSeq experiment)

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

# We get the start and stop positions of the Repliseq and use these intervals to annotate the SV position
idx = pd.IntervalIndex.from_arrays(repli['start'], repli['stop'], closed='left')
repli.index=idx
svs['S50']=repli.loc[svs.POS,'S50'].values

newrow = len(svs)
print('new rows:', newrow)

# Checking if we lose any of the rows in this process. If the new and old number of rows match, we save the annotated file
if ogrow == newrow:
	svs.to_csv('/home/nboev/projects/def-sushant/nboev/preprocess/'+project+'/'+loc+'/RepliSeq/'+chr+'/' +filename+ '.withRepliSeq.csv', sep='\t', index=False)
else:
	print('poor merged')

print('END TIME:', datetime.datetime.now(timezone('EST')))
