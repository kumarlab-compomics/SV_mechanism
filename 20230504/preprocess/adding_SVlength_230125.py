#!/usr/bin/env python3

import pandas as pd
import glob
import numpy as np
from sys import argv
import datetime
from pytz import timezone

print('START TIME:', datetime.datetime.now(timezone('EST')))

# Getting the vcf and opening with pandas
print(argv[1])
svs = pd.read_csv(str(argv[1]), comment='#', sep='\t')
print('PRE ANNOTATION')
print(svs.head())
filename = str(argv[1]).split('/')[9]

# pulling SV_Type column ---> WE WILL NEED TO ALTER THIS!! NEEDS TO MATCH OTHERS, BUT WE NEED TO KNOW WHAT THE OPTIONS ARE (EX. INS WILL BECOME insertion)
svs['SV_Type'] = svs['INFO'].str.split(';', expand=True)[2]
svs['SV_Type'] = svs['SV_Type'].str.split('=', expand=True)[1]
print(svs.SV_Type.unique())

# pulling SVlen column
svs['SVlen'] = svs['INFO'].str.split(';', expand=True)[3]
svs['SVlen'] = svs['SVlen'].str.split('=', expand=True)[1]
svs['SVlen'] = abs(svs['SVlen'].apply(pd.to_numeric))


for idx, row in svs.iterrows():
	if row.SVlen >= 50:
		svs.loc[idx, 'SV_logic'] = True
	else:
		svs.loc[idx, 'SV_logic'] = False

for idx, row in svs.iterrows():
	if row.SV_Type == 'INS':
		svs.loc[idx, 'SV_Type'] = 'insertion'
	elif row.SV_Type == 'DEL':
		svs.loc[idx, 'SV_Type'] = 'deletion'

# post annotation head
print('POST ANNOTATION')
print(svs.head())

svs.to_csv('/home/nboev/projects/def-sushant/nboev/preprocess/HGSVC2_v2_integrated_callset/sv/SVlen/' +filename+ '.SVlength.csv', sep='\t', index=False)

print('end of script')
print('END TIME:', datetime.datetime.now(timezone('EST')))
