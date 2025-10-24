import glob
import numpy as np
import pandas as pd
from sys import argv
import datetime
from pytz import timezone
from scipy.stats import percentileofscore

print('\n **********************')
print('START TIME:', datetime.datetime.now(timezone('EST')))

# Getting the merged feature matrix
print(argv[1])
chromo = str(argv[2])
type = str(argv[3])
project = str(argv[4])
loc = str(argv[5])

# this holds the identities of each sv
result = pd.read_csv(str(argv[1]), comment='#', sep='\t')
print(result.columns)
print(result.head())
print(len(result))

result = result[result['Sim'] == False]
result = result[result.n_zscore==101]
print(result.head())
print(len(result))

ids = list(result.ID)

result[['chr','POS', 'type', 'len']] = result['ID'].str.split('-',expand=True)
result["len"] = pd.to_numeric(result["len"])
result = result.drop(['ID','chr', 'POS', 'type'] , axis=1)

print(result.columns)
forced = result

for idx, row in forced.iterrows():
	if (row.iter1_length>100) & (row.iter1_pident>90)&(row['iter1_prebrkpt_loose']=='True')&(row['iter1_postbrkpt_loose']=='True'):
		forced.loc[idx, 'mechID_homo'] = 'HLH'
	elif (row.iter1_length==0)| (row['iter1_prebrkpt_loose']=='False') | (row['iter1_postbrkpt_loose']=='False')| (row['iter1_prebrkpt_loose']=='0') | (row['iter1_postbrkpt_loose']=='0')| ((row.iter1_pident<80)&(row['iter1_prebrkpt_loose']=='True')&(row['iter1_postbrkpt_loose']=='True')):
		forced.loc[idx, 'mechID_homo'] = 'NLH'
	elif ((row.iter1_length>20) & (row.iter1_length<=100))& (row.iter1_pident>80) &(row['iter1_prebrkpt_loose']=='True')&(row['iter1_postbrkpt_loose']=='True'):
		forced.loc[idx, 'mechID_homo'] = 'ILH'
	else:
		forced.loc[idx, 'mechID_homo'] = 'Undefined'

forced['ID'] = ids

print(forced.head())

forced.to_csv('/home/nboev/projects/def-sushant/nboev/analysis/'+project+'/'+loc+'/IDmechsvSIM/20250616/'+chromo+'.'+type+'_ID0616.tsv', sep='\t', index=False)

print('end of script')
print('END TIME:', datetime.datetime.now(timezone('EST')))

