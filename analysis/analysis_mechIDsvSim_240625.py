import glob
import numpy as np
import pandas as pd
from sys import argv
import datetime
from pytz import timezone
from scipy.stats import percentileofscore
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns

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
#print(result[result.ID == 'chrX-268206-DEL-103'].columns)
#print(result[result.ID == 'chrX-268206-DEL-103'][['Sim_ID']])

# Recall = we are going to use only the real SVs here and only using the ones where we had 101 sims
#result = result[result['Sim_ID'].str.contains('chr')==False]
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

# IN THIS VERSION, WE'RE GOING TO GET ALL THE LABELS SEPARATELY
# For labelling purposes, we will NOT be using the percentiles from the columns, just the raw value

for idx, row in forced.iterrows():
	if (row.prepost_Blastcoverage>0.25) & (row.prepost_pident>90):
		forced.loc[idx, 'mechID_homo'] = 'HR'
	elif (row.prepost_Blastcoverage==0):
		forced.loc[idx, 'mechID_homo'] = 'NHEJ'
	elif ((row.prepost_Blastcoverage>0.0025) & (row.prepost_Blastcoverage<0.2)) & (row.prepost_pident>80) :
		forced.loc[idx, 'mechID_homo'] = 'SSAaEJ'
	else:
		forced.loc[idx, 'mechID_homo'] = 'Undefined'
'''
for idx, row in forced.iterrows():
	if (row.STR_sum >1 ) & (row.Simple_repeat_sum > 1 ) & (row.Simple_repeat_logic == 1) :
		forced.loc[idx, 'mechID_slip'] = 'HSlip'
	else:
		forced.loc[idx, 'mechID_slip'] = 'LSlip'


for idx, row in forced.iterrows():
	if (row.RLoop >3 ) & (row.G4_sum > 1 ) :
		forced.loc[idx, 'mechID_rloop'] = 'HRLoop'
	else:
		forced.loc[idx, 'mechID_rloop'] = 'LRLoop'


forced['mechID'] = forced["mechID_homo"] + "_" + forced["mechID_slip"]+"_" + forced["mechID_rloop"]
'''

forced['ID'] = ids

print(forced.head())

forced.to_csv('/home/nboev/projects/def-sushant/nboev/analysis/'+project+'/'+loc+'/IDmechsvSIM/20240625/'+chromo+'.'+type+'_ID0625.tsv', sep='\t', index=False)

print('end of script')
print('END TIME:', datetime.datetime.now(timezone('EST')))

