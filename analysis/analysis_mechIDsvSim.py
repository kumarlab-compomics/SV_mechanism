import glob
import numpy as np
import pandas as pd
from sys import argv
import datetime
from pytz import timezone
from scipy.stats import percentileofscore

print('\n **********************')
print('START TIME:', datetime.datetime.now(timezone('EST')))

# In this script, we take in a vcf/csv batch file that we have annotated and label the real SVs with their local homology label
# Recall, the file we are passing in contains both real and simulated SVs. Therefore, some filtering will be required.

# Data inputs sent in from the execution files
chromo = str(argv[2])
type = str(argv[3])
project = str(argv[4])
loc = str(argv[5])
result = pd.read_csv(str(argv[1]), comment='#', sep='\t')

# Here is the filtering. We require the Sim binary column to be False. We also require SURVIVOR was successfully able to generate 100 simulations for this real SV
result = result[result['Sim'] == False]
result = result[result.n_zscore==101]
ids = list(result.ID)
result[['chr','POS', 'type', 'len']] = result['ID'].str.split('-',expand=True)
result["len"] = pd.to_numeric(result["len"])
result = result.drop(['ID','chr', 'POS', 'type'] , axis=1)

forced = result

# We use the iter1_pident, iter1_length and iter1_prebrkpt_loose columns from the Blast annotations (processing to calculate the coverage)
# Based on the thresholds, we label SVs with either a : HLH, ILH, NLH, or Undefined label. 

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

# Save the final labelled file
forced.to_csv('/home/nboev/projects/def-sushant/nboev/analysis/'+project+'/'+loc+'/IDmechsvSIM/20250616/'+chromo+'.'+type+'_ID0616.tsv', sep='\t', index=False)

print('end of script')
print('END TIME:', datetime.datetime.now(timezone('EST')))

