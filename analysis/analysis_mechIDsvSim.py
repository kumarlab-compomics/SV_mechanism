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

# We use the prepost_Blastcoverage and prepost_pident columns from the Blast annotations (processing to calculate the coverage)
# Based on the thresholds, we label SVs with either a : HR, NHEJ, SSAaEJ or Undefined label. These labels are akin to the HLH, NLH, ILH and Undefined labels in the manuscript

for idx, row in forced.iterrows():
	if (row.prepost_Blastcoverage>0.25) & (row.prepost_pident>90):
		forced.loc[idx, 'mechID_homo'] = 'HR'
	elif (row.prepost_Blastcoverage==0):
		forced.loc[idx, 'mechID_homo'] = 'NHEJ'
	elif ((row.prepost_Blastcoverage>0.0025) & (row.prepost_Blastcoverage<0.2)) & (row.prepost_pident>80) :
		forced.loc[idx, 'mechID_homo'] = 'SSAaEJ'
	else:
		forced.loc[idx, 'mechID_homo'] = 'Undefined'

forced['ID'] = ids

# Save the final labelled file
forced.to_csv('/home/nboev/projects/def-sushant/nboev/analysis/'+project+'/'+loc+'/IDmechsvSIM/20240625/'+chromo+'.'+type+'_ID0625.tsv', sep='\t', index=False)

print('end of script')
print('END TIME:', datetime.datetime.now(timezone('EST')))

