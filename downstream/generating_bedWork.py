#!/usr/bin/env python3

import pandas as pd
import glob
import os
from sys import argv
import datetime
from pytz import timezone
from Bio import SeqIO

print('\n ********************************')
print('START TIME:', datetime.datetime.now(timezone('EST')))

# In this script, we will generate bed files for each homology-based label from the real SVs. This will be used later for comparison against simulated SVs to calculate enrichment scores (ES). 

# Data inputs sent in from the execution files
svs = pd.read_csv(str(argv[1]), comment='#', sep='\t')
project = str(argv[2])
loc = str(argv[3])
chromo = str(argv[4])
typer = str(argv[5])

# Splitting the ID column
svs[['chr','POS', 'type', 'len']] = svs['ID'].str.split('-',expand=True)
svs["len"] = pd.to_numeric(svs["len"])
svs["POS"] = pd.to_numeric(svs["POS"])

# We loop through all the rows in the svs dataframe. 
# For deletions: We can use the SV's position and length to understand the deleted region (this is what we used for ES), OR, we can use the up and downstream breakpoints (this is what we used for ChromHMM analysis)
# For insertions: We only used the SV's position
for idx, row in svs.iterrows():
	if typer == 'deletion':
		svs.loc[idx, 'start']= int(row.POS)
		svs.loc[idx, 'stop']= int(row.POS + row.len)

		svs.loc[idx, 'start_up']= int(row.POS)
		svs.loc[idx, 'stop_up']= int(row.POS) +1
		svs.loc[idx, 'start_down']= int(row.POS + row.len)-1
		svs.loc[idx, 'stop_down']= int(row.POS + row.len)

	elif typer == 'insertion':
		svs.loc[idx, 'start']= int(row.POS)
		svs.loc[idx, 'stop']= int(row.POS)+1

svs['start'] = svs['start'].astype('Int64')
svs['stop'] = svs['stop'].astype('Int64')

'''
The function, savingdf, takes the svs from above and filters it for the requested homology-based label and saves the results in a format akin to a bed file
The input: 
	homostr (options): 'HLH', 'ILH', 'NLH'
'''

def savingdf (homostr):
	filt = svs[svs.mechID_homo==homostr]
	filt = filt[['chr', 'start', 'stop']]
	filt.to_csv('/home/nboev/projects/def-sushant/nboev/analysis/'+project+'/'+loc+'/IDmechsvSIM/20240625/splitsChromo/bedWork/'+chromo+'.'+typer+'.HOMO.'+homostr+'.bed', sep='\t', index=False, header=False)
	return

# Calling the function for all three labels
savingdf(HLH)
savingdf(ILH)
savingdf(NLH)


# The following steps are similar to those above, however, in this case we save the up and downstream bed files separatedly for deletions
work=['HLH', 'ILH', 'NLH']
if typer == 'deletion':
	svs['start_up'] = svs['start_up'].astype('Int64')
	svs['stop_up'] = svs['stop_up'].astype('Int64')
	svs['start_down'] = svs['start_down'].astype('Int64')
	svs['stop_down'] = svs['stop_down'].astype('Int64')

	for i in work:
		svshr = svs[svs.mechID_homo==i]
		svshrup = svshr[['chr', 'start_up', 'stop_up']]
		svshrdown = svshr[['chr', 'start_down', 'stop_down']]

		svshrup.to_csv('/home/nboev/projects/def-sushant/nboev/analysis/'+project+'/'+loc+'/IDmechsvSIM/20240625/splitsChromo/bedWork/UP'+chromo+'.'+typer+'.HOMO.'+i+'.bed', sep='\t', index=False, header=False)
		svshrdown.to_csv('/home/nboev/projects/def-sushant/nboev/analysis/'+project+'/'+loc+'/IDmechsvSIM/20240625/splitsChromo/bedWork/DOWN'+chromo+'.'+typer+'.HOMO.'+i+'.bed', sep='\t', index=False, header=False)

