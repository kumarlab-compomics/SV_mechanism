import glob
import numpy as np
import pandas as pd
from sys import argv
import datetime
from pytz import timezone
from scipy.stats import fisher_exact
import os

print('\n **********************')
print('START TIME:', datetime.datetime.now(timezone('EST')))

# In this script, we will generate matching bed files for the simulated SVs corresponding to real SVs based on their homology-based label. This will be used later for comparison against simulated SVs to calculate enrichment scores (ES). 

# Data inputs sent in from the execution files. This includes the real SVs, along with the z-score merged file which holds the real and simulated SVs
realspreds = pd.read_csv(str(argv[1]), comment='#', sep='\t')
zscores = pd.read_csv(str(argv[2]), comment='#', sep='\t')
project = str(argv[3])
loc = str(argv[4])
chromo = str(argv[5])
typer = str(argv[6])
pred = str(argv[7])

# We filter for SVs which have 101 entries (ie. 100 sims, 1 real), along with the "Sim" binary flag. 
# Note, this script contains a few important print statements. You want the size of the simulated files to be 100x the size of the real files. 
zscores_sims = zscores[(zscores.n_zscore == 101)&(zscores['Sim']==True)]
print('total number of sims', len(zscores_sims))

# Filtering for the correct homology-based label (ie. HLH, ILH, NLH), as per the input from the bash script. 
# Here, we keep the IDs of the real SVs which we will use to match to the simulated SVs
realspreds_homo = realspreds[realspreds.mechID_homo==pred]
IDs_homo = list(realspreds_homo.ID)
print('number of real svs that map to correct pred', len(realspreds_homo))

# Matching the real SVs' ID to the simulations and double checking the expected size with the if statement
zscores_sims_homo = zscores_sims[zscores_sims['compID'].isin(IDs_homo)]
print('number of sim svs that map to correct real svs within the correct pred', len(zscores_sims_homo))

if (len(IDs_homo)*100) == (len(zscores_sims_homo)) :
	print('no issues with sims and reals')

# Sorting the simulations so all the IDs are grouped together, essential for the part below
zscores_sims_homo = zscores_sims_homo.sort_values(by=['compID'], ascending=False)

# The following is how 1 bed file is generated for each iteration of simulations. 
# We loop through the 100 iterations, and use this to grab the correct row corresponding to a unique real SV ID (ie. for iteration 10, grab the 10th simulated SV across all real SVs)
# We then save this file in a bed file style and save it for later use. 
# Note the way we pull the start/stop locations for insertions/ deletions matches the way we pull them for the real SVs (as per generating_bedWork.py)
df_list = []
for i in range(100):
	df_chunk = zscores_sims_homo.iloc[i::100].reset_index(drop=True)
	df_list.append(df_chunk)

	df_chunk[['chr','POS', 'type', 'len']] = df_chunk['ID'].str.split('-',expand=True)
	df_chunk["len"] = pd.to_numeric(df_chunk["len"])
	df_chunk["POS"] = pd.to_numeric(df_chunk["POS"])

	for idx, row in df_chunk.iterrows():
		if typer == 'deletion':
			df_chunk.loc[idx, 'start']= int(row.POS)
			df_chunk.loc[idx, 'stop']= int(row.POS + row.len)
		else:
			df_chunk.loc[idx, 'start']= int(row.POS)
			df_chunk.loc[idx, 'stop']= int(row.POS)+1

	df_chunk['start'] = df_chunk['start'].astype('Int64')
	df_chunk['stop'] = df_chunk['stop'].astype('Int64')
	df_chunkbed = df_chunk[['chr', 'start', 'stop']]
	df_chunkbed.to_csv('/home/nboev/projects/def-sushant/nboev/analysis/'+project+'/'+loc+'/IDmechsvSIM/20240625/splitsChromo/bedWork100s/'+chromo+'.'+typer+'.HOMOpreds.'+pred+'SIMS'+str(i)+'.bed', sep='\t', index=False, header=False)
