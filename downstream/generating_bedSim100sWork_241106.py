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

# grabbing all the different bed files we need

realspreds = pd.read_csv(str(argv[1]), comment='#', sep='\t')
zscores = pd.read_csv(str(argv[2]), comment='#', sep='\t')

project = str(argv[3])
loc = str(argv[4])
chromo = str(argv[5])
typer = str(argv[6])
pred = str(argv[7])

print(project, loc, chromo, typer)

zscores_sims = zscores[(zscores.n_zscore == 101)&(zscores['Sim']==True)]
print('total number of sims', len(zscores_sims))

# start with fuzzy0
realspreds_fuzzy0 = realspreds[realspreds.mechID_homo==pred]
IDs_fuzzy0 = list(realspreds_fuzzy0.ID)
print('number of real svs that map to correct pred', len(realspreds_fuzzy0))
print(len(IDs_fuzzy0))

# finding the corresponding simulations
zscores_sims_fuzzy0 = zscores_sims[zscores_sims['compID'].isin(IDs_fuzzy0)]
print('number of sim svs that map to correct real svs within the correct pred', len(zscores_sims_fuzzy0))

if (len(IDs_fuzzy0)*100) == (len(zscores_sims_fuzzy0)) :
	print('no issues with sims and reals')


# I need to sort these so the IDs are together
zscores_sims_fuzzy0 = zscores_sims_fuzzy0.sort_values(by=['compID'], ascending=False)

df_list = []
for i in range(100):
	df_chunk = zscores_sims_fuzzy0.iloc[i::100].reset_index(drop=True)
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
	print(df_chunkbed.head())

	df_chunkbed.to_csv('/home/nboev/projects/def-sushant/nboev/analysis/'+project+'/'+loc+'/IDmechsvSIM/20240625/splitsChromo/bedWork100s/'+chromo+'.'+typer+'.HOMOpreds.'+pred+'SIMS'+str(i)+'.bed', sep='\t', index=False, header=False)

print(len(df_list))
print(df_list[0].head())
print(df_list[0][['compID', 'ID']])
print(len(df_list[0]['compID'].unique()))
print(len(df_list[0]['ID'].unique()))
