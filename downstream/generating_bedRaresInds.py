#!/usr/bin/env python3

import glob
import numpy as np
import pandas as pd
from sys import argv
import datetime
from pytz import timezone
import os

print('\n **********************')
print('START TIME:', datetime.datetime.now(timezone('EST')))

# In this script, we take the cluster-based labels, along with genotyping information to create bed-files to represent rare SVs for each individual
# Importantly, we also append ancestry-specific information, allowing us to consider population later

'''
The function, genBed, takes one dataframe and creates a bed-file like file. This will be used at the end for each individual's SVs
The input: 
	df: The df that holds the individuals SVs

Procedure: 
1. We split up the ID column to retrieve the length and position of the deletions
2. We identify the start/stop positions by considering the SV length, making sure the columns are numerics
3. We split these SVs based on if they map to a clusterI or cluster II SV. These two are indepdently returned by the function
'''

def genBed(df):
	df[['chr','POS', 'type', 'len']] = df['ID'].str.split('-',expand=True)
	df["len"] = pd.to_numeric(df["len"])
	df["POS"] = pd.to_numeric(df["POS"])

	for idx, row in df.iterrows():
		df.loc[idx, 'start']= int(row.POS)
		df.loc[idx, 'stop']= int(row.POS + row.len)

	df0 = df[df.ClusterBasedLabel=='clusterI']
	df0 = df0[['chr', 'start', 'stop']]
	df1 = df[df.ClusterBasedLabel=='clusterII']
	df1 = df1[['chr', 'start', 'stop']]

	return[df0,df1]

# The following grabs the homology and clustering-based results after annotation for the 1KG dataset. 
# We grab the tsv files, process them to append chromosome and sv type columns, then append them all together into one large dataframe, 'big'
path = '/home/nboev/projects/def-sushant/nboev/analysis/20220422_3202_phased_SNV_INDEL_SV_bychrom/SVTrue_typedeletion_resTrue/IDmechsvSIM/20240625/splitsChromo/models/results/'
all_files = glob.glob(os.path.join(path , "chr*.tsv"))
li1 = []

for filename in all_files:
	df = pd.read_csv(filename, index_col=None, sep='\t')
	chromo = filename.split('results/')[1].split('.')[0]
	svtype = filename.split('results/')[1].split('.')[1].split('_')[0]
	df['chromo'] = chromo
	df['svtype'] = svtype	
	li1.append(df)
	
big = pd.concat(li1, axis=0, ignore_index=True)

# Next, we grab the original vcf-like file, which holds the SVs, along with genotyping information and the "INFO" column. 
# Since the original file, may or may not, contain an "ID" column, we add the column and merge it to the annotated dataframe from above
path = str(argv[1])
svs_1000 = pd.read_csv(path, index_col=None, sep='\t')

for idx, row in svs_1000.iterrows():
	svs_1000.loc[idx, 'ID'] = row.CHROM +'-'+ str(row.POS+1) +'-'+ 'DEL' +'-'+ str(int(float(row.SVlen)))
merged_1000g = pd.merge(big, svs_1000, how='inner', on=['ID'])

# Since we are interested in rare variants (ie. those with AFs<1%), we use the AC column to filter for rare SVs
merged_1000g['AC']= merged_1000g['INFO'].str.split('AC=',expand=True)[1].str.split(';',expand=True)[0]
merged_1000g['AC'] = merged_1000g['AC'].astype(float)
likelyrare = merged_1000g[merged_1000g.AC<64]

# We need to pull the complete list of individuals in the 1KG dataset
ind_1000g = merged_1000g.loc[:,merged_1000g.columns.str.contains("HG0|NA")].columns.to_list()
ind_1000g.pop(0)
ind_1000g.pop(0)

# We will want to append ancestral information for each individual's SV list. Therefore, we input the metadata sheet
anc = pd.read_csv(str(argv[2]), sep = '\t')
anc.rename({'Sample name': 'ind'}, axis=1, inplace=True)
anc = anc[anc['ind'].isin(ind_1000g)]

# To create individual-level bed files, we move the individual column to the index, to create a dictionary where the individual is the key.
# Next, we loop through these individuals, and use the filtered rare SV dataframe to identify SVs which are present in the individual (ie. looking for the following genotypes: 1|0, 1|1, 0|1)
# Finally, we convert this dictionary into a dataframe for later use
anc.set_index('ind', inplace=True)
mydict = anc.to_dict(orient="ind")

for key, value in mydict.items():
	toadd_rare = likelyrare[((likelyrare[key] == '1|0')|(likelyrare[key] == '1|1')|(likelyrare[key] == '0|1'))]
	toadd_rare['ind'] = key
	mydict[key].update({'on_rare':toadd_rare})

grab_rares = []
for key, value in mydict.items():
	grab_rares.append(mydict[key]['on_rare'])

# This is a subset of rare SVs, which we will now merge with the metadata sheet so we have access to ancestry-specific information
rares = pd.concat(grab_rares, axis=0, ignore_index=True)
rares['role'] = 'rares'
rares_1000g = pd.merge(rares, anc, how='left', on=['ind'])

# For each individual, we send this individual's filtered SVs to the genBed function to generate cluster I and II bed files. 
# These bed files are saved for later use
for i in ind_1000g:
	df0,df1 = genBed(mydict[i]['on_rare'])
	df0.to_csv('/home/nboev/projects/def-sushant/nboev/analysis/20220422_3202_phased_SNV_INDEL_SV_bychrom/SVTrue_typedeletion_resTrue/IDmechsvSIM/20240625/splitsChromo/bedRareInd/'+mydict[i]['Population code']+'.'+i+'.HOMOpreds.clusterI.bed', sep='\t', index=False, header=False)
	df1.to_csv('/home/nboev/projects/def-sushant/nboev/analysis/20220422_3202_phased_SNV_INDEL_SV_bychrom/SVTrue_typedeletion_resTrue/IDmechsvSIM/20240625/splitsChromo/bedRareInd/'+mydict[i]['Population code']+'.'+i+'.HOMOpreds.clusterII.bed', sep='\t', index=False, header=False)
