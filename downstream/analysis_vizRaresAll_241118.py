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

def genBed(df):
	df[['chr','POS', 'type', 'len']] = df['ID'].str.split('-',expand=True)
	df["len"] = pd.to_numeric(df["len"])
	df["POS"] = pd.to_numeric(df["POS"])

	for idx, row in df.iterrows():
		df.loc[idx, 'start']= int(row.POS)
		df.loc[idx, 'stop']= int(row.POS + row.len)

	df['start'] = df['start'].astype('Int64')
	df['stop'] = df['stop'].astype('Int64')

	df0 = df[df.fuzzypred=='fuzzy0']
	df0 = df0[['chr', 'start', 'stop']]

	df1 = df[df.fuzzypred=='fuzzy1']
	df1 = df1[['chr', 'start', 'stop']]

	return[df0,df1]


path = '/home/nboev/projects/def-sushant/nboev/analysis/20220422_3202_phased_SNV_INDEL_SV_bychrom/SVTrue_typedeletion_resTrue/IDmechsvSIM/20240625/splitsChromo/models/results/'
all_files = glob.glob(os.path.join(path , "chr*.tsv"))
li1 = []

for filename in all_files:
	df = pd.read_csv(filename, index_col=None, sep='\t')
	chromo = filename.split('results/')[1].split('.')[0]
	svtype = filename.split('results/')[1].split('.')[1].split('_')[0]
	df['chromo'] = chromo
	df['svtype'] = svtype
	df['dataset'] = '1KG'

	for idx, row in df.iterrows():
		if row.fuzzypred == 'fuzzy0':
			df.loc[idx, 'fuzzypred'] = 'fuzzy1'
		else:
			df.loc[idx, 'fuzzypred'] = 'fuzzy0'

	li1.append(df)
big = pd.concat(li1, axis=0, ignore_index=True)
print(big.head())

path = '/home/nboev/projects/def-sushant/nboev/preprocess/20220422_3202_phased_SNV_INDEL_SV_bychrom/SVTrue_typedeletion_resTrue/SVlen/SVTrue_typedeletion_resTrue.csv'
svs_1000 = pd.read_csv(path, index_col=None, sep='\t')
print(svs_1000)

for idx, row in svs_1000.iterrows():
	svs_1000.loc[idx, 'ID'] = row.CHROM +'-'+ str(row.POS+1) +'-'+ 'DEL' +'-'+ str(int(float(row.SVlen)))
merged_1000g = pd.merge(big, svs_1000, how='inner', on=['ID'])
merged_1000g['dataset'] = '1000G'

merged_1000g['AC']= merged_1000g['INFO'].str.split('AC=',expand=True)[1].str.split(';',expand=True)[0]
merged_1000g['AC'] = merged_1000g['AC'].astype(float)
merged_1000g['AF']= merged_1000g['INFO'].str.split('AF=',expand=True)[1].str.split(';',expand=True)[0]
merged_1000g['AF'] = merged_1000g['AF'].astype(float)

likelyrare = merged_1000g[merged_1000g.AC<64]
likelycommon = merged_1000g[merged_1000g.AC>320]

ind_1000g = merged_1000g.loc[:,merged_1000g.columns.str.contains("HG0|NA")].columns.to_list()
ind_1000g.pop(0)
ind_1000g.pop(0)

anc = pd.read_csv('/home/nboev/projects/def-sushant/nboev/data/phased_SNV_INDEL_SV_20220422_3202/igsr_samples.tsv', sep = '\t')
anc.rename({'Sample name': 'ind'}, axis=1, inplace=True)
anc = anc[anc['ind'].isin(ind_1000g)]
anc.set_index('ind', inplace=True)

mydict = anc.to_dict(orient="ind")

for key, value in mydict.items():
	toadd_rare = likelyrare[((likelyrare[key] == '1|0')|(likelyrare[key] == '1|1')|(likelyrare[key] == '0|1'))]
	toadd_rare['ind'] = key
	mydict[key].update({'on_rare':toadd_rare})

	toadd_common = likelycommon[((likelycommon[key] == '1|0')|(likelycommon[key] == '1|1')|(likelycommon[key] == '0|1'))]
	toadd_common['ind'] = key
	mydict[key].update({'on_common':toadd_common})

grab_rares = []
grab_commons = []
for key, value in mydict.items():
	grab_rares.append(mydict[key]['on_rare'])
	grab_commons.append(mydict[key]['on_common'])

# we need to save these so we can figure out where these SVs are hitting..
rares = pd.concat(grab_rares, axis=0, ignore_index=True)
rares['role'] = 'rares'
rares_1000g = pd.merge(rares, anc, how='left', on=['ind'])

for i in ind_1000g:
	df0,df1 = genBed(mydict[i]['on_rare'])
	df0.to_csv('/home/nboev/projects/def-sushant/nboev/analysis/20220422_3202_phased_SNV_INDEL_SV_bychrom/SVTrue_typedeletion_resTrue/IDmechsvSIM/20240625/splitsChromo/bedRareInd/'+mydict[i]['Population code']+'.'+i+'.HOMOpreds.fuzzy0.bed', sep='\t', index=False, header=False)
	df1.to_csv('/home/nboev/projects/def-sushant/nboev/analysis/20220422_3202_phased_SNV_INDEL_SV_bychrom/SVTrue_typedeletion_resTrue/IDmechsvSIM/20240625/splitsChromo/bedRareInd/'+mydict[i]['Population code']+'.'+i+'.HOMOpreds.fuzzy1.bed', sep='\t', index=False, header=False)


'''
count = rares_1000g.groupby([ 'fuzzypred', 'ind', 'Superpopulation code', 'Population code']).size().unstack(level=0).fillna(0)
#count.to_csv("./rares_1000g_fuzzypredrawALL.csv")
count['sum'] = count['fuzzy0']+count['fuzzy1']
print(count.head())
count_flat = count.reset_index()
count_flat.columns = ['_'.join(filter(None, col)).strip() for col in count_flat.columns.values]

nonAFR_high = list(count_flat[(count_flat['S_u_p_e_r_p_o_p_u_l_a_t_i_o_n_ _c_o_d_e']!='AFR')&(count_flat['s_u_m']>100)]['i_n_d'].unique())
nonAFR_high_mydict = {key: mydict[key]['on_rare'] for key in nonAFR_high if key in mydict}

for i in nonAFR_high:
	df0,df1 = genBed(nonAFR_high_mydict[i])
	df0.to_csv('/home/nboev/projects/def-sushant/nboev/analysis/20220422_3202_phased_SNV_INDEL_SV_bychrom/SVTrue_typedeletion_resTrue/IDmechsvSIM/20240625/splitsChromo/bedRareInd/nonAFR_high/nonAFR_high.'+i+'.HOMOpreds.fuzzy0.bed', sep='\t', index=False, header=False)
	df1.to_csv('/home/nboev/projects/def-sushant/nboev/analysis/20220422_3202_phased_SNV_INDEL_SV_bychrom/SVTrue_typedeletion_resTrue/IDmechsvSIM/20240625/splitsChromo/bedRareInd/nonAFR_high/nonAFR_high.'+i+'.HOMOpreds.fuzzy1.bed', sep='\t', index=False, header=False)


nonAFR_low = list(count_flat[(count_flat['S_u_p_e_r_p_o_p_u_l_a_t_i_o_n_ _c_o_d_e']!='AFR')&(count_flat['s_u_m']<=100)]['i_n_d'].unique())
nonAFR_low_mydict = {key: mydict[key]['on_rare'] for key in nonAFR_low if key in mydict}

for i in nonAFR_low:
	df0,df1 = genBed(nonAFR_low_mydict[i])
	df0.to_csv('/home/nboev/projects/def-sushant/nboev/analysis/20220422_3202_phased_SNV_INDEL_SV_bychrom/SVTrue_typedeletion_resTrue/IDmechsvSIM/20240625/splitsChromo/bedRareInd/nonAFR_low/nonAFR_low.'+i+'.HOMOpreds.fuzzy0.bed', sep='\t', index=False, header=False)
	df1.to_csv('/home/nboev/projects/def-sushant/nboev/analysis/20220422_3202_phased_SNV_INDEL_SV_bychrom/SVTrue_typedeletion_resTrue/IDmechsvSIM/20240625/splitsChromo/bedRareInd/nonAFR_low/nonAFR_low.'+i+'.HOMOpreds.fuzzy1.bed', sep='\t', index=False, header=False)

AFR_high = list(count_flat[(count_flat['S_u_p_e_r_p_o_p_u_l_a_t_i_o_n_ _c_o_d_e']=='AFR')&(count_flat['s_u_m']>250)]['i_n_d'].unique())
AFR_high_mydict = {key: mydict[key]['on_rare'] for key in AFR_high if key in mydict}

for i in AFR_high:
	df0,df1 = genBed(AFR_high_mydict[i])
	df0.to_csv('/home/nboev/projects/def-sushant/nboev/analysis/20220422_3202_phased_SNV_INDEL_SV_bychrom/SVTrue_typedeletion_resTrue/IDmechsvSIM/20240625/splitsChromo/bedRareInd/AFR_high/AFR_high.'+i+'.HOMOpreds.fuzzy0.bed', sep='\t', index=False, header=False)
	df1.to_csv('/home/nboev/projects/def-sushant/nboev/analysis/20220422_3202_phased_SNV_INDEL_SV_bychrom/SVTrue_typedeletion_resTrue/IDmechsvSIM/20240625/splitsChromo/bedRareInd/AFR_high/AFR_high.'+i+'.HOMOpreds.fuzzy1.bed', sep='\t', index=False, header=False)

AFR_low = list(count_flat[(count_flat['S_u_p_e_r_p_o_p_u_l_a_t_i_o_n_ _c_o_d_e']=='AFR')&(count_flat['s_u_m']<=250)]['i_n_d'].unique())
AFR_low_mydict = {key: mydict[key]['on_rare'] for key in AFR_low if key in mydict}

for i in AFR_low:
	df0,df1 = genBed(AFR_low_mydict[i])
	df0.to_csv('/home/nboev/projects/def-sushant/nboev/analysis/20220422_3202_phased_SNV_INDEL_SV_bychrom/SVTrue_typedeletion_resTrue/IDmechsvSIM/20240625/splitsChromo/bedRareInd/AFR_low/AFR_low.'+i+'.HOMOpreds.fuzzy0.bed', sep='\t', index=False, header=False)
	df1.to_csv('/home/nboev/projects/def-sushant/nboev/analysis/20220422_3202_phased_SNV_INDEL_SV_bychrom/SVTrue_typedeletion_resTrue/IDmechsvSIM/20240625/splitsChromo/bedRareInd/AFR_low/AFR_low.'+i+'.HOMOpreds.fuzzy1.bed', sep='\t', index=False, header=False)


count = count.apply(lambda x: x/x.sum(), axis=1).reset_index()
count.to_csv("./rares_1000g_fuzzypredALL.csv")


commons = pd.concat(grab_commons, axis=0, ignore_index=True)
commons['role'] = 'commons'
commons_1000g = pd.merge(commons, anc, how='left', on=['ind'])
count = commons_1000g.groupby([ 'fuzzypred', 'ind', 'Superpopulation code', 'Population code']).size().unstack(level=0).fillna(0)
count.to_csv("./commons_1000g_fuzzypredrawALL.csv")
count = count.apply(lambda x: x/x.sum(), axis=1).reset_index()
count.to_csv("./commons_1000g_fuzzypredALL.csv")
'''






