#!/usr/bin/env python3

import pandas as pd
import glob
from sys import argv
import datetime
from pytz import timezone
import numpy as np
from scipy.stats import zscore
from io import StringIO

print('START TIME:', datetime.datetime.now(timezone('EST')))

# In this script, we take the merged and annotated csv from the real and simulated SVs and calculate z-scores across all numeric columns

#  Data inputs sent in from the execution files
# In the case of the real SVs, we are adding the compID and Sim column which were not originally present 
# The compID column represents the "comparison ID" which will be used to link real-simulated SVs. 
	# Recall, in the simulated vcfs, each SV has its own ID, derived from the simulation's length and location. 
reals = pd.read_csv(str(argv[1]), comment='#', sep='\t', low_memory=False)
project = str(argv[2])
loc = str(argv[3])
filename = str(argv[4])
chr = str(argv[5])
reals['compID'] = reals['ID']
reals['Sim'] = False
print(len(reals))

# In the case of the simulated SVs, we check the length of the simulations are 100x the length of the reals
sims = pd.read_csv(str(argv[6]), comment='#', sep='\t', low_memory=False)
print(len(sims))
sims['compID'] = sims['Sim_ID']

# Concatenate the two dataframes
df_merged = pd.concat([reals, sims], ignore_index=True)

# Counting up the number of simulations which were successfully generated for a given real SV. Recall, this depends on SURVIVOR's success
df_merged['n_zscore'] = df_merged.groupby('compID')['compID'].transform('size')
df_merged = df_merged.replace(np.nan, 0)

'''
The function, z_score, calculates z-scores across a column
The input: 
	x: A list of values to use for calculating a z-score. 
In this case, we use a the nan-omit policy when applying scipy.stats's zscore function
'''

def z_score(x):
	z = zscore(x, nan_policy='omit')
	return z

# The list of columns we can appropriately calculate a z-score (ie. not binary/ logic columns)
cols = ['var_gc', 'var_comp', 'var_flex', 'var_stab', 'pre_flankgc_2000','post_flankgc_2000', 'pre_flankcomp_2000', 'post_flankcomp_2000','pre_flankflex_2000', 'post_flankflex_2000', 'pre_flankstab_2000','post_flankstab_2000', 'LINE_sum','Low_complexity_sum', 'LTR_sum', 'Satellite_sum', 'Simple_repeat_sum','SINE_sum', 'Arep_sum', 'DR_sum', 'G4_sum', 'MR_sum', 'STR_sum','ZDNA_sum', 'avgepi_CTCF', 'stdepi_CTCF','avgepi_H2AFZ', 'stdepi_H2AFZ', 'avgepi_H3K36me3','stdepi_H3K36me3', 'avgepi_H3K9ac','stdepi_H3K9ac', 'avgepi_H3K9me3','stdepi_H3K9me3', 'avgepi_H4K20me1','stdepi_H4K20me1', 'avgepi_H3K4me1','stdepi_H3K4me1', 'avgepi_H3K79me2','stdepi_H3K79me2', 'avgepi_H3K27ac','stdepi_H3K27ac', 'avgepi_H3K27me3','stdepi_H3K27me3', 'avgepi_DNase-seq','stdepi_DNase-seq', 'avgepi_WGB-Seq+','stdepi_WGB-Seq+', 'avgepi_WGB-Seq-','stdepi_WGB-Seq-', 'avg_pre_flankCTCF_2000','avg_post_flankCTCF_2000', 'std_pre_flankCTCF_2000','std_post_flankCTCF_2000', 'avg_pre_flankH2AFZ_2000','avg_post_flankH2AFZ_2000','std_pre_flankH2AFZ_2000', 'std_post_flankH2AFZ_2000','avg_pre_flankH3K36me3_2000', 'avg_post_flankH3K36me3_2000','std_pre_flankH3K36me3_2000', 'std_post_flankH3K36me3_2000','avg_pre_flankH3K9ac_2000', 'avg_post_flankH3K9ac_2000','std_pre_flankH3K9ac_2000', 'std_post_flankH3K9ac_2000','avg_pre_flankH3K9me3_2000', 'avg_post_flankH3K9me3_2000','std_pre_flankH3K9me3_2000', 'std_post_flankH3K9me3_2000','avg_pre_flankH4K20me1_2000', 'avg_post_flankH4K20me1_2000','std_pre_flankH4K20me1_2000', 'std_post_flankH4K20me1_2000','avg_pre_flankH3K4me1_2000', 'avg_post_flankH3K4me1_2000','std_pre_flankH3K4me1_2000', 'std_post_flankH3K4me1_2000','avg_pre_flankH3K79me2_2000', 'avg_post_flankH3K79me2_2000','std_pre_flankH3K79me2_2000', 'std_post_flankH3K79me2_2000', 'avg_pre_flankH3K27ac_2000', 'avg_post_flankH3K27ac_2000','std_pre_flankH3K27ac_2000', 'std_post_flankH3K27ac_2000','avg_pre_flankH3K27me3_2000', 'avg_post_flankH3K27me3_2000','std_pre_flankH3K27me3_2000', 'std_post_flankH3K27me3_2000','avg_pre_flankDNase-seq_2000', 'avg_post_flankDNase-seq_2000','std_pre_flankDNase-seq_2000', 'std_post_flankDNase-seq_2000','avg_pre_flankWGB-Seq+_2000', 'avg_post_flankWGB-Seq+_2000','std_pre_flankWGB-Seq+_2000', 'std_post_flankWGB-Seq+_2000','avg_pre_flankWGB-Seq-_2000', 'avg_post_flankWGB-Seq-_2000','std_pre_flankWGB-Seq-_2000', 'std_post_flankWGB-Seq-_2000','PRE_EP_sd', 'POST_MGW_mean', 'POST_MGW_sd', 'POST_HelT_mean','POST_HelT_sd', 'POST_ProT_mean', 'POST_ProT_sd', 'POST_Roll_mean','POST_Roll_sd', 'POST_EP_mean', 'POST_EP_sd', 'S50','SV_MGW_mean','SV_MGW_sd','SV_HelT_mean','SV_HelT_sd', 'SV_ProT_mean', 'SV_ProT_sd','SV_Roll_mean','SV_Roll_sd', 'SV_EP_mean','SV_EP_sd','PRE_MGW_mean', 'PRE_MGW_sd','PRE_HelT_mean','PRE_HelT_sd', 'PRE_ProT_mean',	'PRE_ProT_sd','PRE_Roll_mean', 'PRE_Roll_sd','PRE_EP_mean', 'prepost_pident', 'prepost_Blastcoverage', 'presv_pident', 'presv_Blastcoverage', 'postsv_pident', 'postsv_Blastcoverage']

# We loop through this list of columns, then check for their presence (recall insertions do no have all these features)
# If a column is present, call the z_score function above, for each unique compID
for i in cols:
	if i in list(df_merged.columns):
		df_merged[i] = pd.to_numeric(df_merged[i], errors='coerce')
		df_merged = df_merged[df_merged[i].notna()]
		df_merged[i+'_zscore'] = df_merged.groupby('compID')[i].transform(z_score)

df_merged = df_merged.replace(np.nan, 0)

# Saving the z-score annotated dataframe
df_merged.to_csv('/home/nboev/projects/def-sushant/nboev/preprocess/'+project+'/'+loc+'/zscore/'+chr+'.'+filename+ 'zscore.csv', sep='\t', index=False)
