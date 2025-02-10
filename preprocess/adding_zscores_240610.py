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

#loading in the chromosome band data, we don't need this yet... our goal is just to see if we can load this in...
# loading in the csv (which holds all the chromo)
reals = pd.read_csv(str(argv[1]), comment='#', sep='\t', low_memory=False)
project = str(argv[2])
loc = str(argv[3])
filename = str(argv[4])
chr = str(argv[5])

reals['compID'] = reals['ID']
reals['Sim'] = False
print(reals.head())

# recall, I think theres a little issue with var_gc (changed biopython version), will need to *100 on sims
# I don't think we'll need this moving forward
#reals['var_gc'] = reals['var_gc'].apply(lambda x: x/100)
#reals['pre_flankgc_2000'] = reals['pre_flankgc_2000'].apply(lambda x: x/100)
#reals['post_flankgc_2000'] = reals['post_flankgc_2000'].apply(lambda x: x/100)

#with open(str(argv[6]), 'r') as f:
#	header = f.readline().strip()
#	filtered_lines = [line.strip() for line in f if 'Sim_ID' not in line]

#filtered_file = StringIO('\n'.join([header] + filtered_lines))
#sims = pd.read_csv(filtered_file, sep='\t')

sims = pd.read_csv(str(argv[6]), comment='#', sep='\t', low_memory=False)
print(len(sims))
sims['compID'] = sims['Sim_ID']
print(sims.head())

df_merged = pd.concat([reals, sims], ignore_index=True)
print(len(df_merged))
print(df_merged.head())
print(df_merged.columns)

def z_score(x):
	z = zscore(x, nan_policy='omit')
	return z

# columns we'll calculate the z score for..
df_merged['n_zscore'] = df_merged.groupby('compID')['compID'].transform('size')
print(df_merged['n_zscore'])

#cols = ['var_gc', 'var_comp', 'var_flex', 'var_stab', 'pre_flankgc_2000', 'post_flankgc_2000', ]
cols = ['var_gc', 'var_comp', 'var_flex', 'var_stab', 'pre_flankgc_2000','post_flankgc_2000', 'pre_flankcomp_2000', 'post_flankcomp_2000','pre_flankflex_2000', 'post_flankflex_2000', 'pre_flankstab_2000','post_flankstab_2000', 'LINE_sum','Low_complexity_sum', 'LTR_sum', 'Satellite_sum', 'Simple_repeat_sum','SINE_sum', 'Arep_sum', 'DR_sum', 'G4_sum', 'MR_sum', 'STR_sum','ZDNA_sum', 'avgepi_CTCF', 'stdepi_CTCF','avgepi_H2AFZ', 'stdepi_H2AFZ', 'avgepi_H3K36me3','stdepi_H3K36me3', 'avgepi_H3K9ac','stdepi_H3K9ac', 'avgepi_H3K9me3','stdepi_H3K9me3', 'avgepi_H4K20me1','stdepi_H4K20me1', 'avgepi_H3K4me1','stdepi_H3K4me1', 'avgepi_H3K79me2','stdepi_H3K79me2', 'avgepi_H3K27ac','stdepi_H3K27ac', 'avgepi_H3K27me3','stdepi_H3K27me3', 'avgepi_DNase-seq','stdepi_DNase-seq', 'avgepi_WGB-Seq+','stdepi_WGB-Seq+', 'avgepi_WGB-Seq-','stdepi_WGB-Seq-', 'avg_pre_flankCTCF_2000','avg_post_flankCTCF_2000', 'std_pre_flankCTCF_2000','std_post_flankCTCF_2000', 'avg_pre_flankH2AFZ_2000','avg_post_flankH2AFZ_2000','std_pre_flankH2AFZ_2000', 'std_post_flankH2AFZ_2000','avg_pre_flankH3K36me3_2000', 'avg_post_flankH3K36me3_2000','std_pre_flankH3K36me3_2000', 'std_post_flankH3K36me3_2000','avg_pre_flankH3K9ac_2000', 'avg_post_flankH3K9ac_2000','std_pre_flankH3K9ac_2000', 'std_post_flankH3K9ac_2000','avg_pre_flankH3K9me3_2000', 'avg_post_flankH3K9me3_2000','std_pre_flankH3K9me3_2000', 'std_post_flankH3K9me3_2000','avg_pre_flankH4K20me1_2000', 'avg_post_flankH4K20me1_2000','std_pre_flankH4K20me1_2000', 'std_post_flankH4K20me1_2000','avg_pre_flankH3K4me1_2000', 'avg_post_flankH3K4me1_2000','std_pre_flankH3K4me1_2000', 'std_post_flankH3K4me1_2000','avg_pre_flankH3K79me2_2000', 'avg_post_flankH3K79me2_2000','std_pre_flankH3K79me2_2000', 'std_post_flankH3K79me2_2000', 'avg_pre_flankH3K27ac_2000', 'avg_post_flankH3K27ac_2000','std_pre_flankH3K27ac_2000', 'std_post_flankH3K27ac_2000','avg_pre_flankH3K27me3_2000', 'avg_post_flankH3K27me3_2000','std_pre_flankH3K27me3_2000', 'std_post_flankH3K27me3_2000','avg_pre_flankDNase-seq_2000', 'avg_post_flankDNase-seq_2000','std_pre_flankDNase-seq_2000', 'std_post_flankDNase-seq_2000','avg_pre_flankWGB-Seq+_2000', 'avg_post_flankWGB-Seq+_2000','std_pre_flankWGB-Seq+_2000', 'std_post_flankWGB-Seq+_2000','avg_pre_flankWGB-Seq-_2000', 'avg_post_flankWGB-Seq-_2000','std_pre_flankWGB-Seq-_2000', 'std_post_flankWGB-Seq-_2000','PRE_EP_sd', 'POST_MGW_mean', 'POST_MGW_sd', 'POST_HelT_mean','POST_HelT_sd', 'POST_ProT_mean', 'POST_ProT_sd', 'POST_Roll_mean','POST_Roll_sd', 'POST_EP_mean', 'POST_EP_sd', 'S50','SV_MGW_mean','SV_MGW_sd','SV_HelT_mean','SV_HelT_sd', 'SV_ProT_mean', 'SV_ProT_sd','SV_Roll_mean','SV_Roll_sd', 'SV_EP_mean','SV_EP_sd','PRE_MGW_mean', 'PRE_MGW_sd','PRE_HelT_mean','PRE_HelT_sd', 'PRE_ProT_mean',	'PRE_ProT_sd','PRE_Roll_mean', 'PRE_Roll_sd','PRE_EP_mean', 'prepost_pident', 'prepost_Blastcoverage', 'presv_pident', 'presv_Blastcoverage', 'postsv_pident', 'postsv_Blastcoverage']

# try this here...
df_merged = df_merged.replace(np.nan, 0)

for i in cols:
	if i in list(df_merged.columns):
#		print(df_merged[[i, 'compID']])

# note... we had an issue where there is NO variation (ie. the denom is 0), we can't calcualte the zscore, so it returns a nan value... when this happens, let's turn these values into 0s
		df_merged[i] = pd.to_numeric(df_merged[i], errors='coerce')

		df_merged = df_merged[df_merged[i].notna()]
#		print(df_merged.groupby('compID').size())

		df_merged[i+'_zscore'] = df_merged.groupby('compID')[i].transform(z_score)
#		print(df_merged[[i, i+'_zscore', 'compID', 'Sim']])

df_merged = df_merged.replace(np.nan, 0)
print(len(df_merged))
df_merged.to_csv('/home/nboev/projects/def-sushant/nboev/preprocess/'+project+'/'+loc+'/zscore/'+chr+'.'+filename+ 'zscore.csv', sep='\t', index=False)
