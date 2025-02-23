#!/usr/bin/env python3
import numpy as np
import pandas as pd
import math
import seaborn as sns
import matplotlib.patches as mpatches
import glob
import itertools
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import hdbscan
import scipy.stats as stats
import statsmodels.api as sm
from statsmodels.formula.api import ols
from scipy import stats
from statsmodels.stats.multitest import multipletests
import hdbscan
from sklearn.metrics import silhouette_score, silhouette_samples
from scipy.spatial import distance_matrix
from scipy.stats import fisher_exact
import os
import datetime
from pytz import timezone
from sys import argv
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
from math import sqrt
from sklearn.metrics import f1_score
from sklearn.metrics import classification_report
from sklearn.neighbors import KNeighborsClassifier
from sklearn.preprocessing import StandardScaler
import joblib

print('\n **********************')
print('START TIME:', datetime.datetime.now(timezone('EST')))

# In this script, we take a dataset of standardized features and apply the clustering + KNN models.

# Data inputs sent in from the execution files. This includes the location of the models we require
df = pd.read_csv(str(argv[1]), index_col=None, sep='\t')
project = str(argv[2])
loc = str(argv[3])
typer = str(argv[4])
chromo = str(argv[5])
st_x = joblib.load(str(argv[6]))
pca = joblib.load(str(argv[7]))
knn = joblib.load(str(argv[8]))

# The name of the columns we ended up using for training. These vary in their lengths and types, as per the column filtering we did
if typer == 'insertion':
	cols = ['SV_MGW_mean_zscore', 'SV_Roll_sd_zscore', 'SV_EP_mean_zscore', 'std_pre_flankH3K27me3_2000_zscore', 'var_stab_zscore', 'avg_pre_flankCTCF_2000_zscore', 'avg_post_flankH3K36me3_2000_zscore', 'POST_ProT_mean_zscore', 'std_post_flankH3K4me1_2000_zscore', 'std_pre_flankH4K20me1_2000_zscore', 'PRE_MGW_sd_zscore', 'std_post_flankH3K36me3_2000_zscore', 'var_gc_zscore', 'postsv_pident_zscore', 'std_post_flankCTCF_2000_zscore', 'SV_HelT_sd_zscore', 'avg_pre_flankDNase-seq_2000_zscore', 'LINE_sum_zscore', 'avg_pre_flankH3K27me3_2000_zscore', 'post_flankcomp_2000_zscore', 'POST_Roll_mean_zscore', 'POST_ProT_sd_zscore', 'Arep_sum_zscore', 'pre_flankflex_2000_zscore', 'avg_post_flankH3K27ac_2000_zscore', 'PRE_HelT_mean_zscore', 'avg_post_flankH3K9ac_2000_zscore', 'POST_MGW_mean_zscore', 'Low_complexity_sum_zscore', 'PRE_Roll_sd_zscore', 'avg_post_flankWGB-Seq-_2000_zscore', 'POST_EP_sd_zscore', 'std_post_flankH3K79me2_2000_zscore', 'std_post_flankDNase-seq_2000_zscore', 'std_post_flankH3K9me3_2000_zscore', 'STR_sum_zscore', 'SV_HelT_mean_zscore', 'avg_post_flankCTCF_2000_zscore', 'avg_post_flankDNase-seq_2000_zscore', 'std_pre_flankH3K36me3_2000_zscore', 'PRE_HelT_sd_zscore', 'avg_post_flankH2AFZ_2000_zscore', 'std_post_flankH4K20me1_2000_zscore', 'avg_pre_flankH3K27ac_2000_zscore', 'var_flex_zscore', 'avg_post_flankH4K20me1_2000_zscore', 'Satellite_sum_zscore', 'avg_pre_flankH3K9me3_2000_zscore', 'avg_pre_flankH2AFZ_2000_zscore', 'avg_post_flankH3K9me3_2000_zscore', 'POST_MGW_sd_zscore', 'std_post_flankH2AFZ_2000_zscore', 'avg_post_flankH3K79me2_2000_zscore', 'SV_Roll_mean_zscore', 'avg_post_flankWGB-Seq+_2000_zscore', 'std_pre_flankH3K9me3_2000_zscore', 'PRE_ProT_mean_zscore', 'LTR_sum_zscore', 'PRE_MGW_mean_zscore', 'std_pre_flankH3K9ac_2000_zscore', 'DR_sum_zscore', 'MR_sum_zscore', 'pre_flankgc_2000_zscore', 'avg_pre_flankH3K36me3_2000_zscore', 'pre_flankcomp_2000_zscore', 'std_post_flankH3K9ac_2000_zscore', 'SINE_sum_zscore', 'POST_HelT_mean_zscore', 'POST_HelT_sd_zscore', 'std_post_flankH3K27ac_2000_zscore', 'post_flankgc_2000_zscore', 'ZDNA_sum_zscore', 'avg_pre_flankWGB-Seq-_2000_zscore', 'prepost_pident_zscore', 'avg_post_flankH3K4me1_2000_zscore', 'G4_sum_zscore', 'std_pre_flankDNase-seq_2000_zscore', 'PRE_EP_sd_zscore', 'avg_pre_flankH3K4me1_2000_zscore', 'SV_MGW_sd_zscore', 'POST_Roll_sd_zscore', 'prepost_Blastcoverage_zscore', 'var_comp_zscore', 'Simple_repeat_sum_zscore', 'std_post_flankH3K27me3_2000_zscore', 'avg_post_flankH3K27me3_2000_zscore', 'PRE_ProT_sd_zscore', 'SV_EP_sd_zscore', 'avg_pre_flankH3K79me2_2000_zscore', 'avg_pre_flankH4K20me1_2000_zscore', 'S50_zscore', 'std_pre_flankH3K27ac_2000_zscore', 'presv_pident_zscore', 'SV_ProT_mean_zscore', 'post_flankflex_2000_zscore', 'avg_pre_flankWGB-Seq+_2000_zscore', 'avg_pre_flankH3K9ac_2000_zscore', 'SV_ProT_sd_zscore', 'PRE_Roll_mean_zscore']
elif typer == 'deletion':
	cols = ['stdepi_H3K4me1_zscore', 'pre_flankflex_2000_zscore', 'avg_post_flankH3K4me1_2000_zscore', 'post_flankflex_2000_zscore', 'SV_Roll_sd_zscore', 'PRE_HelT_mean_zscore', 'LINE_sum_zscore', 'PRE_MGW_sd_zscore', 'avg_pre_flankH3K9ac_2000_zscore', 'SV_EP_sd_zscore', 'SV_HelT_sd_zscore', 'avgepi_H3K9ac_zscore', 'std_post_flankH3K36me3_2000_zscore', 'std_pre_flankH3K27ac_2000_zscore', 'std_post_flankDNase-seq_2000_zscore', 'std_post_flankH3K4me1_2000_zscore', 'avgepi_H3K9me3_zscore', 'avg_pre_flankDNase-seq_2000_zscore', 'stdepi_H4K20me1_zscore', 'var_stab_zscore', 'STR_sum_zscore', 'std_post_flankH3K79me2_2000_zscore', 'SV_Roll_mean_zscore', 'avg_post_flankH4K20me1_2000_zscore', 'avgepi_H3K27me3_zscore', 'std_post_flankH3K9me3_2000_zscore', 'std_post_flankH4K20me1_2000_zscore', 'avgepi_WGB-Seq+_zscore', 'PRE_MGW_mean_zscore', 'avg_pre_flankWGB-Seq+_2000_zscore', 'avg_post_flankCTCF_2000_zscore', 'POST_Roll_mean_zscore', 'std_post_flankWGB-Seq+_2000_zscore', 'avg_post_flankH3K27me3_2000_zscore', 'stdepi_WGB-Seq+_zscore', 'std_post_flankH3K9ac_2000_zscore', 'std_pre_flankH3K9me3_2000_zscore', 'Arep_sum_zscore', 'avg_pre_flankH4K20me1_2000_zscore', 'LTR_sum_zscore', 'std_post_flankH3K27ac_2000_zscore', 'prepost_pident_zscore', 'std_pre_flankH3K36me3_2000_zscore', 'POST_ProT_sd_zscore', 'SV_EP_mean_zscore', 'avg_pre_flankH2AFZ_2000_zscore', 'DR_sum_zscore', 'avg_post_flankH3K79me2_2000_zscore', 'Satellite_sum_zscore', 'avg_pre_flankH3K27ac_2000_zscore', 'S50_zscore', 'SV_MGW_mean_zscore', 'stdepi_WGB-Seq-_zscore', 'var_flex_zscore', 'SV_ProT_sd_zscore', 'var_gc_zscore', 'avgepi_H3K79me2_zscore', 'PRE_ProT_mean_zscore', 'presv_Blastcoverage_zscore', 'avgepi_H3K36me3_zscore', 'POST_ProT_mean_zscore', 'avg_pre_flankH3K9me3_2000_zscore', 'PRE_HelT_sd_zscore', 'post_flankcomp_2000_zscore', 'ZDNA_sum_zscore', 'avg_post_flankH3K27ac_2000_zscore', 'POST_MGW_sd_zscore', 'pre_flankgc_2000_zscore', 'stdepi_H3K9me3_zscore', 'postsv_Blastcoverage_zscore', 'std_post_flankH2AFZ_2000_zscore', 'postsv_pident_zscore', 'stdepi_H3K27me3_zscore', 'var_comp_zscore', 'avg_pre_flankH3K79me2_2000_zscore', 'post_flankgc_2000_zscore', 'avgepi_DNase-seq_zscore', 'Simple_repeat_sum_zscore', 'POST_MGW_mean_zscore', 'avgepi_WGB-Seq-_zscore', 'std_pre_flankH4K20me1_2000_zscore', 'SV_HelT_mean_zscore', 'stdepi_H3K27ac_zscore', 'PRE_ProT_sd_zscore', 'stdepi_H3K36me3_zscore', 'MR_sum_zscore', 'stdepi_H3K79me2_zscore', 'stdepi_DNase-seq_zscore', 'avg_post_flankWGB-Seq+_2000_zscore', 'avg_post_flankH3K36me3_2000_zscore', 'avgepi_H2AFZ_zscore', 'avg_pre_flankCTCF_2000_zscore', 'avg_pre_flankH3K36me3_2000_zscore', 'PRE_Roll_mean_zscore', 'prepost_Blastcoverage_zscore', 'std_post_flankCTCF_2000_zscore', 'Low_complexity_sum_zscore', 'SV_MGW_sd_zscore', 'avg_post_flankH3K9ac_2000_zscore', 'std_pre_flankH3K9ac_2000_zscore', 'G4_sum_zscore', 'SV_ProT_mean_zscore', 'avg_pre_flankH3K4me1_2000_zscore', 'stdepi_H2AFZ_zscore', 'PRE_EP_sd_zscore', 'std_post_flankH3K27me3_2000_zscore', 'avg_post_flankWGB-Seq-_2000_zscore', 'avg_post_flankDNase-seq_2000_zscore', 'avg_post_flankH3K9me3_2000_zscore', 'std_pre_flankH3K4me1_2000_zscore', 'POST_HelT_sd_zscore', 'avgepi_H4K20me1_zscore', 'presv_pident_zscore', 'stdepi_CTCF_zscore', 'stdepi_H3K9ac_zscore', 'avgepi_H3K27ac_zscore', 'avg_pre_flankWGB-Seq-_2000_zscore', 'avg_post_flankH2AFZ_2000_zscore', 'pre_flankcomp_2000_zscore', 'POST_HelT_mean_zscore', 'std_pre_flankDNase-seq_2000_zscore', 'avgepi_CTCF_zscore', 'POST_EP_sd_zscore', 'std_pre_flankH3K27me3_2000_zscore', 'SINE_sum_zscore', 'avgepi_H3K4me1_zscore']

df_filt = df[cols]

# Applying the scaling, PCA and KNN models to the matrix
df_filt = st_x.transform(df_filt)
df_filt = pca.transform(df_filt)
y_pred = knn.predict(df_filt)

# Making clustering-based predictions, and printing out the consistency between the homology-based label and clustering
df['fuzzypred'] = y_pred
print(df.groupby('fuzzypred').size())
df_m = df.groupby(['mechID_homo', 'fuzzypred']).size().unstack(level=0)
df_m_normalized = df_m.T.div(df_m.T.sum(axis=1), axis=0)
print(df_m)
print(df_m_normalized)

# Saving the results
df.to_csv('/home/nboev/projects/def-sushant/nboev/analysis/'+project+'/'+loc+'/IDmechsvSIM/20240625/splitsChromo/models/results/'+chromo+'.'+typer+'.HOMOpreds.tsv', sep='\t', index=False)

