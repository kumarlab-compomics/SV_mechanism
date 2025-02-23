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

# In this script, we take train and save the optimal clustering and KNN models. Very similar in nature to analysis_hdbscanID.py

# Data inputs sent in from the execution files. This includes the # of components to use, the minimum cluster size parameter and the distance metric
path = str(argv[1])
searching = str(argv[2])
searching2 = str(argv[3])
pcs_param = int(argv[4])
min_cluster_size_param = int(argv[5])
metric_param = str(argv[6])
project = str(argv[7])
loc = str(argv[8])
typer = str(argv[9])

print('PARAMS;')
print('PCs', pcs_param)
print('min_cluster_size', min_cluster_size_param)
print('metric', metric_param)

'''
The function, top_two_probs_diff, takes in the clusters' probability of assignment and calculates their difference
The input: 
	probs: For a given sample, we calculate the greatest and 2nd greatest cluster probability.
'''

def top_two_probs_diff(probs):
	sorted_probs = np.sort(probs)
	return sorted_probs[-1] - sorted_probs[-2]

all_files = glob.glob(os.path.join(path , '*'+searching))
li = []

# Looping through all the annotated real csvs and appending them into one large dataframe
for filename in all_files:
	df = pd.read_csv(filename, index_col=None, sep='\t')
	chromo = filename.split(searching2)[1].split('.')[0]
	svtype = filename.split(searching2)[1].split('.')[1].split('_')[0]
	df['chromo'] = chromo
	df['svtype'] = svtype

	li.append(df)
big = pd.concat(li, axis=0, ignore_index=True)

# Here we are holding out test chromosomes
testchromo = ['chr3', 'chr5', 'chr14', 'chr17', 'chr20', 'chr22']
big_train = big[~big.chromo.isin(testchromo)]
big_train = big_train.reset_index(drop=True)


'''
The following two steps were done for dimension reduction (DR), specifically to eliminate redundant features
DR1. We applied ANOVA stats tests across the four categories (NLH, ILH, HLH, Undefined). We kept the features with a p-value was less than 0.05.
	Note: We did not apply post-hoc testing to identify WHICH pairwise test was significant. We are only interested in the features that could be helpful upon discrimination
DR2. We calculated the absolute pairwise correlations across all the features (as per .corr() function). We identified feature pairs where corr>0.95. We then removed the 1/2 of the feature pairs

This leaves us with a filtered dataframe, which will be propogated into the PCA + HDBSCAN steps
'''

# DR1. ANOVA procedure
pval_grab = []
cols_grab = []

for i in list(big_train.filter(like='zscore').columns):
	lol = stats.f_oneway(big_train[i][big_train['mechID_homo'] == 'NLH'],big_train[i][big_train['mechID_homo'] == 'ILH'],big_train[i][big_train['mechID_homo'] == 'NLH'], big_train[i][big_train['mechID_homo'] == 'Undefined'], )
	pval_grab.append(lol[1])
	cols_grab.append(i)

qval_grab = np.array(multipletests(pval_grab)[1])
index_grab = np.where(qval_grab<0.05)[0]
col_want = [cols_grab[i] for i in index_grab]
col_want.append('mechID_homo')
col_want.append('chromo')
big_train_filt = big_train[col_want]

# DR2. Correlation procedure
corr_matrix = big_train_filt.filter(like='zscore').corr().abs()
upper = corr_matrix.where(np.triu(np.ones(corr_matrix.shape), k=1).astype(bool))
to_drop = [column for column in upper.columns if any(upper[column] > 0.95)]
columns_to_keep = list(set(big_train_filt.columns) - set(to_drop))
big_train_filt_reduced = big_train_filt[columns_to_keep]
split1_train = big_train_filt_reduced

'''
The following step uses the # dimension reduction parameter provided and undergoes PCA. 
Recall: We only do this columns with a z-score calculated. This implies that some of the features we annotated will be ignored during this procedure
'''

x_scaled = StandardScaler().fit_transform(split1_train.filter(like='zscore'))
pca = PCA(n_components=pcs_param)
pca_features = pca.fit_transform(x_scaled)

pcs_naming = []
for i in range(1,pcs_param+1):
	take = 'PCA'+str(i)
	pcs_naming.append(take)

# We also save the reduced vectors
pca_df = pd.DataFrame(data=pca_features,columns=pcs_naming)
pca_df['target'] = list(split1_train['mechID_homo'])

'''
The following step applies HDBSCAN itself
1. We pass PCA-reduced matrix to HDBSCAN, specifying the provided minimum size of cluster and distance metric
	Note: We also pass prediction_data=True. This parameter will allow for future use of this clusterer. For this script, this is important since we will need to obtain the membership vector (holds the probabilities)
2. We obtain the membership vector, which holds an SV's membership probability (ie. how confident SV1 belongs to cluster1 vs cluster2)
	Note: We acknowledge, in some cases, no or only one cluster is generated. Therefore, in this case, we exit the script
3. We create the fuzzy and sumprobs columns. These represent the SV's membership as per its greatest membership probability AND the sum of its membership probabilities
4. We then call the top_two_probs_diff function using the top two highest-probability clusters. This is kept in the diffs column
5. We then loop through all the rows and add the Forced column. This column assesses two factors:
	5.1. Diffs>0.25. This means the SV's top two clusters' probability are greater than 0.25. This implies that an SV's membership has confidently assigned to one cluster over the other
 	5.2. sumprobs>0.5. This means the sum of the membership probabilities is greater than 0.5. This implies the SV's was confidently contributing to the clusters INSTEAD of noise
  	Together, these two conditions must be met in order for an SV to be considered "high quality". In this case, the Forced == True
6. We then use a groupby to quantify how many SVs are high confidence, as per their cluster assignment. 
'''

# Applying the hdbscan itself
split1_train_att = pd.concat([split1_train, pca_df], axis=1)
clusterer = hdbscan.HDBSCAN(min_cluster_size=min_cluster_size_param, metric = metric_param, prediction_data=True).fit(pca_features)
soft_clusters = hdbscan.all_points_membership_vectors(clusterer)

# Saving the HDBSCAN model
joblib.dump(clusterer, '/home/nboev/projects/def-sushant/nboev/analysis/'+project+'/'+loc+'/IDmechsvSIM/20240625/splitsChromo/models/HDBSCANtraining'+typer+'.HOMO.sav')

if soft_clusters.shape[1] <1:
	print('These params do not produce suff clusters')
	sys.exit("These params do not produce suff clusters")

print('number of clusters identified', soft_clusters.shape[1])
for j in range(0, soft_clusters.shape[1]):
	split1_train_att['fuzzy'+str(j)] = list(soft_clusters[:,j])

fuzzy_columns = split1_train_att.loc[:,split1_train_att.columns.str.startswith('fuzzy')].columns
split1_train_att['fuzzy'] = split1_train_att.loc[:,fuzzy_columns].idxmax(axis=1)
split1_train_att['sumprobs'] = split1_train_att.loc[:, fuzzy_columns].sum(axis=1)
diffs = np.array([top_two_probs_diff(x) for x in soft_clusters])
split1_train_att['diffs'] = list(diffs)

for idx, row in split1_train_att.iterrows():
	if (row.diffs > 0.25) & (row.sumprobs > 0.5):
		split1_train_att.loc[idx, 'Forced'] = True
	else:
		split1_train_att.loc[idx, 'Forced'] = False

print(split1_train_att.groupby(['fuzzy', 'Forced']).size())

'''
The following steps were used to display the results from these optimal parameters (OP).
OP1. Calculating Silhouette scores using all the points (OP1.1) and the high quality points (OP1.2)
OP2. Identifying clusters with enrichment in a specific homology-based label. We then take NLH and HLH enriched clusters and apply Fisher Exact tests to obtain the odds ratio (OR)
	We consider the homology-label's proportion in the cluster. The greatest proportion is awarded the pseudo-assignment. We do this with high quality cluster assignments only. 
OP3. Quantifying the representation for the homology-labels. In this case, we compare the proportion of homology-based labels that are considered "high quality". 
	We set 50% as the minimum threshold. 
	Note: This does not address cluster enrichment, only label representation.
'''

# OP1.1. Calculating Sihouette scores with all the SVs
print('Calculating silhouette score')
print('including all the entries')
clusterlabel_hdbscan_fuzzy = list(split1_train_att['fuzzy'])
print(silhouette_score(pca_features, clusterlabel_hdbscan_fuzzy))

# OP1.2. Calculating Sihouette scores after filtering out for high quality cluster assignments
forced_index = list(split1_train_att[split1_train_att.Forced == True].index)
forced_df = split1_train_att.iloc[forced_index]
forced_pca = pca_features[forced_index]
print('including forced entries')
clusterlabel_hdbscan_fuzzy = list(forced_df['fuzzy'])
print(silhouette_score(forced_pca, clusterlabel_hdbscan_fuzzy))

# OP2. Identifying and quantifying cluster enrichment, as per homology/threshold based labels
clusterids = []
labels = []

for j in range(0, soft_clusters.shape[1]):
	a = len(split1_train_att[(split1_train_att.Forced ==True) & (split1_train_att.fuzzy=='fuzzy'+str(j))& (split1_train_att.target == 'HLH')])
	b = len(split1_train_att[(split1_train_att.Forced ==True) & (split1_train_att.fuzzy=='fuzzy'+str(j))& (split1_train_att.target == 'ILH')])
	c = len(split1_train_att[(split1_train_att.Forced ==True) & (split1_train_att.fuzzy=='fuzzy'+str(j))& (split1_train_att.target == 'NLH')])
	d = len(split1_train_att[(split1_train_att.Forced ==True) & (split1_train_att.fuzzy=='fuzzy'+str(j))& (split1_train_att.target == 'Undefined')])
	clusterids.append([a,b,c,d])
	labels.append('fuzzy'+str(j))
table = np.array(clusterids)
print(table)

tabledf = pd.DataFrame(table.T, columns = labels, index = ['HLH', 'ILH', 'NLH', 'Undefined'])
tabledf_prop = (tabledf.T/tabledf.sum(axis=1)).T
print(tabledf_prop)

# Looping through the soft-clusters and identifying enrichments in NLH/HLH
nhejlike = []
hrlike = []
for j in range(0, soft_clusters.shape[1]):
	if tabledf_prop.idxmax(axis=0)['fuzzy'+str(j)] == 'NLH':
		print('pseudo identity of fuzzy'+str(j))
		nhejlike.append('fuzzy'+str(j))
	if tabledf_prop.idxmax(axis=0)['fuzzy'+str(j)] == 'HLH':
		hrlike.append('fuzzy'+str(j))


tabledf['NLH_like'] = tabledf[nhejlike].sum(axis=1)
tabledf['HLH_like'] = tabledf[hrlike].sum(axis=1)

a = tabledf.loc['NLH', 'NLH_like']
b = tabledf.loc['HLH', 'NLH_like']
c = tabledf.loc['NLH', 'HLH_like']
d = tabledf.loc['HLH', 'HLH_like']

print([[a,b],[c,d]])
res = fisher_exact([[a,b],[c,d]])
print('OR' ,res[0])

grouped = split1_train_att.groupby(['target', 'Forced']).size().reset_index(name='count')
grouped['proportion'] = grouped.groupby('target')['count'].transform(lambda x: x / x.sum())
print(grouped)

# Saving result from the optimal clusters
big_train.to_csv('/home/nboev/projects/def-sushant/nboev/analysis/'+project+'/'+loc+'/IDmechsvSIM/20240625/splitsChromo/models/results/trainingALL'+typer+'.HOMO.tsv', sep='\t', index=False)
split1_train_att.to_csv('/home/nboev/projects/def-sushant/nboev/analysis/'+project+'/'+loc+'/IDmechsvSIM/20240625/splitsChromo/models/results/trainingSHORT'+typer+'.HOMO.tsv', sep='\t', index=False)

'''
The following steps were used to generate the KNN model
KNN1. Making 6-fold splits across different subsets of chromosomes. We chose this method since we wanted to subject the model to chromosome-specific patterns.
	Essentially, among the training set, we hold-out one of the "splits" here and train. 
	For this dataset, we calculate, (int(np.sqrt(len(tocluster))), which represents the number of "high quality" points in the dataset after clustering. 
 	For each fold, we try a k # of neighbours as per this value, and the +/-50 values, with 10 as a step size. For each fold we apply the KNN and calculate the mean and standard deviation in accuracy, as per the cluster assignments. 
KNN2. We found that mean accuracy, regardless of the split remains high. We did find variation could differ. In this step, we identififed the 'k' which varied the least.
'''

# KNN1. Identifying the ideal number of 'k', as per cross validation
tocluster = split1_train_att[split1_train_att.Forced == True]
st_x= StandardScaler()

# The random 6-splits we used
split1 = ['chr1', 'chr8', 'chr19']
split2 = ['chr6', 'chr16', 'chrX']
split3 = ['chr2', 'chr12', 'chrY']
split4 = ['chr10', 'chr11', 'chr21']
split5 = ['chr7', 'chr13', 'chr18']
split6 = ['chr4', 'chr9', 'chr15']

splitting = [split1, split2, split3, split4, split5, split6]
scoring_avg = []
scoring_std = []
scoring_raw = []

# Looping through the different k-values
k_values = [i for i in range (int(np.sqrt(len(tocluster)))-50,int(np.sqrt(len(tocluster)))+50,10)]
print(k_values)
for i in k_values:
	print(i)
	scoring = []
	for j in splitting:
		print(j)
		train = tocluster[tocluster['chromo'].isin(j)]
		print(len(train))
		test = tocluster[~tocluster['chromo'].isin(j)]
		print(len(test))

		pca = PCA(n_components=pcs_param)
		X_train = train.filter(like='zscore')
		X_train= st_x.fit_transform(X_train)
		X_train = pca.fit_transform(X_train)
		y_train = train.fuzzy
		X_test = test.filter(like='zscore')
		X_test= st_x.transform(X_test)
		X_test = pca.transform(X_test)
		y_test = test.fuzzy
		knn = KNeighborsClassifier(n_neighbors=i, metric = 'euclidean')
		knn.fit(X_train, y_train)
		y_pred = knn.predict(X_test)
		scoring.append(knn.score(X_test, y_test))
		print(knn.score(X_test, y_test))
		pca_df_test = pd.DataFrame(data=X_test, columns=['PC1', 'PC2',  'PC3', 'PC4', 'PC5',] )
		pca_df_test['preds'] = list(y_pred)
	scoring_raw.append(scoring)
	scoring_avg.append(np.mean(scoring))
	scoring_std.append(np.std(scoring))


print('raw scoring')
print(scoring_raw)
print('average scoring')
print(scoring_avg)
print('std scoring')
print(scoring_std)

# KNN2. Optimal # of neighbours to use during the KNN training
print('the k were going to use', k_values[np.argmin(scoring_std)])
ktouse = k_values[np.argmin(scoring_std)]

# Next, we applied these metrics to the whole training dataset and saved the models
train = tocluster
X_train = train.filter(like='zscore')
X_train= st_x.fit_transform(X_train)
pca = PCA(n_components=5)
X_train = pca.fit_transform(X_train)

# Freezing and saving the scaling, PCA and KNN models  
joblib.dump(st_x, '/home/nboev/projects/def-sushant/nboev/analysis/'+project+'/'+loc+'/IDmechsvSIM/20240625/splitsChromo/models/SCALERtraining'+typer+'.HOMO.sav')
joblib.dump(pca, '/home/nboev/projects/def-sushant/nboev/analysis/'+project+'/'+loc+'/IDmechsvSIM/20240625/splitsChromo/models/PCAtraining'+typer+'.HOMO.sav')
y_train = train.fuzzy
knn = KNeighborsClassifier(n_neighbors=ktouse, metric = 'euclidean')
knn.fit(X_train, y_train)
joblib.dump(knn, '/home/nboev/projects/def-sushant/nboev/analysis/'+project+'/'+loc+'/IDmechsvSIM/20240625/splitsChromo/models/KNNtraining'+typer+'.HOMO.sav')

# The final accuracy of the KNN model
print(knn.score(X_train, y_train))

# Applying the model to the held out dataset
colstokeep = tocluster.filter(like='zscore').columns
print(*colstokeep, sep='\n')
big_test = big[big.chromo.isin(testchromo)]
big_test = big_test.reset_index(drop=True)
bigtestX = big_test[colstokeep]
X_test = st_x.transform(bigtestX)
X_test = pca.transform(X_test)
y_pred = knn.predict(X_test)
big_test['fuzzypred'] = y_pred

print(big_test.groupby('fuzzypred').size())






