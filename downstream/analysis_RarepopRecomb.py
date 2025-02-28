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

# In this script, we agglomerate all the recombination-rate annotated bed files (for each individual, and for both cluster I and II), to calculate a mean recombination rate. 

# Grabbing all the bed files in the popRates directory, whereby these are already annotated with recombination rates
# We loop through all these files and split the filename to identify the population code, and cluster, then concatenate across
path='/home/nboev/projects/def-sushant/nboev/analysis/20220422_3202_phased_SNV_INDEL_SV_bychrom/SVTrue_typedeletion_resTrue/IDmechsvSIM/20240625/splitsChromo/bedRareInd/popRates/'
all_files = glob.glob(os.path.join(path , "*.bed"))
li1 = []

for filename in all_files:
	df = pd.read_csv(filename,  sep = '\t', header=None)
	df['Population code'] = filename.split('popRates/')[1].split('.')[0]
	df['ind'] = filename.split('popRates/')[1].split('.')[1]
	df['ClusterBasedLabel'] = filename.split('popRates/')[1].split('.')[3]

	li1.append(df)
lol = pd.concat(li1, axis=0, ignore_index=True)

# We use a groupby the individual, population and cluster to calculate the mean recombination rate. This summarized dataset is saved as meanRareRecombRates.tsv
means = lol.groupby(['ind', 'Population code', 'ClusterBasedLabel'])[6].mean().reset_index()
means.to_csv("/home/nboev/projects/def-sushant/nboev/analysis/20220422_3202_phased_SNV_INDEL_SV_bychrom/SVTrue_typedeletion_resTrue/IDmechsvSIM/20240625/splitsChromo/bedRareInd/popRates/meanRareRecombRates.tsv" , sep='\t', header=True, index=False, )
