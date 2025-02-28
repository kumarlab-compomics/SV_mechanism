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

path='/home/nboev/projects/def-sushant/nboev/analysis/20220422_3202_phased_SNV_INDEL_SV_bychrom/SVTrue_typedeletion_resTrue/IDmechsvSIM/20240625/splitsChromo/bedRareInd/popRates/'
all_files = glob.glob(os.path.join(path , "*.bed"))
print(all_files)
li1 = []

for filename in all_files:
	df = pd.read_csv(filename,  sep = '\t', header=None)
	df['Superpopulation code'] = filename.split('popRates/')[1].split('.')[0]
	df['ind'] = filename.split('popRates/')[1].split('.')[1]
	df['fuzzy'] = filename.split('popRates/')[1].split('.')[3]

	li1.append(df)
lol = pd.concat(li1, axis=0, ignore_index=True)
means = lol.groupby(['ind', 'Superpopulation code', 'fuzzy'])[6].mean().reset_index()

means.to_csv("/home/nboev/projects/def-sushant/nboev/analysis/20220422_3202_phased_SNV_INDEL_SV_bychrom/SVTrue_typedeletion_resTrue/IDmechsvSIM/20240625/splitsChromo/bedRareInd/popRates/meanRareRecombRates.tsv" , sep='\t', header=True, index=False, )
