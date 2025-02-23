import glob
import numpy as np
import pandas as pd
from sys import argv
import datetime
from pytz import timezone
from scipy.stats import fisher_exact
from scipy.stats import norm
import os

print('\n **********************')
print('START TIME:', datetime.datetime.now(timezone('EST')))

# grabbing all the different bed files we need

print(str(argv[5]))
print(str(argv[6]))
print(str(argv[7]))
print(str(argv[8]))
print(str(argv[9]))

if os.path.getsize(str(argv[1])) == 0:
	a = 0
else:
	fuzzy0_in = pd.read_csv(str(argv[1]), comment='#', sep='\t', header=None, )
	a = len(fuzzy0_in)


if os.path.getsize(str(argv[2])) == 0:
	b = 0
else:
	fuzzy1_in = pd.read_csv(str(argv[2]), comment='#', sep='\t', header=None, )
	b = len(fuzzy1_in)


fuzzy0_all = pd.read_csv(str(argv[3]), comment='#', sep='\t', header=None, )
fuzzy1_all = pd.read_csv(str(argv[4]), comment='#', sep='\t', header=None, )

fracobs = a/len(fuzzy0_all)
fracexp = b/len(fuzzy1_all)
enrich = (fracobs-fracexp)/fracexp


lst = [str(argv[8]), str(argv[9]), str(argv[6]), enrich, str(argv[7]), a, len(fuzzy0_all), b, len(fuzzy1_all), str(argv[2])]
df = pd.DataFrame([lst])
df.columns =['Dataset', 'feature', 'svtype', 'Enrich', 'Comp', 'fracobs_num', 'fracobs_denom', 'fracexp_num', 'fracexp_denom', 'file']

df.to_csv('/home/nboev/projects/def-sushant/nboev/analysis/HGSVC2_v2_integrated_callset/sv/IDmechsvSIM/20240625/splitsChromo/bedWork100s/CDS/'+str(argv[9])+'analysis_hdbscanknnEnrichment.csv', mode='a', header=False)

