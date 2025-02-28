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

# In this script, we take a series of bed files, for real and simulated SVs, to calculated enrichment score (ES), for a given genomic element.
# The 1st and 2nd input argument are bed files which map to the intersections, for the real and simulated SVs, respectively. 
# The 3rd and 4th input arguments are the original bed files, for the real and simualted SVs, respectively.
# The remaining input arguments just help create traceable files names

# We check the sizes of the files and pulling out their lengths
if os.path.getsize(str(argv[1])) == 0:
	a = 0
else:
	real = pd.read_csv(str(argv[1]), comment='#', sep='\t', header=None, )
	a = len(real)

if os.path.getsize(str(argv[2])) == 0:
	b = 0
else:
	sim = pd.read_csv(str(argv[2]), comment='#', sep='\t', header=None, )
	b = len(sim)

# Pulling in the "entire" files
real_all = pd.read_csv(str(argv[3]), comment='#', sep='\t', header=None, )
sim_all = pd.read_csv(str(argv[4]), comment='#', sep='\t', header=None, )

# Calculating the fraction observed, and expected, by comparing the intersections to the entire sets. We then use this to calculate the enrichment score for this particular iteration
fracobs = a/len(real_all)
fracexp = b/len(sim_all)
enrich = (fracobs-fracexp)/fracexp

# Creating and saving a dataframe to hold data for all the iterations
lst = [str(argv[8]), str(argv[9]), str(argv[6]), enrich, str(argv[7]), a, len(real_all), b, len(sim_all), str(argv[2])]
df = pd.DataFrame([lst])
df.columns =['Dataset', 'feature', 'svtype', 'Enrich', 'Comp', 'fracobs_num', 'fracobs_denom', 'fracexp_num', 'fracexp_denom', 'file']
df.to_csv('/home/nboev/projects/def-sushant/nboev/analysis/HGSVC2_v2_integrated_callset/sv/IDmechsvSIM/20240625/splitsChromo/bedWork100s/CDS/'+str(argv[9])+'analysis_hdbscanknnEnrichment.csv', mode='a', header=False)
