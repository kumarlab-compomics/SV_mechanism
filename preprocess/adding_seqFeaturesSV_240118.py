#!/usr/bin/env python3

import glob
import numpy as np
import pandas as pd
from sys import argv
import datetime
from pytz import timezone
from Bio.SeqUtils import gc_fraction
import matplotlib.pyplot as plt
import matplotlib as mpl

print('\n **********************')
print('START TIME:', datetime.datetime.now(timezone('EST')))

# in this version --> we're using the actual inserted or deleted sequence in order to calculate one value for each metric

# Getting the vcf and opening with pandas
print(argv[1])
svs = pd.read_csv(str(argv[1]), comment='#', sep='\t')
print('PRE ANNOTATION')
print(svs.head())

# added to generalize = 2023/10/30
project = str(argv[2])
loc = str(argv[3])
filename = str(argv[4])
chr = str(argv[5])

# adding the functions we need (provided by Sushant)
flextab = {"AA":7.6 , "CA":10.9, "GA":8.8 , "TA":12.5, "AC":14.6, "CC":7.2 , "GC":11.1, "TC":8.8, "AG":8.2 , "CG":8.9 , "GG":7.2 , "TG":10.9, "AT":25.0, "CT":8.2 , "GT":14.6, "TT":7.6}
helixtab = {"AA":1.9,  "CA":1.9,  "GA":1.6,  "TA":0.9, "AC":1.3,  "CC":3.1,  "GC":3.1,  "TC":1.6, "AG":1.6,  "CG":3.6,  "GG":3.1,  "TG":1.9, "AT":1.5,  "CT":1.6,  "GT":1.3,  "TT":1.9}

def calc_fragility(sequence, dinucleotides):
	frag = 0
	total = 0
	for i in range(len(sequence)-1):
		dint = str(sequence[i:i+2]).upper()

		if dint in dinucleotides:
			frag = frag + dinucleotides[dint]
			total = total + 1
	return frag/total if total > 0 else 0

def calc_flexibility(sequence):
	return calc_fragility(sequence, flextab)

def calc_stability(sequence):
	return calc_fragility(sequence, helixtab)

# were going to calculate mean Shannon Entropy instead of simple complexxity = 2024/01/18
def shannon_entropy(sequence, window_size):
	entropy_values = []
	for i in range(len(sequence) - window_size + 1):
		window = sequence[i:i+window_size]
		entropy = -sum((window.count(base)/window_size) * (np.log2(window.count(base)/window_size) or 0) for base in set(window))
	entropy_values.append(entropy)
	return np.mean(entropy_values)

for idx, row in svs.iterrows():
	if (row.SV_logic == True) & ((row.SV_Type == 'insertion')|(row.SV_Type == 'INS'))  :
		print(row.ALT)
		svs.loc[idx, 'var_gc'] = gc_fraction(str(row["ALT"]))
		svs.loc[idx, 'var_comp'] = shannon_entropy(str(row["ALT"]), 10)
		svs.loc[idx, 'var_flex'] = calc_flexibility(str(row["ALT"]))
		svs.loc[idx, 'var_stab'] = calc_stability(str(row["ALT"]))


	elif (row.SV_logic == True) & ((row.SV_Type == 'deletion') |(row.SV_Type == 'DEL')) :
		print(row.REF)
		svs.loc[idx, 'var_gc'] = gc_fraction(str(row["REF"]))
		svs.loc[idx, 'var_comp'] = shannon_entropy(str(row["REF"]), 10)
		svs.loc[idx, 'var_flex'] = calc_flexibility(str(row["REF"]))
		svs.loc[idx, 'var_stab'] = calc_stability(str(row["REF"]))


print(svs.head())


# saving the csv, need to generalize
svs.to_csv('/home/nboev/projects/def-sushant/nboev/preprocess/'+project+'/'+loc+'/seqFeaturesSV/' +chr+ '/' +filename+ '.withSVSeqFeat.csv', sep='\t', index=False)

print('end of script')
print('END TIME:', datetime.datetime.now(timezone('EST')))
