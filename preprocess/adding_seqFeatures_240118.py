import glob
import numpy as np
import pandas as pd
from sys import argv
import datetime
from pytz import timezone
from Bio.SeqUtils import gc_fraction
from Bio.SeqUtils import lcc
import matplotlib.pyplot as plt
import matplotlib as mpl

print('\n **********************')
print('START TIME:', datetime.datetime.now(timezone('EST')))

# Getting the vcf and opening with pandas
print(argv[1])
svs = pd.read_csv(str(argv[1]), comment='#', sep='\t')
print('PRE ANNOTATION')
print(svs.head())

project = str(argv[2])
loc = str(argv[3])
filename = str(argv[4])
chr = str(argv[5])
size = str(argv[6])

# calculating GC content
# uses the built in GC function from biopython--> calculates across the whole 2000 flank
def addGC(df, flank):
	pre_gc = "pre_flankgc_"+flank
	post_gc = "post_flankgc_"+flank

	for idx, row in df.iterrows():
		if (row.SV_logic == True) :
			df.loc[idx, pre_gc] = gc_fraction(str(row["pre_flank_seq_"+flank]))
			df.loc[idx, post_gc] = gc_fraction(str(row["post_flank_seq_"+flank]))
	return df

# calling the addGC function
svs = addGC(svs, size)
print('after adding GC content')
print(svs.head())

def shannon_entropy(sequence, window_size):
	entropy_values = []
	for i in range(len(sequence) - window_size + 1):
		window = sequence[i:i+window_size]
		entropy = -sum((window.count(base)/window_size) * (np.log2(window.count(base)/window_size) or 0) for base in set(window))
	entropy_values.append(entropy)
	return np.mean(entropy_values)

# uses the built in GC function from biopython--> calculates across the whole 2000 flank
def addcomp(df, flank):
	pre_comp = "pre_flankcomp_"+flank
	post_comp = "post_flankcomp_"+flank

	for idx, row in df.iterrows():
		if (row.SV_logic == True) :
			df.loc[idx, pre_comp] = shannon_entropy(str(row["pre_flank_seq_"+flank]),10)
			df.loc[idx, post_comp] = shannon_entropy(str(row["post_flank_seq_"+flank]),10)
	return df

# calling the addcomp function
svs = addcomp(svs, size)
print('after adding complexity')
print(svs.head())

# adding functions built by Sushant
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


def addFlexStab(df, flank):
	pre_flex = "pre_flankflex_"+flank
	post_flex = "post_flankflex_"+flank
	pre_stab = "pre_flankstab_"+flank
	post_stab = "post_flankstab_"+flank

	for idx, row in df.iterrows():
		if (row.SV_logic == True) :
			df.loc[idx, pre_flex] = calc_flexibility(str(row["pre_flank_seq_"+flank]))
			df.loc[idx, post_flex] = calc_flexibility(str(row["post_flank_seq_"+flank]))
			df.loc[idx, pre_stab] = calc_stability(str(row["pre_flank_seq_"+flank]))
			df.loc[idx, post_stab] = calc_stability(str(row["post_flank_seq_"+flank]))
	return df

# calling the flexstab function
svs = addFlexStab(svs, size)
print('after adding complexity')
print(svs.head())

# saving the csv, increasing generalization
#svs.to_csv('/home/nboev/projects/def-sushant/nboev/preprocess/HGSVC2_v2_integrated_callset/sv/seqFeatures/' +filename+ '.with' +str(argv[2])+ 'SeqFeat.csv', sep='\t', index=False)
#svs.to_csv('/home/nboev/projects/def-sushant/nboev/preprocess/'+project+'/'+loc+'/seqFeatures/' +filename+ '.with' +size+ 'SeqFeat.csv', sep='\t', index=False)
svs.to_csv('/home/nboev/projects/def-sushant/nboev/preprocess/'+project+'/'+loc+'/seqFeatures/' +chr+ '/' +filename+ '.with'+size+'SeqFeat.csv', sep='\t', index=False)


print('end of script')
print('END TIME:', datetime.datetime.now(timezone('EST')))

