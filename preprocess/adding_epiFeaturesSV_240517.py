import glob
import numpy as np
import pandas as pd
from sys import argv
import datetime
from pytz import timezone
import matplotlib.pyplot as plt
import matplotlib as mpl
import pyBigWig
from scipy.stats import linregress

print('\n **********************')
print('START TIME:', datetime.datetime.now(timezone('EST')))

# Writing this so we can annotate all the files at once
svs = pd.read_csv(str(argv[1]), comment='#', sep='\t')
project = str(argv[2])
loc = str(argv[3])
filename = str(argv[4])
chr = str(argv[5])
episite = str(argv[6])
celline = str(argv[7])

svs = svs[svs['CHROM'] == chr]

# we're going to make it a fct bc we're going to force all the epis

def add_avgepi(epi_file, epi, df):

	to_open = episite+epi_file
	bb = pyBigWig.open(to_open)

# we can only calculate these values for the deletions
	for idx, row in df[(df.SV_logic == True)&((df.SV_Type == 'deletion') | (df.SV_Type == 'DEL'))] .iterrows():
		try:
			avg = bb.stats(row.CHROM, int(row.POS), int(row.POS+row.SVlen))
			df.loc[idx, 'avgepi_'+epi] = avg

			xs = list(range(1, int(row.SVlen)+1))
			ys = bb.values(row.CHROM, int(row.POS), int(row.POS+row.SVlen))
			slope = linregress(xs, ys)[0]
			df.loc[idx, 'slopeepi_'+epi] = slope
			max = bb.stats(row.CHROM, int(row.POS), int(row.POS+row.SVlen), type="max" )
			df.loc[idx, 'maxepi_'+epi] = max
			min = bb.stats(row.CHROM, int(row.POS), int(row.POS+row.SVlen), type="min" )
			df.loc[idx, 'minepi_'+epi] = min
			std = bb.stats(row.CHROM, int(row.POS), int(row.POS+row.SVlen), type="std" )
			df.loc[idx, 'stdepi_'+epi] = std
			cov = bb.stats(row.CHROM, int(row.POS), int(row.POS+row.SVlen), type="coverage" )
			df.loc[idx, 'covepi_'+epi] = cov
		except RuntimeError:
			df.loc[idx, 'avgepi_'+epi] = 0
			df.loc[idx, 'slopeepi_'+epi] =0
			df.loc[idx, 'maxepi_'+epi] = 0
			df.loc[idx, 'minepi_'+epi] = 0
			df.loc[idx, 'stdepi_'+epi] = 0
			df.loc[idx, 'covepi_'+epi] = 0

	return df

if celline == 'H1':
# using the diff epi files
	newsvs = add_avgepi('ENCBS718AAA/CTCF/ENCFF269OPL.bigWig','CTCF', svs)
	newsvs = add_avgepi('ENCBS718AAA/H2AFZ/ENCFF758YFI.bigWig','H2AFZ', newsvs)
	newsvs = add_avgepi('ENCBS718AAA/H3K36me3/ENCFF687LJF.bigWig','H3K36me3', newsvs)
	newsvs = add_avgepi('ENCBS718AAA/H3K9ac/ENCFF890MIB.bigWig','H3K9ac', newsvs)
	newsvs = add_avgepi('ENCBS718AAA/H3K9me3/ENCFF358AWN.bigWig','H3K9me3', newsvs)
	newsvs = add_avgepi('ENCBS718AAA/H4K20me1/ENCFF453OCL.bigWig','H4K20me1', newsvs)
	newsvs = add_avgepi('ENCBS718AAA/H3K4me1/ENCFF335ZGP.bigWig','H3K4me1', newsvs)
	newsvs = add_avgepi('ENCBS718AAA/H3K79me2/ENCFF640WRD.bigWig','H3K79me2', newsvs)
	newsvs = add_avgepi('ENCBS718AAA/H3K27ac/ENCFF771GNB.bigWig','H3K27ac', newsvs)
	newsvs = add_avgepi('ENCBS718AAA/H3K27me3/ENCFF277UCT.bigWig','H3K27me3', newsvs)
	newsvs = add_avgepi('ENCBS111ENC/DNase-seq/ENCFF573NKX.bigWig','DNase-seq', newsvs)
	newsvs = add_avgepi('ENCBS111ENC/WGB-Seq/ENCFF737DXC.bigWig','WGB-Seq+', newsvs)
	newsvs = add_avgepi('ENCBS111ENC/WGB-Seq/ENCFF380UXR.bigWig','WGB-Seq-', newsvs)
elif celline == 'HelaS3':
	newsvs = add_avgepi('ENCSR043WJN/CTCF/ENCFF152FJU.bigWig','CTCF', svs)
	newsvs = add_avgepi('ENCSR043WJN/H2AFZ/ENCFF883SLT.bigWig','H2AFZ', newsvs)
	newsvs = add_avgepi('ENCSR043WJN/H3K36me3/ENCFF533LRJ.bigWig','H3K36me3', newsvs)
	newsvs = add_avgepi('ENCSR043WJN/H3K9ac/ENCFF947CKC.bigWig','H3K9ac', newsvs)
	newsvs = add_avgepi('ENCSR043WJN/H3K9me3/ENCFF395VTK.bigWig','H3K9me3', newsvs)
	newsvs = add_avgepi('ENCSR043WJN/H4K20me1/ENCFF512ZJR.bigWig','H4K20me1', newsvs)
	newsvs = add_avgepi('ENCSR043WJN/H3K4me1/ENCFF884DQE.bigWig','H3K4me1', newsvs)
	newsvs = add_avgepi('ENCSR043WJN/H3K79me2/ENCFF385NSG.bigWig','H3K79me2', newsvs)
	newsvs = add_avgepi('ENCSR043WJN/H3K27ac/ENCFF023UNN.bigWig','H3K27ac', newsvs)
	newsvs = add_avgepi('ENCSR043WJN/H3K27me3/ENCFF175SKU.bigWig','H3K27me3', newsvs)
	newsvs = add_avgepi('ENCSR043WJN/DNase-seq/ENCFF194NJI.bigWig','DNase-seq', newsvs)
	newsvs = add_avgepi('ENCSR043WJN/WGB-Seq/ENCFF193DAN.bigWig','WGB-Seq+', newsvs)
	newsvs = add_avgepi('ENCSR043WJN/WGB-Seq/ENCFF347UWH.bigWig','WGB-Seq-', newsvs)

print(newsvs.head())

newsvs.to_csv('/home/nboev/projects/def-sushant/nboev/preprocess/'+project+'/'+loc+'/ENCODE/'+celline+'/'+chr+'/'+filename+'.EpiFeatSVs.csv', sep='\t', index=False)

print('end of script')
print('END TIME:', datetime.datetime.now(timezone('EST')))

