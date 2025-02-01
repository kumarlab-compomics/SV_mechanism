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

# In this script, we take in a vcf/csv batch file and annotate the SVs themselves with their epigenomic features based on ENCODE pybigwigs
# This annotation requires the use of the pyBigWig package

# Data inputs sent in from the execution files
svs = pd.read_csv(str(argv[1]), comment='#', sep='\t')
project = str(argv[2])
loc = str(argv[3])
filename = str(argv[4])
chr = str(argv[5])
episite = str(argv[6])
celline = str(argv[7])
svs = svs[svs['CHROM'] == chr]

'''
The function, add_avgepi, does the annotations itself
The input: 
	epi_file: Location of the pybigwig files
	epi: The name of the epigenomic feature
	df: The file to annotate
'''
def add_avgepi(epi_file, epi, df):
	to_open = episite+epi_file
	bb = pyBigWig.open(to_open)

	# We must filter out for deletions only since these insertions are not measured in the cell sequenced
	for idx, row in df[(df.SV_logic == True)&((df.SV_Type == 'deletion') | (df.SV_Type == 'DEL'))] .iterrows():
		try:
			# Appending the average, slope of the region, max, min, standard deviation and coverage from the bigwig
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

# Calling the add_avgepi function by cell line and bigwig file
if celline == 'H1':
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

# Save the final appended file
newsvs.to_csv('/home/nboev/projects/def-sushant/nboev/preprocess/'+project+'/'+loc+'/ENCODE/'+celline+'/'+chr+'/'+filename+'.EpiFeatSVs.csv', sep='\t', index=False)

print('end of script')
print('END TIME:', datetime.datetime.now(timezone('EST')))

