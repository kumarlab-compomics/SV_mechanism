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
flanksize = str(argv[7])
celline = str(argv[8])

svs = svs[svs['CHROM'] == chr]


def add_epislope(epi_file, epi, df, flank):
	to_open = episite+epi_file
	bb = pyBigWig.open(to_open)

	pre_name = "pre_flank" +epi+ "_" +flank
	post_name = "post_flank" +epi+ "_" +flank

	for idx, row in df[(df.SV_logic == True)].iterrows():
		if ((row.SV_Type == 'deletion') | (row.SV_Type == 'DEL')):
			avg_pre = bb.stats(row.CHROM, int(row.POS) - int(flank), int(row.POS))
			avg_post = bb.stats(row.CHROM, int(row.POS)+int(row.SVlen), int(row.POS) + int(row.SVlen) +  int(flank) )
			df.loc[idx, 'avg_'+pre_name] = avg_pre
			df.loc[idx, 'avg_'+post_name] = avg_post

			xs = list(range(1, int(flank)+1))
			ys_pre = bb.values(row.CHROM, int(row.POS) - int(flank), int(row.POS))
			ys_post = bb.values(row.CHROM, int(row.POS)+int(row.SVlen), int(row.POS) + int(row.SVlen) +  int(flank) )
			m_pre = linregress(xs, ys_pre)[0]
			m_post = linregress(xs, ys_post)[0]
			df.loc[idx, 'slope_'+pre_name] = m_pre
			df.loc[idx, 'slope_'+post_name] = m_post

			max_pre = bb.stats(row.CHROM, int(row.POS) - int(flank), int(row.POS), type="max")
			max_post = bb.stats(row.CHROM, int(row.POS)+int(row.SVlen), int(row.POS) + int(row.SVlen) +  int(flank), type="max" )
			df.loc[idx, 'max_'+pre_name] = max_pre
			df.loc[idx, 'max_'+post_name] = max_post

			min_pre = bb.stats(row.CHROM, int(row.POS) - int(flank), int(row.POS), type="min")
			min_post = bb.stats(row.CHROM, int(row.POS)+int(row.SVlen), int(row.POS) + int(row.SVlen) +  int(flank), type="min" )
			df.loc[idx, 'min_'+pre_name] = min_pre
			df.loc[idx, 'min_'+post_name] = min_post

			std_pre = bb.stats(row.CHROM, int(row.POS) - int(flank), int(row.POS), type="std")
			std_post = bb.stats(row.CHROM, int(row.POS)+int(row.SVlen), int(row.POS) + int(row.SVlen) +  int(flank), type="std" )
			df.loc[idx, 'std_'+pre_name] = std_pre
			df.loc[idx, 'std_'+post_name] = std_post

			cov_pre = bb.stats(row.CHROM, int(row.POS) - int(flank), int(row.POS), type="coverage")
			cov_post = bb.stats(row.CHROM, int(row.POS)+int(row.SVlen), int(row.POS) + int(row.SVlen) +  int(flank), type="coverage" )
			df.loc[idx, 'cov_'+pre_name] = cov_pre
			df.loc[idx, 'cov_'+post_name] = cov_post


		elif ((row.SV_Type == 'insertion') | (row.SV_Type == 'INS')):

			avg_pre = bb.stats(row.CHROM, int(row.POS) - int(flank), int(row.POS))
			avg_post = bb.stats(row.CHROM, int(row.POS), int(row.POS) + int(flank) )
			df.loc[idx, 'avg_'+pre_name] = avg_pre
			df.loc[idx, 'avg_'+post_name] = avg_post

			xs = list(range(1, int(flank)+1))
			ys_pre = bb.values(row.CHROM, int(row.POS) - int(flank), int(row.POS))
			ys_post = bb.values(row.CHROM, int(row.POS), int(row.POS) + int(flank) )
			m_pre = linregress(xs, ys_pre)[0]
			m_post = linregress(xs, ys_post)[0]
			df.loc[idx, 'slope_'+pre_name] = m_pre
			df.loc[idx, 'slope_'+post_name] = m_post


			max_pre = bb.stats(row.CHROM, int(row.POS) - int(flank), int(row.POS), type="max")
			max_post = bb.stats(row.CHROM, int(row.POS), int(row.POS) + int(flank) , type="max" )
			df.loc[idx, 'max_'+pre_name] = max_pre
			df.loc[idx, 'max_'+post_name] = max_post

			min_pre = bb.stats(row.CHROM, int(row.POS) - int(flank), int(row.POS), type="min")
			min_post = bb.stats(row.CHROM, int(row.POS), int(row.POS) + int(flank) , type="min" )
			df.loc[idx, 'min_'+pre_name] = min_pre
			df.loc[idx, 'min_'+post_name] = min_post

			std_pre = bb.stats(row.CHROM, int(row.POS) - int(flank), int(row.POS), type="std")
			std_post = bb.stats(row.CHROM, int(row.POS), int(row.POS) + int(flank) , type="std" )
			df.loc[idx, 'std_'+pre_name] = std_pre
			df.loc[idx, 'std_'+post_name] = std_post

			cov_pre = bb.stats(row.CHROM, int(row.POS) - int(flank), int(row.POS), type="coverage")
			cov_post = bb.stats(row.CHROM, int(row.POS), int(row.POS) + int(flank) , type="coverage" )
			df.loc[idx, 'cov_'+pre_name] = cov_pre
			df.loc[idx, 'cov_'+post_name] = cov_post

	return df

if celline == 'H1':
	newsvs = add_epislope('ENCBS718AAA/CTCF/ENCFF269OPL.bigWig','CTCF', svs, flanksize)
	newsvs = add_epislope('ENCBS718AAA/H2AFZ/ENCFF758YFI.bigWig','H2AFZ', newsvs, flanksize)
	newsvs = add_epislope('ENCBS718AAA/H3K36me3/ENCFF687LJF.bigWig','H3K36me3', newsvs, flanksize)
	newsvs = add_epislope('ENCBS718AAA/H3K9ac/ENCFF890MIB.bigWig','H3K9ac', newsvs, flanksize)
	newsvs = add_epislope('ENCBS718AAA/H3K9me3/ENCFF358AWN.bigWig','H3K9me3', newsvs, flanksize)
	newsvs = add_epislope('ENCBS718AAA/H4K20me1/ENCFF453OCL.bigWig','H4K20me1', newsvs, flanksize)
	newsvs = add_epislope('ENCBS718AAA/H3K4me1/ENCFF335ZGP.bigWig','H3K4me1', newsvs, flanksize)
	newsvs = add_epislope('ENCBS718AAA/H3K79me2/ENCFF640WRD.bigWig','H3K79me2', newsvs, flanksize)
	newsvs = add_epislope('ENCBS718AAA/H3K27ac/ENCFF771GNB.bigWig','H3K27ac', newsvs, flanksize)
	newsvs = add_epislope('ENCBS718AAA/H3K27me3/ENCFF277UCT.bigWig','H3K27me3', newsvs, flanksize)
	newsvs = add_epislope('ENCBS111ENC/DNase-seq/ENCFF573NKX.bigWig','DNase-seq', newsvs, flanksize)
	newsvs = add_epislope('ENCBS111ENC/WGB-Seq/ENCFF737DXC.bigWig','WGB-Seq+', newsvs, flanksize)
	newsvs = add_epislope('ENCBS111ENC/WGB-Seq/ENCFF380UXR.bigWig','WGB-Seq-', newsvs, flanksize)

	print(newsvs.head())

# UPDATE THIS SECTION ONCE ALL THE FILES ARE DOWNLOADED
elif celline == 'HelaS3':
	newsvs = add_epislope('ENCSR043WJN/CTCF/ENCFF152FJU.bigWig','CTCF', svs, flanksize)
	newsvs = add_epislope('ENCSR043WJN/H2AFZ/ENCFF883SLT.bigWig','H2AFZ', newsvs, flanksize)
	newsvs = add_epislope('ENCSR043WJN/H3K36me3/ENCFF533LRJ.bigWig','H3K36me3', newsvs, flanksize)
	newsvs = add_epislope('ENCSR043WJN/H3K9ac/ENCFF947CKC.bigWig','H3K9ac', newsvs, flanksize)
	newsvs = add_epislope('ENCSR043WJN/H3K9me3/ENCFF395VTK.bigWig','H3K9me3', newsvs, flanksize)
	newsvs = add_epislope('ENCSR043WJN/H4K20me1/ENCFF512ZJR.bigWig','H4K20me1', newsvs, flanksize)
	newsvs = add_epislope('ENCSR043WJN/H3K4me1/ENCFF884DQE.bigWig','H3K4me1', newsvs, flanksize)
	newsvs = add_epislope('ENCSR043WJN/H3K79me2/ENCFF385NSG.bigWig','H3K79me2', newsvs, flanksize)
	newsvs = add_epislope('ENCSR043WJN/H3K27ac/ENCFF023UNN.bigWig','H3K27ac', newsvs, flanksize)
	newsvs = add_epislope('ENCSR043WJN/H3K27me3/ENCFF175SKU.bigWig','H3K27me3', newsvs, flanksize)
	newsvs = add_epislope('ENCSR043WJN/DNase-seq/ENCFF194NJI.bigWig','DNase-seq', newsvs, flanksize)
	newsvs = add_epislope('ENCSR043WJN/WGB-Seq/ENCFF193DAN.bigWig','WGB-Seq+', newsvs, flanksize)
	newsvs = add_epislope('ENCSR043WJN/WGB-Seq/ENCFF347UWH.bigWig','WGB-Seq-', newsvs, flanksize)

newsvs.to_csv('/home/nboev/projects/def-sushant/nboev/preprocess/'+project+'/'+loc+'/ENCODE/'+celline+'/'+chr+'/'+filename+'.'+flanksize+'.EpiFeat.csv', sep='\t', index=False)


print('end of script')
print('END TIME:', datetime.datetime.now(timezone('EST')))




