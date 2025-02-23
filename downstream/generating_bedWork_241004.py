#!/usr/bin/env python3

import pandas as pd
import glob
import os
from sys import argv
import datetime
from pytz import timezone
from Bio import SeqIO

print('\n ********************************')
print('START TIME:', datetime.datetime.now(timezone('EST')))

# loading in the csv (which holds all the chromo)
# we're only doing deletions for now

svs = pd.read_csv(str(argv[1]), comment='#', sep='\t')
project = str(argv[2])
loc = str(argv[3])
chromo = str(argv[4])
typer = str(argv[5])

# we need to find the POS from the ID
svs[['chr','POS', 'type', 'len']] = svs['ID'].str.split('-',expand=True)
svs["len"] = pd.to_numeric(svs["len"])
svs["POS"] = pd.to_numeric(svs["POS"])

for idx, row in svs.iterrows():
	if typer == 'deletion':
#		svs.loc[idx, 'start']= int(row.POS) - 2000
#		svs.loc[idx, 'stop']= int(row.POS + row.len) + 2000
		svs.loc[idx, 'start']= int(row.POS)
		svs.loc[idx, 'stop']= int(row.POS + row.len)

# were also going to do the version where we have one base for each breakpoint
		svs.loc[idx, 'start_up']= int(row.POS)
		svs.loc[idx, 'stop_up']= int(row.POS) +1
		svs.loc[idx, 'start_down']= int(row.POS + row.len)-1
		svs.loc[idx, 'stop_down']= int(row.POS + row.len)

	elif typer == 'insertion':
#		svs.loc[idx, 'start']= int(row.POS) - 2000
#		svs.loc[idx, 'stop']= int(row.POS) + 2000
		svs.loc[idx, 'start']= int(row.POS)
		svs.loc[idx, 'stop']= int(row.POS)+1

svs['start'] = svs['start'].astype('Int64')
svs['stop'] = svs['stop'].astype('Int64')


# HR related
svshr = svs[svs.mechID_homo=='HR']
svshr = svshr[['chr', 'start', 'stop']]
print(svshr.head())

svsssa = svs[svs.mechID_homo=='SSAaEJ']
svsssa = svsssa[['chr', 'start', 'stop']]
print(svsssa.head())

svsnhej = svs[svs.mechID_homo=='NHEJ']
svsnhej = svsnhej[['chr', 'start', 'stop']]
print(svsnhej.head())

svshr.to_csv('/home/nboev/projects/def-sushant/nboev/analysis/'+project+'/'+loc+'/IDmechsvSIM/20240625/splitsChromo/bedWork/'+chromo+'.'+typer+'.HOMO.HR.bed', sep='\t', index=False, header=False)
svsssa.to_csv('/home/nboev/projects/def-sushant/nboev/analysis/'+project+'/'+loc+'/IDmechsvSIM/20240625/splitsChromo/bedWork/'+chromo+'.'+typer+'.HOMO.SSAaEJ.bed', sep='\t', index=False, header=False)
svsnhej.to_csv('/home/nboev/projects/def-sushant/nboev/analysis/'+project+'/'+loc+'/IDmechsvSIM/20240625/splitsChromo/bedWork/'+chromo+'.'+typer+'.HOMO.NHEJ.bed', sep='\t', index=False, header=False)

work=['HR', 'SSAaEJ', 'NHEJ']

if typer == 'deletion':
	svs['start_up'] = svs['start_up'].astype('Int64')
	svs['stop_up'] = svs['stop_up'].astype('Int64')
	svs['start_down'] = svs['start_down'].astype('Int64')
	svs['stop_down'] = svs['stop_down'].astype('Int64')

	for i in work:
		svshr = svs[svs.mechID_homo==i]
		svshrup = svshr[['chr', 'start_up', 'stop_up']]
		svshrdown = svshr[['chr', 'start_down', 'stop_down']]

		svshrup.to_csv('/home/nboev/projects/def-sushant/nboev/analysis/'+project+'/'+loc+'/IDmechsvSIM/20240625/splitsChromo/bedWork/UP'+chromo+'.'+typer+'.HOMO.'+i+'.bed', sep='\t', index=False, header=False)
		svshrdown.to_csv('/home/nboev/projects/def-sushant/nboev/analysis/'+project+'/'+loc+'/IDmechsvSIM/20240625/splitsChromo/bedWork/DOWN'+chromo+'.'+typer+'.HOMO.'+i+'.bed', sep='\t', index=False, header=False)

