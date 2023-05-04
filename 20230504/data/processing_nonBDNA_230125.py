#!/usr/bin/env python3

import pandas as pd
import glob
import numpy as np
from sys import argv
import datetime
from pytz import timezone

print('\n **************************')
print('START TIME:', datetime.datetime.now(timezone('EST')))

# start by input of nonBDNA data --> We need to rename the file
print(argv[1])
chromo = str(argv[1]).split('/')[10].split('_')[0]
print('chr:', chromo)

gquad = pd.read_csv(str(argv[1]), delimiter = "\t")
print(argv[2])
# the annotation we're doing
annotation = str(argv[2])
gquad.rename(columns={'Type': annotation}, inplace=True)
# sort so order is smallest to greatest
gquad = gquad.sort_values(by=['Start'])
print(gquad.head())

# pulling the folder location by the annotation
loc = str(argv[1]).split('/')[8]

short = gquad[['Sequence_name', 'Start', 'Stop',annotation]]
short.to_csv('/home/nboev/projects/def-sushant/nboev/data/nonBDNA/' +loc+ '/processed/' +chromo+ '.' +annotation+  'forbedtools.bed', sep='\t', header = False, index=False)


# we need to pull the data of where our "max" chromosome position is
print(argv[3])
band = pd.read_csv(str(argv[3]), delimiter = "\t")
chromoband = band[band['#chrom'] == chromo]
print(chromoband.head())
max_end = max(chromoband.chromEnd)

# getting the combos of the nonBDNA we have listed in file
listing = gquad[['Start', 'Stop']].values.tolist()

for i in range(len(listing)):
	if i == 0:
		nexter = listing[i][0] - 1
		new_row = {'Sequence_name':chromo, annotation:'Not_Motif', 'Start':0, 'Stop': nexter}
		gquad = gquad.append(new_row, ignore_index=True)
	elif i == len(listing)-1:
		zero = listing[i-1][1]+1
		nexter = max_end
		new_row = {'Sequence_name':chromo, annotation:'Not_Motif', 'Start':zero, 'Stop': nexter}
		gquad = gquad.append(new_row, ignore_index=True)
	else:
		zero = listing[i-1][1]+1
		nexter = listing[i][0] - 1
		new_row = {'Sequence_name': chromo, annotation:'Not_Motif', 'Start':zero, 'Stop': nexter}
		gquad = gquad.append(new_row, ignore_index=True)


gquad = gquad.sort_values('Start')
print(gquad.head())
print(gquad.tail())

gquad.to_csv('/home/nboev/projects/def-sushant/nboev/data/nonBDNA/' +loc+ '/processed/' +chromo+ '.' +annotation+  '.csv', sep='\t', index=False)


print('END TIME:', datetime.datetime.now(timezone('EST')))

