#!/usr/bin/env python3

import pandas as pd
import glob
import numpy as np
from sys import argv
import os
import datetime
from pytz import timezone

print('\n**********************')
print('START TIME:', datetime.datetime.now(timezone('EST')))

project = str(argv[1])
loc = str(argv[2])
filename = str(argv[3])
chr = str(argv[4])
typer = str(argv[5])

# Seq features - variant
feat = pd.read_csv(str(argv[6]), comment='#', sep='\t')

types = feat.SV_Type.unique()
if ('deletion' in types) | ('insertion' in types) :
	if typer == 'DEL':
		type = 'deletion'
	elif typer == 'INS':
		type = 'insertion'
elif ('DEL' in types) | ('INS' in types):
	if typer == 'DEL':
		type = 'DEL'
	elif typer == 'INS':
		type = 'INS'

feat = feat[(feat.CHROM == chr) & (feat.SV_Type == type)]
print(feat.head())

# were missing some columns!! need to add the way a way to add ID here!!!
for idx, row in feat.iterrows():
	if (row.SV_Type == 'deletion') | (row.SV_Type == 'DEL') :
		feat.loc[idx, 'ID'] = row.CHROM +'-'+ str(row.POS+1) +'-'+ 'DEL' +'-'+ str(int(float(row.SVlen)))
	elif (row.SV_Type == 'insertion') | (row.SV_Type == 'INS') :
		feat.loc[idx, 'ID'] = row.CHROM +'-'+ str(row.POS+1) +'-'+ 'INS' +'-'+ str(int(float(row.SVlen)))

# REACLL!!! THIS WILL NEED TO CHANGE!!! WE NEED TO SAY.. CONTAINS!
if ('SIM' in loc ) | ('sim' in loc ) | ('Sim' in loc ):
	feat_filt = feat[['ID', 'var_gc', 'var_comp', 'var_flex', 'var_stab', 'Sim', 'Sim_ID']]
else:
	feat_filt = feat[['ID', 'var_gc', 'var_comp', 'var_flex', 'var_stab']]

feat_filt = feat_filt.replace(np. nan,0)
feat_filt = feat_filt.reset_index(drop=True)

# Getting the correct number of rows we should end with
ogrows = len(feat_filt)
print('number of rows we expect', ogrows)

print(feat_filt.head())
print('number of lines from feat_filt', len(feat_filt.ID.unique()))

# Seq features - flanks
flankfeat = pd.read_csv(str(argv[7]), comment='#', sep='\t')
flankfeat = flankfeat[(flankfeat.CHROM == chr) & (flankfeat.SV_Type == type)]

# were missing some columns
for idx, row in flankfeat.iterrows():
	if (row.SV_Type == 'deletion') | (row.SV_Type == 'DEL') :
		flankfeat.loc[idx, 'ID'] = row.CHROM +'-'+ str(row.POS+1) +'-'+ 'DEL' +'-'+ str(int(float(row.SVlen)))
	elif (row.SV_Type == 'insertion') | (row.SV_Type == 'INS') :
		flankfeat.loc[idx, 'ID'] = row.CHROM +'-'+ str(row.POS+1) +'-'+ 'INS' +'-'+ str(int(float(row.SVlen)))

flankfeat = flankfeat.loc[:, flankfeat.columns.str.contains('ID|pre|post')]
flankfeat = flankfeat.loc[:, ~flankfeat.columns.str.contains('_seq_')]
flankfeat = flankfeat.replace(np. nan,0)
flankfeat = flankfeat.reset_index(drop=True)

print(flankfeat.head())
print('number of lines from flankfeat', len(flankfeat.ID.unique()))


# Repeat Masker - variant
repSV = pd.read_csv(str(argv[8]), comment='#', sep='\t')
repSV = repSV[(repSV.CHROM == chr) & (repSV.SV_Type == type)]

# ALTERNG THIS TO INCLUDE SUBCLASSES!!
for idx, row in repSV.iterrows():
	if (row.SV_Type == 'deletion') | (row.SV_Type == 'DEL') :
		repSV.loc[idx, 'ID'] = row.CHROM +'-'+ str(row.POS+1) +'-'+ 'DEL' +'-'+ str(int(float(row.SVlen)))
	elif (row.SV_Type == 'insertion') | (row.SV_Type == 'INS') :
		repSV.loc[idx, 'ID'] = row.CHROM +'-'+ str(row.POS+1) +'-'+ 'INS' +'-'+ str(int(float(row.SVlen)))

	if 'Simple_repeat' in str(row.RepMasker_repClass):
		repSV.loc[idx, 'Simple_repeat_logic'] = 1
	else:
		repSV.loc[idx, 'Simple_repeat_logic'] = 0

	if 'Low_complexity' in str(row.RepMasker_repClass):
		repSV.loc[idx, 'Low_complexity_logic'] = 1
	else:
		repSV.loc[idx, 'Low_complexity_logic'] = 0

# differing between l1 vs l2
	if 'LINE/L1' in str(row.RepMasker_repClass):
		repSV.loc[idx, 'LINE1_logic'] = 1
	else:
		repSV.loc[idx, 'LINE1_logic'] = 0
	if 'LINE/L2' in str(row.RepMasker_repClass):
		repSV.loc[idx, 'LINE2_logic'] = 1
	else:
		repSV.loc[idx, 'LINE2_logic'] = 0
	if 'LINE/CR1' in str(row.RepMasker_repClass):
		repSV.loc[idx, 'CR1_logic'] = 1
	else:
		repSV.loc[idx, 'CR1_logic'] = 0

# differing between SINEs
	if 'SINE/Alu' in str(row.RepMasker_repClass):
		repSV.loc[idx, 'Alu_logic'] = 1
	else:
		repSV.loc[idx, 'Alu_logic'] = 0
	if 'SINE/MIR' in str(row.RepMasker_repClass):
		repSV.loc[idx, 'MIR_logic'] = 1
	else:
		repSV.loc[idx, 'MIR_logic'] = 0

# adding ltr erv elements
	if 'LTR/ERV' in str(row.RepMasker_repClass):
		repSV.loc[idx, 'ERV_logic'] = 1
	else:
		repSV.loc[idx, 'ERV_logic'] = 0
	if 'LTR/Gypsy' in str(row.RepMasker_repClass):
		repSV.loc[idx, 'Gypsy_logic'] = 1
	else:
		repSV.loc[idx, 'Gypsy_logic'] = 0

	if 'Satellite' in str(row.RepMasker_repClass):
		repSV.loc[idx, 'Satellite_logic'] = 1
	else:
		repSV.loc[idx, 'Satellite_logic'] = 0

# adding SVA elements
	if 'Retroposon/SVA' in str(row.RepMasker_repClass):
		repSV.loc[idx, 'SVA_logic'] = 1
	else:
		repSV.loc[idx, 'SVA_logic'] = 0


repSV = repSV[['ID', 'Simple_repeat_logic', 'Low_complexity_logic','LINE1_logic', 'LINE2_logic', 'CR1_logic', 'Alu_logic', 'MIR_logic','ERV_logic', 'Gypsy_logic', 'Satellite_logic', 'SVA_logic' ]]
repSV = repSV.replace(np. nan,0)
repSV = repSV.reset_index(drop=True)

print(repSV.head())
print('number of lines from repSV', len(repSV.ID.unique()))

# Repeat Masker - flanks, now its separate files for each annotation
repLINE = pd.read_csv(str(argv[9]), comment='#', sep='\t')
repLINE = repLINE[(repLINE.CHROM == chr) & (repLINE.SV_Type == type)]
for idx, row in repLINE.iterrows():
	if (row.SV_Type == 'deletion') | (row.SV_Type == 'DEL') :
		repLINE.loc[idx, 'ID'] = row.CHROM +'-'+ str(row.POS+1) +'-'+ 'DEL' +'-'+ str(int(float(row.SVlen)))
	elif (row.SV_Type == 'insertion') | (row.SV_Type == 'INS') :
		repLINE.loc[idx, 'ID'] = row.CHROM +'-'+ str(row.POS+1) +'-'+ 'INS' +'-'+ str(int(float(row.SVlen)))
repLINE = repLINE[['ID', 'LINE_sum']]
print(repLINE.head())
print('number of lines from repLINE', len(repLINE.ID.unique()))

repLOW = pd.read_csv(str(argv[10]), comment='#', sep='\t')
repLOW = repLOW[(repLOW.CHROM == chr) & (repLOW.SV_Type == type)]
for idx, row in repLOW.iterrows():
	if (row.SV_Type == 'deletion') | (row.SV_Type == 'DEL') :
		repLOW.loc[idx, 'ID'] = row.CHROM +'-'+ str(row.POS+1) +'-'+ 'DEL' +'-'+ str(int(float(row.SVlen)))
	elif (row.SV_Type == 'insertion') | (row.SV_Type == 'INS') :
		repLOW.loc[idx, 'ID'] = row.CHROM +'-'+ str(row.POS+1) +'-'+ 'INS' +'-'+ str(int(float(row.SVlen)))
repLOW = repLOW[['ID', 'Low_complexity_sum']]
print(repLOW.head())
print('number of lines from repLOW', len(repLOW.ID.unique()))

repLTR = pd.read_csv(str(argv[11]), comment='#', sep='\t')
repLTR = repLTR[(repLTR.CHROM == chr) & (repLTR.SV_Type == type)]
for idx, row in repLTR.iterrows():
	if (row.SV_Type == 'deletion') | (row.SV_Type == 'DEL') :
		repLTR.loc[idx, 'ID'] = row.CHROM +'-'+ str(row.POS+1) +'-'+ 'DEL' +'-'+ str(int(float(row.SVlen)))
	elif (row.SV_Type == 'insertion') | (row.SV_Type == 'INS') :
		repLTR.loc[idx, 'ID'] = row.CHROM +'-'+ str(row.POS+1) +'-'+ 'INS' +'-'+ str(int(float(row.SVlen)))
repLTR = repLTR[['ID', 'LTR_sum']]
print(repLTR.head())
print('number of lines from repLTR', len(repLTR.ID.unique()))

repSAT = pd.read_csv(str(argv[12]), comment='#', sep='\t')
repSAT = repSAT[(repSAT.CHROM == chr) & (repSAT.SV_Type == type)]
for idx, row in repSAT.iterrows():
	if (row.SV_Type == 'deletion') | (row.SV_Type == 'DEL') :
		repSAT.loc[idx, 'ID'] = row.CHROM +'-'+ str(row.POS+1) +'-'+ 'DEL' +'-'+ str(int(float(row.SVlen)))
	elif (row.SV_Type == 'insertion') | (row.SV_Type == 'INS') :
		repSAT.loc[idx, 'ID'] = row.CHROM +'-'+ str(row.POS+1) +'-'+ 'INS' +'-'+ str(int(float(row.SVlen)))
repSAT = repSAT[['ID', 'Satellite_sum']]
print(repSAT.head())
print('number of lines from rep', len(repSAT.ID.unique()))

repSIM = pd.read_csv(str(argv[13]), comment='#', sep='\t')
repSIM = repSIM[(repSIM.CHROM == chr) & (repSIM.SV_Type == type)]
for idx, row in repSIM.iterrows():
	if (row.SV_Type == 'deletion') | (row.SV_Type == 'DEL') :
		repSIM.loc[idx, 'ID'] = row.CHROM +'-'+ str(row.POS+1) +'-'+ 'DEL' +'-'+ str(int(float(row.SVlen)))
	elif (row.SV_Type == 'insertion') | (row.SV_Type == 'INS') :
		repSIM.loc[idx, 'ID'] = row.CHROM +'-'+ str(row.POS+1) +'-'+ 'INS' +'-'+ str(int(float(row.SVlen)))
repSIM = repSIM[['ID', 'Simple_repeat_sum']]
print(repSIM.head())
print('number of lines from rep', len(repSIM.ID.unique()))

repSINE = pd.read_csv(str(argv[14]), comment='#', sep='\t')
repSINE = repSINE[(repSINE.CHROM == chr) & (repSINE.SV_Type == type)]
for idx, row in repSINE.iterrows():
	if (row.SV_Type == 'deletion') | (row.SV_Type == 'DEL') :
		repSINE.loc[idx, 'ID'] = row.CHROM +'-'+ str(row.POS+1) +'-'+ 'DEL' +'-'+ str(int(float(row.SVlen)))
	elif (row.SV_Type == 'insertion') | (row.SV_Type == 'INS') :
		repSINE.loc[idx, 'ID'] = row.CHROM +'-'+ str(row.POS+1) +'-'+ 'INS' +'-'+ str(int(float(row.SVlen)))
repSINE = repSINE[['ID', 'SINE_sum']]
print(repSINE.head())
print('number of lines from rep', len(repSINE.ID.unique()))

nonbAPR = pd.read_csv(str(argv[15]), comment='#', sep='\t')
nonbAPR = nonbAPR[(nonbAPR.CHROM == chr) & (nonbAPR.SV_Type == type)]
for idx, row in nonbAPR.iterrows():
	if (row.SV_Type == 'deletion') | (row.SV_Type == 'DEL') :
		nonbAPR.loc[idx, 'ID'] = row.CHROM +'-'+ str(row.POS+1) +'-'+ 'DEL' +'-'+ str(int(float(row.SVlen)))
	elif (row.SV_Type == 'insertion') | (row.SV_Type == 'INS') :
		nonbAPR.loc[idx, 'ID'] = row.CHROM +'-'+ str(row.POS+1) +'-'+ 'INS' +'-'+ str(int(float(row.SVlen)))
nonbAPR = nonbAPR[['ID', 'Arep_sum']]
print(nonbAPR.head())
print('number of lines from nonb', len(nonbAPR.ID.unique()))

nonbDR = pd.read_csv(str(argv[16]), comment='#', sep='\t')
nonbDR = nonbDR[(nonbDR.CHROM == chr) & (nonbDR.SV_Type == type)]
for idx, row in nonbDR.iterrows():
	if (row.SV_Type == 'deletion') | (row.SV_Type == 'DEL') :
		nonbDR.loc[idx, 'ID'] = row.CHROM +'-'+ str(row.POS+1) +'-'+ 'DEL' +'-'+ str(int(float(row.SVlen)))
	elif (row.SV_Type == 'insertion') | (row.SV_Type == 'INS') :
		nonbDR.loc[idx, 'ID'] = row.CHROM +'-'+ str(row.POS+1) +'-'+ 'INS' +'-'+ str(int(float(row.SVlen)))
nonbDR = nonbDR[['ID', 'DR_sum']]
print(nonbDR.head())
print('number of lines from nonb', len(nonbDR.ID.unique()))

nonbGQ = pd.read_csv(str(argv[17]), comment='#', sep='\t')
nonbGQ = nonbGQ[(nonbGQ.CHROM == chr) & (nonbGQ.SV_Type == type)]
for idx, row in nonbGQ.iterrows():
	if (row.SV_Type == 'deletion') | (row.SV_Type == 'DEL') :
		nonbGQ.loc[idx, 'ID'] = row.CHROM +'-'+ str(row.POS+1) +'-'+ 'DEL' +'-'+ str(int(float(row.SVlen)))
	elif (row.SV_Type == 'insertion') | (row.SV_Type == 'INS') :
		nonbGQ.loc[idx, 'ID'] = row.CHROM +'-'+ str(row.POS+1) +'-'+ 'INS' +'-'+ str(int(float(row.SVlen)))
nonbGQ = nonbGQ[['ID', 'G4_sum']]
print(nonbGQ.head())
print('number of lines from nonb', len(nonbGQ.ID.unique()))

nonbMR = pd.read_csv(str(argv[18]), comment='#', sep='\t')
nonbMR = nonbMR[(nonbMR.CHROM == chr) & (nonbMR.SV_Type == type)]
for idx, row in nonbMR.iterrows():
	if (row.SV_Type == 'deletion') | (row.SV_Type == 'DEL') :
		nonbMR.loc[idx, 'ID'] = row.CHROM +'-'+ str(row.POS+1) +'-'+ 'DEL' +'-'+ str(int(float(row.SVlen)))
	elif (row.SV_Type == 'insertion') | (row.SV_Type == 'INS') :
		nonbMR.loc[idx, 'ID'] = row.CHROM +'-'+ str(row.POS+1) +'-'+ 'INS' +'-'+ str(int(float(row.SVlen)))
nonbMR = nonbMR[['ID', 'MR_sum']]
print(nonbMR.head())
print('number of lines from nonb', len(nonbMR.ID.unique()))

nonbSTR = pd.read_csv(str(argv[19]), comment='#', sep='\t')
nonbSTR = nonbSTR[(nonbSTR.CHROM == chr) & (nonbSTR.SV_Type == type)]
for idx, row in nonbSTR.iterrows():
	if (row.SV_Type == 'deletion') | (row.SV_Type == 'DEL') :
		nonbSTR.loc[idx, 'ID'] = row.CHROM +'-'+ str(row.POS+1) +'-'+ 'DEL' +'-'+ str(int(float(row.SVlen)))
	elif (row.SV_Type == 'insertion') | (row.SV_Type == 'INS') :
		nonbSTR.loc[idx, 'ID'] = row.CHROM +'-'+ str(row.POS+1) +'-'+ 'INS' +'-'+ str(int(float(row.SVlen)))
nonbSTR = nonbSTR[['ID', 'STR_sum']]
print(nonbSTR.head())
print('number of lines from nonb', len(nonbSTR.ID.unique()))

nonbZ = pd.read_csv(str(argv[20]), comment='#', sep='\t')
nonbZ = nonbZ[(nonbZ.CHROM == chr) & (nonbZ.SV_Type == type)]
for idx, row in nonbZ.iterrows():
	if (row.SV_Type == 'deletion') | (row.SV_Type == 'DEL') :
		nonbZ.loc[idx, 'ID'] = row.CHROM +'-'+ str(row.POS+1) +'-'+ 'DEL' +'-'+ str(int(float(row.SVlen)))
	elif (row.SV_Type == 'insertion') | (row.SV_Type == 'INS') :
		nonbZ.loc[idx, 'ID'] = row.CHROM +'-'+ str(row.POS+1) +'-'+ 'INS' +'-'+ str(int(float(row.SVlen)))
nonbZ = nonbZ[['ID', 'ZDNA_sum']]
print(nonbZ.head())
print('number of lines from nonb', len(nonbZ.ID.unique()))


# Chromoband of breakpoint
chromoband = pd.read_csv(str(argv[21]), comment='#', sep='\t')
chromoband = chromoband[(chromoband.CHROM == chr) & (chromoband.SV_Type == type)]
for idx, row in chromoband.iterrows():
	if (row.SV_Type == 'deletion') | (row.SV_Type == 'DEL') :
		chromoband.loc[idx, 'ID'] = row.CHROM +'-'+ str(row.POS+1) +'-'+ 'DEL' +'-'+ str(int(float(row.SVlen)))
	elif (row.SV_Type == 'insertion') | (row.SV_Type == 'INS') :
		chromoband.loc[idx, 'ID'] = row.CHROM +'-'+ str(row.POS+1) +'-'+ 'INS' +'-'+ str(int(float(row.SVlen)))
chromoband = chromoband[['ID', 'gieStain']]
chromoband_dum = pd.get_dummies(chromoband, columns = ['gieStain'])
chromoband_dum = chromoband_dum.reset_index(drop=True)
print(chromoband_dum.head())
print('number of lines from chromoband_dum', len(chromoband_dum.ID.unique()))

# Blast results  - variants and flanks
blast = pd.read_csv(str(argv[22]), comment='#', sep='\t')
print(blast.head())

blast = blast[(blast.CHROM == chr) & (blast.SV_Type == type)]
blast[['chr', 'pos', 'len', 'type']] = blast['ID'].str.split('_',  expand=True)
blast[["pos", "len"]] = blast[["pos", "len"]].apply(pd.to_numeric)

for idx, row in blast.iterrows():
	blast.loc[idx, 'prepost_Blastcoverage'] = (row.prepost_length)/2000
	blast.loc[idx, 'presv_Blastcoverage'] = (row.presv_qend - row.presv_qstart-row.presv_gapopen)/row.SVlen
	blast.loc[idx, 'postsv_Blastcoverage'] = (row.postsv_qend - row.postsv_qstart-row.postsv_gapopen)/row.SVlen

	if (row.SV_Type == 'deletion') | (row.SV_Type == 'DEL') :
		blast.loc[idx, 'ID'] = row.CHROM +'-'+ str(row.POS+1) +'-'+ 'DEL' +'-'+ str(int(float(row.SVlen)))
	elif (row.SV_Type == 'insertion') | (row.SV_Type == 'INS') :
		blast.loc[idx, 'ID'] = row.CHROM +'-'+ str(row.POS+1) +'-'+ 'INS' +'-'+ str(int(float(row.SVlen)))

blast = blast[['ID', 'prepost_pident', 'prepost_Blastcoverage', 'presv_pident', 'presv_Blastcoverage', 'postsv_pident', 'postsv_Blastcoverage']]
blast = blast.fillna(0)
blast = blast.reset_index(drop=True)

print(blast.head())
print('number of lines from blasting', len(blast.ID.unique()))

# Epigenetics results - depends on type of variant we're doing
epifeat = pd.read_csv(str(argv[23]), comment='#', sep='\t')
epifeat = epifeat[(epifeat.CHROM == chr) & (epifeat.SV_Type == type)]
for idx, row in epifeat.iterrows():
	if (row.SV_Type == 'deletion') | (row.SV_Type == 'DEL') :
		epifeat.loc[idx, 'ID'] = row.CHROM +'-'+ str(row.POS+1) +'-'+ 'DEL' +'-'+ str(int(float(row.SVlen)))
	elif (row.SV_Type == 'insertion') | (row.SV_Type == 'INS') :
		epifeat.loc[idx, 'ID'] = row.CHROM +'-'+ str(row.POS+1) +'-'+ 'INS' +'-'+ str(int(float(row.SVlen)))

epifeat = epifeat.loc[:, epifeat.columns.str.contains('ID|avg|slope|std')]
epifeat = epifeat.fillna(0)
epifeat = epifeat.reset_index(drop=True)

print(epifeat.head())
print('number of lines from epifeat', len(epifeat.ID.unique()))

# Adding DNA shape features
shape = pd.read_csv(str(argv[24]), sep='\t')
shape[['chr', 'pos', 'len', 'type']] = shape['ID'].str.split('_',  expand=True)
shape[["pos", "len"]] = shape[["pos", "len"]].apply(pd.to_numeric)
shape = shape[(shape.chr == chr) & (shape.type == type)]

for idx, row in shape.iterrows():
	if (row.type == 'deletion') | (row.type == 'DEL') :
		shape.loc[idx, 'ID'] = row.chr +'-'+ str(row.pos+1) +'-'+ 'DEL' +'-'+ str(int(float(row.len)))
	elif (row.type == 'insertion') | (row.type == 'INS') :
		shape.loc[idx, 'ID'] = row.chr +'-'+ str(row.pos+1) +'-'+ 'INS' +'-'+ str(int(float(row.len)))

shape = shape.loc[:, shape.columns.str.contains('ID|SV|PRE|POST')]
shape = shape.replace(np. nan,0)
shape = shape.reset_index(drop=True)
print(shape.head())
print('number of lines from shape', len(shape.ID.unique()))


# Repliseq value at breakpoint
repli = pd.read_csv(str(argv[25]), sep='\t')
repli = repli[(repli.CHROM == chr) & (repli.SV_Type == type)]
for idx, row in repli.iterrows():
	if (row.SV_Type == 'deletion') | (row.SV_Type == 'DEL') :
		repli.loc[idx, 'ID'] = row.CHROM +'-'+ str(row.POS+1) +'-'+ 'DEL' +'-'+ str(int(float(row.SVlen)))
	elif (row.SV_Type == 'insertion') | (row.SV_Type == 'INS') :
		repli.loc[idx, 'ID'] = row.CHROM +'-'+ str(row.POS+1) +'-'+ 'INS' +'-'+ str(int(float(row.SVlen)))
repli = repli[['ID', 'S50']]
print(repli.head())
print('number of lines from repli', len(repli.ID.unique()))

# Rloop flank data
rloop = pd.read_csv(str(argv[26]), sep='\t')
rloop = rloop.rename(columns={'unique_id': 'ID'})
rloop[['chr', 'pos', 'len', 'type']] = rloop['ID'].str.split('_',  expand=True)
rloop[["pos", "len"]] = rloop[["pos", "len"]].apply(pd.to_numeric)
rloop = rloop[(rloop.chr == chr) & (rloop.type == type)]

for idx, row in rloop.iterrows():
	if (row.type == 'deletion') | (row.type == 'DEL') :
		rloop.loc[idx, 'ID'] = row.chr +'-'+ str(row.pos+1) +'-'+ 'DEL' +'-'+ str(int(float(row.len)))
	elif (row.type == 'insertion') | (row.type == 'INS') :
		rloop.loc[idx, 'ID'] = row.chr +'-'+ str(row.pos+1) +'-'+ 'INS' +'-'+ str(int(float(row.len)))

rloop = rloop[['ID', 'RLoop']]
rloop = rloop.replace(np. nan,0)
rloop = rloop.reset_index(drop=True)
print(rloop.head())
print('number of lines from Rloop', len(rloop.ID.unique()))

if (type == 'insertion') | (type=='INS'):
	dataframes = [feat_filt, flankfeat, repSV, repLINE, repLOW, repLTR, repSAT, repSIM, repSINE, nonbAPR, nonbDR, nonbGQ, nonbMR, nonbSTR, nonbZ, chromoband_dum, blast, epifeat, shape, repli, rloop]
	results = dataframes[0]
	for df in dataframes[1:]:
		results = results.drop_duplicates(subset=['ID'])
		df = df.drop_duplicates(subset=['ID'])

		results = results.merge(df, on='ID', how='inner', suffixes=('', '_delme'))
		results = results[[c for c in results.columns if not c.endswith('_delme')]]

	results = results.replace(np.nan, 0)
	print(results.head())
	newrows = len(results)
	print('number of rows we got', newrows)

elif (type == 'deletion') | (type=='DEL'):
	epi = pd.read_csv(str(argv[27]), comment='#', sep='\t')
	epi = epi[(epi.CHROM == chr) & (epi.SV_Type == type)]
	for idx, row in epi.iterrows():
		if (row.SV_Type == 'deletion') | (row.SV_Type == 'DEL') :
			epi.loc[idx, 'ID'] = row.CHROM +'-'+ str(row.POS+1) +'-'+ 'DEL' +'-'+ str(int(float(row.SVlen)))
		elif (row.SV_Type == 'insertion') | (row.SV_Type == 'INS') :
			epi.loc[idx, 'ID'] = row.CHROM +'-'+ str(row.POS+1) +'-'+ 'INS' +'-'+ str(int(float(row.SVlen)))

	epi = epi.loc[:, epi.columns.str.contains('ID|avg|slope|std')]
	epi = epi.replace(np. nan,0)
	epi = epi.reset_index(drop=True)

# I don't think we need this!
#	epi = epi.drop(columns=['Sim_ID'])

	print(epi.head())
	print('number of lines from epi', len(epi.ID.unique()))

	dataframes = [feat_filt, flankfeat, repSV, repLINE, repLOW, repLTR, repSAT, repSIM, repSINE, nonbAPR, nonbDR, nonbGQ, nonbMR, nonbSTR, nonbZ, chromoband_dum, blast, epifeat, shape, repli, rloop, epi]
	results = dataframes[0]
	for df in dataframes[1:]:
		results = results.drop_duplicates(subset=['ID'])
		df = df.drop_duplicates(subset=['ID'])

		results = results.merge(df, on='ID', how='inner', suffixes=('', '_delme'))
		results = results[[c for c in results.columns if not c.endswith('_delme')]]

	results = results.replace(np.nan, 0)
	print(results.head())
	newrows = len(results)
	print('number of rows we got', newrows)

results.to_csv('/home/nboev/projects/def-sushant/nboev/preprocess/'+project+'/'+loc+'/merge_FeatMatrix/'+chr+'/' +filename+'.'+typer+ '.FeatMatrix.csv', sep='\t', index=False)

print('END TIME:', datetime.datetime.now(timezone('EST')))

