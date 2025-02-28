#!/bin/bash
#SBATCH --nodes=2
#SBATCH --mem=60000
#SBATCH --time=04:30:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nb.boev@mail.utoronto.ca
#SBATCH --job-name=bedWork
#SBATCH --output=bedWork.o
#SBATCH --error=bedWork.e
#SBATCH --qos=privilege

# Loading in virtual environment and modules
source $HOME/virtenv/bin/activate
module load samtools
module load bedtools

# The input: 
# 1: Project name/ location
# 2: sv
# Note: We need to point to csvs which hold the homology based results along with the z-score merged files. We are generating bed files to represent these strata by SV type

# Creating a directory to store the bed files
mkdir /home/nboev/projects/def-sushant/nboev/analysis/$1/$2/IDmechsvSIM/20240625/splitsChromo/bedWork100s

# We loop through SV types, homology-based label, and each chromosome to push the generating_bedSim100sWork.py script. 
# We then concatenate these together and save the larger file (ie. regardless of chromosome).
for k in insertion deletion ;
	do
	for j in HLH ILH NLH ;
		do
		for i in X Y 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 ;
			do
			python generating_bedSim100sWork.py \
/home/nboev/projects/def-sushant/nboev/analysis/$1/$2/IDmechsvSIM/20240625/splitsChromo/models/results/chr${i}.${k}.HOMOpreds.tsv \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/zscore/chr${i}.20240625.chr${i}.${k}_FeatMatrixzscore.csv \
$1 $2 \
chr${i} ${k} \
${j} \
>> /home/nboev/projects/def-sushant/nboev/analysis/$1/$2/generating_bedSim100sWork_241106.py.txt

			for m in {0..100} ;
				do
				cat /home/nboev/projects/def-sushant/nboev/analysis/$1/$2/IDmechsvSIM/20240625/splitsChromo/bedWork100s/chr${i}.${k}.HOMOpreds.${j}SIMS${m}.bed \
>> /home/nboev/projects/def-sushant/nboev/analysis/$1/$2/IDmechsvSIM/20240625/splitsChromo/bedWork100s/${k}.HOMOpreds.${j}SIMS${m}.bed

				rm /home/nboev/projects/def-sushant/nboev/analysis/$1/$2/IDmechsvSIM/20240625/splitsChromo/bedWork100s/chr${i}.${k}.HOMOpreds.${j}SIMS${m}.bed
				done
			done
		done
	done


#How to run:
# sbatch execute_genbedSims100Work_241106.sh HGSVC2_v2_integrated_callset sv


