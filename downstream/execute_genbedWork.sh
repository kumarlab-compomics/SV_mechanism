#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=30000
#SBATCH --time=00:05:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nb.boev@mail.utoronto.ca
#SBATCH --job-name=bed
#SBATCH --output=bed.o
#SBATCH --error=bed.e
#SBATCH --qos=privilege

# Loading in virtual environment and modules
source $HOME/virtenv/bin/activate
module load samtools
module load bedtools

# The input: 
# 1: Project name and location
# 2: The type of SVs to annotate (ie. sv or svSim). This is a subdirectory that holds the vcf/csv file
# Note: We need to point to csvs which hold the homology based results. We are generating bed files to represent these strata by SV type

# Creating a directory to store the bed files
mkdir /home/nboev/projects/def-sushant/nboev/analysis/$1/$2/IDmechsvSIM/20240625/splitsChromo/bedWork

# We loop through SV types and each chromosome to push the generating_bedWork.py script
for j in insertion deletion ;
	do
	for i in X Y 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 ;
		do
		python generating_bedWork.py \
/home/nboev/projects/def-sushant/nboev/analysis/$1/$2/IDmechsvSIM/20240625/splitsChromo/models/results/chr${i}.${j}.HOMOpreds.tsv \
$1 $2 \
chr${i} ${j} \
>> /home/nboev/projects/def-sushant/nboev/analysis/$1/$2/generating_bedWork_241004.py.txt
		done
	done


#How to run:
# sbatch execute_genbedWork_241004.sh HGSVC2_v2_integrated_callset sv
