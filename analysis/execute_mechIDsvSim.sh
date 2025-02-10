#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=10G
#SBATCH --time=02:30:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nb.boev@mail.utoronto.ca
#SBATCH --job-name=mechID
#SBATCH --output=mechID.o
#SBATCH --error=mechID.e
#SBATCH --qos=privileged

# Loading in virtual environment and modules
module --force purge ; module load StdEnv/2020
source $HOME/virtenv/bin/activate
module load gcc/9.3.0 r/4.0.2
module load scipy-stack

# The input: 
# 1: Project name and location
# 2: The z-score subdir annotated file to annotate

# Creating directories to hold the labelled SVs
mkdir /home/nboev/projects/def-sushant/nboev/analysis/$1
mkdir /home/nboev/projects/def-sushant/nboev/analysis/$1/$2
mkdir /home/nboev/projects/def-sushant/nboev/analysis/$1/$2/IDmechsvSIM
mkdir /home/nboev/projects/def-sushant/nboev/analysis/$1/$2/IDmechsvSIM/20240625

# We loop through the chromosomes and the SV types to send the z-score annotated file to the python script, analysis_mechIDsvSim_240625.py which will annotate the SVs by their local homology
for i in X Y 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 ;
	do
	for j in deletion insertion;
		do
		python analysis_mechIDsvSim_240625.py \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/zscore/chr${i}.20240625.chr${i}.${j}_FeatMatrixzscore.csv \
chr${i} ${j} \
$1 $2 \
>> /home/nboev/projects/def-sushant/nboev/analysis/$1/$2/analysis_mechIDsvSim_240625.py.txt
		done
	done

# example:
# sbatch execute_mechIDsvSim_240625.sh ONT_1000Gs_100FREEZE sv
