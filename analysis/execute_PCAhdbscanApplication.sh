#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=30000
#SBATCH --time=00:15:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nb.boev@mail.utoronto.ca
#SBATCH --job-name=hdbscan_new
#SBATCH --output=hdbscan_new.o
#SBATCH --error=hdbscan_new.e
#SBATCH --qos=privileged

# Loading in virtual environment and modules
module --force purge ; module load StdEnv/2020
source $HOME/virtenv/bin/activate
module load gcc/9.3.0 r/4.0.2
module load scipy-stack

# In this script, we are applying the trained scaler, PCA and KNN models to other datasets.
# We generated an insertion and deletion models separately, therefore we point to the appropriate SVs and saved models for their use. 
# Recall, the models are trained on HGSVC2 SVs

# The input: 
# 1: Project name and location
# 2: The svs to annotate (ie. sv)

# Creating directories to save the outputs
mkdir /home/nboev/projects/def-sushant/nboev/analysis/$1/$2/IDmechsvSIM/20240625/splitsChromo
mkdir /home/nboev/projects/def-sushant/nboev/analysis/$1/$2/IDmechsvSIM/20240625/splitsChromo/models
mkdir /home/nboev/projects/def-sushant/nboev/analysis/$1/$2/IDmechsvSIM/20240625/splitsChromo/models/results

# Looping through ins/dels, across each chromosome separately
for k in insertion deletion ;
	do
	for i in X Y 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 ;
		do
		python analysis_hdbscanIDApplicationHOMO.py \
/home/nboev/projects/def-sushant/nboev/analysis/$1/$2/IDmechsvSIM/20240625/chr${i}.${k}_ID0625.tsv \
$1 $2 \
${k} \
chr${i} \
/home/nboev/projects/def-sushant/nboev/analysis/HGSVC2_v2_integrated_callset/sv/IDmechsvSIM/20240625/splitsChromo/models/SCALERtraining${k}.HOMO.sav \
/home/nboev/projects/def-sushant/nboev/analysis/HGSVC2_v2_integrated_callset/sv/IDmechsvSIM/20240625/splitsChromo/models/PCAtraining${k}.HOMO.sav \
/home/nboev/projects/def-sushant/nboev/analysis/HGSVC2_v2_integrated_callset/sv/IDmechsvSIM/20240625/splitsChromo/models/KNNtraining${k}.HOMO.sav \
>> /home/nboev/projects/def-sushant/nboev/analysis/$1/$2/analysis_hdbscanIDApplicationHOMO_240708.py.txt
		done
	done

# How to run:
# sbatch execute_PCAhdbscanApplication_240709.sh ONT_1000Gs_100FREEZE sv
