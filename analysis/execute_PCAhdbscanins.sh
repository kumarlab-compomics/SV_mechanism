#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=30000
#SBATCH --time=09:00:00
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

# In this script, we looped through different parameters for HDBSCAN clustering. This includes:
# 1. Building a model for insertions and deletions separately
# 2. Reducing to different n number of dimensions with PCA (2 5 10 20)
# 3. Requiring different minimum cluster sizes (300 500 700 900 1100 1300)
# 4. Testing different distance metrics (braycurtis manhattan canberra euclidean chebyshev)
# Using the above parameters, we push the python script, analysis_hdbscanID.py which will calculate the values we'll use for selecting the optimal parameters
# Recall: This is exclusively done on the HGSVC2 SVs, whereby in the python script we remove the testing chromosomes

for a in insertion deletion ;
	do
	for b in 2 5 10 20 ;
		do
		for c in 300 500 700 900 1100 1300 ;
			do
			for d in braycurtis manhattan canberra euclidean chebyshev ;
				do
				python analysis_hdbscanID.py \
/home/nboev/projects/def-sushant/nboev/analysis/HGSVC2_v2_integrated_callset/sv/IDmechsvSIM/20240625/ \
${a}_ID0625.tsv \
240625/ \
${b} \
${c} \
${d} \
>> ./sv/IDmechsvSIM/20240625/splitsChromo/models/${a}/analysis_hdbscanID_240703.${a}.${b}.${c}.${d}.py.txt
				done
			done
		done
	done

# how to run: order : sv type, number of pcs, min size of cluster, distance metric
# sbatch execute_PCAhdbscan_240627.sh
