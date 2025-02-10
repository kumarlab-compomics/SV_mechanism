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

module --force purge ; module load StdEnv/2020
source $HOME/virtenv/bin/activate
module load gcc/9.3.0 r/4.0.2
module load scipy-stack

# 2024/06/28
# we're going to start to do experiments in order to create our pseudo labelled data
# were trying to optimize the features

# Different version where we're trying to optimize for nhej instead of hr + using combined clusters

for a in insertion ;
	do
	for b in 2 5 10 20 ;
		do
		for c in 300 500 700 900 1100 1300 ;
			do
			for d in braycurtis manhattan canberra euclidean chebyshev ;
				do
				python analysis_hdbscanID_240703.py \
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
# sbatch execute_PCAhdbscan_240627.sh deletion 2 500 braycurtis
# sbatch execute_PCAhdbscan_240627.sh deletion 5 500 braycurtis
# sbatch execute_PCAhdbscan_240627.sh deletion 10 500 braycurtis
# sbatch execute_PCAhdbscan_240627.sh deletion 20 500 braycurtis

# sbatch execute_PCAhdbscan_240627.sh deletion 2 100 braycurtis
# sbatch execute_PCAhdbscan_240627.sh deletion 5 100 braycurtis
# sbatch execute_PCAhdbscan_240627.sh deletion 10 100 braycurtis
# sbatch execute_PCAhdbscan_240627.sh deletion 20 100 braycurtis

# sbatch execute_PCAhdbscan_240627.sh deletion 2 1000 braycurtis
# sbatch execute_PCAhdbscan_240627.sh deletion 5 1000 braycurtis
# sbatch execute_PCAhdbscan_240627.sh deletion 10 1000 braycurtis
# sbatch execute_PCAhdbscan_240627.sh deletion 20 1000 braycurtis


# other distance metrics
# sbatch execute_PCAhdbscan_240627.sh deletion 2 500 euclidean
# sbatch execute_PCAhdbscan_240627.sh deletion 2 500 chebyshev
# sbatch execute_PCAhdbscan_240627.sh deletion 2 500 manhattan
# sbatch execute_PCAhdbscan_240627.sh deletion 2 500 minkowski
# sbatch execute_PCAhdbscan_240627.sh deletion 2 500 jaccard
# sbatch execute_PCAhdbscan_240627.sh deletion 2 500 canberra

# sbatch execute_PCAhdbscan_240627.sh insertion 2 500 braycurtis

# didnt work...
# sbatch execute_PCAhdbscan_240627.sh deletion 2 500 mahalanobis
