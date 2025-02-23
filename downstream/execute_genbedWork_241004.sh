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

source $HOME/virtenv/bin/activate
module load samtools
module load bedtools

# Comparing predictions (testing) to other genomic features
mkdir /home/nboev/projects/def-sushant/nboev/analysis/$1/$2/IDmechsvSIM/20240625/splitsChromo/bedWork

for j in insertion deletion ;
	do
	for i in X Y 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 ;
		do
		python generating_bedWork_241004.py \
/home/nboev/projects/def-sushant/nboev/analysis/$1/$2/IDmechsvSIM/20240625/splitsChromo/models/results/chr${i}.${j}.HOMOpreds.tsv \
$1 $2 \
chr${i} ${j} \
>> /home/nboev/projects/def-sushant/nboev/analysis/$1/$2/generating_bedWork_241004.py.txt
		done
	done


#How to run:
# sbatch execute_genbedWork_241004.sh HGSVC2_v2_integrated_callset sv
# sbatch execute_genbedWork_241004.sh 20220422_3202_phased_SNV_INDEL_SV_bychrom SVTrue_typedeletion_resTrue
# sbatch execute_genbedWork_241004.sh ONT_1000Gs_100FREEZE sv

