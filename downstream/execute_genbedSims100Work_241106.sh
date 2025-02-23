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

source $HOME/virtenv/bin/activate
module load samtools
module load bedtools

# Comparing predictions (testing) to other genomic features

mkdir /home/nboev/projects/def-sushant/nboev/analysis/$1/$2/IDmechsvSIM/20240625/splitsChromo/bedWork100s

for k in insertion deletion ;
	do
	for j in HR SSAaEJ NHEJ ;
		do
		for i in X Y 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 ;
			do
			python generating_bedSim100sWork_241106.py \
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
# sbatch execute_genbedSims100Work_241106.sh 20220422_3202_phased_SNV_INDEL_SV_bychrom SVTrue_typedeletion_resTrue
# sbatch execute_genbedSims100Work_241106.sh ONT_1000Gs_100FREEZE sv



