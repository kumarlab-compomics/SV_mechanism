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

module --force purge ; module load StdEnv/2020
source $HOME/virtenv/bin/activate
module load gcc/9.3.0 r/4.0.2
module load scipy-stack

# Using the optimized parameters = 2024/07/09 (see the progress presentation on 2024/07/05)

mkdir /home/nboev/projects/def-sushant/nboev/analysis/$1/$2/IDmechsvSIM/20240625/splitsChromo
mkdir /home/nboev/projects/def-sushant/nboev/analysis/$1/$2/IDmechsvSIM/20240625/splitsChromo/models
mkdir /home/nboev/projects/def-sushant/nboev/analysis/$1/$2/IDmechsvSIM/20240625/splitsChromo/models/results

# For insertions: Do for all chromosomes!! + save the outputs!!!

for k in insertion deletion ;
	do
#	for i in X Y 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 ;
	for i in 4 ;
		do
		python analysis_hdbscanIDApplicationHOMO_240708.py \
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
# execute_PCAhdbscanApplication_240709.sh HGSVC2_v2_integrated_callset sv

# sbatch execute_PCAhdbscanApplication_240709.sh 20220422_3202_phased_SNV_INDEL_SV_bychrom SVTrue_typedeletion_resTrue


# Running on the CPTAC pilot val = 2024/07/19
# sbatch execute_PCAhdbscanApplication_240709.sh CPTAC3 bf584e4c-f981-473e-b40a-c73d7f3699d2.wgs.sanger_raw_pindel.raw_somatic_mutation
# sbatch execute_PCAhdbscanApplication_240709.sh CPTAC3 fcc8a59b-19b0-4b78-9ca1-f5b9e580f183.wgs.sanger_raw_pindel.raw_somatic_mutation


# Running on the Belyeu2021 set = 2024/0729
# sbatch execute_PCAhdbscanApplication_240709.sh Belyeu2021 sv

# Running on HGSVC3 = 2024/08/21
# execute_PCAhdbscanApplication_240709.sh allHGSVC3_20240415_Freeze4_GRCh38 sv


# Running on ONT = 2024/09/11
# sbatch execute_PCAhdbscanApplication_240709.sh ONT_1000Gs_100FREEZE sv
