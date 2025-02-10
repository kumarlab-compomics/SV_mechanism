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

module --force purge ; module load StdEnv/2020
source $HOME/virtenv/bin/activate
module load gcc/9.3.0 r/4.0.2
module load scipy-stack

mkdir /home/nboev/projects/def-sushant/nboev/analysis/$1
mkdir /home/nboev/projects/def-sushant/nboev/analysis/$1/$2
mkdir /home/nboev/projects/def-sushant/nboev/analysis/$1/$2/IDmechsvSIM
mkdir /home/nboev/projects/def-sushant/nboev/analysis/$1/$2/IDmechsvSIM/20240625

# Running the sex chromosomes
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


# how to run:
# sbatch execute_mechIDsvSim_240625.sh HGSVC2_v2_integrated_callset sv

# sbatch execute_mechIDsvSim_240625.sh 20220422_3202_phased_SNV_INDEL_SV_bychrom SVTrue_typedeletion_resTrue

# Running on CPTAC val = 2024/07/19
# sbatch execute_mechIDsvSim_240625.sh CPTAC3 bf584e4c-f981-473e-b40a-c73d7f3699d2.wgs.sanger_raw_pindel.raw_somatic_mutation
# sbatch execute_mechIDsvSim_240625.sh CPTAC3 fcc8a59b-19b0-4b78-9ca1-f5b9e580f183.wgs.sanger_raw_pindel.raw_somatic_mutation


# Running on the Belyeu2021 dataset = 2024/07/29
# sbatch execute_mechIDsvSim_240625.sh Belyeu2021 sv

# Running on the HGSVC3 dataset = 2024/08/20
# sbatch execute_mechIDsvSim_240625.sh allHGSVC3_20240415_Freeze4_GRCh38 sv

# Running on ONT dataset = 2024/09/11
# sbatch execute_mechIDsvSim_240625.sh ONT_1000Gs_100FREEZE sv

# Running on DepMap = 2024/10/21
# sbatch execute_mechIDsvSim_240625.sh DepMap CA922
# sbatch execute_mechIDsvSim_240625.sh DepMap DOTC24510
# sbatch execute_mechIDsvSim_240625.sh DepMap FU97
# sbatch execute_mechIDsvSim_240625.sh DepMap HCC1428
# sbatch execute_mechIDsvSim_240625.sh DepMap HSKTC
# sbatch execute_mechIDsvSim_240625.sh DepMap KMCH1
# sbatch execute_mechIDsvSim_240625.sh DepMap KP363T
# sbatch execute_mechIDsvSim_240625.sh DepMap KYSE30
# sbatch execute_mechIDsvSim_240625.sh DepMap KYSE410
# sbatch execute_mechIDsvSim_240625.sh DepMap KYSE450
# sbatch execute_mechIDsvSim_240625.sh DepMap NCIH748
# sbatch execute_mechIDsvSim_240625.sh DepMap SH10TC
# sbatch execute_mechIDsvSim_240625.sh DepMap SNU119
# sbatch execute_mechIDsvSim_240625.sh DepMap SNU251
# sbatch execute_mechIDsvSim_240625.sh DepMap SUDHL1
# sbatch execute_mechIDsvSim_240625.sh DepMap UPCISCC074
# sbatch execute_mechIDsvSim_240625.sh DepMap VMRCMELG


