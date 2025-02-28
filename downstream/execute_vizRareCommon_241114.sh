#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=30000
#SBATCH --time=03:30:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nb.boev@mail.utoronto.ca
#SBATCH --job-name=viz
#SBATCH --output=viz.o
#SBATCH --error=viz.e
#SBATCH --qos=privilege

# Loading in virtual environment and modules
module --force purge ; module load StdEnv/2020
source $HOME/virtenv/bin/activate
module load bedtools

# The input: 
# 1: The location of the original SV .csv (holds SVlength AND genotype information)
# 2: The location of the igsr_samples.tsv, which holds the ancestral information for each individual

# Making a directory to hold outputs
mkdir /home/nboev/projects/def-sushant/nboev/analysis/20220422_3202_phased_SNV_INDEL_SV_bychrom/SVTrue_typedeletion_resTrue/IDmechsvSIM/20240625/splitsChromo/bedRareInd/popRates

# will need to delete these two!!!
#$1 = '/home/nboev/projects/def-sushant/nboev/preprocess/20220422_3202_phased_SNV_INDEL_SV_bychrom/SVTrue_typedeletion_resTrue/SVlen/SVTrue_typedeletion_resTrue.csv'
#$2 = /home/nboev/projects/def-sushant/nboev/data/phased_SNV_INDEL_SV_20220422_3202/igsr_samples.tsv'

# Pushing the python script, generating_bedRaresInds.py. This script generates individual-level SVs for both clusters (clusterI and clusterII), among rare SVs. 
python generating_bedRaresInds.py \
$1 \
$2 \
>> ./SVTrue_typedeletion_resTrue/analysis_vizRaresAll_241118.py.txt

# We loop through all the populuations in 1KG (note: We exluded the individual who is a known admixture)
# For each population, we concatenate the population rates. 
# Then, for the two bed files (ClusterI and clusterII), we find the individual-level bed files.
# Finally, we use bedtools to intersect the appropriate recombination rate bed file to the SV bed file. 
for i in ACB ASW BEB CDX CEU CHB CHS CLM ESN FIN GBR GIH GWD IBS ITU JPT KHV LWK MSL MXL PEL PJL PUR STU TSI YRI;
	do
	cat /home/nboev/projects/def-sushant/nboev/data/popRecombRate/hg38/${i}/${i}_recombination_map_hg38_chr_{1..22}.bed |uniq \
> /home/nboev/projects/def-sushant/nboev/data/popRecombRate/hg38/${i}/${i}_recombination_map_hg38.bed

	for k in clusterI clusterII ;
		do
		for filename in /home/nboev/projects/def-sushant/nboev/analysis/20220422_3202_phased_SNV_INDEL_SV_bychrom/SVTrue_typedeletion_resTrue/IDmechsvSIM/20240625/splitsChromo/bedRareInd/${i}.*.HOMOpreds.${k}.bed ;
			do
			file=$(basename "$filename" .bed)

			bedtools intersect -wa -wb \
-a ${filename} \
-b /home/nboev/projects/def-sushant/nboev/data/popRecombRate/hg38/${i}/${i}_recombination_map_hg38.bed \
> /home/nboev/projects/def-sushant/nboev/analysis/20220422_3202_phased_SNV_INDEL_SV_bychrom/SVTrue_typedeletion_resTrue/IDmechsvSIM/20240625/splitsChromo/bedRareInd/popRates/${file}.bed

			done
		done
	done

# ** work on this part!!
python analysis_RarepopRecomb.py \
>> ./SVTrue_typedeletion_resTrue/analysis_RarepopRecomb_241127.py.txt
