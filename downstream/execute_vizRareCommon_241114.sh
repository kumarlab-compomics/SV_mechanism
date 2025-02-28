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

module --force purge ; module load StdEnv/2020
source $HOME/virtenv/bin/activate
module load bedtools

# This script pulls out the common, rares and denovos amonng the 602 trios but removes all duplicates
# We also care about inheritence patterns in this case (ie. not all 3000 inds are represented)
#python analysis_vizRareCommon_241114.py \
#>> ./SVTrue_typedeletion_resTrue/analysis_vizRareCommon_241114.py.txt


# In this case, we've used all the 3000 inds and pulled out the rare vars, which we will use to calculate % Fuzzys per subpops
#mkdir /home/nboev/projects/def-sushant/nboev/analysis/20220422_3202_phased_SNV_INDEL_SV_bychrom/SVTrue_typedeletion_resTrue/IDmechsvSIM/20240625/splitsChromo/bedRareInd
#mkdir /home/nboev/projects/def-sushant/nboev/analysis/20220422_3202_phased_SNV_INDEL_SV_bychrom/SVTrue_typedeletion_resTrue/IDmechsvSIM/20240625/splitsChromo/bedRareInd/nonAFR_high
#mkdir /home/nboev/projects/def-sushant/nboev/analysis/20220422_3202_phased_SNV_INDEL_SV_bychrom/SVTrue_typedeletion_resTrue/IDmechsvSIM/20240625/splitsChromo/bedRareInd/nonAFR_low
#mkdir /home/nboev/projects/def-sushant/nboev/analysis/20220422_3202_phased_SNV_INDEL_SV_bychrom/SVTrue_typedeletion_resTrue/IDmechsvSIM/20240625/splitsChromo/bedRareInd/AFR_high
#mkdir /home/nboev/projects/def-sushant/nboev/analysis/20220422_3202_phased_SNV_INDEL_SV_bychrom/SVTrue_typedeletion_resTrue/IDmechsvSIM/20240625/splitsChromo/bedRareInd/AFR_low
mkdir /home/nboev/projects/def-sushant/nboev/analysis/20220422_3202_phased_SNV_INDEL_SV_bychrom/SVTrue_typedeletion_resTrue/IDmechsvSIM/20240625/splitsChromo/bedRareInd/popRates

#python analysis_vizRaresAll_241118.py \
#>> ./SVTrue_typedeletion_resTrue/analysis_vizRaresAll_241118.py.txt

# Usisng all the rare ones, we're now going to see which genes have been hit
#for k in fuzzy0 fuzzy1 ;
#	do
#	for j in AFR_high AFR_low nonAFR_high nonAFR_low ;
#		do
#		for filename in /home/nboev/projects/def-sushant/nboev/analysis/20220422_3202_phased_SNV_INDEL_SV_bychrom/SVTrue_typedeletion_resTrue/IDmechsvSIM/20240625/splitsChromo/bedRareInd/${j}/*bed
#			do
#			file=$(basename "$filename")
#			hg_string=$(echo "$file" | awk -F. '{print $2}')
#
#			bedtools intersect -wa \
#-a /home/nboev/projects/def-sushant/nboev/analysis/20220422_3202_phased_SNV_INDEL_SV_bychrom/SVTrue_typedeletion_resTrue/IDmechsvSIM/20240625/splitsChromo/bedRareInd/${j}/${j}.${hg_string}.HOMOpreds.${k}.bed \
#-b /home/nboev/projects/def-sushant/nboev/data/RefSeq/hg38.ncbiRefSeq_onlyCDS.bed | \
#awk '!seen[$0]++' \
#> /home/nboev/projects/def-sushant/nboev/analysis/20220422_3202_phased_SNV_INDEL_SV_bychrom/SVTrue_typedeletion_resTrue/IDmechsvSIM/20240625/splitsChromo/bedRareInd/${j}/${j}.${hg_string}.HOMOpreds.${k}clean.bed

#			bedtools intersect -wa \
#-a /home/nboev/projects/def-sushant/nboev/data/RefSeq/hg38.ncbiRefSeq_onlyCDS.bed \
#-b /home/nboev/projects/def-sushant/nboev/analysis/20220422_3202_phased_SNV_INDEL_SV_bychrom/SVTrue_typedeletion_resTrue/IDmechsvSIM/20240625/splitsChromo/bedRareInd/${j}/${j}.${hg_string}.HOMOpreds.${k}clean.bed | \
#awk '{ print $NF }' | sort | uniq \
#> /home/nboev/projects/def-sushant/nboev/analysis/20220422_3202_phased_SNV_INDEL_SV_bychrom/SVTrue_typedeletion_resTrue/IDmechsvSIM/20240625/splitsChromo/bedRareInd/${j}/${j}.${hg_string}.HOMOpreds.${k}cleanGENES.bed

#			rm /home/nboev/projects/def-sushant/nboev/analysis/20220422_3202_phased_SNV_INDEL_SV_bychrom/SVTrue_typedeletion_resTrue/IDmechsvSIM/20240625/splitsChromo/bedRareInd/${j}/${j}.${hg_string}.HOMOpreds.${k}clean.bed
#			done
#		done
#	done

# we need to figure out the ind level recombination rates:
for i in ACB ASW BEB CDX CEU CHB CHS CLM ESN FIN GBR GIH GWD IBS ITU JPT KHV LWK MSL MXL PEL PJL PUR STU TSI YRI;
	do
	cat /home/nboev/projects/def-sushant/nboev/data/popRecombRate/hg38/${i}/${i}_recombination_map_hg38_chr_{1..22}.bed |uniq \
> /home/nboev/projects/def-sushant/nboev/data/popRecombRate/hg38/${i}/${i}_recombination_map_hg38.bed

	for k in fuzzy0 fuzzy1 ;
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


python analysis_RarepopRecomb_241127.py \
>> ./SVTrue_typedeletion_resTrue/analysis_RarepopRecomb_241127.py.txt
