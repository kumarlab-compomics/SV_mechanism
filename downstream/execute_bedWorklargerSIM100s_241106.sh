#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=30000
#SBATCH --time=01:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nb.boev@mail.utoronto.ca
#SBATCH --job-name=bedFisher
#SBATCH --output=bedFisher.o
#SBATCH --error=bedFisher.e
#SBATCH --qos=privilege

source $HOME/virtenv/bin/activate
module load samtools
module load bedtools


#for j in insertion ;
for j in insertion deletion ;
       do
	for i in HR SSAaEJ NHEJ ;
		do
		mkdir /home/nboev/projects/def-sushant/nboev/analysis/$3/$4/IDmechsvSIM/20240625/splitsChromo/bedWork100s/$2

		bedtools intersect -wa \
-a /home/nboev/projects/def-sushant/nboev/analysis/$3/$4/IDmechsvSIM/20240625/splitsChromo/bedWork/${j}.HOMO.${i}.bed \
-b $1 | \
awk '!seen[$0]++' \
> /home/nboev/projects/def-sushant/nboev/analysis/$3/$4/IDmechsvSIM/20240625/splitsChromo/bedWork100s/$2/${j}.HOMOpreds.${i}.bed_$2clean.bed

		for m in {0..100} ;
			do
			bedtools intersect -wa \
-a /home/nboev/projects/def-sushant/nboev/analysis/$3/$4/IDmechsvSIM/20240625/splitsChromo/bedWork100s/${j}.HOMOpreds.${i}SIMS${m}.bed \
-b $1 | \
awk '!seen[$0]++' \
> /home/nboev/projects/def-sushant/nboev/analysis/$3/$4/IDmechsvSIM/20240625/splitsChromo/bedWork100s/$2/${j}.HOMOpreds.${i}SIMS${m}.bed_$2clean.bed

			python analysis_hdbscanknnEnrichment.py \
/home/nboev/projects/def-sushant/nboev/analysis/$3/$4/IDmechsvSIM/20240625/splitsChromo/bedWork100s/$2/${j}.HOMOpreds.${i}.bed_$2clean.bed \
/home/nboev/projects/def-sushant/nboev/analysis/$3/$4/IDmechsvSIM/20240625/splitsChromo/bedWork100s/$2/${j}.HOMOpreds.${i}SIMS${m}.bed_$2clean.bed \
/home/nboev/projects/def-sushant/nboev/analysis/$3/$4/IDmechsvSIM/20240625/splitsChromo/bedWork/${j}.HOMO.${i}.bed \
/home/nboev/projects/def-sushant/nboev/analysis/$3/$4/IDmechsvSIM/20240625/splitsChromo/bedWork100s/${j}.HOMOpreds.${i}SIMS${m}.bed \
chr ${j} realvssimWork${i} \
$3 $2 \
>> /home/nboev/projects/def-sushant/nboev/analysis/$3/$4/analysis_hdbscanknnEnrichment.py.txt

			rm /home/nboev/projects/def-sushant/nboev/analysis/$3/$4/IDmechsvSIM/20240625/splitsChromo/bedWork100s/$2/${j}.HOMOpreds.${i}SIMS${m}.bed_CDSclean.bed

			done
		rm /home/nboev/projects/def-sushant/nboev/analysis/$3/$4/IDmechsvSIM/20240625/splitsChromo/bedWork100s/$2/${j}.HOMOpreds.${i}.bed_$2clean.bed
		done
	done

# sbatch execute_bedWorklargerSIM100s_241106.sh /home/nboev/projects/def-sushant/nboev/data/RefSeq/hg38.ncbiRefSeq_onlyCDS.bed CDS HGSVC2_v2_integrated_callset sv
# sbatch execute_bedWorklargerSIM100s_241106.sh /home/nboev/projects/def-sushant/nboev/data/SCREEN/GRCh38-cCREs_dELS.bed ALLdELS HGSVC2_v2_integrated_callset sv
# sbatch execute_bedWorklargerSIM100s_241106.sh /home/nboev/projects/def-sushant/nboev/data/SCREEN/GRCh38-cCREs_pELS.bed ALLpELS HGSVC2_v2_integrated_callset sv
# sbatch execute_bedWorklargerSIM100s_241106.sh /home/nboev/projects/def-sushant/nboev/data/FragileSites/fragile_site_bed/fragile_site.bed fragile HGSVC2_v2_integrated_callset sv
# sbatch execute_bedWorklargerSIM100s_241106.sh /home/nboev/projects/def-sushant/nboev/data/SegmentalDups/genomicSuperDups.txt.min0.95.bed segdup HGSVC2_v2_integrated_callset sv

# sbatch execute_bedWorklargerSIM100s_241106.sh /home/nboev/projects/def-sushant/nboev/data/RefSeq/hg38.ncbiRefSeq_onlyCDS.bed CDS ONT_1000Gs_100FREEZE sv


# sbatch execute_bedWorklargerSIM100s_241106.sh /home/nboev/projects/def-sushant/nboev/data/RefSeq/hg38.ncbiRefSeq_onlyCDS.bed CDS 20220422_3202_phased_SNV_INDEL_SV_bychrom SVTrue_typedeletion_resTrue


# sbatch execute_bedWorklargerSIM100s_241106.sh /home/nboev/projects/def-sushant/nboev/data/RefSeq/hg38.ncbiRefSeq_subCDSfromtrans.bed subCDSfromtrans HGSVC2_v2_integrated_callset sv

# sbatch execute_bedWorklargerSIM100s_241106.sh /home/nboev/projects/def-sushant/nboev/data/TAD/hg38/H1-ESC_Dixon_2015-raw_TADs.txt TAD HGSVC2_v2_integrated_callset sv
# sbatch execute_bedWorklargerSIM100s_241106.sh /home/nboev/projects/def-sushant/nboev/data/TAD/hg38/H1-ESC_Dixon_2015-raw_TADsCOMP.bed TADCOMP HGSVC2_v2_integrated_callset sv

# sbatch execute_bedWorklargerSIM100s_241106.sh /home/nboev/projects/def-sushant/nboev/data/LAD/13059_2023_2849_MOESM5_ESM.xlsxtad1.bed LAD1 HGSVC2_v2_integrated_callset sv
# sbatch execute_bedWorklargerSIM100s_241106.sh /home/nboev/projects/def-sushant/nboev/data/LAD/13059_2023_2849_MOESM5_ESM.xlsxtad2.bed LAD2 HGSVC2_v2_integrated_callset sv

# sbatch execute_bedWorklargerSIM100s_241106.sh /home/nboev/projects/def-sushant/nboev/data/FragileSites/FragileRare.bed FragileRare HGSVC2_v2_integrated_callset sv
# sbatch execute_bedWorklargerSIM100s_241106.sh /home/nboev/projects/def-sushant/nboev/data/FragileSites/FragileCommon.bed FragileCommon HGSVC2_v2_integrated_callset sv
