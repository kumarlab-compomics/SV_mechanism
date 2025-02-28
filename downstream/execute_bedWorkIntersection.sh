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

# Loading in virtual environment and modules
source $HOME/virtenv/bin/activate
module load samtools
module load bedtools

# The input: 
# 1: bed file to use for annotation/ intersection
# 2: The shorthand name for the above bed file/ genomie element
# 3: Project name/ location
# 4: sv

# Here, we loop through insertions and deletions separately, along with for each homology-based label
# We annotate the bed files generated (either real or simulated), using the bed file provided. We clean up the file so there are no duplicates
	# Note: For the simulations, we loop through all 100 iterations separately
 # For each simulated iteration, we send the intersected and orginal files to the analysis_ES.py script. 
 # We clean up files as we go. 

for j in insertion deletion ;
       do
	for i in HLH ILH NLH ;
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

			python analysis_ES.py \
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


# sbatch execute_bedWorklargerSIM100s_241106.sh /home/nboev/projects/def-sushant/nboev/data/SCREEN/GRCh38-cCREs_dELS.bed ALLdELS HGSVC2_v2_integrated_callset sv
