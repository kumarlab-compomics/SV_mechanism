#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=300000
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nb.boev@mail.utoronto.ca
#SBATCH --job-name=nonBDNAproc
#SBATCH --output=nonBDNAproc.o
#SBATCH --error=nonBDNAproc.e
#SBATCH --qos=privileged

source $HOME/virtenv/bin/activate
module load scipy-stack

# moving the files around... only need to do this once..
#mv /home/nboev/projects/def-sushant/nboev/data/nonBDNA/g-quadruplex_forming_repeats/*.tsv /home/nboev/projects/def-sushant/nboev/data/nonBDNA/g-quadruplex_forming_repeats/raw
#mv /home/nboev/projects/def-sushant/nboev/data/nonBDNA/z-dna_motifs/*.tsv /home/nboev/projects/def-sushant/nboev/data/nonBDNA/z-dna_motifs/raw


# running on 2023/03/15 = after adding a version built for bedtools
#python processing_nonBDNA_230125.py /home/nboev/projects/def-sushant/nboev/data/nonBDNA/g-quadruplex_forming_repeats/raw/chrX_GQ.tsv G_Quadruplex_Motif \
#/home/nboev/projects/def-sushant/nboev/data/chromoBand/ChromosomeBand_hg38.txt \
#>> ./processing_nonBDNA_230125.py.txt
#python processing_nonBDNA_230125.py /home/nboev/projects/def-sushant/nboev/data/nonBDNA/z-dna_motifs/raw/chrX_Z.tsv Z_DNA_Motif \
#/home/nboev/projects/def-sushant/nboev/data/chromoBand/ChromosomeBand_hg38.txt \
#>> ./processing_nonBDNA_230125.py.txt
#echo "Processing of Chromosome X was successful"

#python processing_nonBDNA_230125.py /home/nboev/projects/def-sushant/nboev/data/nonBDNA/g-quadruplex_forming_repeats/raw/chrY_GQ.tsv G_Quadruplex_Motif \
#/home/nboev/projects/def-sushant/nboev/data/chromoBand/ChromosomeBand_hg38.txt \
#>> ./processing_nonBDNA_230125.py.txt
#python processing_nonBDNA_230125.py /home/nboev/projects/def-sushant/nboev/data/nonBDNA/z-dna_motifs/raw/chrY_Z.tsv Z_DNA_Motif \
#/home/nboev/projects/def-sushant/nboev/data/chromoBand/ChromosomeBand_hg38.txt \
#>> ./processing_nonBDNA_230125.py.txt
#echo "Processing of Chromosome Y was successful"


# for all the autosomes
#for i in {1..22};
#do
#	python processing_nonBDNA_230125.py /home/nboev/projects/def-sushant/nboev/data/nonBDNA/g-quadruplex_forming_repeats/raw/chr${i}_GQ.tsv G_Quadruplex_Motif \
#/home/nboev/projects/def-sushant/nboev/data/chromoBand/ChromosomeBand_hg38.txt \
#>> ./processing_nonBDNA_230125.py.txt
#	python processing_nonBDNA_230125.py /home/nboev/projects/def-sushant/nboev/data/nonBDNA/z-dna_motifs/raw/chr${i}_Z.tsv Z_DNA_Motif \
#/home/nboev/projects/def-sushant/nboev/data/chromoBand/ChromosomeBand_hg38.txt \
#>> ./processing_nonBDNA_230125.py.txt
#	echo "Processing of Chromosome "$i" was successful"
#done


# ADDING THE OTHERS I ADDED =2023/04/26

# looping through the nonB for X, Y and autosomes
#array=('a-phased_repeats' 'direct_repeats' 'inverted_repeats' 'mirror_repeats' 'short_tandem_repeats' )
#array2=('APR' 'DR' 'IR' 'MR' 'STR')

# we had a timeout issue --> I believe IR? files are big? = need to rerun with new list
#array=('inverted_repeats' 'mirror_repeats' 'short_tandem_repeats' )
#array2=('IR' 'MR' 'STR')

# inverteds are not working? must be super big files? skip them for now = 2023/05/04
array=('mirror_repeats' 'short_tandem_repeats' )
array2=('MR' 'STR')

for index in ${!array[*]};
do
	python processing_nonBDNA_230125.py /home/nboev/projects/def-sushant/nboev/data/nonBDNA/${array[$index]}/raw/chrX_${array2[$index]}.tsv ${array[$index]} \
/home/nboev/projects/def-sushant/nboev/data/chromoBand/ChromosomeBand_hg38.txt \
>> ./processing_nonBDNA_230125.py.txt

	python processing_nonBDNA_230125.py /home/nboev/projects/def-sushant/nboev/data/nonBDNA/${array[$index]}/raw/chrY_${array2[$index]}.tsv ${array[$index]} \
/home/nboev/projects/def-sushant/nboev/data/chromoBand/ChromosomeBand_hg38.txt \
>> ./processing_nonBDNA_230125.py.txt

	for i in {1..22};
	do
		python processing_nonBDNA_230125.py /home/nboev/projects/def-sushant/nboev/data/nonBDNA/${array[$index]}/raw/chr${i}_${array2[$index]}.tsv ${array[$index]} \
/home/nboev/projects/def-sushant/nboev/data/chromoBand/ChromosomeBand_hg38.txt \
>> ./processing_nonBDNA_230125.py.txt
	done
done

