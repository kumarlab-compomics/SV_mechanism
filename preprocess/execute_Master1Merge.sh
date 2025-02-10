#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=3G
#SBATCH --time=03:30:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nb.boev@mail.utoronto.ca
#SBATCH --job-name=MastersMerge
#SBATCH --output=MasterMerge_%j.o
#SBATCH --error=MasterMerge_%j.e
#SBATCH --qos=privilege
#SBATCH --array=0-24

# Loading in virtual environment and modules
module --force purge ; module load StdEnv/2020
source $HOME/virtenv/bin/activate
module load gcc/9.3.0 r/4.0.2
module load scipy-stack

# To run 24 jobs in array, based on this list of chromosomes
names=($(cat ./chromos.txt))
files=${names[${SLURM_ARRAY_TASK_ID}]}

# The input: 
# 1: Project name and location
# 2: The type of SVs to annotate (ie. sv or svSim). This is a subdirectory that holds the vcf/csv file
# 3: Cell line name to use for epigenomic annotation (must lead to the location of the ENCODE files)

# Creating a directory which will hold these chromosomes files
mkdir /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/merge_FeatMatrix
mkdir /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/merge_FeatMatrix/${files}

# The following:
# Loops through all the splits for a particular chromosome
# For insertions and deletions, we handle merging separately, since epigenetic features differ (ie. we cannot annotate the epigenetic state as per the inserted sequence)
# We push the merging python script, merge_FeatMatrix.py, which systematically merges the SVs as per their unique ID
# The processing we do for each file differs sligtly, therefore we "manually" control each file input

cat /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/SVlen/${files}/splitnames.txt | while read line
	do
	loc=$(basename "$line")
	locnoext=$(basename "${loc%.*}")

	python merge_FeatMatrix.py \
$1 $2 ${locnoext} \
${files} DEL \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/seqFeaturesSV/${files}/${locnoext}.withSVSeqFeat.csv \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/seqFeatures/${files}/${locnoext}.with2000SeqFeat.csv \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/searchRepeatMasker/${files}/${locnoext}/${locnoext}.searchRepMasker.csv \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/RepeatMasker/${files}/postmerge/LINE.${locnoext}.flanks.csv \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/RepeatMasker/${files}/postmerge/Low_complexity.${locnoext}.flanks.csv \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/RepeatMasker/${files}/postmerge/LTR.${locnoext}.flanks.csv \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/RepeatMasker/${files}/postmerge/Satellite.${locnoext}.flanks.csv \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/RepeatMasker/${files}/postmerge/Simple_repeat.${locnoext}.flanks.csv \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/RepeatMasker/${files}/postmerge/SINE.${locnoext}.flanks.csv \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/nonBDNA/${files}/postmerge/a-phased_repeats.${locnoext}.flanks.csv \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/nonBDNA/${files}/postmerge/direct_repeats.${locnoext}.flanks.csv \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/nonBDNA/${files}/postmerge/g-quadruplex_forming_repeats.${locnoext}.flanks.csv \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/nonBDNA/${files}/postmerge/mirror_repeats.${locnoext}.flanks.csv \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/nonBDNA/${files}/postmerge/short_tandem_repeats.${locnoext}.flanks.csv \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/nonBDNA/${files}/postmerge/z-dna_motifs.${locnoext}.flanks.csv \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/chromoBand/${files}/${locnoext}.SVlength.csv \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/Blast/${files}/postmerge/${locnoext}.Blastmergesungap.csv \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/ENCODE/$3/${files}/${locnoext}.2000.EpiFeat.csv \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/DNAshape/${files}/postmerge/${locnoext}.DNAshape.csv \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/RepliSeq/${files}/${locnoext}.withRepliSeq.csv \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/RLoopForming/${files}/postmerge/${locnoext}.RLoopflanks_merged.csv \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/ENCODE/$3/${files}/${locnoext}.EpiFeatSVs.csv \
>> /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/merge_FeatMatrix_240530.py.txt

        python merge_FeatMatrix.py \
$1 $2 ${locnoext} \
${files} INS \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/seqFeaturesSV/${files}/${locnoext}.withSVSeqFeat.csv \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/seqFeatures/${files}/${locnoext}.with2000SeqFeat.csv \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/searchRepeatMasker/${files}/${locnoext}/${locnoext}.searchRepMasker.csv \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/RepeatMasker/${files}/postmerge/LINE.${locnoext}.flanks.csv \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/RepeatMasker/${files}/postmerge/Low_complexity.${locnoext}.flanks.csv \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/RepeatMasker/${files}/postmerge/LTR.${locnoext}.flanks.csv \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/RepeatMasker/${files}/postmerge/Satellite.${locnoext}.flanks.csv \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/RepeatMasker/${files}/postmerge/Simple_repeat.${locnoext}.flanks.csv \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/RepeatMasker/${files}/postmerge/SINE.${locnoext}.flanks.csv \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/nonBDNA/${files}/postmerge/a-phased_repeats.${locnoext}.flanks.csv \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/nonBDNA/${files}/postmerge/direct_repeats.${locnoext}.flanks.csv \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/nonBDNA/${files}/postmerge/g-quadruplex_forming_repeats.${locnoext}.flanks.csv \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/nonBDNA/${files}/postmerge/mirror_repeats.${locnoext}.flanks.csv \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/nonBDNA/${files}/postmerge/short_tandem_repeats.${locnoext}.flanks.csv \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/nonBDNA/${files}/postmerge/z-dna_motifs.${locnoext}.flanks.csv \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/chromoBand/${files}/${locnoext}.SVlength.csv \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/Blast/${files}/postmerge/${locnoext}.Blastmergesungap.csv \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/ENCODE/$3/${files}/${locnoext}.2000.EpiFeat.csv \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/DNAshape/${files}/postmerge/${locnoext}.DNAshape.csv \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/RepliSeq/${files}/${locnoext}.withRepliSeq.csv \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/RLoopForming/${files}/postmerge/${locnoext}.RLoopflanks_merged.csv \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/ENCODE/$3/${files}/${locnoext}.EpiFeatSVs.csv \
>> /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/merge_FeatMatrix_240530.py.txt
	done

# We merge these files all together 
for ext in DEL INS ;
	do
	cat /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/merge_FeatMatrix/${files}/*.${ext}.FeatMatrix.csv \
> /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/merge_FeatMatrix/${files}.full.csv.${ext}.FeatMatrix.csv
	done


# example :
# sbatch execute_Master1Merge.sh HGSVC2_v2_integrated_callset sv H1
# sbatch execute_Master1Merge.sh HGSVC2_v2_integrated_callset svSim H1
