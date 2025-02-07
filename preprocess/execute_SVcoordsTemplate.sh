#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=300000
#SBATCH --time=02:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nb.boev@mail.utoronto.ca
#SBATCH --job-name=SVcoordsCHROM
#SBATCH --output=/home/nboev/projects/def-sushant/nboev/preprocess/PROJECT/LOCATION/SVcoords/log/SVcoordsCHROM_%j.o
#SBATCH --error=/home/nboev/projects/def-sushant/nboev/preprocess/PROJECT/LOCATION/SVcoords/log/SVcoordsCHROM_%j.e
#SBATCH --array=0-0

# Loading in virtual environment and modules
module --force purge ; module load StdEnv/2020
source $HOME/virtenv/bin/activate
module load scipy-stack
module load samtools
module load bedtools

# For this chromosome, we have a job for each split
names=($(cat /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/SVlen/CHROM/splitnames.txt))
files=${names[${SLURM_ARRAY_TASK_ID}]}
filename=$(basename "$files")
filenamenoext=$(basename "${filename%.*}")

# The input (comes from execute_Master1.sh)
# 1: Project name and location
# 2: The type of SVs to annotate (ie. sv or svSim). This is a subdirectory that holds the vcf/csv file

# We push this python script to generate bed files which represent the 2000bp flanking regions, these will be used for intersections with bedtools
python adding_SVcoords_230321.py \
${files} \
$1 $2 \
${filenamenoext} \
CHROM \
2000 \
>> /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/adding_SVcoords_230321.py.txt

# Annotating Rloops
# Rloops are sites where DNA:RNA hybrids which can cause transcriptional stress and breaks
# Using bed files from here: https://rloopbase.nju.edu.cn/download.jsp, we annotate the beds created above by different "Rloop levels"
	# A level refers to how many pieces of experimental evidence showed Rloop formation at a particular site
 # We annotate using counts for each level
for i in {1..9};
	do
	bedtools annotate -counts -i /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/SVcoords/CHROM/${filenamenoext}.SVcoords2000flank_bedtools.bed \
-files /home/nboev/projects/def-sushant/nboev/data/RLoopForming/processed/level${i}_R-loop_zones.processed.bed \
> /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/RLoopForming/CHROM/premerge/${filenamenoext}.Level${i}.flanks.bed
	done

# We push this python script which takes in these counts and merges + defines the greatest level for a region
python adding_flankRLoop_240521.py \
${files} \
$1 $2 \
${filenamenoext} \
CHROM \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/RLoopForming/CHROM/premerge/${filenamenoext}.Level9.flanks.bed \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/RLoopForming/CHROM/premerge/${filenamenoext}.Level8.flanks.bed \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/RLoopForming/CHROM/premerge/${filenamenoext}.Level7.flanks.bed \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/RLoopForming/CHROM/premerge/${filenamenoext}.Level6.flanks.bed \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/RLoopForming/CHROM/premerge/${filenamenoext}.Level5.flanks.bed \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/RLoopForming/CHROM/premerge/${filenamenoext}.Level4.flanks.bed \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/RLoopForming/CHROM/premerge/${filenamenoext}.Level3.flanks.bed \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/RLoopForming/CHROM/premerge/${filenamenoext}.Level2.flanks.bed \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/RLoopForming/CHROM/premerge/${filenamenoext}.Level1.flanks.bed \
>> /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/adding_flankRLoop_240521.py.txt

rm /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/RLoopForming/CHROM/premerge/${filenamenoext}.Level9.flanks.bed
rm /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/RLoopForming/CHROM/premerge/${filenamenoext}.Level8.flanks.bed
rm /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/RLoopForming/CHROM/premerge/${filenamenoext}.Level7.flanks.bed
rm /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/RLoopForming/CHROM/premerge/${filenamenoext}.Level6.flanks.bed
rm /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/RLoopForming/CHROM/premerge/${filenamenoext}.Level5.flanks.bed
rm /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/RLoopForming/CHROM/premerge/${filenamenoext}.Level4.flanks.bed
rm /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/RLoopForming/CHROM/premerge/${filenamenoext}.Level3.flanks.bed
rm /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/RLoopForming/CHROM/premerge/${filenamenoext}.Level2.flanks.bed
rm /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/RLoopForming/CHROM/premerge/${filenamenoext}.Level1.flanks.bed


# Annotating Repeat Masker elements
# Elements like satellites and LINES can leverage different repair mechanisms for their introduction or expansion
# Bed files were originally downloaded here: https://www.repeatmasker.org/genomes/hg38/RepeatMasker-rm405-db20140131/hg38.fa.out.gz
# Whereby specific motifs were identified, and "collapsed" across a regions. This means there are no overlaps for specific motif at a specific region
	# The motifs we used include: LINE, Low_complexity, LTR, Satellite, Simple_repeat, SINE
# Similar to the above, we merge the counts in the flanks using the py script
for i in 'LINE' 'Low_complexity' 'LTR' 'Satellite' 'Simple_repeat' 'SINE' ;
	do
	bedtools annotate -counts -i /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/SVcoords/CHROM/${filenamenoext}.SVcoords2000flank_bedtools.bed \
-files /home/nboev/projects/def-sushant/nboev/data/RepeatMasker/postprocessed/${i}.post_forbedtools.bed \
> /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/RepeatMasker/CHROM/premerge/${filenamenoext}.${i}.flanks.bed

	python adding_flankRepeatMasker_240521.py \
${files} \
$1 $2 \
${filenamenoext} \
CHROM \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/RepeatMasker/CHROM/premerge/${filenamenoext}.${i}.flanks.bed \
${i} \
>> /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/adding_flankRepeatMasker_240521.py.txt

	rm /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/RepeatMasker/CHROM/premerge/${filenamenoext}.${i}.flanks.bed
	done


# Annotating nonB DNA structures
# More evidence shows that certain types of repeats + G4s can induce genomic instability
# Bed files were originally downloaded here: https://ncifrederick.cancer.gov/bacs/ftp/?nonb
# We similarly process this files
	# The motifs we used include: G4s, ZDNA, a-phased repeats, DRs, MRs, STRs
# Similar to the above, we merge the counts in the flanks using the py script
array=('g-quadruplex_forming_repeats' 'z-dna_motifs' 'a-phased_repeats' 'direct_repeats' 'mirror_repeats' 'short_tandem_repeats' )

for index in ${!array[*]};
	do
	bedtools annotate -counts -i /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/SVcoords/CHROM/${filenamenoext}.SVcoords2000flank_bedtools.bed \
-files /home/nboev/projects/def-sushant/nboev/data/nonBDNA/${array[$index]}/processed/CHROM.${array[$index]}forbedtools.bed \
> /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/nonBDNA/CHROM/premerge/${filenamenoext}.${array[$index]}.flanks.bed

	python adding_flanknonBDNA_240521.py \
${files} \
$1 $2 \
${filenamenoext} \
CHROM \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/nonBDNA/CHROM/premerge/${filenamenoext}.${array[$index]}.flanks.bed \
${array[$index]} \
>> /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/adding_flanknonBDNA_240521.py.txt

	rm /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/nonBDNA/CHROM/premerge/${filenamenoext}.${array[$index]}.flanks.bed

	done
