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

module --force purge ; module load StdEnv/2020
source $HOME/virtenv/bin/activate
module load scipy-stack
module load samtools
module load bedtools

# looping through files from the split
names=($(cat /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/SVlen/CHROM/splitnames.txt))
files=${names[${SLURM_ARRAY_TASK_ID}]}
echo "${files}"
filename=$(basename "$files")
filenamenoext=$(basename "${filename%.*}")

echo "${filename}"
echo "${filenamenoext}"


# I think we need to add something in the python file so they don't overwrite!
python adding_SVcoords_230321.py \
${files} \
$1 $2 \
${filenamenoext} \
CHROM \
2000 \
>> /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/adding_SVcoords_230321.py.txt

# ANNOTATION OF R LOOPS - good!!
for i in {1..9};
	do
	bedtools annotate -counts -i /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/SVcoords/CHROM/${filenamenoext}.SVcoords2000flank_bedtools.bed \
-files /home/nboev/projects/def-sushant/nboev/data/RLoopForming/processed/level${i}_R-loop_zones.processed.bed \
> /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/RLoopForming/CHROM/premerge/${filenamenoext}.Level${i}.flanks.bed
	done

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


# ANNOTATION OF REPEAT MASKER
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


# ANNOTATION OF NONB !!! need to alter this!!!

array=('g-quadruplex_forming_repeats' 'z-dna_motifs' 'a-phased_repeats' 'direct_repeats' 'mirror_repeats' 'short_tandem_repeats' )
#array2=('G_Quadruplex_Motif' 'Z_DNA_Motif' 'a-phased_repeats' 'direct_repeats' 'mirror_repeats' 'short_tandem_repeats' )

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
