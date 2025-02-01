#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=30000
#SBATCH --time=00:10:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nb.boev@mail.utoronto.ca
#SBATCH --job-name=Masters
#SBATCH --output=Master_%j.o
#SBATCH --error=Master_%j.e
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
echo "${files}"

# The input: 
# 1: Project name and location
# 2: The type of SVs to annotate (ie. sv or svSim). This is a subdirectory that holds the vcf/csv file
# 3: The name of the file to annotate 
# 4: Cell line name to use for epigenomic annotation (must lead to the location of the ENCODE files)

# Creating a directory which will hold these chromosomes files and filtering the original vcf for this chromosome
mkdir /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/SVlen/${files}
awk '$1 == "'${files}'" { print }' /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/SVlen/$3 \
> /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/SVlen/${files}/${files}.full.csv

# We take the filtered file and make 10K splits for this chromosomes
split -l 10000 -d --additional-suffix=.csv \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/SVlen/${files}/${files}.full.csv \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/SVlen/${files}/${files}.full.csv.split

# Get the header so we can append to the batched files
header=$(head -n 1 /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/SVlen/$3)
for f in /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/SVlen/${files}/${files}.full.csv.split*;
	do
	echo "$header" > "${f}.csv"
	tail -n +1 "$f" >> "${f}.csv"
	rm "$f"
	done

splits=$(ls -1 /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/SVlen/${files}/${files}.full.csv.split* | wc -l)
echo "${splits}"

# Getting a list of the number of splits we created for this chromosome
ls /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/SVlen/${files}/${files}.full.csv.split* \
> /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/SVlen/${files}/splitnames.txt
rm /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/SVlen/${files}/${files}.full.csv

# For the following annotations (chromoBand seqFeaturesSV flankSeq SVcoords RepliSeq searchRepeatMasker ENCODE):
# We make a copy of the "template" script and use the number of splits to push job arrays for this annotation time
for script in chromoBand seqFeaturesSV flankSeq SVcoords RepliSeq searchRepeatMasker ;
	do
	mkdir /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/${script}
	cp execute_${script}Template.sh /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/${script}/execute_${script}_${files}.sh
	sed -i "s/array=0-0/array=0-${splits}/g" /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/${script}/execute_${script}_${files}.sh
	sed -i "s/CHROM/${files}/g" /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/${script}/execute_${script}_${files}.sh
	sed -i "s/PROJECT/$1/g" /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/${script}/execute_${script}_${files}.sh
	sed -i "s/LOCATION/$2/g" /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/${script}/execute_${script}_${files}.sh
	mkdir /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/${script}/${files}
	sbatch /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/${script}/execute_${script}_${files}.sh $1 $2
	done

# In this case, we specify the cell type to use for epigenomic annotation
mkdir /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/ENCODE
mkdir /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/ENCODE/$4
cp execute_ENCODETemplate.sh /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/ENCODE/execute_ENCODE_${files}.sh
sed -i "s/array=0-0/array=0-${splits}/g" /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/ENCODE/execute_ENCODE_${files}.sh
sed -i "s/CHROM/${files}/g" /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/ENCODE/execute_ENCODE_${files}.sh
sed -i "s/PROJECT/$1/g" /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/ENCODE/execute_ENCODE_${files}.sh
sed -i "s/LOCATION/$2/g" /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/ENCODE/execute_ENCODE_${files}.sh
mkdir /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/ENCODE/$4/${files}
sbatch /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/ENCODE/execute_ENCODE_${files}.sh $1 $2 $4

# We require the following directories to hold the pre and post processed files
for script in RLoopForming RepeatMasker nonBDNA Blast DNAshape ;
	do
	mkdir /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/${script}
	mkdir /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/${script}/${files}
	mkdir /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/${script}/${files}/premerge
	mkdir /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/${script}/${files}/postmerge
	done

for script in seqFeatures ;
	do
	mkdir /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/${script}
	mkdir /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/${script}/${files}
	done


# example:
# sbatch execute_Master1_240502.sh allHGSVC3_20240415_Freeze4_GRCh38 sv freeze4_GRCh38_sv_insdel_alt.vcf.SVlength.csv H1
# sbatch execute_Master1_240502.sh allHGSVC3_20240415_Freeze4_GRCh38 svSim variants_freeze4_sv_insdel_alt.vcf.SVlength.Simulations.csv H1

