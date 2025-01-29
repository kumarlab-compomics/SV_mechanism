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

# in this script, we need to push each chromosome through job arrays
module --force purge ; module load StdEnv/2020
source $HOME/virtenv/bin/activate
module load gcc/9.3.0 r/4.0.2
module load scipy-stack

# for now, we'll use the chromo list in the simulation directory
names=($(cat ./chromos.txt))

files=${names[${SLURM_ARRAY_TASK_ID}]}
echo "${files}"


# I need to split up the file by chromosome and make splits we can actually annotate
# In this case, maybe we can make row splits of 10K per chromosome

# we might need to use an argument for the executes, for example which location we'd actually be placing these files
# ie. HGSVC + svSim vs sv etc..
mkdir /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/SVlen/${files}

# need to make the directory below include arguments!
# want to keep the header of the file
awk '$1 == "'${files}'" { print }' /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/SVlen/$3 \
> /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/SVlen/${files}/${files}.full.csv

# making splits of 10K, might need to increase later
split -l 10000 -d --additional-suffix=.csv \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/SVlen/${files}/${files}.full.csv \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/SVlen/${files}/${files}.full.csv.split

# getting header
header=$(head -n 1 /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/SVlen/$3)

for f in /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/SVlen/${files}/${files}.full.csv.split*;
	do
	echo "$header" > "${f}.csv"
	tail -n +1 "$f" >> "${f}.csv"
	rm "$f"
	done

splits=$(ls -1 /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/SVlen/${files}/${files}.full.csv.split* | wc -l)
echo "${splits}"

ls /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/SVlen/${files}/${files}.full.csv.split* \
> /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/SVlen/${files}/splitnames.txt

rm /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/SVlen/${files}/${files}.full.csv

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

# we need this for the ENCODE files, we need to add the specificatio of the cell line
mkdir /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/ENCODE
mkdir /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/ENCODE/$4
cp execute_ENCODETemplate.sh /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/ENCODE/execute_ENCODE_${files}.sh
sed -i "s/array=0-0/array=0-${splits}/g" /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/ENCODE/execute_ENCODE_${files}.sh
sed -i "s/CHROM/${files}/g" /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/ENCODE/execute_ENCODE_${files}.sh
sed -i "s/PROJECT/$1/g" /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/ENCODE/execute_ENCODE_${files}.sh
sed -i "s/LOCATION/$2/g" /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/ENCODE/execute_ENCODE_${files}.sh
mkdir /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/ENCODE/$4/${files}
sbatch /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/ENCODE/execute_ENCODE_${files}.sh $1 $2 $4

# we need to add this for the flank annotations:
for script in RLoopForming RepeatMasker nonBDNA Blast DNAshape ;
	do
	mkdir /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/${script}
	mkdir /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/${script}/${files}
	mkdir /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/${script}/${files}/premerge
	mkdir /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/${script}/${files}/postmerge
	done

# we need this for flank seq annotations:
for script in seqFeatures ;
	do
	mkdir /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/${script}
	mkdir /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/${script}/${files}
	done



# how to generate: variants_freeze4_sv_insdel_alt.vcf.SVlength.Simulations.csv within SVlen: See /home/nboev/projects/def-sushant/nboev/data/SimulatedSVs/HGSVC2_v2_integrated_callset/execute_PostProcessing_240502.sh


# we need to include arguments to actually use: Running = 2024/06/03

# Doing HGSVC3 heere instead of Narval = 2024/08/
# sbatch execute_Master1_240502.sh allHGSVC3_20240415_Freeze4_GRCh38 sv freeze4_GRCh38_sv_insdel_alt.vcf.SVlength.csv H1
# sbatch execute_Master1_240502.sh allHGSVC3_20240415_Freeze4_GRCh38 svSim variants_freeze4_sv_insdel_alt.vcf.SVlength.Simulations.csv H1

