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
#SBATCH --array=0-25

# remember to change the array number back

# in this script, we need to push each chromosome through job arrays
module --force purge ; module load StdEnv/2020
source $HOME/virtenv/bin/activate
module load gcc/9.3.0 r/4.0.2
module load scipy-stack

# We need to actually merge all the files together now
# Generating the files we need to merge, we do this in parallel within one chromosome for all the splits

# for now, we'll use the chromo list in the simulation directory
#names=($(cat ./chromos14.txt))
names=($(cat ./chromos.txt))
#names=($(cat ./chromosontthird.txt))

files=${names[${SLURM_ARRAY_TASK_ID}]}
echo "${files}"

mkdir /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/merge_FeatMatrix
mkdir /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/merge_FeatMatrix/${files}

cat /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/SVlen/${files}/splitnames.txt | while read line
	do
	loc=$(basename "$line")
	locnoext=$(basename "${loc%.*}")
	echo "${locnoext}"

#	head /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/seqFeaturesSV/${files}/${locnoext}.withSVSeqFeat.csv
#	head /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/seqFeatures/${files}/${locnoext}.with2000SeqFeat.csv

	python merge_FeatMatrix_240530.py \
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

        python merge_FeatMatrix_240530.py \
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

# we need to merge the files together
#for ext in DEL deletion INS insertion ;
for ext in DEL INS ;
	do
	cat /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/merge_FeatMatrix/${files}/*.${ext}.FeatMatrix.csv \
> /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/merge_FeatMatrix/${files}.full.csv.${ext}.FeatMatrix.csv
	done



# NOTE = 30MINS IS NOT ENOUGH FOR 20 SPLITS

# how to generate: variants_freeze4_sv_insdel_alt.vcf.SVlength.Simulations.csv within SVlen: See /home/nboev/projects/def-sushant/nboev/data/SimulatedSVs/HGSVC2_v2_integrated_callset/execute_PostProcessing_240502.sh

# Running:
# sbatch execute_Master1Merge_240530.sh HGSVC2_v2_integrated_callset sv H1
# sbatch execute_Master1Merge_240530.sh HGSVC2_v2_integrated_callset svSim H1

# Running for 1000Gs
# sbatch execute_Master1Merge_240530.sh 20220422_3202_phased_SNV_INDEL_SV_bychrom SVTrue_typedeletion_resTrueIND
# Running for 1000Gs simulations
# sbatch execute_Master1Merge_240530.sh 20220422_3202_phased_SNV_INDEL_SV_bychrom SVTrue_typedeletion_resTrueINDsvSIM

# Running for 1000Gs updated file
# sbatch execute_Master1Merge_240530.sh 20220422_3202_phased_SNV_INDEL_SV_bychrom SVTrue_typedeletion_resTrue H1
# sbatch execute_Master1Merge_240530.sh 20220422_3202_phased_SNV_INDEL_SV_bychrom SVTrue_typedeletion_resTrueSIM H1

# Running on the CPTAC validation set: AFTER ADDING THE epigenetic part
# sbatch execute_Master1Merge_240530.sh CPTAC3 bf584e4c-f981-473e-b40a-c73d7f3699d2.wgs.sanger_raw_pindel.raw_somatic_mutation HelaS3
# sbatch execute_Master1Merge_240530.sh CPTAC3 bf584e4c-f981-473e-b40a-c73d7f3699d2.wgs.sanger_raw_pindel.SIM HelaS3
# sbatch execute_Master1Merge_240530.sh CPTAC3 fcc8a59b-19b0-4b78-9ca1-f5b9e580f183.wgs.sanger_raw_pindel.raw_somatic_mutation HelaS3
# sbatch execute_Master1Merge_240530.sh CPTAC3 fcc8a59b-19b0-4b78-9ca1-f5b9e580f183.wgs.sanger_raw_pindel.SIM HelaS3

# Running on Belyeu2021 = 2024/07/25
# sbatch execute_Master1Merge_240530.sh Belyeu2021 sv H1
# sbatch execute_Master1Merge_240530.sh Belyeu2021 svSim H1

# Running for HGSVC3 = 2024/08/19
# sbatch execute_Master1Merge_240530.sh allHGSVC3_20240415_Freeze4_GRCh38 sv H1
# sbatch execute_Master1Merge_240530.sh allHGSVC3_20240415_Freeze4_GRCh38 svSim H1

# Running for ONT data = 2024/09/10
# sbatch execute_Master1Merge_240530.sh ONT_1000Gs_100FREEZE sv H1
# sbatch execute_Master1Merge_240530.sh ONT_1000Gs_100FREEZE svSim H1

# Running for DepMap = 2024/10/20
# sbatch execute_Master1Merge_240530.sh DepMap CA922 H1
