#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=5G
#SBATCH --time=5:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nb.boev@mail.utoronto.ca
#SBATCH --job-name=Survivor_CHROM
#SBATCH --output=/home/nboev/projects/def-sushant/nboev/data/SimulatedSVs/PROJECT/LOCATION/SurvivorCHROM_%j.o
#SBATCH --error=/home/nboev/projects/def-sushant/nboev/data/SimulatedSVs/PROJECT/LOCATION/SurvivorCHROM_%j.e
#SBATCH --array=0-0

# changed the time +memory req, was 15G and 3.5 days

module --force purge ; module load StdEnv/2020
module load samtools
module load survivor
source $HOME/virtenv/bin/activate

# looping through files from the split
names=($(cat /home/nboev/scratch/data/SimulatedSVs/$1/$2/rawIND/CHROM/splitnames.txt))
files=${names[${SLURM_ARRAY_TASK_ID}]}
echo "${files}"

split=$(basename ${files} .txt)
mkdir /home/nboev/scratch/data/SimulatedSVs/$1/$2/rawIND/CHROM/${split}

# Add something to capture different referece versions!! = instead of hg38_seq, allow user to pick!, we don't need to do this here!
cat ${files} | while read line;
        do
	echo $line
        echo /home/nboev/scratch/data/SimulatedSVs/$1/$2/parameterIND/CHROM/$line

        SURVIVOR simSV \
/home/nboev/projects/def-sushant/nboev/data/Genome/hg38_seq/chroms/CHROM.fa \
/home/nboev/scratch/data/SimulatedSVs/$1/$2/parameterIND/CHROM/$line \
0.1 0 \
/home/nboev/scratch/data/SimulatedSVs/$1/$2/rawIND/CHROM/${split}/$line

        rm /home/nboev/scratch/data/SimulatedSVs/$1/$2/rawIND/CHROM/${split}/*txt.fasta
        rm /home/nboev/scratch/data/SimulatedSVs/$1/$2/rawIND/CHROM/${split}/*txt.vcf

        if [[ ${line} == *"INS"* ]]; then
                echo "insertion"

                python /home/nboev/projects/def-sushant/nboev/data/SimulatedSVs/20220422_3202_phased_SNV_INDEL_SV_bychrom/generating_vcfINS_240430.py \
/home/nboev/scratch/data/SimulatedSVs/$1/$2/rawIND/CHROM/${split}/${line}.insertions.fa \
/home/nboev/scratch/data/SimulatedSVs/$1/$2/rawIND/CHROM/${split}/${line}.bed \
CHROM ${line} \
100 \
>> /home/nboev/scratch/data/SimulatedSVs/$1/$2/rawIND/CHROM/${split}/generating_vcfINS_240430.py.txt

                rm /home/nboev/scratch/data/SimulatedSVs/$1/$2/rawIND/CHROM/${split}/${line}.insertions.fa
                rm /home/nboev/scratch/data/SimulatedSVs/$1/$2/rawIND/CHROM/${split}/${line}.bed
                rm /home/nboev/scratch/data/SimulatedSVs/$1/$2/parameterIND/CHROM/$line
        else
            	echo "deletion"
                rm /home/nboev/scratch/data/SimulatedSVs/$1/$2/rawIND/CHROM/${split}/${line}.insertions.fa

                python /home/nboev/projects/def-sushant/nboev/data/SimulatedSVs/20220422_3202_phased_SNV_INDEL_SV_bychrom/generating_vcfDEL_240430.py \
/home/nboev/scratch/data/SimulatedSVs/$1/$2/rawIND/CHROM/${split}/${line}.bed \
CHROM ${line} \
100 \
$3 \
>> /home/nboev/scratch/data/SimulatedSVs/$1/$2/rawIND/CHROM/${split}/generating_vcfDEL_240430.py.txt

                rm /home/nboev/scratch/data/SimulatedSVs/$1/$2/rawIND/CHROM/${split}/${line}.bed
                rm /home/nboev/scratch/data/SimulatedSVs/$1/$2/parameterIND/CHROM/$line

        fi
        done

awk '(NR == 1) || (FNR > 1)' \
/home/nboev/scratch/data/SimulatedSVs/$1/$2/rawIND/CHROM/${split}/*.csv \
> /home/nboev/projects/def-sushant/nboev/data/SimulatedSVs/$1/$2/processedIND/CHROM/${split}.csv

rm /home/nboev/scratch/data/SimulatedSVs/$1/$2/rawIND/CHROM/${split}/*.csv



