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

# Loading in virtual environment and modules
module --force purge ; module load StdEnv/2020
module load samtools
module load survivor
source $HOME/virtenv/bin/activate

# We run the correct number of jobs, as per the number of SVs were present, and therefore "splits" created. 
names=($(cat /home/nboev/scratch/data/SimulatedSVs/$1/$2/rawIND/CHROM/splitnames.txt))
files=${names[${SLURM_ARRAY_TASK_ID}]}
split=$(basename ${files} .txt)
mkdir /home/nboev/scratch/data/SimulatedSVs/$1/$2/rawIND/CHROM/${split}

# The input: (sent from execute_svSimsMaster.sh)
# 1: Project name and location
# 2: The type of SVs to annotate (ie. sv or svSim). This is a subdirectory that holds the vcf/csv file

# The following: 
# Loops through the configuration files
# For each configuration file, we create simulations using the simSV function from SURVIVOR. We supplement the correct reference fasta by providing the chromosome-specific fasta
# For insertions, we require the fasta file produced by SURVIVOR. The bed and the fasta are sent to the generating_vcfINS.py script
# For deletions, we only need the bed file. The REF sequences are taken from the reference genome using faidx. This is done in generating_vcfDEL.py
# Since we are only generating 100 simulations/ real SV, we use 100 as the input argument for these python functions. 
# We delete unnecessary files as we go to mitigate overwriting and space issues

cat ${files} | while read line;
        do
        SURVIVOR simSV \
/home/nboev/projects/def-sushant/nboev/data/Genome/hg38_seq/chroms/CHROM.fa \
/home/nboev/scratch/data/SimulatedSVs/$1/$2/parameterIND/CHROM/$line \
0.1 0 \
/home/nboev/scratch/data/SimulatedSVs/$1/$2/rawIND/CHROM/${split}/$line

        rm /home/nboev/scratch/data/SimulatedSVs/$1/$2/rawIND/CHROM/${split}/*txt.fasta
        rm /home/nboev/scratch/data/SimulatedSVs/$1/$2/rawIND/CHROM/${split}/*txt.vcf

        if [[ ${line} == *"INS"* ]]; then
                echo "insertion"

                python generating_vcfINS.py \
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

                python generating_vcfDEL.py \
/home/nboev/scratch/data/SimulatedSVs/$1/$2/rawIND/CHROM/${split}/${line}.bed \
CHROM ${line} \
100 \
$3 \
>> /home/nboev/scratch/data/SimulatedSVs/$1/$2/rawIND/CHROM/${split}/generating_vcfDEL_240430.py.txt

                rm /home/nboev/scratch/data/SimulatedSVs/$1/$2/rawIND/CHROM/${split}/${line}.bed
                rm /home/nboev/scratch/data/SimulatedSVs/$1/$2/parameterIND/CHROM/$line

        fi
        done

# We merge across all the files for this split
awk '(NR == 1) || (FNR > 1)' \
/home/nboev/scratch/data/SimulatedSVs/$1/$2/rawIND/CHROM/${split}/*.csv \
> /home/nboev/projects/def-sushant/nboev/data/SimulatedSVs/$1/$2/processedIND/CHROM/${split}.csv

rm /home/nboev/scratch/data/SimulatedSVs/$1/$2/rawIND/CHROM/${split}/*.csv



