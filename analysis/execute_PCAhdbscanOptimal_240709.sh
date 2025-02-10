#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=30000
#SBATCH --time=00:15:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nb.boev@mail.utoronto.ca
#SBATCH --job-name=hdbscan_new
#SBATCH --output=hdbscan_new.o
#SBATCH --error=hdbscan_new.e
#SBATCH --qos=privileged

module --force purge ; module load StdEnv/2020
source $HOME/virtenv/bin/activate
module load gcc/9.3.0 r/4.0.2
module load scipy-stack

# Using the optimized parameters = 2024/07/09 (see the progress presentation on 2024/07/05)

# For insertions: 2024/08/20
python analysis_hdbscanIDOptimalHOMO_240708.py \
/home/nboev/projects/def-sushant/nboev/analysis/HGSVC2_v2_integrated_callset/sv/IDmechsvSIM/20240625/ \
insertion_ID0625.tsv \
240625/ \
5 \
1300 \
braycurtis \
HGSVC2_v2_integrated_callset \
sv \
insertion \
>> ./sv/analysis_hdbscanIDOptimalHOMO_240708.py.txt

# For deletions: 2024/08/20
python analysis_hdbscanIDOptimalHOMO_240708.py \
/home/nboev/projects/def-sushant/nboev/analysis/HGSVC2_v2_integrated_callset/sv/IDmechsvSIM/20240625/ \
deletion_ID0625.tsv \
240625/ \
5 \
700 \
braycurtis \
HGSVC2_v2_integrated_callset \
sv \
deletion \
>> ./sv/analysis_hdbscanIDOptimalHOMO_240708.py.txt

