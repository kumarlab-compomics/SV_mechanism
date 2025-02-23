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

# Loading in virtual environment and modules
module --force purge ; module load StdEnv/2020
source $HOME/virtenv/bin/activate
module load gcc/9.3.0 r/4.0.2
module load scipy-stack

# No inputs required.
# Based on the results from execute_PCAhdbscanins.sh, we established the optimal models:
  # Insertions: PCA using 5 components, a minimum of 1300 SVs required for a cluster, Bray-Curtis distance metric
  # Deletions: PCA using 5 components, a minimum of 700 SVs required for a cluster, Bray-Curtis distance metric

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
