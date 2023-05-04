#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=300000
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nb.boev@mail.utoronto.ca
#SBATCH --job-name=SVlen
#SBATCH --output=SVlen.o
#SBATCH --error=SVlen.e
#SBATCH --qos=privileged

source $HOME/virtenv/bin/activate
module load scipy-stack

#python adding_SVlength_230125.py /home/nboev/projects/def-sushant/nboev/data/HGSVC2_v2_integrated_callset/sv/variants_freeze4_sv_insdel_alt.vcf >> ./sv/adding_SVlength_230125.py.txt
#echo "Adding length annotation for alt file  was successful"

#python adding_SVlength_230125.py /home/nboev/projects/def-sushant/nboev/data/HGSVC2_v2_integrated_callset/sv/variants_freeze4_sv_insdel_sym.vcf >> ./sv/adding_SVlength_230125.py.txt
#echo "Adding length annotation for sym file  was successful"


# adding short svs
#python adding_SVlength_230125.py /home/nboev/projects/def-sushant/nboev/data/HGSVC2_v2_integrated_callset/indel/variants_freeze4_indel_insdel_alt.vcf >> ./indel/adding_SVlength_230125.py.txt
#echo "Adding length annotation for indel file  was successful"

# need to move the file
mv /home/nboev/projects/def-sushant/nboev/preprocess/HGSVC2_v2_integrated_callset/sv/SVlen/variants_freeze4_indel_insdel_alt.vcf.SVlength.csv /home/nboev/projects/def-sushant/nboev/preprocess/HGSVC2_v2_integrated_callset/indel/SVlen
