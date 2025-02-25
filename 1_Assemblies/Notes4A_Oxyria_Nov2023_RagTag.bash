###################################
# Oxyria Svalbard genome RagTag mapping
# Oct 2023
###############################
# download the other genome
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_029168935.1/
 
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/029/168/935/GCA_029168935.1_ASM2916893v1/GCA_029168935.1_ASM2916893v1_genomic.fna.gz
gunzip GCA_029168935.1_ASM2916893v1_genomic.fna.gz

#---------------------------------------------
# RagTag
tmux new-session -s ragtag
tmux attach-session -t ragtag

salloc -c32 --time 2:55:00 --mem 120000m --account def-rieseber

source ~/scratch/Annotation/Poison_Oak/RagTag/RagTag_vir_env/bin/activate

module load StdEnv/2020
module load python/3.10.2
module load minimap2/2.24
module load scipy-stack/2022a

pip install RagTag

cd /home/celphin/scratch/repeats/input_chromosomes/Oxyria

/home/celphin/scratch/Annotation/Poison_Oak/RagTag/RagTag-master/ragtag.py \
scaffold Oxy_dig_1.scaffolds_FINAL.final.review.fasta GCA_029168935.1_ASM2916893v1_genomic.fna \
-t 28 -o ./Oxyria_ragtag_output/

#-------------------
# Check assemblies
module load StdEnv/2020 
module load gcc/9.3.0
module load r-bundle-bioconductor/3.14
cd ~/assembly-stats/

./assembly-stats ~/scratch/repeats/input_chromosomes/Oxyria/Oxyria_ragtag_output/ragtag.scaffold.fasta

stats for /home/celphin/scratch/repeats/input_chromosomes/Oxyria/Oxyria_ragtag_output/ragtag.scaffold.fasta
sum = 561229550, n = 333, ave = 1685374.02, largest = 86582034
N50 = 78410798, n = 4
N60 = 76064323, n = 5
N70 = 76064323, n = 5
N80 = 73303751, n = 6
N90 = 72361354, n = 7
N100 = 17254, n = 333
N_count = 2400
Gaps = 24

./assembly-stats ~/scratch/repeats/input_chromosomes/Oxyria/Oxy_dig_1.scaffolds_FINAL.final.review.fasta

stats for /home/celphin/scratch/repeats/input_chromosomes/Oxyria/Oxy_dig_1.scaffolds_FINAL.final.review.fasta
sum = 590427020, n = 1560, ave = 378478.86, largest = 88146246
N50 = 76036264, n = 4
N60 = 73842794, n = 5
N70 = 73289246, n = 6
N80 = 70136240, n = 7
N90 = 70136240, n = 7
N100 = 7, n = 1560
N_count = 87158
Gaps = 250

##########################################