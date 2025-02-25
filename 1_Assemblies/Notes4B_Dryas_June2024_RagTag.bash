######################################
# RagTag Dryas genomes
# June 2024
###################################
# download genomes
cd ~/scratch/Dryas/Dryas_genomes/

# Dryas alaskensis 
# Alaska
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_026122645.1/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/026/122/645/GCA_026122645.1_ASM2612264v1/GCA_026122645.1_ASM2612264v1_genomic.fna.gz
gunzip GCA_026122645.1_ASM2612264v1_genomic.fna.gz
mv GCA_026122645.1_ASM2612264v1_genomic.fna Dry-alask.fasta

# Dryas integrifolia
# Alaska
# Entireleaf Mountain-Avens
#https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_028570905.1
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/028/570/905/GCA_028570905.1_ASM2857090v1/GCA_028570905.1_ASM2857090v1_genomic.fna.gz
gunzip GCA_028570905.1_ASM2857090v1_genomic.fna.gz
mv GCA_028570905.1_ASM2857090v1_genomic.fna Dry-int.fasta

# Dryas octopetala H2 - D. ajanensis
# Alaska
# Eight-petal Mountain-Avens
#https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_025631235.1
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/025/631/235/GCA_025631235.1_ASM2563123v1/GCA_025631235.1_ASM2563123v1_genomic.fna.gz
gunzip GCA_025631235.1_ASM2563123v1_genomic.fna.gz
mv GCA_025631235.1_ASM2563123v1_genomic.fna Dry-octo-H2.fasta

# Dryas drummondii
# unknown location
# Yellow Mountain-avens
#https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_003254865.1
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/254/865/GCA_003254865.1_ASM325486v1/GCA_003254865.1_ASM325486v1_genomic.fna.gz
gunzip GCA_003254865.1_ASM325486v1_genomic.fna.gz
mv GCA_003254865.1_ASM325486v1_genomic.fna Dry-drumm.fasta

# Dryas octopetala H0
# 	ROYAL BOTANIC GARDEN EDINBURGH H1
# Eight-petal Mountain-Avens
#https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_963921425.1
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/963/921/425/GCA_963921425.1_drDryOcto1.1/GCA_963921425.1_drDryOcto1.1_genomic.fna.gz
gunzip GCA_963921425.1_drDryOcto1.1_genomic.fna.gz
mv GCA_963921425.1_drDryOcto1.1_genomic.fna Dry-octo-H0.fasta

# Dryas octopetala H1
# 	ROYAL BOTANIC GARDEN EDINBURGH H2
# Eight-petal Mountain-Avens
#https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_963921435.1
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/963/921/435/GCA_963921435.1_drDryOcto1.1_alternate_haplotype/GCA_963921435.1_drDryOcto1.1_alternate_haplotype_genomic.fna.gz
gunzip GCA_963921435.1_drDryOcto1.1_alternate_haplotype_genomic.fna.gz
mv GCA_963921435.1_drDryOcto1.1_alternate_haplotype_genomic.fna Dry-octo-H1.fasta




############################################################
# Map to the Dryas octopetala assembly

# RagTag - with order and join scaffolds based on reference sequence - will not change the scaffolds though
# https://github.com/malonge/RagTag
# https://github.com/malonge/RagTag/wiki/scaffold 

# Needs
    # Minimap2
    # Python 3 (with the following auto-installed packages)
        # numpy
        # intervaltree
        # pysam
        # networkx

tmux new-session -s RagTag
tmux attach-session -t RagTag

cd ~

module load StdEnv/2020
module load python/3.10.2
module load minimap2/2.24
module load scipy-stack/2022a

python3 -m venv ~/RagTag
source ~/RagTag/bin/activate
# deactivate # to close

module load StdEnv/2020
module load python/3.10.2
module load minimap2/2.24
module load scipy-stack/2022a

pip install RagTag

pip list --local

# test
ragtag.py --version
#v2.1.0

#--------------------------------
cd ~/scratch/RagTag/
wget https://github.com/malonge/RagTag/archive/refs/heads/master.zip
unzip master.zip
chmod -R 755 .

#--------------------------------
salloc -c40 --time 2:55:00 --mem 190000m --account def-cronk

cd /home/celphin/scratch/RagTag/
source ~/RagTag/bin/activate

module load StdEnv/2020
module load python/3.10.2
module load minimap2/2.24
module load scipy-stack/2022a

pip install RagTag

cd ~/scratch/Dryas/Dryas_genomes/

#----------
# Try running
/home/celphin/scratch/RagTag/RagTag-master/ragtag.py \
scaffold Dry-octo-H0.fasta Dry-int.fasta \
-t 39 -u -o ./Dry-int_ragtag_output/

#more ./Dry-int_ragtag_output/ragtag.scaffold.stats
#placed_sequences        placed_bp       unplaced_sequences      unplaced_bp         gap_bp        gap_sequences
#29211                  176 587 019           40 915                46 139 852       2 918 700        29187

#----------
# Try running
/home/celphin/scratch/RagTag/RagTag-master/ragtag.py \
scaffold Dry-octo-H0.fasta Dry-drumm.fasta \
-t 39 -u -o ./Dry-drumm_ragtag_output/

#more ./Dry-drumm_ragtag_output/ragtag.scaffold.stats
#placed_sequences        placed_bp       unplaced_sequences      unplaced_bp        gap_bp         gap_sequences
#1047                  211 560 089       12 310                    13 986 971        103100          1031

#----------
# Try running
/home/celphin/scratch/RagTag/RagTag-master/ragtag.py \
scaffold Dry-octo-H0.fasta Dry-alask.fasta \
-t 39 -u -o ./Dry-alask_ragtag_output/

#more ./Dry-alask_ragtag_output/ragtag.scaffold.stats
#placed_sequences        placed_bp       unplaced_sequences      unplaced_bp       gap_bp          gap_sequences
#33344                  176 172 628       30718                    40 533 169        3332100            33321

#----------
# Try running
/home/celphin/scratch/RagTag/RagTag-master/ragtag.py \
scaffold Dry-octo-H0.fasta Dry-octo-H1.fasta \
-t 39 -u -o ./Dry-octo-H1_ragtag_output/

#more ./Dry-octo-H1_ragtag_output/ragtag.scaffold.stats
#placed_sequences        placed_bp       unplaced_sequences      unplaced_bp     gap_bp       gap_sequences
#3041                    230427780        1089                   25 532 077        303000          3030

#----------
# Try running
/home/celphin/scratch/RagTag/RagTag-master/ragtag.py \
scaffold Dry-octo-H0.fasta Dry-octo-H2.fasta \
-t 39 -u -o ./Dry-octo-H2_ragtag_output/

#more ./Dry-octo-H2_ragtag_output/ragtag.scaffold.stats
#placed_sequences           placed_bp       unplaced_sequences       unplaced_bp       gap_bp      gap_sequences
#27879                    176 802 386       58071                   61 762 932        2 785 700          27857

#----------
# Try running
/home/celphin/scratch/RagTag/RagTag-master/ragtag.py \
scaffold Dry-octo-H0.fasta Dryas_octopetala_H1.supercontigs.fa \
-t 39 -u -o ./Old_Dryas_octopetala_H1_ragtag_output/


######################################
# check stats
# https://github.com/sanger-pathogens/assembly-stats

cd ~
git clone https://github.com/sanger-pathogens/assembly-stats.git
cd assembly-stats
mkdir build
cd build
cmake -DINSTALL_DIR:PATH=/home/celphin/assembly-stats/ ..
make
make test
make install

#-----------------------
module load StdEnv/2020 
module load gcc/9.3.0
module load r-bundle-bioconductor/3.14
cd ~/assembly-stats/build

./assembly-stats ~/scratch/Dryas/Dryas_genomes/Dryas_octopetala_H1.supercontigs.fa
stats for /home/celphin/scratch/Dryas/Dryas_genomes/Dryas_octopetala_H1.supercontigs.fa
sum = 233557465, n = 179, ave = 1304790.31, largest = 19685959
N50 = 12105649, n = 8
N60 = 11145057, n = 10
N70 = 8228389, n = 12
N80 = 4503589, n = 16
N90 = 2104560, n = 24
N100 = 9458, n = 179
N_count = 134679
Gaps = 190

./assembly-stats ~/scratch/Dryas/Dryas_genomes/Old_Dryas_octopetala_H1_ragtag_output/ragtag.scaffold.fasta
sum = 233562465, n = 129, ave = 1810561.74, largest = 29810504
N50 = 26945421, n = 5
N60 = 26945421, n = 5
N70 = 22319059, n = 6
N80 = 20495270, n = 8
N90 = 19685959, n = 9
N100 = 9458, n = 129
N_count = 139679
Gaps = 227



./assembly-stats ~/scratch/Dryas/Dryas_genomes/Dry-int.fasta
stats for /home/celphin/scratch/Dryas/Dryas_genomes/Dry-int.fasta
sum = 222 726 871, n = 70126, ave = 3176.10, largest = 19 977 900
N50 = 15 434, n = 248
N60 = 6 657, n = 2629
N70 = 3 563, n = 7285
N80 = 1 785, n = 16226
N90 = 855, n = 34956
N100 = 500, n = 70126
N_count = 10732131
Gaps = 17052

./assembly-stats ~/scratch/Dryas/Dryas_genomes/Dry-drumm.fasta
stats for /home/celphin/scratch/Dryas/Dryas_genomes/Dry-drumm.fasta
sum = 225 547 060, n = 13357, ave = 16886.06, largest = 4 847 072
N50 = 931 783, n = 63
N60 = 712 713, n = 91
N70 = 504 074, n = 129
N80 = 324 189, n = 182
N90 = 160 079, n = 275
N100 = 200, n = 13357
N_count = 3855360
Gaps = 12620


./assembly-stats ~/scratch/Dryas/Dryas_genomes/Dry-alask.fasta
stats for /home/celphin/scratch/Dryas/Dryas_genomes/Dry-alask.fasta
sum = 216 705 797, n = 64062, ave = 3382.75, largest = 19 176 389
N50 = 13 449, n = 342
N60 = 6 196, n = 2875
N70 = 3 528, n = 7582
N80 = 1 925, n = 15994
N90 = 968, n = 32122
N100 = 492, n = 64062
N_count = 12343629
Gaps = 18414


./assembly-stats ~/scratch/Dryas/Dryas_genomes/Dry-octo-H0.fasta
stats for /home/celphin/scratch/Dryas/Dryas_genomes/Dry-octo-H0.fasta
sum = 230 084 053, n = 28, ave = 8217287.61, largest = 33 681 925
N50 = 26 579 288, n = 4
N60 = 22 578 703, n = 5
N70 = 21 929 824, n = 6
N80 = 21 730 711, n = 7
N90 = 21 219 003, n = 8
N100 = 2095, n = 28
N_count = 87600
Gaps = 438


./assembly-stats ~/scratch/Dryas/Dryas_genomes/Dry-octo-H1.fasta
stats for /home/celphin/scratch/Dryas/Dryas_genomes/Dry-octo-H1.fasta
sum = 255 959 857, n = 4130, ave = 61975.75, largest = 1 259 362
N50 = 134 855, n = 493
N60 = 94595, n = 720
N70 = 58265, n = 1062
N80 = 33343, n = 1655
N90 = 22267, n = 2583
N100 = 5860, n = 4130
N_count = 200
Gaps = 1


./assembly-stats ~/scratch/Dryas/Dryas_genomes/Dry-octo-H2.fasta
stats for /home/celphin/scratch/Dryas/Dryas_genomes/Dry-octo-H2.fasta
sum = 238565318, n = 85950, ave = 2775.63, largest = 20 084 902
N50 = 11 214, n = 761
N60 = 5243, n = 4018
N70 = 2770, n = 10359
N80 = 1386, n = 22777
N90 = 769, n = 46600
N100 = 500, n = 85950
N_count = 10501541
Gaps = 17561


#------------------
./assembly-stats ~/scratch/Dryas/Dryas_genomes/Dry-int_ragtag_output/ragtag.scaffold.fasta
stats for /home/celphin/scratch/Dryas/Dryas_genomes/Dry-int_ragtag_output/ragtag.scaffold.fasta
sum = 225645571, n = 40939, ave = 5511.75, largest = 28381354
N50 = 21599655, n = 5
N60 = 21410073, n = 6
N70 = 18488805, n = 7
N80 = 12375, n = 82
N90 = 1212, n = 8781
N100 = 500, n = 40939
N_count = 13650831
Gaps = 46239


./assembly-stats ~/scratch/Dryas/Dryas_genomes/Dry-drumm_ragtag_output/ragtag.scaffold.fasta
stats for /home/celphin/scratch/Dryas/Dryas_genomes/Dry-drumm_ragtag_output/ragtag.scaffold.fasta
sum = 225650160, n = 12326, ave = 18306.84, largest = 32137890
N50 = 21834029, n = 5
N60 = 20569620, n = 6
N70 = 20106494, n = 7
N80 = 19617937, n = 8
N90 = 19018337, n = 9
N100 = 200, n = 12326
N_count = 3958460
Gaps = 13651


./assembly-stats ~/scratch/Dryas/Dryas_genomes/Dry-alask_ragtag_output/ragtag.scaffold.fasta
stats for /home/celphin/scratch/Dryas/Dryas_genomes/Dry-alask_ragtag_output/ragtag.scaffold.fasta
sum = 220037897, n = 30741, ave = 7157.80, largest = 27876026
N50 = 21923467, n = 5
N60 = 20333848, n = 6
N70 = 18775577, n = 7
N80 = 7463978, n = 9
N90 = 1919, n = 4852
N100 = 492, n = 30741
N_count = 15675729
Gaps = 51735


./assembly-stats ~/scratch/Dryas/Dryas_genomes/Dry-octo-H1_ragtag_output/ragtag.scaffold.fasta
stats for /home/celphin/scratch/Dryas/Dryas_genomes/Dry-octo-H1_ragtag_output/ragtag.scaffold.fasta
sum = 256262857, n = 1100, ave = 232966.23, largest = 34281279
N50 = 23221946, n = 5
N60 = 23098449, n = 6
N70 = 22619986, n = 7
N80 = 22085426, n = 8
N90 = 917941, n = 10
N100 = 7495, n = 1100
N_count = 303200
Gaps = 3031


./assembly-stats ~/scratch/Dryas/Dryas_genomes/Dry-octo-H2_ragtag_output/ragtag.scaffold.fasta
stats for /home/celphin/scratch/Dryas/Dryas_genomes/Dry-octo-H2_ragtag_output/ragtag.scaffold.fasta
sum = 241351018, n = 58093, ave = 4154.56, largest = 28248211
N50 = 21483562, n = 5
N60 = 21008324, n = 6
N70 = 7342160, n = 8
N80 = 2861, n = 2719
N90 = 883, n = 20140
N100 = 500, n = 58093
N_count = 13287241
Gaps = 45418



