###############################################
#January 2024:
    #Running stats on all assemblies
###############################################
#Install assembly stats:
#-------------------------------
# check stats of assemblies
# https://github.com/sanger-pathogens/assembly-stats

cd ~
git clone https://github.com/sanger-pathogens/assembly-stats.git
cd assembly-stats
mkdir build
cd build
cmake -DINSTALL_DIR:PATH=/home/msandler/assembly-stats/ ..
make
make test
make install
#----------------------------
module load StdEnv/2020 
module load gcc/9.3.0
module load r-bundle-bioconductor/3.14
cd ~/assembly-stats/

./assembly-stats /home/msandler/projects/def-rieseber/Dryas_shared_data/Oxyria/Polygonaceae_Genomes_Annotations/OTHER_Fagopyrum_tataricum/F_tataricum_H1.fasta

sum = 525965937, n = 1986, ave = 264836.83, largest = 62825076
N50 = 55820945, n = 5
N60 = 55018663, n = 6
N70 = 52827897, n = 7
N80 = 49984048, n = 8
N90 = 52097, n = 227
N100 = 12664, n = 1986
N_count = 9000
Gaps = 18


./assembly-stats /home/msandler/projects/def-rieseber/Dryas_shared_data/Oxyria/Polygonaceae_Genomes_Annotations/OTHER_Fagopyrum_tataricum/F_tataricum_H2.fasta

sum = 484958728, n = 856, ave = 566540.57, largest = 62576178
N50 = 55909694, n = 5
N60 = 55909694, n = 5
N70 = 54304797, n = 6
N80 = 52475509, n = 7
N90 = 43921284, n = 8
N100 = 12355, n = 856
N_count = 41000
Gaps = 82


./assembly-stats /home/msandler/projects/def-rieseber/Dryas_shared_data/Oxyria/Polygonaceae_Genomes_Annotations/Fagopyrum_escelentum/F_escelentum_H1.fasta

sum = 1263710397, n = 676, ave = 1869394.08, largest = 167747230
N50 = 156067145, n = 4
N60 = 155016467, n = 5
N70 = 153613506, n = 6
N80 = 139285584, n = 7
N90 = 136735094, n = 8
N100 = 4278, n = 676
N_count = 146000
Gaps = 292

./assembly-stats /home/msandler/projects/def-rieseber/Dryas_shared_data/Oxyria/Polygonaceae_Genomes_Annotations/Fagopyrum_escelentum/F_escelentum_H2.fasta

sum = 1209527273, n = 206, ave = 5871491.62, largest = 167865472
N50 = 150259016, n = 4
N60 = 143691915, n = 5
N70 = 142407128, n = 6
N80 = 140199158, n = 7
N90 = 134059295, n = 8
N100 = 5803, n = 206
N_count = 135000
Gaps = 270

./assembly-stats /home/msandler/projects/def-rieseber/Dryas_shared_data/Oxyria/Polygonaceae_Genomes_Annotations/MS_CE_Fagopyrum_tataricum/Fagopyrum_tataricum_Main.fasta

sum = 505882852, n = 7020, ave = 72063.08, largest = 68031765
N50 = 53883329, n = 5
N60 = 52287906, n = 6
N70 = 51545819, n = 7
N80 = 49982843, n = 8
N90 = 231255, n = 19
N100 = 1000, n = 7020
N_count = 16633491
Gaps = 2456

./assembly-stats /home/msandler/projects/def-rieseber/Dryas_shared_data/Oxyria/Polygonaceae_Genomes_Annotations/NCBI_Oxyria_digyna/Oxyria_digyna.fasta

sum = 561227150, n = 357, ave = 1572064.85, largest = 79472951
N50 = 36875860, n = 6
N60 = 34330788, n = 7
N70 = 22272060, n = 9
N80 = 17901664, n = 12
N90 = 10794573, n = 16
N100 = 17254, n = 357
N_count = 0
Gaps = 0

./assembly-stats ~/projects/def-rieseber/Dryas_shared_data/Oxyria/Polygonaceae_Genomes_Annotations/Rheum_nobile/R_nobile.fasta

sum = 1503241085, n = 11, ave = 136658280.45, largest = 162681177
N50 = 137450686, n = 6
N60 = 135503078, n = 7
N70 = 133562642, n = 8
N80 = 130015115, n = 9
N90 = 126994042, n = 10
N100 = 113113601, n = 11
N_count = 208000
Gaps = 416

./assembly-stats ~/projects/def-rieseber/Dryas_shared_data/Oxyria/Polygonaceae_Genomes_Annotations/Rheum_tangaticum/R_tangaticum.fasta
sum = 2710931324, n = 11, ave = 246448302.18, largest = 303179193
N50 = 245343095, n = 6
N60 = 242862938, n = 7
N70 = 240708114, n = 8
N80 = 240013100, n = 9
N90 = 239847201, n = 10
N100 = 194260622, n = 11
N_count = 1139165
Gaps = 17920

./assembly-stats /home/msandler/projects/def-rieseber/Dryas_shared_data/Oxyria/Polygonaceae_Genomes_Annotations/Polygonum_aviculare/Polavi_Main.fasta 

sum = 352071428, n = 88, ave = 4000811.68, largest = 42788684
N50 = 33827450, n = 5
N60 = 32556271, n = 6
N70 = 31959286, n = 7
N80 = 31043210, n = 8
N90 = 28844451, n = 10
N100 = 23458, n = 88
N_count = 3198
Gaps = 54

./assembly-stats ~/projects/def-rieseber/Dryas_shared_data/Oxyria/Polygonaceae_Genomes_Annotations/Rumex_hastalus/Rumex_genome.fasta
sum = 1647029117, n = 17152, ave = 96025.48, largest = 344501033
N50 = 158235527, n = 4
N60 = 150393790, n = 5
N70 = 412064, n = 65
N80 = 121687, n = 960
N90 = 49169, n = 3114
N100 = 2, n = 17152
N_count = 2777700
Gaps = 27777
