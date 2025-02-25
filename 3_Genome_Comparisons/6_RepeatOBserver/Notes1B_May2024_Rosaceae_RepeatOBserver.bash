############################################
# Run the new Dryas genome
# May 2024
# Cassandra Elphinstone
# D. octo (chr) https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_963921425.1/
	# Hap https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_963921435.1/
	# Hap https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_025631235.1/
# D. drumm (sca) https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_003254865.1/
# D. int (sca) https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_028570905.1/
# D. alas (sca) https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_026122645.1/

###########################################

cd /home/celphin/scratch/repeats/auto_script/

# download the Setup_Run_Repeats.sh from github
wget https://raw.githubusercontent.com/celphin/RepeatOBserverV1/main/Setup_Run_Repeats.sh
chmod +x Setup_Run_Repeats.sh
dos2unix Setup_Run_Repeats.sh

# install the R package
# already done here
install.packages("devtools")
library(devtools)
install_github("celphin/RepeatOBserverV1") #to install the package
# Select 1:All to install all the required packages
library(RepeatOBserverV1) # to load the package

#------------------------------
# go to FTP directory for GenBank assembly
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/963/921/425/GCA_963921425.1_drDryOcto1.1/GCA_963921425.1_drDryOcto1.1_genomic.fna.gz
# Unzip
gunzip GCA_963921425.1_drDryOcto1.1_genomic.fna.gz
# rename
mv GCA_963921425.1_drDryOcto1.1_genomic.fna DryOcto.fasta

#---------------------------------
# run
cat << EOF > Auto_DryOcto_H0.sh
#!/bin/bash
#SBATCH --account=def-cronk
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=191000M

module load StdEnv/2020
module load seqkit/2.3.1
module load emboss/6.6.0
module load r/4.1.2

srun Setup_Run_Repeats.sh -i DryOcto -f DryOcto.fasta -h H0 -c 20 -m 191000M -g FALSE

EOF

sbatch Auto_DryOcto_H0.sh

sq

# "summary plots 2 starting"
# Error in checkForRemoteErrors(val) :
  # 9 nodes produced errors; first error: unused argument (filename = base::paste0(ofi, "bp", rangebp[1], "_", rangebp[2], "seq", rangeseq1[1], "_", rangeseq1[2], flaglog, ".pdf"))
# Calls: <Anonymous> ... clusterApply -> staticClusterApply -> checkForRemoteErrors
# Execution halted

#-----------------
# try again
#  load new version, filename and change to name in pdf and try again, try again with width and height


##########################
# Centromere positions

# histograms find most of them
# DryOcto_H0-AT Chr1 2.7e+07 29695000
# DryOcto_H0-AT Chr2 7e+06 29030000
# DryOcto_H0-AT Chr3 2.3e+07 24470000
# DryOcto_H0-AT Chr4 3e+06 22595000
# DryOcto_H0-AT Chr5 3e+06 18595000
# DryOcto_H0-AT Chr6 1.9e+07 17945000
# DryOcto_H0-AT Chr7 1.9e+07 17730000
# DryOcto_H0-AT Chr8 1.9e+07 17220000
# DryOcto_H0-AT Chr9 3e+06 16045000


# DryOcto_H0-AT_Chr1_cent500 26 457 501 3.4e+07 DryOcto_H0-AT Chr1 - matches spectra
# DryOcto_H0-AT_Chr2_cent500 7 282 501 33500000 DryOcto_H0-AT Chr2 - matches spectra
# DryOcto_H0-AT_Chr3_cent500 20 332 501 28500000 DryOcto_H0-AT Chr3 
# DryOcto_H0-AT_Chr4_cent500 332 501 2.7e+07 DryOcto_H0-AT Chr4 - misses centromere in spectra
# DryOcto_H0-AT_Chr5_cent500 17 501 2.3e+07 DryOcto_H0-AT Chr5
# DryOcto_H0-AT_Chr6_cent500 16 387 501 2.2e+07 DryOcto_H0-AT Chr6
# DryOcto_H0-AT_Chr7_cent500 16 942 501 2.2e+07 DryOcto_H0-AT Chr7
# DryOcto_H0-AT_Chr8_cent500 15 922 501 21500000 DryOcto_H0-AT Chr8
# DryOcto_H0-AT_Chr9_cent500 2 412 501 20500000 DryOcto_H0-AT Chr9


#########################
# General plot

cd /home/celphin/scratch/repeats/auto_script/output_chromosomes/DryOcto_H0-AT/Summary_output/

module load StdEnv/2020
module load seqkit/2.3.1
module load emboss/6.6.0
module load r/4.1.2

R

library(RepeatOBserverV1)
inpath="~/scratch/repeats/auto_script/input_chromosomes/DryOcto_H0-AT/chromosome_files/"
fname="DryOcto_H0-AT"
outpath="~/scratch/repeats/auto_script/output_chromosomes"

plot_all_chromosomes(fname=fname, inpath=inpath, outpath=outpath)



###############################
# search for centromeric repeats
cd /home/celphin/scratch/repeats/auto_script/EMBOSS_repeats

# To run this script EMBOSS6.6.0 is required. 
module load nixpkgs/16.09  intel/2018.3  emboss/6.6.0

# download the script
wget https://raw.githubusercontent.com/celphin/RepeatOBserverV1/main/repeat_seq_finder.sh
chmod +x repeat_seq_finder.sh
dos2unix repeat_seq_finder.sh

# uncomment seqkit load
nano repeat_seq_finder.sh

# direct the script to your chromosome files: "/home/celphin/scratch/repeats/auto_script/input_chromosomes/NerLuet_H0-AT/chromosome_files"
# Choose the Chromosome fasta file that you want to look at the repeat in: "NerLuet_H0-AT_Chr1part01.fasta"
# Define the bp range: 13200000 13500000
# Define the repeat length you are looking for: 4 (will work for repeats up to 600bp long)
# Define the number of bp per line in your fasta file: 60 (default for RepeatOBserverV1)

./repeat_seq_finder.sh "/home/celphin/scratch/repeats/auto_script/input_chromosomes/DryOcto_H0-AT/chromosome_files" \
"DryOcto_H0-AT_Chr3part01.fasta" 16000000 16500000 350 60 5000

#342bp repeat
cttaataagtcgagtattccaagattgtttcgacaaattcattgtacggcgtattgatattttttattaaaaagaaagattgaTTGATACGTATATGCATTGACTACAGGGattttaacaaagttattattaagaacaagtatcctttatcacgaattgatgatttgtttaactagttaaaggatgttttagtattttcaaagatttatcTAAGTTTAAGTTATcgttaagtctgagtggtctaaagagacattttgaagactacctttcaaactagatatgacattatgacttgtgttgatgtcgttttgggtttgcTAGTGTgccaacagtgtttatgaa
cttaagaaGTCAAGTATTCCAAGATTAttttgacaaattcattgtggggcgtgttgaactttttcatgaaaaagaaagatggatcgatacgtatgtgtattgattatagggatcttaacaaagttactattaagagcaagtatcctttaccacgaattgatgatttgtttaactagttaaagtatgttttagtattttcaaagattgatctaagattaagttatcattaagtctgagtggtcgagagaaacattttgaagactacctttcaaactagataaggcattatgagttgtgttgatgtcgttttgggttgactagtgtgctgacagtgtttatgag
cttaataagtcgagtattccaagattgtttgaacaaattcattgtgggccgtattgatcttttttatgaaaaagaaagatggatcgacacgtatgtgtattgaatacggggatcttaacaaagttactattaagaacaagtatcttttatcacgaattgatgatttgtttaactagttaaaggatgttttgatattttcaaagattgatccaagattaagttatcattaagtctgagtggtcgagagaaacattttgaagactaccgttaaaactagatatggcattatgagttgtgttgatgttgttttgggtttactagcgtaccgacagtgtttttgaa
cttaataagtcgagtattccaagattatttcgacaaattcattgtacggcgtattgatattttttattaaaaagaaagattgaTTGATACGTatatgtattgactacagggatcttaacaaagttattattaagaacaagtatcctttatcacgaattgatgatttgtttaactagttaaaggatgttttagtattttcaaagattgatctaagattaagttatcattaagtctgagtggtctagagagacattttgaagactacctttcaaactagatatggcattatgacttgtgttgatgtcgttttgggtttactagtgtgctgacagtgtttatgaa
cttaataagtcgagtattccaagattgcttcgacaaattcctTGTGGGGcttattgatcttttttatgaaaaaaaaagatggatcgatacgtatatgtattgactacagggatcttaacgaagttattattaagaacaagtatcctttatcacgaattgatgatttgtttaactagttaaatgatgttttggttttttcaaagatttatcttacattaagttatcattaagtctgagtggtccatagagacattttgaaaactacctttcaaactagatatggcattatgacttgtgttgatgtcattttgggtttactagtgtgccgacagtgtttatgaa


#-----------------------------
./repeat_seq_finder.sh "/home/celphin/scratch/repeats/auto_script/input_chromosomes/DryOcto_H0-AT/chromosome_files" \
"DryOcto_H0-AT_Chr7part01.fasta" 6500000 7200000 500 60 5000

# 452bp - 111 copies
gaggggtgcaacacgaggacttcccaggaggtcacccatcctagtactgctctcgcccaagcacacttaacttcggagttctgatgggatccggtgcattagtgctggtatggtcgcacccaccattcatcgcgcaaaatatttatttaaccctgcttcccgcatcctcgacgctcccattcgccttcctctccatttcaccacgcgatcccgaatcccacgggcaccggagaaaggtaaaacaaaacgtaaaacgcaccgacgtaaaatgaagggtgtaacacgaggactccccaccgttcattgcactttgaatgcacaaatgagtgagtatataatttttatatctcatcttgtcgttatgtgattttgcattaatagtcaatcacatgattttgaattaatgcaccggataagggtaaaacaaaacgtaaaacgcaccgacgcaaaac
gaggggtgcaacacgaggacttcccaggaggtcacccatcctagtactgctctcgcccaagcacacttaacttcggagttctgatgggatccggtgcattagtgctggtatggtcgcacccaccattcatcgcgcaaaatatttatttaaccctgcttcccgcatcctcgacgctcccattcgccttcctctccatttcaccacgcgatcccgaatcccacgggcaccggagaaaggtaaaacaaaacgtaaaacgcaccgacgtaaaatgaagggtgtaacacgaggactccccaccgttcattgcactttgaatgcacaaatgagtgagtatataatttttatatctcatcttgtcgttatgtgattttgcattaatagtcaatcacatgattttgaattaatgcaccggataagggtaaaacaaaacgtaaaacgcaccgacgcaaaac
gaggggtgcaacacgaggacttcccaggaggtcacccatcctagtactgctctcgcccaagcacacttaacttcggagttctgatgggatccggtgcattagtgctggtatggtcgcacccaccattcatcgcgcaaaatatttatttaaccctgcttcccgcatcctcgacgctcccattcgccttcctctccatttcaccacgcgatcccgaatcccacgggcaccggagaaaggtaaaacaaaacgtaaaacgcaccgacgtaaaatgaagggtgtaacacgaggactccccaccgttcattgcactttgaatgcacaaatgagtgagtatataatttttatatctcatcttgtcgttatgtgattttgcattaatagtcaatcacatgattttgaattaatgcaccggataagggtaaaacaaaacgtaaaacgcaccgacgcaaaac
gaggggtgcaacacgaggacttcccaggaggtcacccatcctagtactgctctcgcccaagcacacttaacttcggagttctgatgggatccggtgcattagtgctggtatggtcgcacccaccattcatcgcgcaaaatatttatttaaccctgcttcccgcatcctcgacgctcccattcgccttcctctccatttcaccacgcgatcccgaatcccacgggcaccggagaaaggtaaaacaaaacgtaaaacgcaccgacgtaaaatgaagggtgtaacacgaggactccccaccgttcattgcactttgaatgcacaaatgagtgagtatataatttttatatctcatcttgtcgttatgtgattttgcattaatagtcaatcacatgattttgaattaatgcaccggataagggtaaaacaaaacgtaaaacgcaccgacgcaaaac


######################################
# Running all Dryas genomes in RepeatOBserverV1


#########################################
# Run RepeatOBserver on Dint, Doct, Ddrum, and Dalas

cd /home/celphin/scratch/Dryas/Dryas_genomes/RepeatOBserver/

# download the script
wget https://raw.githubusercontent.com/celphin/RepeatOBserverV1/main/Setup_Run_Repeats.sh
chmod +x Setup_Run_Repeats.sh
dos2unix Setup_Run_Repeats.sh

#-------------------------
# copy over genomes
cd ~/scratch/Dryas/Dryas_genomes/RepeatOBserver
module load StdEnv/2020 seqkit/2.3.1
seqkit seq -m 1000000 ~/scratch/Dryas/Dryas_genomes/Dry-int_ragtag_output/ragtag.scaffold.fasta > Dry-int-chr.fasta
seqkit seq -m 1000000 ~/scratch/Dryas/Dryas_genomes/Dry-drumm_ragtag_output/ragtag.scaffold.fasta > Dry-drumm-chr.fasta
seqkit seq -m 1000000 ~/scratch/Dryas/Dryas_genomes/Dry-alask_ragtag_output/ragtag.scaffold.fasta > Dry-alask-chr.fasta
seqkit seq -m 1000000 ~/scratch/Dryas/Dryas_genomes/Dry-octo-H1_ragtag_output/ragtag.scaffold.fasta > Dry-octo-H1-chr.fasta
seqkit seq -m 1000000 ~/scratch/Dryas/Dryas_genomes/Dry-octo-H2_ragtag_output/ragtag.scaffold.fasta > Dry-octo-H2-chr.fasta
seqkit seq -m 1000000 ~/scratch/Dryas/Dryas_genomes/Dry-octo-H0.fasta > Dry-octo-H0-chr.fasta

#-----------------------------
# copy over the genomes into one folder

cd /home/celphin/scratch/Dryas/Dryas_genomes/RepeatOBserver
# start the script
cat << EOF > Auto_Dryint.sh
#!/bin/bash
#SBATCH --account=def-cronk
#SBATCH --time=5:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=191000M

module load StdEnv/2020
module load seqkit/2.3.1
module load emboss/6.6.0
module load r/4.1.2

srun Setup_Run_Repeats.sh -i Dryint -f Dry-int-chr.fasta -h H0 -c 19 -m 191000M -g FALSE

EOF

sbatch Auto_Dryint.sh

#------------------------- 
cat << EOF > Auto_Dryoct.sh
#!/bin/bash
#SBATCH --account=def-cronk
#SBATCH --time=5:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=191000M

module load StdEnv/2020
module load seqkit/2.3.1
module load emboss/6.6.0
module load r/4.1.2

srun Setup_Run_Repeats.sh -i Dryoct -f Dry-octo-H0-chr.fasta -h H0 -c 19 -m 191000M -g FALSE

EOF

sbatch Auto_Dryoct.sh

#------------------------- 
cat << EOF > Auto_Drydrum.sh
#!/bin/bash
#SBATCH --account=def-cronk
#SBATCH --time=5:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=191000M

module load StdEnv/2020
module load seqkit/2.3.1
module load emboss/6.6.0
module load r/4.1.2

srun Setup_Run_Repeats.sh -i Drydrum -f Dry-drumm-chr.fasta -h H0 -c 19 -m 191000M -g FALSE

EOF

sbatch Auto_Drydrum.sh

#------------------------- 
cat << EOF > Auto_Dryalask.sh
#!/bin/bash
#SBATCH --account=def-cronk
#SBATCH --time=5:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=191000M

module load StdEnv/2020
module load seqkit/2.3.1
module load emboss/6.6.0
module load r/4.1.2

srun Setup_Run_Repeats.sh -i Dryalask -f Dry-alask-chr.fasta -h H0 -c 19 -m 191000M -g FALSE

EOF

sbatch Auto_Dryalask.sh
















###########################################
# run RepeatOBserver on other Rosaceae genomes to compare centromeres

cd /home/celphin/scratch/Dryas/Dryas_genomes/RepeatOBserver

wget https://raw.githubusercontent.com/celphin/RepeatOBserverV1/main/Setup_Run_Repeats.sh
chmod +x Setup_Run_Repeats.sh
dos2unix Setup_Run_Repeats.sh

#------------------------
# download Roseaceae genomes 

# https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_916048215.2/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/916/048/215/GCF_916048215.2_drMalSylv7.2/GCF_916048215.2_drMalSylv7.2_genomic.fna.gz
gunzip GCF_916048215.2_drMalSylv7.2_genomic.fna.gz
mv   GCF_916048215.2_drMalSylv7.2_genomic.fna  Malussylvestris.fasta


# https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000346465.2/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/346/465/GCF_000346465.2_Prunus_persica_NCBIv2/GCF_000346465.2_Prunus_persica_NCBIv2_genomic.fna.gz
gunzip GCF_000346465.2_Prunus_persica_NCBIv2_genomic.fna.gz
mv   GCF_000346465.2_Prunus_persica_NCBIv2_genomic.fna  Prunuspersica.fasta


# https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000184155.1/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/184/155/GCF_000184155.1_FraVesHawaii_1.0/GCF_000184155.1_FraVesHawaii_1.0_genomic.fna.gz
gunzip GCF_000184155.1_FraVesHawaii_1.0_genomic.fna.gz
mv GCF_000184155.1_FraVesHawaii_1.0_genomic.fna Fragariavesca.fasta


# https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_958449725.1/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/958/449/725/GCF_958449725.1_drRosRugo1.1/GCF_958449725.1_drRosRugo1.1_genomic.fna.gz
gunzip GCF_958449725.1_drRosRugo1.1_genomic.fna.gz
mv  GCF_958449725.1_drRosRugo1.1_genomic.fna   Rosarugosa.fasta


# https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_933775445.1/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/933/775/445/GCF_933775445.1_drPotAnse1.1/GCF_933775445.1_drPotAnse1.1_genomic.fna.gz
gunzip GCF_933775445.1_drPotAnse1.1_genomic.fna.gz
mv   GCF_933775445.1_drPotAnse1.1_genomic.fna  Argentinaanserina.fasta

# https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_040183295.1/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/040/183/295/GCA_040183295.1_ASM4018329v1/GCA_040183295.1_ASM4018329v1_genomic.fna.gz
gunzip GCA_040183295.1_ASM4018329v1_genomic.fna.gz
mv GCA_040183295.1_ASM4018329v1_genomic.fna     Rubusargutus.fasta

# https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_030142095.1/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/030/142/095/GCA_030142095.1_RiMJ/GCA_030142095.1_RiMJ_genomic.fna.gz
gunzip GCA_030142095.1_RiMJ_genomic.fna.gz
mv  GCA_030142095.1_RiMJ_genomic.fna Rubusidaeus.fasta

#--------------------------
# start the script
cat << EOF > Auto_Argentinaanserina.sh
#!/bin/bash
#SBATCH --account=def-cronk
#SBATCH --time=5:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=191000M

module load StdEnv/2020
module load seqkit/2.3.1
module load emboss/6.6.0
module load r/4.1.2

srun Setup_Run_Repeats.sh -i Argentinaanserina -f Argentinaanserina.fasta -h H0 -c 19 -m 191000M -g FALSE

EOF

sbatch Auto_Argentinaanserina.sh

# 7 chromosomes

#--------------------------
# start the script
cat << EOF > Auto_Rosarugosa.sh
#!/bin/bash
#SBATCH --account=def-cronk
#SBATCH --time=5:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=191000M

module load StdEnv/2020
module load seqkit/2.3.1
module load emboss/6.6.0
module load r/4.1.2

srun Setup_Run_Repeats.sh -i Rosarugosa -f Rosarugosa.fasta -h H0 -c 19 -m 191000M -g FALSE

EOF

sbatch Auto_Rosarugosa.sh

# 7 chromosomes
#--------------------------
# start the script
cat << EOF > Auto_Fragariavesca.sh
#!/bin/bash
#SBATCH --account=def-cronk
#SBATCH --time=5:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=191000M

module load StdEnv/2020
module load seqkit/2.3.1
module load emboss/6.6.0
module load r/4.1.2

srun Setup_Run_Repeats.sh -i Fragariavesca -f Fragariavesca.fasta -h H0 -c 19 -m 191000M -g FALSE

EOF

sbatch Auto_Fragariavesca.sh

# 7 chromosomes

#--------------------------
# start the script
cat << EOF > Auto_Malussylvestris.sh
#!/bin/bash
#SBATCH --account=def-cronk
#SBATCH --time=5:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=191000M

module load StdEnv/2020
module load seqkit/2.3.1
module load emboss/6.6.0
module load r/4.1.2

srun Setup_Run_Repeats.sh -i Malussylvestris -f Malussylvestris.fasta -h H0 -c 19 -m 191000M -g FALSE

EOF

sbatch Auto_Malussylvestris.sh

# 17 chromosomes

#--------------------------
# start the script
cat << EOF > Auto_Prunuspersica.sh
#!/bin/bash
#SBATCH --account=def-cronk
#SBATCH --time=5:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=191000M

module load StdEnv/2020
module load seqkit/2.3.1
module load emboss/6.6.0
module load r/4.1.2

srun Setup_Run_Repeats.sh -i Prunuspersica -f Prunuspersica.fasta -h H0 -c 19 -m 191000M -g FALSE

EOF

sbatch Auto_Prunuspersica.sh

# 8 chromosomes

#--------------------------------
# start the script
cat << EOF > Auto_Rubusargutus.sh
#!/bin/bash
#SBATCH --account=def-cronk
#SBATCH --time=5:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=191000M

module load StdEnv/2020
module load seqkit/2.3.1
module load emboss/6.6.0
module load r/4.1.2

srun Setup_Run_Repeats.sh -i Rubusargutus -f Rubusargutus.fasta -h H0 -c 19 -m 191000M -g FALSE

EOF

sbatch Auto_Rubusargutus.sh

# 7 chromosomes

#--------------------------
# start the script
cat << EOF > Auto_Rubusidaeus.sh
#!/bin/bash
#SBATCH --account=def-cronk
#SBATCH --time=5:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=191000M

module load StdEnv/2020
module load seqkit/2.3.1
module load emboss/6.6.0
module load r/4.1.2

srun Setup_Run_Repeats.sh -i Rubusidaeus -f Rubusidaeus.fasta -h H0 -c 19 -m 191000M -g FALSE

EOF

sbatch Auto_Rubusidaeus.sh

# 7 chromosomes

#########################################
