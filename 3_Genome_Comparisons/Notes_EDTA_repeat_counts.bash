##############################
# EDTA -- TE counts in genomes
# July 2024
#############################

# in Graham 
# https://github.com/rieseberglab/Genome_assemblies_annotations/blob/main/Oxyria_Assembly_Annotations/2_Annotations/1_Annotation_Maker_Pipeline/Notes4_EDTA_File_Organization

head -n 25 /home/celphin/projects/def-rieseber/Dryas_shared_data/Oxyria/Polygonaceae_Genomes_Annotations/*/EDTA_Files/*.EDTA.TEanno.sum

#####################################

# ==> MS_CE_Fagopyrum_tataricum/EDTA_Files/F_tataricum.fasta.mod.EDTA.TEanno.sum <==
# Repeat Classes
# ==============
# Total Sequences: 1
# Total Length: 66284675 bp
# Class                  Count        bpMasked    %masked
# =====                  =====        ========     =======
# LTR                    --           --           --
    # Copia              378          182511       0.28%
    # Gypsy              416          221042       0.33%
# TIR                    --           --           --
    # CACTA              6            1141         0.00%
    # Mutator            29           7945         0.01%
    # Tc1_Mariner        2            877          0.00%
    # hAT                31           18947        0.03%
    # unknown            14           3824         0.01%
# nonTIR                 --           --           --
    # helitron           33427        19378241     29.23%
                      # ---------------------------------
    # total interspersed 34303        19814528     29.89%

# ---------------------------------------------------------
# Total                  34303        19814528     29.89%

# Repeat Stats
# ============

# ==> MS_CE_Oxyria_digyna/EDTA_Files/Oxyria.fasta.mod.EDTA.TEanno.sum <==
# Repeat Classes
# ==============
# Total Sequences: 1
# Total Length: 88144195 bp
# Class                  Count        bpMasked    %masked
# =====                  =====        ========     =======
# LTR                    --           --           --
    # Copia              6408         7456013      8.46%
    # Gypsy              9350         12686339     14.39%
    # unknown            5982         4230541      4.80%
# TIR                    --           --           --
    # CACTA              2363         1421296      1.61%
    # Mutator            8193         4033821      4.58%
    # PIF_Harbinger      756          288402       0.33%
    # Tc1_Mariner        6210         2042486      2.32%
    # hAT                2579         1274597      1.45%
# nonTIR                 --           --           --
    # helitron           25558        14125944     16.03%
                      # ---------------------------------
    # total interspersed 67399        47559439     53.96%

# ---------------------------------------------------------
# Total                  67399        47559439     53.96%

# Repeat Stats

# ==> NCBI_Oxyria_digyna/EDTA_Files/Oxyria_NCBI.fasta.mod.EDTA.TEanno.sum <==
# Repeat Classes
# ==============
# Total Sequences: 1
# Total Length: 79472951 bp
# Class                  Count        bpMasked    %masked
# =====                  =====        ========     =======
# LTR                    --           --           --
    # Copia              5907         6163333      7.76%
    # Gypsy              9909         12548396     15.79%
    # unknown            3640         2356054      2.96%
# TIR                    --           --           --
    # CACTA              2507         1492839      1.88%
    # Mutator            5123         2463171      3.10%
    # PIF_Harbinger      878          294192       0.37%
    # Tc1_Mariner        5380         2879077      3.62%
    # hAT                2035         908128       1.14%
# nonTIR                 --           --           --
    # helitron           21536        12268232     15.44%
                      # ---------------------------------
    # total interspersed 56915        41373422     52.06%

# ---------------------------------------------------------
# Total                  56915        41373422     52.06%

# Repeat Stats

# ==> Polygonum_aviculare/EDTA_Files/P_aviculare.fasta.mod.EDTA.TEanno.sum <==
# Repeat Classes
# ==============
# Total Sequences: 1
# Total Length: 42788138 bp
# Class                  Count        bpMasked    %masked
# =====                  =====        ========     =======
# LTR                    --           --           --
    # Copia              2581         4203648      9.82%
    # Gypsy              4449         5346770      12.50%
    # unknown            2693         2022869      4.73%
# TIR                    --           --           --
    # CACTA              359          167737       0.39%
    # Mutator            2431         678927       1.59%
    # PIF_Harbinger      697          202154       0.47%
    # Tc1_Mariner        998          396308       0.93%
    # hAT                2323         671231       1.57%
# nonTIR                 --           --           --
    # helitron           8212         2936589      6.86%
                      # ---------------------------------
    # total interspersed 24743        16626233     38.86%

# ---------------------------------------------------------
# Total                  24743        16626233     38.86%


##################################################
# look for TE files in other genomes
# Run EDTA on genomes
# in Beluga
#--------------------------------------------------------------------------------
# EDTA_individual: Look for TE elements on individual Fasta 
#--------------------------------------------------------------------------------
cd /home/celphin/scratch/Oxyria/EDTA
################################
# Polygonaceae

#----------------------------------
#Fagopyrum esculentum : 
#Haploid 1,2
#~/projects/def-rieseber/Dryas_shared_data/Oxyria/Polygonaceae_Genomes_Annotations/Fagopyrum_escelentum/
# https://figshare.com/articles/dataset/Fagopyrum_genome_data/21617562/2

#Haploid 1:
wget https://figshare.com/ndownloader/files/39865477
mv 39865477 F_escelentum_H1.fasta

#Haploid 2:
wget https://figshare.com/ndownloader/files/39866572
mv 39866572 F_escelentum_H2.fasta

#-----------------------------------------
#Rheum nobile: (not chromosome level)
#~/projects/def-rieseber/Dryas_shared_data/Oxyria/Polygonaceae_Genomes_Annotations/Rheum_nobile
#https://figshare.com/articles/dataset/Rheum_nobile_genome/19662933

#Rno.genome.Chr.fa.gz
wget https://figshare.com/ndownloader/files/38250711
mv 38250711 Rno.genome.Chr.fa.gz
gunzip Rno.genome.Chr.fa.gz
mv Rno.genome.Chr.fa R_nobile.fasta


#-------------------------------------------------
#Rheum tangaticum:
#https://figshare.com/articles/dataset/Rheum_tanguticum_Genome/19663062
#Download of proteins also in Notes5_Oxyria_Annotation_Rerun.bash
wget https://figshare.com/ndownloader/files/41810904
mv 41810904 R_tangaticum.fasta.gz
gunzip R_tangaticum.fasta.gz

#################################################
# Rosaceae

# download protein files from other genomes
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_019419815.1/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/019/419/815/GCF_019419815.1_Pyrus_bretschneideri_v1/GCF_019419815.1_Pyrus_bretschneideri_v1_genomic.fna.gz
gunzip GCF_019419815.1_Pyrus_bretschneideri_v1_genomic.fna.gz
mv GCF_019419815.1_Pyrus_bretschneideri_v1_genomic.fna  Pyrus_bretschneideri.fasta


# https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_916048215.2/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/916/048/215/GCF_916048215.2_drMalSylv7.2/GCF_916048215.2_drMalSylv7.2_genomic.fna.gz
gunzip GCF_916048215.2_drMalSylv7.2_genomic.fna.gz
mv   GCF_916048215.2_drMalSylv7.2_genomic.fna  Malus_sylvestris.fasta


# https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000346465.2/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/346/465/GCF_000346465.2_Prunus_persica_NCBIv2/GCF_000346465.2_Prunus_persica_NCBIv2_genomic.fna.gz
gunzip GCF_000346465.2_Prunus_persica_NCBIv2_genomic.fna.gz
mv   GCF_000346465.2_Prunus_persica_NCBIv2_genomic.fna  Prunus_persica.fasta


# https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000184155.1/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/184/155/GCF_000184155.1_FraVesHawaii_1.0/GCF_000184155.1_FraVesHawaii_1.0_genomic.fna.gz
gunzip GCF_000184155.1_FraVesHawaii_1.0_genomic.fna.gz
mv GCF_000184155.1_FraVesHawaii_1.0_genomic.fna Fragaria_vesca.fasta


# https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_958449725.1/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/958/449/725/GCF_958449725.1_drRosRugo1.1/GCF_958449725.1_drRosRugo1.1_genomic.fna.gz
gunzip GCF_958449725.1_drRosRugo1.1_genomic.fna.gz
mv  GCF_958449725.1_drRosRugo1.1_genomic.fna   Rosa_rugosa.fasta


# https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_933775445.1/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/933/775/445/GCF_933775445.1_drPotAnse1.1/GCF_933775445.1_drPotAnse1.1_genomic.fna.gz
gunzip GCF_933775445.1_drPotAnse1.1_genomic.fna.gz
mv   GCF_933775445.1_drPotAnse1.1_genomic.fna  Argentina_anserina.fasta



##############################################
# Brassicaceae


# https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_000733195.1/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/733/195/GCA_000733195.1_A_alpina_V4/GCA_000733195.1_A_alpina_V4_genomic.fna.gz
gunzip GCA_000733195.1_A_alpina_V4_genomic.fna.gz
mv GCA_000733195.1_A_alpina_V4_genomic.fna  Arabis_alpina.fasta


# https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000004255.2/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/004/255/GCF_000004255.2_v.1.0/GCF_000004255.2_v.1.0_genomic.fna.gz
gunzip GCF_000004255.2_v.1.0_genomic.fna.gz
mv GCF_000004255.2_v.1.0_genomic.fna Arabidopsis_lyrata.fasta


# https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001735.4/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.fna.gz
gunzip GCF_000001735.4_TAIR10.1_genomic.fna.gz
mv GCF_000001735.4_TAIR10.1_genomic.fna Arabidopsis_thaliana.fasta



# https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000375325.1/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/375/325/GCF_000375325.1_Caprub1_0/GCF_000375325.1_Caprub1_0_genomic.fna.gz
gunzip GCF_000375325.1_Caprub1_0_genomic.fna.gz
mv GCF_000375325.1_Caprub1_0_genomic.fna Capsella_rubella.fasta


# https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000695525.1/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/695/525/GCF_000695525.1_BOL/GCF_000695525.1_BOL_genomic.fna.gz
gunzip GCF_000695525.1_BOL_genomic.fna.gz
mv GCF_000695525.1_BOL_genomic.fna Brassica_oleracea.fasta



# https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_911865555.2/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/911/865/555/GCA_911865555.2_T_arvense_v2/GCA_911865555.2_T_arvense_v2_genomic.fna.gz 
gunzip GCA_911865555.2_T_arvense_v2_genomic.fna.gz 
mv GCA_911865555.2_T_arvense_v2_genomic.fna Thlaspi_arvense.fasta


mkdir Draba_nivalis; cd Draba_nivalis
# https://datadryad.org/stash/dataset/doi:10.5061/dryad.pg4f4qrm4
# wget https://datadryad.org/stash/downloads/file_stream/403400

#####################################
# Fix repeatmasker dependancy 
cd ~/scratch/Annotation/RepeatMasker/RepeatMasker

module load StdEnv/2020
module load python/3.10.2
module load scipy-stack/2022a
module load gcc/9.3.0
module load cuda/11.7
module load hdf5/1.12.1
module load parasail/2.5 
module load  h5py/3.6.0

perl ./configure

# gives prompts below
# needs to be the full path not using the ~/scratch
The full path including the name for the TRF program.
TRF_PRGM:

/home/celphin/scratch/Annotation/RepeatMasker/TRF/build/src/trf

Add a Search Engine:
   1. Crossmatch: [ Un-configured ]
   2. RMBlast: [ Un-configured ]
   3. HMMER3.1 & DFAM: [ Un-configured ]
   4. ABBlast: [ Un-configured ]

   5. Done


Enter Selection:2

The path to the installation of the RMBLAST sequence alignment program.
RMBLAST_DIR: 

/home/celphin/scratch/Annotation/RepeatMasker/rmblast/bin/


#-------------------------
Add a Search Engine:
   1. Crossmatch: [ Un-configured ]
   2. RMBlast: [ Configured, Default ]
   3. HMMER3.1 & DFAM: [ Un-configured ]
   4. ABBlast: [ Un-configured ]

   5. Done
Enter Selection: 5
Building FASTA version of RepeatMasker.lib .....


###########################
# Update EDTA
git clone https://github.com/oushujun/EDTA.git
perl ./EDTA/EDTA/pl


########################
nano EDTA_script.sh

#!/bin/bash
#SBATCH --account=def-rieseber
#SBATCH --time=0-20:00:00
#SBATCH --cpus-per-task=40
#SBATCH --mem=191000M

SPP_Hap=$1

cd /home/celphin/scratch/Oxyria/EDTA/
module load StdEnv/2020 intel/2020.1.217 metagenome-atlas/2.5.0
set +u
source ~/miniconda2/bin/activate EDTA
set -u

#whole run
 EDTA.pl --genome /home/celphin/scratch/Oxyria/EDTA/${SPP_Hap} --sensitive 0 --anno 1 --evaluate 0 --threads 40 --force 1 \
 --repeatmasker /home/celphin/scratch/Annotation/RepeatMasker/RepeatMasker

# finalize annotation
perl /home/celphin/scratch/Oxyria/EDTA/EDTA/EDTA.pl --genome /home/celphin/scratch/Oxyria/EDTA/${SPP_Hap} \
--sensitive 0 --step anno --evaluate 0 --threads 40 --overwrite 0 \
--repeatmasker /home/celphin/scratch/Annotation/RepeatMasker/RepeatMasker

conda deactivate

#---------------------
sbatch EDTA_script.sh Brassica_oleracea.fasta  
sbatch EDTA_script.sh Arabidopsis_lyrata.fasta    
sbatch EDTA_script.sh Fragaria_vesca.fasta        
sbatch EDTA_script.sh R_nobile.fasta
sbatch EDTA_script.sh Arabidopsis_thaliana.fasta  
sbatch EDTA_script.sh Capsella_rubella.fasta   
sbatch EDTA_script.sh Malus_sylvestris.fasta      
sbatch EDTA_script.sh Rosa_rugosa.fasta
sbatch EDTA_script.sh Arabis_alpina.fasta         
sbatch EDTA_script.sh F_escelentum_H1.fasta    
sbatch EDTA_script.sh Prunus_persica.fasta        
sbatch EDTA_script.sh R_tangaticum.fasta
sbatch EDTA_script.sh Argentina_anserina.fasta    
sbatch EDTA_script.sh F_escelentum_H2.fasta    
sbatch EDTA_script.sh Pyrus_bretschneideri.fasta  
sbatch EDTA_script.sh Thlaspi_arvense.fasta
sbatch EDTA_script.sh Fagopyrum_tataricum_Main.fasta
sbatch EDTA_script.sh Oxyria_digyna.fasta
sbatch EDTA_script.sh Oxyria_Main.fasta
sbatch EDTA_script.sh Draniv.fasta
sbatch EDTA_script.sh DryOcto.fasta 
sbatch EDTA_script.sh Polavi_Main.fasta


grep ">" Cochlearia_groenlandica_h1.fasta |head
awk '/^>/ {print ">Cgro_" ++count} !/^>/ {print}' Cochlearia_groenlandica_h1.fasta > Cochlearia_groenlandica_h1chr.fasta

sbatch EDTA_script.sh Cochlearia_groenlandica_h1chr.fasta


#-------------------
# copy over genomes
cd /home/celphin/scratch/Oxyria/EDTA/Dryas

cp /home/celphin/scratch/Annotation/liftoff/data/Dry-*fasta .

#----------------
# rename long sequnece IDs
grep ">" Dry-alask-chr1.fasta |head
awk '/^>/ {print ">Dalask_" ++count} !/^>/ {print}' Dry-alask-chr.fasta > Dry-alask-chr1.fasta

grep ">" Dry-int-chr1.fasta |head
awk '/^>/ {print ">Dint_" ++count} !/^>/ {print}' Dry-int-chr.fasta > Dry-int-chr1.fasta

grep ">" Dry-octo-H2-chr1.fasta |head
awk '/^>/ {print ">Dajan_" ++count} !/^>/ {print}' Dry-octo-H2-chr.fasta > Dry-octo-H2-chr1.fasta

grep ">" Dry-octo-H1-chr1.fasta |head
awk '/^>/ {print ">Doct_" ++count} !/^>/ {print}' Dry-octo-H1-chr.fasta > Dry-octo-H1-chr1.fasta


# worked
sbatch ../EDTA_script.sh Dryas/Dry-drumm-chr1.fasta
sbatch ../EDTA_script.sh Dryas/Dry-octo-H1-chr1.fasta

# add other Dryas spp
sbatch ../EDTA_script.sh Dryas/DoctH0_Main.fasta
sbatch ../EDTA_script.sh Dryas/Dry-alask-chr1.fasta
sbatch ../EDTA_script.sh Dryas/Dry-int-chr1.fasta
sbatch ../EDTA_script.sh Dryas/Dry-octo-H2-chr1.fasta


cd ..
sbatch EDTA_script.sh Cochlearia_groenlandica_h1chr.fasta
sbatch EDTA_script.sh Oxyria_Main.fasta
sbatch EDTA_script.sh Dryas_octopetala_H1.supercontigs.fa


#-------------------
# New annotate
# test

cd /home/celphin/scratch/Oxyria/EDTA/
module load StdEnv/2020 intel/2020.1.217 metagenome-atlas/2.5.0
set +u
source ~/miniconda2/bin/activate EDTA
set -u

perl /home/celphin/scratch/Oxyria/EDTA/EDTA/EDTA.pl --genome /home/celphin/scratch/Oxyria/EDTA/Oxyria_Main.fasta \
--sensitive 0 --step anno --evaluate 0 --threads 40 --overwrite 0 \
--repeatmasker /home/celphin/scratch/Annotation/RepeatMasker/RepeatMasker

# Sat Nov  2 17:33:09 EDT 2024    Dependency checking:
# Error: AnnoSINE is not found in the AnnoSINE path ./!

#-----------
# Run
sbatch EDTA_anno_script.sh Oxyria_Main.fasta
sbatch ../EDTA_anno_script.sh Dryas/DoctH0_Main.fasta


###########################
# Plot summary output

# https://github.com/oushujun/EDTA/issues/92

head -n 25 *mod.EDTA.TEanno.sum > Repeat_Stats.txt
more Repeat_Stats.txt

# make Arctic_plant_EDTA_results.csv on computer and upload
# make Arctic_Dryas_EDTA_results.csv on computer and upload

#-------------------------------
cd /home/celphin/scratch/Oxyria/EDTA/plots
cd /lustre04/scratch/celphin/Oxyria/EDTA/plots

tmux new-session -s EDTA
tmux attach-session -t EDTA

module load StdEnv/2023
module load r/4.4.0
 
R 

# check if ggplot2 is installed, if so, load it, 
# if not, install and load it
if("ggplot2" %in% rownames(installed.packages())){
    library(ggplot2)
} else {
    install.packages("ggplot2")
    library(ggplot2)
}

library(tidyverse)
#library(statebins)

# Make overall barplot of all species EDTA results

# https://datatricks.co.uk/multiple-bar-charts-in-r

data  <- read.csv("/lustre04/scratch/celphin/Oxyria/EDTA/plots/Arctic_plant_EDTA_results.txt",sep="\t")
Dryas_data  <- read.csv("/lustre04/scratch/celphin/Oxyria/EDTA/plots/Dryas_EDTA_results.txt",sep="\t")

# Reclassify Arctic and alpine into a single group
data$Combined_Arctic <- ifelse(data$Arctic %in% c("Arctic", "alpine"), "Arctic/Alpine", "non-Arctic")

# Ensure your factors are set correctly
data$Repeat <- as.factor(data$Repeat)
data$Combined_Arctic <- as.factor(data$Combined_Arctic)
data <- as.data.frame(data)
data$Arctic <- as.factor(data$Arctic)
data$Family<- as.factor(data$Family)
data$Spp<- as.factor(data$Spp)

head(data)
#------------------
# Check significance

# Running the Kruskal-Wallis test
kruskal_test_result <- kruskal.test(Count ~ Combined_Arctic, data = data)

# Viewing the results
print(kruskal_test_result)

# data:  Count by Combined_Arctic
# Kruskal-Wallis chi-squared = 2.9113, df = 1, p-value = 0.08796

#---------
# check TIR
TIR_data <- data[which(data$Repeat=="TIR"),]

# Running the Kruskal-Wallis test
kruskal_test_result <- kruskal.test(Count ~ Combined_Arctic, data = TIR_data)

# Viewing the results
print(kruskal_test_result)

# data:  Count by Combined_Arctic
# Kruskal-Wallis chi-squared = 20.427, df = 1, p-value = 6.194e-06

#---------
# check CACTA
CACTA_data <- data[which(data$TE_type=="CACTA"),]

# Running the Kruskal-Wallis test
kruskal_test_result <- kruskal.test(Count ~ Combined_Arctic, data = CACTA_data)

# Viewing the results
print(kruskal_test_result)

#Kruskal-Wallis chi-squared = 4.0404, df = 1, p-value = 0.04442

#---------
# check hAT
hAT_data <- data[which(data$TE_type=="hAT"),]

# Running the Kruskal-Wallis test
kruskal_test_result <- kruskal.test(Count ~ Combined_Arctic, data = hAT_data)

# Viewing the results
print(kruskal_test_result)

# Kruskal-Wallis chi-squared = 4.8889, df = 1, p-value = 0.02703

#---------
# check Tc1_Mariner
Tc1_Mariner_data <- data[which(data$TE_type=="Tc1_Mariner"),]

# Running the Kruskal-Wallis test
kruskal_test_result <- kruskal.test(Count ~ Combined_Arctic, data = Tc1_Mariner_data)

# Viewing the results
print(kruskal_test_result)

# Kruskal-Wallis chi-squared = 4.0404, df = 1, p-value = 0.04442

# Mutator: Kruskal-Wallis chi-squared = 3.6465, df = 1, p-value = 0.05619


#------------------
# Plotting
# Custom theme
theme <- theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               axis.line.x = element_line(color="#919191", size = 0.1),
               axis.line.y = element_line(color="#919191", size = 0.1),
               axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1),
               axis.ticks.x = element_blank()) # Optional: remove ticks if needed

# For the total results plot
grDevices::png(file = "/home/celphin/scratch/Oxyria/EDTA/plots/Arctic_plant_EDTA_results_total.png", width = 1000, height = 800)
gg <- ggplot(data) + 
    geom_bar(aes(x = Spp, y = Count, fill = TE_type), position = "stack", stat = "identity") + 
    facet_wrap(~ Arctic, scales = "free_x") + # Free x scales for each facet
    theme +
    scale_x_discrete(labels = function(x) ifelse(x %in% unique(data$Spp[data$Arctic == unique(data$Arctic)]), x, "")) # Only show relevant species
print(gg)
grDevices::dev.off()

# For the by family plot
grDevices::png(file = "/home/celphin/scratch/Oxyria/EDTA/plots/Arctic_plant_EDTA_results_byfamily.png", width = 1000, height = 800)
gg <- ggplot(data) + 
    geom_bar(aes(x = Spp, y = Count, fill = TE_type), position = "stack", stat = "identity") + 
    facet_wrap(~ Family, scales = "free_x") + # Free x scales for each facet
    theme +
    scale_x_discrete(labels = function(x) ifelse(x %in% unique(data$Spp[data$Family == unique(data$Family)]), x, "")) # Only show relevant species
print(gg)
grDevices::dev.off()

# For the by repeat type
# Create the boxplot with points
grDevices::png(file = "/home/celphin/scratch/Oxyria/EDTA/plots/Repeat_abundance_boxplot_with_points.png", width = 1000, height = 800)
gg <- ggplot(data) + 
    geom_boxplot(aes(x = Arctic, y = Count, fill = Arctic), outlier.size = 1, alpha = 0.6) + 
    geom_jitter(aes(x = Arctic, y = Count), color = "black", width = 0.2, size = 1, alpha = 0.7) + 
    facet_wrap(~ Repeat, scales = "free_y") + # Free y scales for each facet
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line.x = element_line(color="#919191", size = 0.1),
          axis.line.y = element_line(color="#919191", size = 0.1),
          axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5)) +
    labs(x = "Habitat Type", y = "Abundance (Count)") +
    scale_fill_brewer(palette = "Set2") # Change color palette as needed

print(gg)
grDevices::dev.off()

# Create the boxplot with points for each specific type of repeat
grDevices::png(file = "/home/celphin/scratch/Oxyria/EDTA/plots/TE_type_abundance_boxplot_with_points.png", width = 1000, height = 800)
gg <- ggplot(data) + 
    geom_boxplot(aes(x = Arctic, y = Count, fill = Arctic), outlier.size = 1, alpha = 0.6) + 
    geom_jitter(aes(x = Arctic, y = Count), color = "black", width = 0.2, size = 1, alpha = 0.7) + 
    facet_wrap(~ TE_type, scales = "free_y") + # Free y scales for each facet
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line.x = element_line(color="#919191", size = 0.1),
          axis.line.y = element_line(color="#919191", size = 0.1),
          axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5)) +
    labs(x = "Habitat Type", y = "Abundance (Count)") +
    scale_fill_brewer(palette = "Set2") # Change color palette as needed

print(gg)
grDevices::dev.off()

# Create the boxplot with points
grDevices::png(file = "/home/celphin/scratch/Oxyria/EDTA/plots/TE_type_abundance_boxplot_combined.png", width = 1000, height = 800)
gg <- ggplot(data) + 
    geom_boxplot(aes(x = Combined_Arctic, y = Count, fill = Combined_Arctic), outlier.size = 1, alpha = 0.6) + 
    geom_jitter(aes(x = Combined_Arctic, y = Count), color = "black", width = 0.2, size = 1, alpha = 0.7) + 
    facet_wrap(~ TE_type, scales = "free_y") + # Free y scales for each facet
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line.x = element_line(color="#919191", size = 0.1),
          axis.line.y = element_line(color="#919191", size = 0.1),
          axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5)) +
    labs(x = "Habitat Type", y = "Abundance (Count)") +
    scale_fill_brewer(palette = "Set2") # Change color palette as needed

print(gg)
grDevices::dev.off()

#----------------------------
# make pdf not png

# For the total results plot
grDevices::pdf(file = "/home/celphin/scratch/Oxyria/EDTA/plots/Arctic_plant_EDTA_results_total.pdf")
gg <- ggplot(data) + 
    geom_bar(aes(x = Spp, y = Count, fill = TE_type), position = "stack", stat = "identity") + 
    facet_wrap(~ Arctic, scales = "free_x") + # Free x scales for each facet
    theme +
    scale_x_discrete(labels = function(x) ifelse(x %in% unique(data$Spp[data$Arctic == unique(data$Arctic)]), x, "")) # Only show relevant species
print(gg)
grDevices::dev.off()

# For the by family plot
grDevices::pdf(file = "/home/celphin/scratch/Oxyria/EDTA/plots/Arctic_plant_EDTA_results_byfamily.pdf")
gg <- ggplot(data) + 
    geom_bar(aes(x = Spp, y = Count, fill = TE_type), position = "stack", stat = "identity") + 
    facet_wrap(~ Family, scales = "free_x") + # Free x scales for each facet
    theme +
    scale_x_discrete(labels = function(x) ifelse(x %in% unique(data$Spp[data$Family == unique(data$Family)]), x, "")) # Only show relevant species
print(gg)
grDevices::dev.off()

# For the by repeat type
# Create the boxplot with points
grDevices::pdf(file = "/home/celphin/scratch/Oxyria/EDTA/plots/Repeat_abundance_boxplot_with_points.pdf")
gg <- ggplot(data) + 
    geom_boxplot(aes(x = Arctic, y = Count, fill = Arctic), outlier.size = 1, alpha = 0.6) + 
    geom_jitter(aes(x = Arctic, y = Count), color = "black", width = 0.2, size = 1, alpha = 0.7) + 
    facet_wrap(~ Repeat, scales = "free_y") + # Free y scales for each facet
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line.x = element_line(color="#919191", size = 0.1),
          axis.line.y = element_line(color="#919191", size = 0.1),
          axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5)) +
    labs(x = "Habitat Type", y = "Abundance (Count)") +
    scale_fill_brewer(palette = "Set2") # Change color palette as needed

print(gg)
grDevices::dev.off()

# Create the boxplot with points for each specific type of repeat
grDevices::pdf(file = "/home/celphin/scratch/Oxyria/EDTA/plots/TE_type_abundance_boxplot_with_points.pdf")
gg <- ggplot(data) + 
    geom_boxplot(aes(x = Arctic, y = Count, fill = Arctic), outlier.size = 1, alpha = 0.6) + 
    geom_jitter(aes(x = Arctic, y = Count), color = "black", width = 0.2, size = 1, alpha = 0.7) + 
    facet_wrap(~ TE_type, scales = "free_y") + # Free y scales for each facet
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line.x = element_line(color="#919191", size = 0.1),
          axis.line.y = element_line(color="#919191", size = 0.1),
          axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5)) +
    labs(x = "Habitat Type", y = "Abundance (Count)") +
    scale_fill_brewer(palette = "Set2") # Change color palette as needed

print(gg)
grDevices::dev.off()

# Create the boxplot with points
grDevices::pdf(file = "/home/celphin/scratch/Oxyria/EDTA/plots/TE_type_abundance_boxplot_combined.pdf")
gg <- ggplot(data) + 
    geom_boxplot(aes(x = Combined_Arctic, y = Count, fill = Combined_Arctic), outlier.size = 1, alpha = 0.6) + 
    geom_jitter(aes(x = Combined_Arctic, y = Count), color = "black", width = 0.2, size = 1, alpha = 0.7) + 
    facet_wrap(~ TE_type, scales = "free_y") + # Free y scales for each facet
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line.x = element_line(color="#919191", size = 0.1),
          axis.line.y = element_line(color="#919191", size = 0.1),
          axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5)) +
    labs(x = "Habitat Type", y = "Abundance (Count)") +
    scale_fill_brewer(palette = "Set2") # Change color palette as needed

print(gg)
grDevices::dev.off()

############################
# Plot gene or repeat densities
# https://genviz.org/module-02-r/0002/03/03/ggplot2_exercises/

tmux new-session -s EDTA
tmux attach-session -t EDTA

cd /home/celphin/scratch/Oxyria/EDTA

# get chromosome lengths

# remove sequneces that are too short
module load StdEnv/2023 seqkit/2.5.1
seqkit seq -m 10000000 ./Dryas/Dry-octo-H2_DoctH0_Main.fasta > DryOcto_chr.fasta
seqkit seq -m 10000000 Oxyria_Main.fasta > Oxyria_digyna_chr.fasta

module load  StdEnv/2020 bioawk/1.0
bioawk -c fastx '{print $name "\t" length($seq)}' DryOcto_chr.fasta  > DryOcto_chr_sizes.txt
bioawk -c fastx '{print $name "\t" length($seq)}' Oxyria_digyna_chr.fasta > Oxyria_digyna_chr_sizes.txt

#-----------------------

module load StdEnv/2023
module load r/4.4.0
 
R

library(dplyr)
library(ggplot2) 
library(tidyverse)
#library(statebins)

# import a text file with gene positions
# Dryas
Dry_genes0 <- read.table("/home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/genomes/Dryas_octopetala/Dryas_octopetala.gff3",sep="\t",header=F)
Dry_chr_sizes0 <- read.table("/home/celphin/scratch/Oxyria/EDTA/DryOcto_chr_sizes.txt",sep="\t",header=F)
Dry_TE_repeats0 <- read.table("/home/celphin/scratch/Oxyria/EDTA/DoctH0_Main.fasta.mod.EDTA.TEanno.gff3",sep="\t",header=F)
Dry_wgd0 <- read.table("/home/celphin/scratch/Oxyria/DupGen_finder/output/7.wgd.pairs",sep="\t",header=T)
Dry_tandem0 <- read.table("/home/celphin/scratch/Oxyria/DupGen_finder/output/7.tandem.pairs",sep="\t",header=T)
Dry_proximal0 <- read.table("/home/celphin/scratch/Oxyria/DupGen_finder/output/7.proximal.pairs",sep="\t",header=T)
Dry_transposed0 <- read.table("/home/celphin/scratch/Oxyria/DupGen_finder/output/7.transposed.pairs",sep="\t",header=T)
Dry_dispersed0 <- read.table("/home/celphin/scratch/Oxyria/DupGen_finder/output/7.dispersed.pairs",sep="\t",header=T)

# Oxyria
Oxy_genes0 <- read.table("/home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/genomes/Oxyria_digyna_H1/Oxyria_digyna_H1.gff3",sep="\t",header=F)
Oxy_chr_sizes0 <- read.table("/home/celphin/scratch/Oxyria/EDTA/Oxyria_digyna_chr_sizes.txt",sep="\t",header=F)
Oxy_TE_repeats0 <- read.table("/home/celphin/scratch/Oxyria/EDTA/Oxyria_digyna.fasta.mod.EDTA.TEanno.gff3",sep="\t",header=F)
Oxy_wgd0 <- read.table("/home/celphin/scratch/Oxyria/DupGen_finder/output/12.wgd.pairs",sep="\t",header=T)
Oxy_tandem0 <- read.table("/home/celphin/scratch/Oxyria/DupGen_finder/output/12.tandem.pairs",sep="\t",header=T)
Oxy_proximal0 <- read.table("/home/celphin/scratch/Oxyria/DupGen_finder/output/12.proximal.pairs",sep="\t",header=T)
Oxy_transposed0 <- read.table("/home/celphin/scratch/Oxyria/DupGen_finder/output/12.transposed.pairs",sep="\t",header=T)
Oxy_dispersed0 <- read.table("/home/celphin/scratch/Oxyria/DupGen_finder/output/12.dispersed.pairs",sep="\t",header=T)

# DupGen seq IDs
SequenceIDs <- read.table("/home/celphin/scratch/Oxyria/DupGen_finder/data/SequenceIDs.txt",sep=":",header=F)

# Interproscan data
Gene_ont_file <- "/home/celphin/scratch/Oxyria/synteny_quantity/Total_interproscan_output_edited3.tsv"
gene_ont <- read.delim(Gene_ont_file, header = TRUE, sep = "\t", na.strings = "-", colClasses = c("character", "character", "character", "character"))

#############################

Spp_genes0 <- Oxy_genes0
Spp_TE_repeats0 <- Oxy_TE_repeats0
Spp="Oxyria_digyna"


Spp_genes0 <- Dry_genes0
Spp_TE_repeats0 <- Dry_TE_repeats0
Spp="Dryas_octopetala_H0"

#-------------------------
# Plot genes
#---------------------------
scaffold_lengths <- Spp_genes0 %>%
  as.data.frame() %>%
  group_by(V1) %>%
  summarise(length = max(V5)) %>%
  ungroup()
  
threshold <- 1e7
long_scaffolds <- scaffold_lengths %>%
  filter(length > threshold) %>%
  pull(V1)

filtered_genes <- Spp_genes0[which(Spp_genes0$V1 %in% long_scaffolds & Spp_genes0$V3 == "gene"),]

# Edit so columns are: chr, position (no end or gene name required)
Spp_genes <- as.data.frame(cbind(filtered_genes$V1, filtered_genes$V4))
colnames(Spp_genes) <- c("chr", "pos")
Spp_genes$pos <- as.numeric(Spp_genes$pos)

# make a histogram plot of genes over the provided chromosomes 
plottedSppGenes <- ggplot(Spp_genes) + 
	geom_histogram(aes(x=pos),binwidth=1000000) + 
	facet_wrap(~chr,ncol=1) + 
	xlab("Genomic position (bins 1 Mb)") + 
	ggplot2::theme_classic() +
	ylab("Number of genes")

# save it to an image
png(paste0("./plots/", Spp, "_gene_density.png"),width=700,height=1500)
print(plottedSppGenes)
dev.off()

pdf(paste0("./plots/", Spp, "_gene_density.pdf"))
print(plottedSppGenes)
dev.off()

#---------------------------------
# run through all the repeat types
#--------------------------------
# Step 1: Filter out repeats on short scaffolds
scaffold_lengths <- Spp_TE_repeats0 %>%
  as.data.frame() %>%
  group_by(V1) %>%
  summarise(length = max(V5)) %>%
  ungroup()

threshold <- 1e7
long_scaffolds <- scaffold_lengths %>%
  filter(length > threshold) %>%
  pull(V1)

# Step 2: Get unique repeat types
unique_repeat_types <- unique(Spp_TE_repeats0$V3)

# Step 3: Loop through each repeat type and create plots
for (repeat_type in unique_repeat_types) {
  
  # Filter for the current repeat type
  filtered_repeats <- Spp_TE_repeats0 %>%
    filter(V1 %in% long_scaffolds & V3 == repeat_type)

  # Check if the filtered data is not empty
  if (nrow(filtered_repeats) > 0) {
    # Prepare the data for plotting
    Spp_repeats <- as.data.frame(cbind(filtered_repeats$V1, filtered_repeats$V4))
    colnames(Spp_repeats) <- c("chr", "pos")
    Spp_repeats$pos <- as.numeric(Spp_repeats$pos)

    # Create the histogram plot
    plottedSpp_repeats <- ggplot(Spp_repeats) + 
      geom_histogram(aes(x = pos), binwidth = 1000000) + 
      facet_wrap(~ chr, ncol = 1) + 
      xlab("Genomic position (bins 1 Mb)") + 
      theme_classic() +
      ylab("Number of repeats") +
      ggtitle(paste("Histogram of", repeat_type))

    # Step 4: Save the plot to an image file
    filename <- paste0("./plots/", Spp, "_repeats_", gsub(" ", "_", repeat_type), "_density.png")
	filename1 <- paste0("./plots/", Spp, "_repeats_", gsub(" ", "_", repeat_type), "_density.pdf")
    png(filename, width = 700, height = 1500)
    print(plottedSpp_repeats)
    dev.off()
	pdf(filename1)
	print(plottedSpp_repeats)
	dev.off()
  } else {
    message(paste("No data on more than 1 chromosome for repeat type:", repeat_type))
  }
}

# Oxyria - No data for repeat type: PIF_Harbinger_TIR_transposon

#----------------------------
# Gene duplicates
#--------------------------
Spp="Dryas_octopetala_H0"
wgddata <- Dry_wgd0
tanddata <- Dry_tandem0
proxdata <- Dry_proximal0
transdata <- Dry_transposed0 
dispdata <- Dry_dispersed0

Spp="Oxyria_digyna"
wgddata <- Oxy_wgd0
tanddata <- Oxy_tandem0
proxdata <- Oxy_proximal0
transdata <- Oxy_transposed0 
dispdata <- Oxy_dispersed0

#---------------------------
#WGD
 # Step 1: Split Location columns to extract chromosome and position
location_data <- wgddata %>%
  select(Location, Location.1) %>%
  pivot_longer(cols = everything(), names_to = "Type", values_to = "Location") %>%
  separate(Location, into = c("chr", "pos"), sep = ":", convert = TRUE)

scaffold_lengths <- location_data %>%
  as.data.frame() %>%
  group_by(chr) %>%
  summarise(length = max(pos)) %>%
  ungroup()

threshold <- 1e7
long_scaffolds <- scaffold_lengths %>%
  filter(length > threshold) %>%
  pull(chr)

  # Filter for the current repeat type
  filtered_data <- location_data %>%
    filter(chr %in% long_scaffolds)

location_data <- filtered_data

# Convert pos to numeric
location_data$pos <- as.numeric(location_data$pos)

# Step 2: Create the histogram plot
plotted_locations <- ggplot(location_data) + 
  geom_histogram(aes(x = pos), binwidth = 1000000) + 
  facet_wrap(~chr, ncol = 1) + 
  xlab("Genomic position (bins 1 Mb)") + 
  theme_classic() +
  ylab("Number of locations") +
  ggtitle("Histogram of Locations")

# Step 3: Save the plot to an image file
png(paste0("./plots/", Spp, "_wgd_location_histogram.png"), width = 700, height = 1500)
print(plotted_locations)
dev.off()
pdf(paste0("./plots/", Spp, "_wgd_location_histogram.pdf"))
print(plotted_locations)
dev.off()
#---------------------------
#Tandem
 # Step 1: Split Location columns to extract chromosome and position
location_data <- tanddata %>%
  select(Location, Location.1) %>%
  pivot_longer(cols = everything(), names_to = "Type", values_to = "Location") %>%
  separate(Location, into = c("chr", "pos"), sep = ":", convert = TRUE)

scaffold_lengths <- location_data %>%
  as.data.frame() %>%
  group_by(chr) %>%
  summarise(length = max(pos)) %>%
  ungroup()

threshold <- 1e7
long_scaffolds <- scaffold_lengths %>%
  filter(length > threshold) %>%
  pull(chr)

  # Filter for the current repeat type
  filtered_data <- location_data %>%
    filter(chr %in% long_scaffolds)

location_data <- filtered_data

# Convert pos to numeric
location_data$pos <- as.numeric(location_data$pos)

# Step 2: Create the histogram plot
plotted_locations <- ggplot(location_data) + 
  geom_histogram(aes(x = pos), binwidth = 1000000) + 
  facet_wrap(~chr, ncol = 1) + 
  xlab("Genomic position (bins 1 Mb)") + 
  theme_classic() +
  ylab("Number of locations") +
  ggtitle("Histogram of Locations")

# Step 3: Save the plot to an image file
png(paste0("./plots/", Spp, "_tandem_location_histogram.png"), width = 700, height = 1500)
print(plotted_locations)
dev.off()
pdf(paste0("./plots/", Spp, "_tandem_location_histogram.pdf"))
print(plotted_locations)
dev.off()
#---------------------------
# Proximal
 # Step 1: Split Location columns to extract chromosome and position
location_data <- proxdata %>%
  select(Location, Location.1) %>%
  pivot_longer(cols = everything(), names_to = "Type", values_to = "Location") %>%
  separate(Location, into = c("chr", "pos"), sep = ":", convert = TRUE)

scaffold_lengths <- location_data %>%
  as.data.frame() %>%
  group_by(chr) %>%
  summarise(length = max(pos)) %>%
  ungroup()

threshold <- 1e7
long_scaffolds <- scaffold_lengths %>%
  filter(length > threshold) %>%
  pull(chr)

  # Filter for the current repeat type
  filtered_data <- location_data %>%
    filter(chr %in% long_scaffolds)

location_data <- filtered_data

# Convert pos to numeric
location_data$pos <- as.numeric(location_data$pos)

# Step 2: Create the histogram plot
plotted_locations <- ggplot(location_data) + 
  geom_histogram(aes(x = pos), binwidth = 1000000) + 
  facet_wrap(~chr, ncol = 1) + 
  xlab("Genomic position (bins 1 Mb)") + 
  theme_classic() +
  ylab("Number of locations") +
  ggtitle("Histogram of Locations")

# Step 3: Save the plot to an image file
png(paste0("./plots/", Spp, "_proximal_location_histogram.png"), width = 700, height = 1500)
print(plotted_locations)
dev.off()
pdf(paste0("./plots/", Spp, "_proximal_location_histogram.pdf"))
print(plotted_locations)
dev.off()
#---------------------------
# Transposed
 # Step 1: Split Location columns to extract chromosome and position
location_data <- transdata %>%
  select(Location, Location.1) %>%
  pivot_longer(cols = everything(), names_to = "Type", values_to = "Location") %>%
  separate(Location, into = c("chr", "pos"), sep = ":", convert = TRUE)

scaffold_lengths <- location_data %>%
  as.data.frame() %>%
  group_by(chr) %>%
  summarise(length = max(pos)) %>%
  ungroup()

threshold <- 1e7
long_scaffolds <- scaffold_lengths %>%
  filter(length > threshold) %>%
  pull(chr)

  # Filter for the current repeat type
  filtered_data <- location_data %>%
    filter(chr %in% long_scaffolds)

location_data <- filtered_data

# Convert pos to numeric
location_data$pos <- as.numeric(location_data$pos)

# Step 2: Create the histogram plot
plotted_locations <- ggplot(location_data) + 
  geom_histogram(aes(x = pos), binwidth = 1000000) + 
  facet_wrap(~chr, ncol = 1) + 
  xlab("Genomic position (bins 1 Mb)") + 
  theme_classic() +
  ylab("Number of locations") +
  ggtitle("Histogram of Locations")

# Step 3: Save the plot to an image file
png(paste0("./plots/", Spp, "_transposed_location_histogram.png"), width = 700, height = 1500)
print(plotted_locations)
dev.off()
pdf(paste0("./plots/", Spp, "_transposed_location_histogram.pdf"))
print(plotted_locations)
dev.off()
#---------------------------
# Dispersed
 # Step 1: Split Location columns to extract chromosome and position
location_data <- dispdata %>%
  select(Location, Location.1) %>%
  pivot_longer(cols = everything(), names_to = "Type", values_to = "Location") %>%
  separate(Location, into = c("chr", "pos"), sep = ":", convert = TRUE)

scaffold_lengths <- location_data %>%
  as.data.frame() %>%
  group_by(chr) %>%
  summarise(length = max(pos)) %>%
  ungroup()

threshold <- 1e7
long_scaffolds <- scaffold_lengths %>%
  filter(length > threshold) %>%
  pull(chr)

  # Filter for the current repeat type
  filtered_data <- location_data %>%
    filter(chr %in% long_scaffolds)

location_data <- filtered_data

# Convert pos to numeric
location_data$pos <- as.numeric(location_data$pos)

# Step 2: Create the histogram plot
plotted_locations <- ggplot(location_data) + 
  geom_histogram(aes(x = pos), binwidth = 1000000) + 
  facet_wrap(~chr, ncol = 1) + 
  xlab("Genomic position (bins 1 Mb)") + 
  theme_classic() +
  ylab("Number of locations") +
  ggtitle("Histogram of Locations")

# Step 3: Save the plot to an image file
png(paste0("./plots/", Spp, "_dispersed_location_histogram.png"), width = 700, height = 1500)
print(plotted_locations)
dev.off()
pdf(paste0("./plots/", Spp, "_dispersed_location_histogram.pdf"))
print(plotted_locations)
dev.off()

#############################
# Join gene duplicates with Sequence IDs
colnames(SequenceIDs) <- c("Duplicate.1", "gene")
# join with Interproscan data
colnames(gene_ont) <- c("spp", "gene", "INTPRO", "descrip", "GOterm")
length(unique(gene_ont$INTPRO))
unique(gene_ont$spp)
# [1] "Arabis_alpina_interproscan_output.tsv"
# [2] "Cochlearia_groenlandica_interproscan_output.tsv"
# [3] "Draba_nivalis_interproscan_output.tsv"
# [4] "Dryas_octopetala_interproscan_output.tsv"
# [5] "Oxyria_digyna_H1_interproscan_output.tsv"
# [6] "Rheum_nobile_H0_interproscan_output.tsv"

#-----------------------

Spp_tand_genes<- dplyr::left_join(tanddata, SequenceIDs, by="Duplicate.1")
gene_ont_Spp <- gene_ont[which(gene_ont$spp=="Dryas_octopetala_interproscan_output.tsv"),]
gene_ont_Spp$gene <- as.factor(gene_ont_Spp$gene)
Spp_tand_genes$gene <- gsub(" ", "", Spp_tand_genes$gene)
Spp_tand_genes_ont <- dplyr::left_join(Spp_tand_genes, gene_ont_Spp, by="gene")

unique(Spp_tand_genes_ont$GOterm)
unique(Spp_tand_genes_ont$descrip)

# Count occurrences of each unique GOterm
go_counts <- table(Spp_tand_genes_ont$descrip)

# Convert to a data frame for easier viewing
go_counts_df <- as.data.frame(go_counts)

# Rename columns
colnames(go_counts_df) <- c("GOterm", "Count")

# Order by Count in descending order and get the top 10
top_go_counts <- go_counts_df %>%
  arrange(desc(Count)) %>%
  head(10)

# Print the top 10
print(top_go_counts)

#----------------------
Spp_wgd_genes<- dplyr::left_join(wgddata, SequenceIDs, by="Duplicate.1")
gene_ont_Spp <- gene_ont[which(gene_ont$spp=="Dryas_octopetala_interproscan_output.tsv"),]
gene_ont_Spp$gene <- as.factor(gene_ont_Spp$gene)
Spp_wgd_genes$gene <- gsub(" ", "", Spp_wgd_genes$gene)
Spp_wgd_genes_ont <- dplyr::left_join(Spp_wgd_genes, gene_ont_Spp, by="gene")

unique(Spp_wgd_genes_ont$GOterm)
unique(Spp_wgd_genes_ont$descrip)

# Count occurrences of each unique GOterm
go_counts <- table(Spp_wgd_genes_ont$descrip)

# Convert to a data frame for easier viewing
go_counts_df <- as.data.frame(go_counts)

# Rename columns
colnames(go_counts_df) <- c("GOterm", "Count")

# Order by Count in descending order and get the top 10
top_go_counts <- go_counts_df %>%
  arrange(desc(Count)) %>%
  head(10)

# Print the top 10
print(top_go_counts)
# DNA integration 
# https://www.ebi.ac.uk/QuickGO/term/GO:0015074
# Integrase,Integrase zinc-binding domain,Retrotransposon gag domain,Reverse transcriptase,Reverse transcriptase domain,Reverse transcriptase/Diguanylate cyclase domain,Ribonuclease H superfamily,Ribonuclease H-like superfamily


####################################
# make a density plot of genes over the provided chromosomes 

library(dplyr)
library(tidyr)
library(ggplot2)
library(grid)

# Prepare gene data
gene_data <- Spp_genes %>%
  mutate(Type = "Gene")

# Prepare repeat data
repeat_data_list <- lapply(unique_repeat_types, function(repeat_type) {
  filtered_repeats <- Spp_TE_repeats0 %>%
    filter(V3 == repeat_type) %>%
    select(chr = V1, pos = V4) %>%
    mutate(Type = repeat_type)
  return(filtered_repeats)
})

# Combine all repeat datasets
repeat_data <- bind_rows(repeat_data_list)

# Define the minimum chromosome length
min_length <- 1e7  # Set your desired minimum length here

# Combine gene and repeat data
combined_data <- bind_rows(
  gene_data,
  repeat_data
)

# Set up a plotting area for each chromosome and type
plot_data <- combined_data %>%
  group_by(chr) %>%
  summarise(Start = min(pos), End = max(pos), .groups = 'drop') %>%
  ungroup() %>%
  mutate(Length = End - Start) %>%
  filter(Length >= min_length)

# Filter combined data based on calculated chromosome lengths
filtered_combined_data <- combined_data %>%
  inner_join(plot_data, by = "chr")

# Create a data frame with y positions dynamically
n_types <- length(unique_types)  # Count unique types

# Generate y positions based on the number of unique types
y_positions <- data.frame(
  Type = unique_types,
  ymin = seq(0, (n_types - 1) * 0.1, by = 0.1),  # Each type spaced by 0.1
  ymax = seq(0.1, n_types * 0.1, by = 0.1)        # Each type spaced by 0.1
)

# Convert chr to a factor or numeric if it's not already
filtered_combined_data$chr <- as.factor(filtered_combined_data$chr)
plot_data$chr <- as.factor(plot_data$chr)

# Create plot
png("./plots/Total_gene_repeat_density_plot.png", width = 1500, height = 900)
print(ggplot() +
  geom_rect(data = plot_data,
            aes(xmin = Start, xmax = End, ymin = as.numeric(chr)-0.5, ymax = as.numeric(chr)+0.5),
            colour="black", fill = NA, alpha = 0.2) +

  # Adjust the y aesthetics to use factor levels for chr
  geom_segment(data = filtered_combined_data %>% filter(Type == "Gene"),
               aes(x = pos, xend = pos, 
                   y = as.numeric(chr)-0.5, yend = as.numeric(chr)-0.5+0.33),
               color = "tomato3", alpha = 0.05, size = 1) +

  geom_segment(data = filtered_combined_data %>% filter(Type == "LTR_retrotransposon"),
               aes(x = pos, xend = pos, 
                   y = as.numeric(chr)-0.5+0.33, yend = as.numeric(chr)-0.5 +0.67),
               color = "blue", alpha = 0.05, size = 1) +

  geom_segment(data = filtered_combined_data %>% filter(Type == "repeat_region"),
               aes(x = pos, xend = pos, 
                   y = as.numeric(chr)-0.5+0.67, yend = as.numeric(chr)+0.5),
               color = "green", alpha = 0.05, size = 1) +

  scale_y_reverse() +
  labs(x = "Genomic Position (Mbp)", y = "Chromosome") +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")
)
dev.off()

# Create plot
pdf("./plots/Total_gene_repeat_density_plot.pdf")
print(ggplot() +
  geom_rect(data = plot_data,
            aes(xmin = Start, xmax = End, ymin = as.numeric(chr)-0.5, ymax = as.numeric(chr)+0.5),
            colour="black", fill = NA, alpha = 0.2) +

  # Adjust the y aesthetics to use factor levels for chr
  geom_segment(data = filtered_combined_data %>% filter(Type == "Gene"),
               aes(x = pos, xend = pos, 
                   y = as.numeric(chr)-0.5, yend = as.numeric(chr)-0.5+0.33),
               color = "tomato3", alpha = 0.05, size = 0.1) +

  geom_segment(data = filtered_combined_data %>% filter(Type == "LTR_retrotransposon"),
               aes(x = pos, xend = pos, 
                   y = as.numeric(chr)-0.5+0.33, yend = as.numeric(chr)-0.5 +0.67),
               color = "blue", alpha = 0.05, size = 0.1) +

  geom_segment(data = filtered_combined_data %>% filter(Type == "repeat_region"),
               aes(x = pos, xend = pos, 
                   y = as.numeric(chr)-0.5+0.67, yend = as.numeric(chr)+0.5),
               color = "green", alpha = 0.05, size = 0.1) +

  scale_y_reverse() +
  labs(x = "Genomic Position (Mbp)", y = "Chromosome") +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")
)
dev.off()

###############################
# Make some similar plots of Interproscan data 

# join with chr and start Position
library(dplyr)
library(stringr)

Spp_genes0 <- Spp_genes0 %>%
  mutate(gene = str_extract(V9, "(?<=ID=)[^;]+"))

gene_ont_Spp0 <- left_join(gene_ont_Spp, Spp_genes0, by="gene")

geneont_data1 <- gene_ont_Spp0 %>%
  filter(grepl("oxidoreductase", descrip))%>%
  select(chr = V1, pos = V4) %>%
  mutate(Type = "oxidoreductase")

geneont_data2 <- gene_ont_Spp0 %>%
  filter(grepl("stress", descrip))%>%
  select(chr = V1, pos = V4) %>%
  mutate(Type = "stress")

geneont_data3 <- gene_ont_Spp0 %>%
  filter(grepl("Aquaporin", descrip))%>%
  select(chr = V1, pos = V4) %>%
  mutate(Type = "water")

# Define the minimum chromosome length
min_length <- 1e7  # Set your desired minimum length here

# Combine gene and repeat data
combined_data <- bind_rows(
  geneont_data1,
  geneont_data2,
  geneont_data3
)

unique_combined_data <- combined_data %>%
  distinct()

# Set up a plotting area for each chromosome and type
plot_data <- unique_combined_data %>%
  group_by(chr) %>%
  summarise(Start = min(pos), End = max(pos), .groups = 'drop') %>%
  ungroup() %>%
  mutate(Length = End - Start) %>%
  filter(Length >= min_length)

# Filter combined data based on calculated chromosome lengths
filtered_combined_data <- unique_combined_data %>%
  inner_join(plot_data, by = "chr")

# Create a data frame with y positions dynamically
n_types <- length(unique_types)  # Count unique types

# Generate y positions based on the number of unique types
y_positions <- data.frame(
  Type = unique_types,
  ymin = seq(0, (n_types - 1) * 0.1, by = 0.1),  # Each type spaced by 0.1
  ymax = seq(0.1, n_types * 0.1, by = 0.1)        # Each type spaced by 0.1
)

# Convert chr to a factor or numeric if it's not already
filtered_combined_data$chr <- as.factor(filtered_combined_data$chr)
plot_data$chr <- as.factor(plot_data$chr)

# Create plot
png("./plots/Water_temp_pos_plot.png", width = 1500, height = 900)
print(ggplot() +
  geom_rect(data = plot_data,
            aes(xmin = Start, xmax = End, ymin = as.numeric(chr)-0.5, ymax = as.numeric(chr)+0.5),
            colour="black", fill = NA, alpha = 0.2) +

  # Adjust the y aesthetics to use factor levels for chr
  geom_segment(data = filtered_combined_data %>% filter(Type == "stress"),
               aes(x = pos, xend = pos, 
                   y = as.numeric(chr)-0.5, yend = as.numeric(chr)-0.5+0.33),
               color = "tomato3", alpha = 0.5, size = 2) +

  geom_segment(data = filtered_combined_data %>% filter(Type == "water"),
               aes(x = pos, xend = pos, 
                   y = as.numeric(chr)-0.5+0.33, yend = as.numeric(chr)-0.5 +0.67),
               color = "blue", alpha = 0.5, size = 2) +

  geom_segment(data = filtered_combined_data %>% filter(Type == "oxidoreductase"),
               aes(x = pos, xend = pos, 
                   y = as.numeric(chr)-0.5+0.67, yend = as.numeric(chr)+0.5),
               color = "green", alpha = 0.5, size = 2) +

  scale_y_reverse() +
  labs(x = "Genomic Position (Mbp)", y = "Chromosome") +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")
)
dev.off()

# Create plot
pdf("./plots/Water_temp_pos_plot.pdf")
print(ggplot() +
  geom_rect(data = plot_data,
            aes(xmin = Start, xmax = End, ymin = as.numeric(chr)-0.5, ymax = as.numeric(chr)+0.5),
            colour="black", fill = NA, alpha = 0.2) +

  # Adjust the y aesthetics to use factor levels for chr
  geom_segment(data = filtered_combined_data %>% filter(Type == "stress"),
               aes(x = pos, xend = pos, 
                   y = as.numeric(chr)-0.5, yend = as.numeric(chr)-0.5+0.33),
               color = "tomato3", alpha = 0.5, size = 1) +

  geom_segment(data = filtered_combined_data %>% filter(Type == "water"),
               aes(x = pos, xend = pos, 
                   y = as.numeric(chr)-0.5+0.33, yend = as.numeric(chr)-0.5 +0.67),
               color = "blue", alpha = 0.5, size = 1) +

  geom_segment(data = filtered_combined_data %>% filter(Type == "oxidoreductase"),
               aes(x = pos, xend = pos, 
                   y = as.numeric(chr)-0.5+0.67, yend = as.numeric(chr)+0.5),
               color = "green", alpha = 0.5, size = 1) +

  scale_y_reverse() +
  labs(x = "Genomic Position (Mbp)", y = "Chromosome") +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")
)
dev.off()

#############################################

# join with chr and start Position
library(dplyr)
library(stringr)

Spp_genes0 <- Spp_genes0 %>%
  mutate(gene = str_extract(V9, "(?<=ID=)[^;]+"))

gene_ont_Spp0 <- left_join(gene_ont_Spp, Spp_genes0, by="gene")

geneont_data1 <- gene_ont_Spp0 %>%
  filter(grepl("methylation", descrip))%>%
  select(chr = V1, pos = V4) %>%
  mutate(Type = "methylation")

geneont_data2 <- gene_ont_Spp0 %>%
  filter(grepl("GO:0048658", GOterm))%>%
  select(chr = V1, pos = V4) %>%
  mutate(Type = "anther")

geneont_data3 <- gene_ont_Spp0 %>%
  filter(grepl("histone", descrip))%>%
  select(chr = V1, pos = V4) %>%
  mutate(Type = "histone")

# Define the minimum chromosome length
min_length <- 1e7  # Set your desired minimum length here

# Combine gene and repeat data
combined_data <- bind_rows(
  geneont_data1,
  geneont_data2,
  geneont_data3
)

unique_combined_data <- combined_data %>%
  distinct()

# Set up a plotting area for each chromosome and type
plot_data <- unique_combined_data %>%
  group_by(chr) %>%
  summarise(Start = 0, End = max(pos), .groups = 'drop') %>%
  ungroup() %>%
  mutate(Length = End - Start) %>%
  filter(Length >= min_length)

# Filter combined data based on calculated chromosome lengths
filtered_combined_data <- unique_combined_data %>%
  inner_join(plot_data, by = "chr")

# Create a data frame with y positions dynamically
n_types <- length(unique_types)  # Count unique types

# Generate y positions based on the number of unique types
y_positions <- data.frame(
  Type = unique_types,
  ymin = seq(0, (n_types - 1) * 0.1, by = 0.1),  # Each type spaced by 0.1
  ymax = seq(0.1, n_types * 0.1, by = 0.1)        # Each type spaced by 0.1
)

# Convert chr to a factor or numeric if it's not already
filtered_combined_data$chr <- as.factor(filtered_combined_data$chr)
plot_data$chr <- as.factor(plot_data$chr)

# Create plot
png("./plots/Epigenetic_plot.png", width = 1500, height = 900)
print(ggplot() +
  geom_rect(data = plot_data,
            aes(xmin = Start, xmax = End, ymin = as.numeric(chr)-0.5, ymax = as.numeric(chr)+0.5),
            colour="black", fill = NA, alpha = 0.2) +

  # Adjust the y aesthetics to use factor levels for chr
  geom_segment(data = filtered_combined_data %>% filter(Type == "methylation"),
               aes(x = pos, xend = pos, 
                   y = as.numeric(chr)-0.5, yend = as.numeric(chr)-0.5+0.33),
               color = "tomato3", alpha = 0.5, size = 2) +

  geom_segment(data = filtered_combined_data %>% filter(Type == "histone"),
               aes(x = pos, xend = pos, 
                   y = as.numeric(chr)-0.5+0.33, yend = as.numeric(chr)-0.5 +0.67),
               color = "blue", alpha = 0.5, size = 2) +

  geom_segment(data = filtered_combined_data %>% filter(Type == "anther"),
               aes(x = pos, xend = pos, 
                   y = as.numeric(chr)-0.5+0.67, yend = as.numeric(chr)+0.5),
               color = "green", alpha = 0.5, size = 2) +

  scale_y_reverse() +
  labs(x = "Genomic Position (Mbp)", y = "Chromosome") +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")
)
dev.off()

# Create plot
pdf("./plots/Water_temp_pos_plot.pdf")
print(ggplot() +
  geom_rect(data = plot_data,
            aes(xmin = Start, xmax = End, ymin = as.numeric(chr)-0.5, ymax = as.numeric(chr)+0.5),
            colour="black", fill = NA, alpha = 0.2) +

  # Adjust the y aesthetics to use factor levels for chr
  geom_segment(data = filtered_combined_data %>% filter(Type == "stress"),
               aes(x = pos, xend = pos, 
                   y = as.numeric(chr)-0.5, yend = as.numeric(chr)-0.5+0.33),
               color = "tomato3", alpha = 0.5, size = 1) +

  geom_segment(data = filtered_combined_data %>% filter(Type == "water"),
               aes(x = pos, xend = pos, 
                   y = as.numeric(chr)-0.5+0.33, yend = as.numeric(chr)-0.5 +0.67),
               color = "blue", alpha = 0.5, size = 1) +

  geom_segment(data = filtered_combined_data %>% filter(Type == "oxidoreductase"),
               aes(x = pos, xend = pos, 
                   y = as.numeric(chr)-0.5+0.67, yend = as.numeric(chr)+0.5),
               color = "green", alpha = 0.5, size = 1) +

  scale_y_reverse() +
  labs(x = "Genomic Position (Mbp)", y = "Chromosome") +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")
)
dev.off()






##########################################
# plotting Kimura distance for repeats

RepeatMasker -pa 2 -s -a -inv -dir ./RepMask -no_is -norna -xsmall -nolow -div 40 -lib EDTA.TElib.fa -cutoff 225 genome.fasta

calcDivergenceFromAlign.pl -s genome.divsum genome.fasta.align


#--------------------
# in R again

# install.packages("reshape")
# install.packages("ggplot2")
# install.packages("viridis")
# install.packages("hrbrthemes")
# install.packages("tidyverse")
# install.packages("gridExtra")

library(reshape)
library(ggplot2)
library(viridis)
library(hrbrthemes)
library(tidyverse)
library(gridExtra)

KimuraDistance <- read.csv("/home/celphin/scratch/Oxyria/EDTA/",sep=" ")

#add here the genome size in bp
genomes_size=230000000

kd_melt = melt(KimuraDistance,id="Div")
kd_melt$norm = kd_melt$value/genomes_size * 100

ggplot(kd_melt, aes(fill=variable, y=norm, x=Div)) + 
  geom_bar(position="stack", stat="identity",color="black") +
  scale_fill_viridis(discrete = T) +
  theme_classic() +
  xlab("Kimura substitution level") +
  ylab("Percent of the genome") + 
  labs(fill = "") +
  coord_cartesian(xlim = c(0, 55)) +
  theme(axis.text=element_text(size=11),axis.title =element_text(size=12))
  