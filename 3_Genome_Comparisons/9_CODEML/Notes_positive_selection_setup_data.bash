#########################
# Prep orthogroups data for CODEML or RELAX
# Jan 2026
############################

# follow instructions here
# https://github.com/siribi/CODEML_WORKFLOW
# https://github.com/siribi/RELAX-workflow

# other tutorial: https://bioinformaticsworkbook.org/dataAnalysis/ComparativeGenomics/Finding_Positive_Selection_With_Codeml.html#gsc.tab=0

###############################
# To run either program we need:

# Make a giant fasta file containing all proteins from every species
# Make a giant fasta file containing all transcripts from every species

# Alignments of genes (MAFFT) and then codon based alignment (pal2nal)

# A set of files with gene identifiers made from  the Orthogroups.txt file

# Extract the newick tree and gene trees
	# https://timetree.org/ 
	# OR use STAG tree


######################
# Transfer orthofinder data to Narval from Beluga

# Untar files
tmux new-session -s total
tmux attach-session -t total
cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/
tar -xzvf Total_Results_Aug18.tar.gz

# tmux new-session -s Brass
# tmux attach-session -t Brass
# cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Brassicaceae_genomes/orthofinder
# tar -xzvf Brassicaceae_Results_Oct09.tar.gz

# tmux new-session -s Rose
# tmux attach-session -t Rose
# cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Rosaceae_genomes/orthofinder
# tar -xzvf Rosaceae_Results_Aug16.tar.gz

# tmux new-session -s Poly
# tmux attach-session -t Poly
# cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Polygonaceae_genomes/orthofinder
# tar -xzvf Polygonaceae_Results_Oct07.tar.gz

#######################
# Check total data

cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Species_Tree

more SpeciesTree_rooted_node_labels.txt

# (((Fagopyrum_escelentum_H2:0.0377569,Fagopyrum_tataricum_H1:0.0354292)N3:0.10263
# 9,(Polygunum_aviculare_H0:0.156221,((Rheum_nobile_H0:0.035189,Rheum_tangaticum_H
# 0:0.0440871)N11:0.0634952,Oxyria_digyna_H1:0.136143)N7:0.0338668)N4:0.0351737)N1
# :0.093076,(((Rosa_rugosa:0.0478462,(Fragaria_vesca:0.0593422,Argentina_anserina:
# 0.076546)N12:0.0183635)N8:0.0789374,(Dryas_octopetala:0.142376,(Prunus_persica:0
# .0787626,(Malus_sylvestris:0.0169615,Pyrus_bretschneideri:0.0188564)N16:0.080969
# 9)N13:0.0316866)N9:0.0295128)N5:0.115343,(Cochlearia_groenlandica:0.121011,((((A
# rabidopsis_lyrata:0.022788,Arabidopsis_thaliana:0.0793045)N19:0.0265357,Capsella
# _rubella:0.0514998)N17:0.0309554,(Arabis_alpina:0.0831942,Draba_nivalis:0.085671
# 7)N18:0.0322681)N14:0.0134594,(Thlaspi_arvense:0.102948,Brassica_oleracea:0.0829
# 899)N15:0.0188365)N10:0.0265485)N6:0.207906)N2:0.093076)N0;


#####################
# Load modules
module load StdEnv/2020  gcc/9.3.0
module load cufflinks/2.2.1 
module load cdbfasta/0.99
module load diamond/2.0.4 
module load paml/4.9j
module load mafft/7.471

#module load orthofinder/2.5.2 
#module load clustalo/1.2.4
#module load pal2nal.v14 

#--------------
# After running orthofinder
# Make a giant fasta file containing all proteins from every species
cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/

cat 01_CleanProteins/*fasta >AllCleanProteins.fasta
ml cdbfasta; cdbfasta AllCleanProteins.fasta

#----------------
# need to download transcripts for all species
# Arabidopsis_lyrata       Dryas_octopetala         Prunus_persica
# Arabidopsis_thaliana     Fagopyrum_escelentum_H2  Pyrus_bretschneideri
# Arabis_alpina            Fagopyrum_tataricum_H1   Rheum_nobile_H0
# Argentina_anserina       Fragaria_vesca           Rheum_tangaticum_H0
# Brassica_oleracea        Malus_sylvestris         Rosa_rugosa
# Capsella_rubella         other                    Thlaspi_arvense
# Cochlearia_groenlandica  Oxyria_digyna_H1
# Draba_nivalis            Polygunum_aviculare_H0

mkdir transcripts

#----------------------
# Downloaded Arctic genomes

mkdir Cochlearia_groenlandica; cd Cochlearia_groenlandica
# https://springernature.figshare.com/collections/Whole-genome_sequencing_of_13_Arctic_plants_and_draft_genomes_of_Oxyria_digyna_and_Cochlearia_groenlandica/6965802/1
# transfer by globus
mv ../cg.h1.cds.fa Cochlearia_groenlandica.fna
cd ..

#----------------------------
# copy over Dryas annotation we made
mkdir Dryas_octopetala; cd Dryas_octopetala
cp /home/celphin/scratch/Annotation/CAP_Snakemake/DoctH0/Summary_data/DoctH0.AED_0.6_protein.fasta \
/home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/genomes/Dryas_octopetala/Dryas_octopetala.fa
cd ..
# GeneID: DoctH0_Chr100000001 

module load StdEnv/2020 gffread/0.12.3
gffread Dryas_octopetala.gff3 \
  -g DoctH0_Main.fasta \
  -x Dryas_octopetala.fna

module load StdEnv/2020 gffread/0.12.3
gffread Dryas_octopetala.gff3 \
  -g DoctH0_Main.fasta \
  -y Dryas_octopetala.fa

#-----------------------------
mkdir Draba_nivalis; cd Draba_nivalis
# https://datadryad.org/stash/dataset/doi:10.5061/dryad.pg4f4qrm4
wget Dniv87_Chicago_SSPACE_LINKS14_1kb_ChromonomerRun4_integrated_20July2018_unwrap.all.Run5.maker.genes.gff
wget Dniv87_Chicago_SSPACE_LINKS14_1kb_ChromonomerRun4_integrated_20July2018_unwrap.all.maker.transcripts.Run5.fasta
wget https://datadryad.org/stash/downloads/file_stream/403395
mv Dniv87_Chicago_SSPACE_LINKS14_1kb_ChromonomerRun4_integrated_20July2018_unwrap.all.Run5.maker.genes.gff Draba_nivalis.gff
mv Dniv87_Chicago_SSPACE_LINKS14_1kb_ChromonomerRun4_integrated_20July2018_unwrap.all.maker.transcripts.Run5.fasta Draba_nivalis.fa
cd ..

# make into protein file

cd /home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/genomes/Draba_nivalis

module load StdEnv/2020 python/3.10.2
source ~/gff3_env/bin/activate
gff3_to_fasta -g  Draba_nivalis.gff -f Draba_nivalis_genome.fasta -st pep -d simple -o Draba_nivalis



#-----------------------------
# Polygonaceae
cd cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/transcripts

cp /home/celphin/scratch/Oxyria/Polygonaceae_Genomes_Annotations/MS_CE_Fagopyrum_tataricum/Fagopyrum.AED_0.6.genes.fasta .
cp /home/celphin/scratch/Oxyria/Polygonaceae_Genomes_Annotations/MS_CE_Oxyria_digyna/Oxyria.AED_0.6.genes.fasta .
cp /home/celphin/scratch/Oxyria/Polygonaceae_Genomes_Annotations/Fagopyrum_escelentum/F_escelentum_H1.genes.fasta .
cp /home/celphin/scratch/Oxyria/Polygonaceae_Genomes_Annotations/Fagopyrum_escelentum/F_escelentum_H2.genes.fasta .
cp /home/celphin/scratch/Oxyria/Polygonaceae_Genomes_Annotations/OTHER_Fagopyrum_tataricum/F_tataricum_H1.genes.fasta .
cp /home/celphin/scratch/Oxyria/Polygonaceae_Genomes_Annotations/OTHER_Fagopyrum_tataricum/F_tataricum_H2.genes.fasta .
cp /home/celphin/scratch/Oxyria/Polygonaceae_Genomes_Annotations/Rheum_tangaticum/R_tangaticum.genes.fasta .

#-------------------------
# rename files

mkdir Rheum_tangaticum_H0; cd Rheum_tangaticum_H0
mv ../R_tangaticum.genes.cds.fasta  Rheum_tangaticum_H0.fna
cd ..

mkdir Fagopyrum_escelentum_H2; cd Fagopyrum_escelentum_H2
mv ../F_escelentum_H2.genes.cds.fasta  Fagopyrum_escelentum_H2.fna
cd ..

mkdir Fagopyrum_tataricum_H1; cd Fagopyrum_tataricum_H1
mv ../F_tataricum_H1.genes.cds.fasta  Fagopyrum_tataricum_H1.fna
cd ..

#---------------
# Build CDS files from genome
module load StdEnv/2020 gffread/0.12.3

cp /home/celphin/scratch/Oxyria/Polygonaceae_Genomes_Annotations/Polygonum_aviculare/Polavi_Main.fasta .
cp /home/celphin/scratch/Oxyria/Polygonaceae_Genomes_Annotations/Rheum_nobile/R_nobile.genes.fasta .
cp /home/celphin/scratch/Oxyria/Polygonaceae_Genomes_Annotations/NCBI_Oxyria_digyna/Oxyria_NCBI.AED_0.6.genes.fasta .

# Polygunum_aviculare_H0
mkdir Polygunum_aviculare_H0; cd Polygunum_aviculare_H0
mv Polavi_Main.fasta  Polygunum_aviculare_H0.fasta
gffread Polygunum_aviculare_H0.gff3 \
  -g Polygunum_aviculare_H0.fasta \
  -x Polygunum_aviculare_H0.fna
cd ..

# Rheum_nobile_H0
# cds>RnoChr01:263056-263398
# protein>RnoG0000001.1
# gff3 > RnoChr01 RnoG0000001.1

mkdir Rheum_nobile_H0; cd Rheum_nobile_H0
mv R_nobile.fasta  Rheum_nobile_H0.fasta
gffread Rheum_nobile_H0.gff3 \
  -g Rheum_nobile_H0.fasta \
  -x Rheum_nobile_H0.fna
cd ..

#-------------
# Oxyria_digyna_H1
# gff3 GeneID=Oxyria_NCBI_Chr100000001-RA  
# protein Oxyria_NCBI_Chr100000005
# Interpro: Oxyria_NCBI_Chr600003726

# gff3 Chr=Oxyrt-1-86582034
# Looking for Oxyria_ncbi.fasta
# /home/celphin/scratch/repeats/input_chromosomes/Oxyria/Oxyria_ragtag_output/ragtag.scaffold.fasta

#Reformat chromosomes: cannot find original file
module load StdEnv/2020 bioawk/1.0
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}'  Oxyria_ragtag.fasta  |\
awk -F '\t' '{printf("%d\t%s\n",length($2),$0);}' |\
sort -k1,1nr | cut -f 2- | tr "\t" "\n" > Oxyria_ragtag.scaffold_sorted.fasta
bioawk -c fastx '{ print ">Oxyrt-" ++i  "-" length($seq) "\n" $seq}'  < Oxyria_ragtag.scaffold_sorted.fasta > Oxyria_ncbi.fasta
awk 'BEGIN{FS="[> ]"} /^>/{val=$2;next}  {print val,length($0)}'  Oxyria_ncbi.fasta | head -n 23
# Oxyrt-1-86582034 86582034
# Oxyrt-2-79714091 79714091
# Oxyrt-3-79472951 79472951
# Oxyrt-4-78410798 78410798
# Oxyrt-5-76064323 76064323
# Oxyrt-6-73303751 73303751
# Oxyrt-7-72361354 72361354
# Oxyrt-8-755699 755699
# Oxyrt-9-538353 538353
# Oxyrt-10-287654 287654
# Oxyrt-11-284167 284167
# Oxyrt-12-166496 166496
# Oxyrt-13-127425 127425
# Oxyrt-14-124062 124062
# Oxyrt-15-123133 123133
# Oxyrt-16-117240 117240
# Oxyrt-17-111818 111818
# Oxyrt-18-110261 110261
# Oxyrt-19-109555 109555
# Oxyrt-20-105686 105686
# Oxyrt-21-99206 99206
# Oxyrt-22-90264 90264
# Oxyrt-23-85086 85086

grep ">" Oxyria_ncbi.fasta | head -n 20
# >Oxyrt-1-86582034
# >Oxyrt-2-79714091
# >Oxyrt-3-79472951
# >Oxyrt-4-78410798
# >Oxyrt-5-76064323
# >Oxyrt-6-73303751
# >Oxyrt-7-72361354
# >Oxyrt-8-755699
# >Oxyrt-9-538353
# >Oxyrt-10-287654
# >Oxyrt-11-284167
# >Oxyrt-12-166496
# >Oxyrt-13-127425
# >Oxyrt-14-124062
# >Oxyrt-15-123133
# >Oxyrt-16-117240
# >Oxyrt-17-111818
# >Oxyrt-18-110261
# >Oxyrt-19-109555
# >Oxyrt-20-105686


mv Oxyria_ncbi.fasta  Oxyria_digyna_H1.fasta
gffread Oxyria_digyna_H1.gff3 \
  -g Oxyria_digyna_H1.fasta \
  -x Oxyria_digyna_H1.fna

# now CDS: >Oxyria_NCBI_Chr100000001-RA

#-------------------
# Rosaceae

mkdir Pyrus_bretschneideri; cd Pyrus_bretschneideri
# download protein files from other genomes
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_019419815.1/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/019/419/815/GCF_019419815.1_Pyrus_bretschneideri_v1/GCF_019419815.1_Pyrus_bretschneideri_v1_cds_from_genomic.fna.gz
gunzip GCF_019419815.1_Pyrus_bretschneideri_v1_cds_from_genomic.fna.gz
mv GCF_019419815.1_Pyrus_bretschneideri_v1_cds_from_genomic.fna  Pyrus_bretschneideri.fna
cd ..

mkdir Malus_sylvestris; cd Malus_sylvestris
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_916048215.2/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/916/048/215/GCF_916048215.2_drMalSylv7.2/GCF_916048215.2_drMalSylv7.2_cds_from_genomic.fna.gz
gunzip GCF_916048215.2_drMalSylv7.2_cds_from_genomic.fna.gz
mv   GCF_916048215.2_drMalSylv7.2_cds_from_genomic.fna  Malus_sylvestris.fna
cd ..

mkdir Prunus_persica; cd Prunus_persica
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000346465.2/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/346/465/GCF_000346465.2_Prunus_persica_NCBIv2/GCF_000346465.2_Prunus_persica_NCBIv2_cds_from_genomic.fna.gz
gunzip GCF_000346465.2_Prunus_persica_NCBIv2_cds_from_genomic.fna.gz
mv   GCF_000346465.2_Prunus_persica_NCBIv2_cds_from_genomic.fna  Prunus_persica.fna
cd ..

mkdir Fragaria_vesca; cd Fragaria_vesca
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000184155.1/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/184/155/GCF_000184155.1_FraVesHawaii_1.0/GCF_000184155.1_FraVesHawaii_1.0_cds_from_genomic.fna.gz
gunzip GCF_000184155.1_FraVesHawaii_1.0_cds_from_genomic.fna.gz
mv GCF_000184155.1_FraVesHawaii_1.0_cds_from_genomic.fna Fragaria_vesca.fna
cd ..

mkdir Rosa_rugosa; cd Rosa_rugosa
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_958449725.1/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/958/449/725/GCF_958449725.1_drRosRugo1.1/GCF_958449725.1_drRosRugo1.1_cds_from_genomic.fna.gz
gunzip GCF_958449725.1_drRosRugo1.1_cds_from_genomic.fna.gz
mv  GCF_958449725.1_drRosRugo1.1_cds_from_genomic.fna   Rosa_rugosa.fna
cd ..

mkdir Argentina_anserina; cd Argentina_anserina
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_933775445.1/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/933/775/445/GCF_933775445.1_drPotAnse1.1/GCF_933775445.1_drPotAnse1.1_cds_from_genomic.fna.gz
gunzip GCF_933775445.1_drPotAnse1.1_cds_from_genomic.fna.gz
mv   GCF_933775445.1_drPotAnse1.1_cds_from_genomic.fna  Argentina_anserina.fna
cd ..

mkdir Rubus_argutus; cd Rubus_argutus
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_040183295.1/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/040/183/295/GCA_040183295.1_ASM4018329v1/GCA_040183295.1_ASM4018329v1_cds_from_genomic.fna.gz
gunzip GCA_040183295.1_ASM4018329v1_cds_from_genomic.fna.gz
mv GCA_040183295.1_ASM4018329v1_cds_from_genomic.fna     Rubus_argutus.fna
cd ..


#-------------------------
# Brassicaceae
mkdir Arabis_alpina; cd Arabis_alpina
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_000733195.1/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/733/195/GCA_000733195.1_A_alpina_V4/GCA_000733195.1_A_alpina_V4_cds_from_genomic.fna.gz
gunzip GCA_000733195.1_A_alpina_V4_cds_from_genomic.fna.gz
mv  GCA_000733195.1_A_alpina_V4_cds_from_genomic.fna Arabis_alpina.fna
cd ..

mkdir Arabidopsis_lyrata; cd Arabidopsis_lyrata # arctic species too
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000004255.2/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/004/255/GCF_000004255.2_v.1.0/GCF_000004255.2_v.1.0_cds_from_genomic.fna.gz
gunzip GCF_000004255.2_v.1.0_cds_from_genomic.fna.gz
mv GCF_000004255.2_v.1.0_cds_from_genomic.fna Arabidopsis_lyrata.fna
cd ..

mkdir Arabidopsis_thaliana; cd Arabidopsis_thaliana
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001735.4/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_cds_from_genomic.fna.gz
gunzip GCF_000001735.4_TAIR10.1_cds_from_genomic.fna.gz
mv GCF_000001735.4_TAIR10.1_cds_from_genomic.fna  Arabidopsis_thaliana.fna
cd ..

mkdir Capsella_rubella; cd Capsella_rubella
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000375325.1/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/375/325/GCF_000375325.1_Caprub1_0/GCF_000375325.1_Caprub1_0_cds_from_genomic.fna.gz
gunzip GCF_000375325.1_Caprub1_0_cds_from_genomic.fna.gz
mv GCF_000375325.1_Caprub1_0_cds_from_genomic.fna Capsella_rubella.fna
cd ..

mkdir Brassica_oleracea; cd Brassica_oleracea
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000695525.1/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/695/525/GCF_000695525.1_BOL/GCF_000695525.1_BOL_cds_from_genomic.fna.gz
gunzip GCF_000695525.1_BOL_cds_from_genomic.fna.gz
mv GCF_000695525.1_BOL_cds_from_genomic.fna Brassica_oleracea.fna
cd ..

mkdir Thlaspi_arvense; cd Thlaspi_arvense
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_911865555.2/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/911/865/555/GCA_911865555.2_T_arvense_v2/GCA_911865555.2_T_arvense_v2_cds_from_genomic.fna.gz
gunzip GCA_911865555.2_T_arvense_v2_cds_from_genomic.fna.gz
mv GCA_911865555.2_T_arvense_v2_cds_from_genomic.fna Thlaspi_arvense.fna
cd ..


######################
# rename gene names to match proteins
# missing Rubus_argutus, Draba_nivalis, Dryas_octopetala

# See R script - Notes_CDS_file_parsing.

# 21 genomes - and 4 Arctic

mv peptide/ CDS/

################################
# Load modules
module load StdEnv/2020  gcc/9.3.0
module load cufflinks/2.2.1 
module load cdbfasta/0.99
module load diamond/2.0.4 
module load paml/4.9j
module load mafft/7.471

#module load orthofinder/2.5.2 
#module load clustalo/1.2.4
#module load pal2nal.v14 

#----------------------
# Make a giant fasta file containing proteins from every species
cat /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/proteins/*fa >AllCleanProteins.fasta

# Make a giant fasta file containing all transcripts from every species
cat /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/CDS/*fa >AllCleanTranscripts.fasta

#---
# index files
cdbfasta AllCleanTranscripts.fasta
cdbfasta AllCleanProteins.fasta

#-------------
# check for duplicated names
grep "^>" AllCleanTranscripts.fasta | sort | uniq -d
grep "^>" AllCleanProteins.fasta | sort | uniq -d

# >accD
# >atpA
# >atpB
# >atpE
# >atpF
# >atpH
# >atpI
# >ccsA
# >cemA
# >clpP
# >matK
# >ndhA
# >ndhB
# >ndhC
# >ndhD
# >ndhE
# >ndhF
# >ndhG
# >ndhH
# >ndhI
# >ndhJ
# >ndhK
# >petA
# >petB
# >petD
# >petG
# >petL
# >petN
# >psaA
# >psaB
# >psaC
# >psaI
# >psaJ
# >psbA
# >psbB
# >psbC
# >psbD
# >psbE
# >psbF
# >psbH
# >psbI
# >psbJ
# >psbK
# >psbL
# >psbM
# >psbN
# >psbT
# >psbZ
# >rbcL
# >rpl14
# >rpl16
# >rpl2
# >rpl20
# >rpl22
# >rpl23
# >rpl32
# >rpl33
# >rpl36
# >rpoA
# >rpoB
# >rpoC1
# >rpoC2
# >rps11
# >rps12
# >rps14
# >rps15
# >rps16
# >rps18
# >rps19
# >rps2
# >rps3
# >rps4
# >rps7
# >rps8
# >ycf1
# >ycf2
# >ycf3
# >ycf4

# check for Dryas and Draba gene names
cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/
grep "DoctH0" AllCleanTranscripts.fasta | wc -l
grep "DoctH0" AllCleanProteins.fasta | wc -l
grep "maker-lg" AllCleanTranscripts.fasta | wc -l
grep "maker-lg" AllCleanProteins.fasta | wc -l
# 39696
# 39696
# 26653
# 26653

#-----------------------
# Explore tree files - maybe duplicate names are excluded
# Continue for now

grep ycf4 /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees/OG002*_tree.txt

#####################
# to clear and rerun 
cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees/

rm *_pal2nal
rm *_Proteins_alignment.fasta
rm *_Proteins.fasta
rm *_Transcripts.fasta

#--------------------
# narval2
tmux new-session -s total
tmux attach-session -t total

cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/

# Extract orthologue protein and transcript names by 
# looping through the tree files
for f in /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees/*tree.txt; do sed 's/(//g' $f |sed 's/,/\n/g' |sed 's/:.*//g' |sort|uniq|sed 's/Proteins_/\t/1'|cut -f 2  |cdbyank AllCleanProteins.fasta.cidx >${f}_Proteins.fasta;done

# check results 
cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees/
ls *_Proteins.fasta | wc -l
# 23145
ls *_tree.txt | wc -l
# 23145

grep "Doct" *_Proteins.fasta | wc -l
35790
grep "maker-lg" *_Proteins.fasta | wc -l
23626

#------------------
# Format to get gene names
cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/

for f in /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees/*tree.txt; do sed 's/(//g' $f |sed 's/,/\n/g' |sed 's/:.*//g' |sort|uniq|sed 's/Proteins_/\t/1'|cut -f 2  |cdbyank AllCleanTranscripts.fasta.cidx >${f}_Transcripts.fasta;done

cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees/
ls *_tree.txt_Transcripts.fasta | wc -l
#23145

#---------------
grep "Doct" *_tree.txt_Transcripts.fasta | wc -l
0

grep "Doct" AllCleanProteins.fasta |head
# >DoctH0_Chr100000001
# >DoctH0_Chr100000002
# >DoctH0_Chr100000003
# >DoctH0_Chr100000004
# >DoctH0_Chr100000005
# >DoctH0_Chr100000006
# >DoctH0_Chr100000009
# >DoctH0_Chr100000011
# >DoctH0_Chr100000012
# >DoctH0_Chr100000013
grep "Doct" AllCleanTranscripts.fasta | head
# >DoctH0_Chr100000001-RA
# >DoctH0_Chr100000002-RA
# >DoctH0_Chr100000003-RA
# >DoctH0_Chr100000004-RA
# >DoctH0_Chr100000005-RA
# >DoctH0_Chr100000006-RA
# >DoctH0_Chr100000009-RA
# >DoctH0_Chr100000011-RA
# >DoctH0_Chr100000012-RA
# >DoctH0_Chr100000013-RA

# Remove -RA from end of transcript names and rejoin
cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/CDS/
sed -i 's/\-RA//g' Dryas_octopetala.fa
sed -i 's/\-RA//g' Oxyria_digyna_H1.fa
sed -i 's/\-RA//g' Polygunum_aviculare_H0.fa

# rejoin and index
cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/
cat /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/CDS/*fa >AllCleanTranscripts.fasta
cdbfasta AllCleanTranscripts.fasta

# Check
grep "Doct" AllCleanTranscripts.fasta | head
grep "Polavi" AllCleanTranscripts.fasta | head
grep "Oxyria" AllCleanTranscripts.fasta | head

# rerun
cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/
for f in /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees/*tree.txt; do sed 's/(//g' $f |sed 's/,/\n/g' |sed 's/:.*//g' |sort|uniq|sed 's/Proteins_/\t/1'|cut -f 2  |cdbyank AllCleanTranscripts.fasta.cidx >${f}_Transcripts.fasta;done

# check again
cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees/
grep "Doct" *_tree.txt_Transcripts.fasta 
grep "Oxyria" *_tree.txt_Transcripts.fasta 
grep "Polavi" *_tree.txt_Transcripts.fasta 

#------------------------
# Generate protein alignments with MAFFT
tmux attach-session -t total

cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees/

for f in *_Proteins.fasta; do
  echo "mafft \"$f\" > \"${f%.*}_alignment.fasta\"" >> generateMAFFTAlignment.sh
done > generateMAFFTAlignment.sh
# add #!/bin/bash

chmod +x generateMAFFTAlignment.sh
./generateMAFFTAlignment.sh
# there is a --thread option if too slow

#------------------------
# Should we do the protein alignments with PRANK?
# More phylogenetically correct but more variable too
# https://gensoft.pasteur.fr/docs/prank/170427/
# https://ariloytynoja.github.io/prank-msa/

module load StdEnv/2020  gcc/9.3.0  prank/170427

for f in *_Proteins.fasta; do
  echo "prank -protein -d=\"$f\" -o=\"${f%.*}_alignment.fasta\""
done > generatePRANKAlignment.sh

#---------------------------
# Run pal2nal – convert alignments into paml format
# The suggestion is to use pal2nal with “-nogap”, though codeml has a way to deal with these.
# However, if you do not have sequence in your pal2nal output, likely you had gaps
# along the entire sequence across all of the species analyzed. 

cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/

wget http://www.bork.embl.de/pal2nal/distribution/pal2nal.v14.tar.gz
tar -zxvf pal2nal.v14.tar.gz

cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees

for f in *_Proteins_alignment.fasta; do
  base_name="${f%_Proteins_alignment.fasta}"  # Remove the "_Proteins_alignment.fasta" part
  echo "/home/celphin/scratch/Oxyria_Positive_Selection_Test/pal2nal.v14/pal2nal.pl \
  \"$f\" \"${base_name}_Transcripts.fasta\" -output paml -nogap > \"${base_name}_pal2nal\""
done > pal2nal.sh

# add #!/bin/bash

chmod +x pal2nal.sh

./pal2nal.sh

# ERROR: number of input seqs differ (aa: 1688;  nuc: 1381)!!
# ERROR: number of input seqs differ (aa: 1688;  nuc: 1248)!!
# Missing Dryas, Oxyria and Polyav nucleotides
# FIXED by removing the -RA in these gene names

# Also
#---  ERROR: inconsistency between the following pep and nuc seqs  ---#
# >RnoG0010531.1
# FSSRPLSCLLRRLDHALHHSALLASRLPSSPPSTAPRCLPPVCLQACLRDRCPSTTPSAF
# ETAAFDDSLLLRACLARVFSWKVTMMLRCLLYFMTRNLNWPLDLLSHHQKNWSWHEVLNQ
# REMQRTHYFQSQVMEVNPLLVGEFLSHLQEEEGYMLQNQVPLIRDLDEQLWRDIVDKKLM
# RTCLWMALLNMQNKNQRQRWNHVLYRILLQMGRIHIVRTMILLELLHLMKQFHERNDKKF
# IPNGLMTTNFMRLSSPSLIRRSLCSLV
# >RnoG0010531.1
# TTCTCCTCGCGTCCTCTATCTTGTCTCCTTCGGCGGCTCGATCACGCCCTCCACCACTCC
# GCGTTGCTCGCCTCCCGTTTGCCTTCGAGCCCGCCCTCCACCGCTCCGCGTTGCTTGCCT
# CCCGTCTGCCTTCAAGCCTGCCTTCGAGACCGCTGCCCTTCAACGACTCCCTCCGCCTTC
# GAGACCGCTGCCTTCGACGACTCCCTCCTCCTTCGAGCCTGCCTGGCCAGAGTTTTCAGC
# TGGAAAGTGACAATGATGTTGTGAAGATGTTTATTGTACTTCATGACAAGAAATTGATTA
# AATTGGCCATTAGATCTGCTGAGTCACCATCAAAAAAATTGGAGTTGGCACGAGGTTTGA
# TTGAATCAAAGAGAAATGCAGAGAACTCATTACTTCCAAAGTCAAGTCATGGAAGTCAAC
# CCACTTCTAGTAGGAGAGTTCCTGTCACATCTGCAAGAAGAAGAGGGATACATGCTCCAA
# AATCAGTAAGTACCTCTCATTCGATAGGACTTAGACGAACAATTATGGAGAGACATTGTA
# TGAGACAAGAAGTTGATGCGAACATGCCTTTGGATGGCATTACTGAACATGCAAAACAAA
# AACTAGCAAAGGCAGAGATGGAATCACGTTCTTTACAGAATTCTGTTGCAAATGGGAAGA
# ATCCATATAGTAAGAACAATGATTCTTCTGGAGTTGTTACACCTCATGAAACAGTTTCAC
# TAGGAAAGGAATTGATAGGACAAAAAGTTTATACCAAATGGCTTGATGACAACAAATTTT
# ATGAGGCTGTCATCACCCAGTTTAATCCGCAGATCATTATGTTCATTGGTTGA
# /usr/bin/which: no bl2seq in (/cvmfs/soft.computecanada.ca/easybuild/software/20
# 20/avx2/Core/mafft/7.471/bin:/cvmfs/soft.computecanada.ca/easybuild/software/20$
# 0/avx2/Compiler/gcc9/paml/4.9j/bin:/cvmfs/soft.computecanada.ca/easybuild/softwa
# re/2020/avx2/Core/diamond/2.0.4/bin:/cvmfs/soft.computecanada.ca/easybuild/softw
# are/2020/avx2/Core/cdbfasta/0.99/bin:/cvmfs/soft.computecanada.ca/easybuild/soft
# Run bl2seq (-p tblastn) or GeneWise to see the inconsistency.


#------------------------
# narval2
tmux new-session -s total
tmux attach-session -t total

cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees

ls | head -n 50

OG0000000_tree.txt
OG0000000_tree.txt_pal2nal
OG0000000_tree.txt_Proteins_alignment.fasta
OG0000000_tree.txt_Proteins.fasta
OG0000000_tree.txt_Transcripts.fasta

# How to match genes: https://academic.oup.com/plphys/article/180/1/404/6117721?login=false











