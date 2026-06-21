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

#----------------
# First need to download transcripts for all species
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
#35790
grep "maker-lg" *_Proteins.fasta | wc -l
#23626

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

###############################
# Guidance and PRANK
# PRANK more phylogenetically correct but more variable too 
# - needs bootstrapping in guidance

# https://gensoft.pasteur.fr/docs/prank/170427/
# https://ariloytynoja.github.io/prank-msa/

cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees/Arctic_trees

ls *_tree.txt_unique.nxh.ABSREL.json | \
sed 's/_tree.txt_unique.nxh.ABSREL.json/_tree.txt_Proteins.fasta/g' > prank_list.txt

head prank_list.txt

wc -l prank_list.txt
# 16328 

#----------------
# Use GUIDANCE2 to filter the prank alignments and 
# only keep those with scoresa above 0.93
# https://github.com/anzaika/guidance

git clone https://github.com/anzaika/guidance.git

cd guidance
make
# Running `make' takes a while

module load StdEnv/2020 gcc/9.3.0 seqkit/2.3.1
module load prank/170427 perl/5.30.2 ruby/2.7.1 bioperl/1.7.8

cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees/Arctic_trees

seqkit seq -m 30 OG0022767_tree.txt_Proteins.fasta -o OG0022767_tree.txt_Proteins.filtered.fasta

# test for one
perl /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees/Arctic_trees/guidance/www/Guidance/guidance.pl \
--seqFile OG0022767_tree.txt_Proteins.filtered.fasta \
--msaProgram PRANK \
--seqType aa \
--userTree OG0022767_tree.txt \
--outDir /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees/Arctic_trees/OG0022767_guidance \
--bootstraps 20 \
--seqCutoff 0.90 \
--colCutoff 0.90

cp ./OG0022767_guidance/Seqs.Orig.fas.FIXED.Without_low_SP_Seq.With_Names OG0022767.MSA.guidanceFiltered.fasta

################
# Run script below until disk quota error

cat << 'EOF' > prank_guidance_array.sh
#!/bin/bash
#SBATCH --account=def-henryg
#SBATCH --time=0-7:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=200G
#SBATCH --output=logs/prank_guidance_%A_%a.out
#SBATCH --error=logs/prank_guidance_%A_%a.err

# Load modules
module load StdEnv/2020 gcc/9.3.0 seqkit/2.3.1
module load prank/170427 perl/5.30.2 ruby/2.7.1 bioperl/1.7.8

# Working directory
cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees/Arctic_trees

# Number of files per SLURM array task
CHUNK=20
NCPU=$SLURM_CPUS_PER_TASK  # 20 CPUs

# Determine which lines this array task should process
START=$(( (SLURM_ARRAY_TASK_ID - 1) * CHUNK + 1 ))
END=$(( START + CHUNK - 1 ))

TOTAL=$(wc -l < prank_remaining.txt)
if [ $START -gt $TOTAL ]; then
    echo "No files to process for this task."
    exit 0
fi
if [ $END -gt $TOTAL ]; then
    END=$TOTAL
fi

echo "Processing lines $START to $END"

# Export needed vars for parallel jobs
export NCPU

# Run up to NCPU GUIDANCE jobs in parallel
sed -n "${START},${END}p" prank_remaining.txt | \
xargs -n 1 -P $NCPU -I {} bash -c '
FILE="{}"
PREFIX=$(basename "$FILE" _tree.txt_Proteins.fasta)
OUT="${PREFIX}_guidance"
FINAL="${PREFIX}_MSA.guidanceFiltered.fasta"

# Skip if output already exists
if [[ -s "$FINAL" ]]; then
    echo "Final alignment exists for $PREFIX — skipping"
    exit 0
fi

# Pre-filter very short proteins
FILTERED_FILE="${PREFIX}_tree.txt_Proteins.filtered.fasta"
seqkit seq -m 30 "$FILE" -o "$FILTERED_FILE"

# Check that filtered file exists
if [[ ! -s "$FILTERED_FILE" ]]; then
    echo "No sequences remaining after filtering for $PREFIX — skipping"
    exit 0
fi

echo "Running GUIDANCE for $PREFIX"

perl /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees/Arctic_trees/guidance/www/Guidance/guidance.pl \
    --seqFile "$FILTERED_FILE" \
    --msaProgram PRANK \
    --seqType aa \
    --outDir "$OUT" \
    --bootstraps 20 \
    --seqCutoff 0.9 \
    --colCutoff 0.9

echo "Finished $PREFIX"

# Rename filtered MSA to include OG name
if [[ -f "${OUT}/Seqs.Orig.fas.FIXED.Without_low_SP_Seq.With_Names" ]]; then
    cp "${OUT}/Seqs.Orig.fas.FIXED.Without_low_SP_Seq.With_Names" "${PREFIX}_MSA.guidanceFiltered.fasta"
    cp "${OUT}/MSA.PRANK.Without_low_SP_Col.With_Names" "${PREFIX}_MSA.PRANK.Without_low_SP_Col.With_Names"
    cp "${OUT}/MSA.PRANK.Guidance2_res_pair_seq.scr_with_Names" "${PREFIX}_MSA.PRANK.Guidance2_res_pair_seq.scr_with_Names"
    cp "${OUT}/MSA.PRANK.aln.With_Names" "${PREFIX}_MSA.PRANK.aln.With_Names"
    rm -rf "$OUT"
fi

'

echo "All files in this array task finished."

EOF


chmod +x prank_guidance_array.sh
dos2unix prank_guidance_array.sh

CHUNK=20
TOTAL=$(wc -l < prank_remaining.txt)
JOBS=$(( (TOTAL + CHUNK - 1) / CHUNK ))
echo "Submitting $JOBS array jobs"
# Submitting 743 array jobs

sbatch --array=1-${JOBS}%100 prank_guidance_array.sh
# Submitted batch job 57729659

ls *_MSA.guidanceFiltered.fasta | wc -l 
# 1476

##########################
# Check how many finished

cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees/Arctic_trees

# Some failed for time limit
grep "TIME LIMIT" logs/prank_guidance_*.err | wc -l
# 99

# Some failed for Disk quota but still say finished
grep "Disk quota exceeded" logs/prank_guidance_*.out | wc -l
# 261

# Say finished if disk quota error too
grep "Finished" logs/prank_guidance_*.out | wc -l
# 3440

# 1620 have files
# Maybe 1537 before copying 

# How many are small or incomplete?
find . -name "*_MSA.guidanceFiltered.fasta" -size -2k |wc -l
# 144 are small
find . -name "*_MSA.guidanceFiltered.fasta" -size -2k > failed_small.txt

# less than 2 sequences
for f in *_MSA.guidanceFiltered.fasta; do
    if [[ $(grep -c "^>" "$f") -lt 2 ]]; then
        echo "$f"
    fi
done > failed_fewseq.txt

wc -l failed_fewseq.txt

#--------------
# Combine filed lists

cat failed_* | sort -u > all_failed.txt
wc -l all_failed.txt
#201

sed 's|^\./||' all_failed.txt | sort -u > all_failed_clean.txt
wc -l all_failed_clean.txt
#144

# Delete failed files
xargs rm < all_failed_clean.txt

# Verify
ls *_MSA.guidanceFiltered.fasta | wc -l
# 1476


##########################
# tar files to make more room

tmux new-session -s total1
tmux attach-session -t total1

cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees

tar -cf prank_no_guidance.tar prank_no_guidance

rm -r prank_no_guidance/

#---------------------
# Clean up guidance directories so far

tmux new-session -s total
tmux attach-session -t total

cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees/Arctic_trees

mkdir guidance_directories
mv OG*_guidance guidance_directories

tar -cf guidance_directories.tar guidance_directories

rm -r guidance_directories

###############################
# Run for all files still missing 

while read f; do
    PREFIX=$(basename "$f" _tree.txt_Proteins.fasta)
    if [[ ! -s ${PREFIX}_MSA.guidanceFiltered.fasta ]]; then
        echo "$f"
    fi
done < prank_list.txt > prank_remaining.txt

wc -l prank_list.txt
wc -l prank_remaining.txt

# 16328 prank_list.txt
# 14852 prank_remaining.txt

#---------------------
# Pre-filter

module load StdEnv/2023 seqkit/2.5.1

INPUT_LIST=prank_remaining.txt
CLEAN_LIST=prank_remaining_clean.txt

> too_few_sequences.txt
> too_short_sequences.txt
> extreme_length_variation.txt

while read FILE; do

    # Count sequences
    SEQCOUNT=$(grep -c "^>" "$FILE")

    # Skip if too few sequences
    if [ "$SEQCOUNT" -lt 3 ]; then
        echo "$FILE" >> too_few_sequences.txt
        continue
    fi

    # Get min and max lengths
    read MIN MAX <<< $(seqkit stats -T "$FILE" | awk 'NR==2 {print $6, $8}')

    # Skip if too short
    if [ "$MIN" -lt 30 ]; then
        echo "$FILE" >> too_short_sequences.txt
        continue
    fi

    # Skip if extreme variation
    if [ "$MIN" -gt 0 ]; then
        RATIO=$(awk "BEGIN {print $MAX/$MIN}")
        if (( $(awk "BEGIN {print ($RATIO > 10)}") )); then
            echo "$FILE" >> extreme_length_variation.txt
            continue
        fi
    fi

    # If passed all filters → keep it
    echo "$FILE" >> "$CLEAN_LIST"

done < "$INPUT_LIST"

wc -l prank_remaining_clean.txt
# 14471 prank_remaining_clean.txt
# 16328 prank_list.txt
# 14852 prank_remaining.txt

# Removed 381

###############################

cat << 'EOF' > prank_guidance_array.sh
#!/bin/bash
#SBATCH --account=def-henryg
#SBATCH --time=0-2:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=logs/guidance_%A_%a.out
#SBATCH --error=logs/guidance_%A_%a.err

set -euo pipefail

# -----------------------------
# Modules
# -----------------------------
# Load modules
module load StdEnv/2020 gcc/9.3.0 seqkit/2.3.1
module load prank/170427 perl/5.30.2 ruby/2.7.1 bioperl/1.7.8

# -----------------------------
# Work directory (Lustre)
# -----------------------------
cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees/Arctic_trees

INPUT_LIST=prank_remaining_clean.txt

FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$INPUT_LIST")

if [[ -z "$FILE" ]]; then
    exit 0
fi

PREFIX=$(basename "$FILE" _tree.txt_Proteins.fasta)
FINAL="${PREFIX}_MSA.guidanceFiltered.fasta"

# Skip if already done (restart safe)
if [[ -s "$FINAL" ]]; then
    echo "Skipping $PREFIX (already exists)"
    exit 0
fi

# -----------------------------
# Use NODE-LOCAL TEMP SPACE
# -----------------------------
WORKDIR="${SLURM_TMPDIR}/${PREFIX}"
mkdir -p "$WORKDIR"
export TMPDIR="$WORKDIR"

# -----------------------------
# Light preprocessing in temp
# -----------------------------
cp "$FILE" "$WORKDIR/input.fasta"

# Remove very short sequences
seqkit seq -m 30 "$WORKDIR/input.fasta" -o "$WORKDIR/filtered.fasta"

SEQCOUNT=$(grep -c "^>" "$WORKDIR/filtered.fasta" || true)

if [[ "$SEQCOUNT" -lt 3 ]]; then
    echo "Skipping $PREFIX (too few sequences)"
    exit 0
fi

# -----------------------------
# Run GUIDANCE (reduced output load)
# -----------------------------
perl /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees/Arctic_trees/guidance/www/Guidance/guidance.pl \
    --seqFile "$WORKDIR/filtered.fasta" \
    --msaProgram PRANK \
    --seqType aa \
    --outDir "$WORKDIR/guidance_out" \
    --bootstraps 10 \
    --seqCutoff 0.9 \
    --colCutoff 0.9

# -----------------------------
# Copy only FINAL alignment back
# -----------------------------
if [[ -f "$WORKDIR/guidance_out/Seqs.Orig.fas.FIXED.Without_low_SP_Seq.With_Names" ]]; then
    cp "$WORKDIR/guidance_out/Seqs.Orig.fas.FIXED.Without_low_SP_Seq.With_Names" \
       /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees/Arctic_trees/"$FINAL"
    echo "Finished $PREFIX"
else
    echo "$PREFIX failed" >> guidance_failures.txt
fi

# Cleanup node temp (automatic anyway)
rm -rf "$WORKDIR"

EOF

chmod +x prank_guidance_array.sh
dos2unix prank_guidance_array.sh

#------------------
# Submit in batches

TOTAL=$(wc -l < prank_remaining_clean.txt)
BATCH=900
CONCURRENCY=50

echo "Total genes: $TOTAL"
echo "Batch size: $BATCH"


for ((START=1; START<=TOTAL; START+=BATCH)); do

    END=$((START + BATCH - 1))
    if [ $END -gt $TOTAL ]; then
        END=$TOTAL
    fi

    echo sbatch --array=${START}-${END}%${CONCURRENCY} prank_guidance_array.sh

done

sbatch --array=1-900%25 prank_guidance_array.sh
# Submitted batch job 57751254 - start 4pm

# update to 100 allowed
scontrol update JobId=57751254 ArrayTaskThrottle=100

# Cancel PD jobs and update script to 4 hours of time
scancel -t PD 57751254

#############################
# For next batches

# Narval1
tmux new-session -s total
tmux attach-session -t total

INPUT_LIST="prank_remaining_clean.txt"
SCRIPT="prank_guidance_array.sh"

TOTAL=$(wc -l < $INPUT_LIST)

BATCH=900
CONCURRENCY=100
START=103

while [ $START -le $TOTAL ]; do

    CURRENT=$(squeue -h -u $USER -r | wc -l)

    if [ $CURRENT -lt 50 ]; then

        END=$((START + BATCH - 1))
        if [ $END -gt $TOTAL ]; then
            END=$TOTAL
        fi

        echo "Submitting $START-$END"

        sbatch --array=${START}-${END}%${CONCURRENCY} $SCRIPT

        START=$((END + 1))

    else

        echo "Queue full ($CURRENT jobs). Waiting..."
        sleep 1000

    fi

done

# Queue full (896 jobs). Waiting...
# Queue full (70 jobs). Waiting...

# Working but slow for 15k files

#########################################
# Restart it to run faster in chunks

# Cancel PD jobs and update script to 4 hours of time
scancel -t PD 57760128

cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees/Arctic_trees

# Get completed OG IDs - files larger than 1k
find . -name "*_MSA.guidanceFiltered.fasta" -size +1k \
| sed 's|./||' \
| sed 's/_MSA.guidanceFiltered.fasta//' \
> completed_OGs.txt

# 1550 completed_OGs.txt without file size
# 1526 completed_OGs.txt with files size filter

# Remove them from remaining list
grep -v -F -f completed_OGs.txt prank_remaining_clean.txt > prank_remaining_restart.txt

# Replace list
mv prank_remaining_restart.txt prank_remaining_clean.txt

# Check
wc -l prank_remaining_clean.txt
# 14397 prank_remaining_clean.txt


######################
# Rerun with SLURM_TMP dir for output and prank tmp dir

cat << 'EOF' > prank_guidance_array.sh
#!/bin/bash
#SBATCH --account=def-henryg
#SBATCH --time=0-7:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=80G
#SBATCH --output=logs/guidance_%A_%a.out
#SBATCH --error=logs/guidance_%A_%a.err

set -euo pipefail

# ----------------------------
# Load modules
# ----------------------------
module load StdEnv/2020
module load gcc/9.3.0
module load seqkit/2.3.1
module load prank/170427
module load perl/5.30.2
module load bioperl/1.7.8

# ----------------------------
# Working directory
# ----------------------------
WORKDIR=/home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees/Arctic_trees
cd $WORKDIR

INPUT_LIST="prank_remaining_clean.txt"

# ----------------------------
# Parallel settings
# ----------------------------
CHUNK=20
NCPU=$SLURM_CPUS_PER_TASK

TOTAL=$(wc -l < $INPUT_LIST)

START=$(( (SLURM_ARRAY_TASK_ID - 1) * CHUNK + 1 ))
END=$(( START + CHUNK - 1 ))

if [ $START -gt $TOTAL ]; then
    exit 0
fi

if [ $END -gt $TOTAL ]; then
    END=$TOTAL
fi

# Ensure tmp is set
export TMPDIR=$SLURM_TMPDIR

# ----------------------------
# GUIDANCE function
# ----------------------------

run_guidance() {

    FILE="$1"

    PREFIX=$(basename "$FILE" _tree.txt_Proteins.fasta)
    FINAL="${PREFIX}_MSA.guidanceFiltered.fasta"

    if [[ -s "$WORKDIR/$FINAL" ]]; then
        echo "Skipping $PREFIX"
        return
    fi

    # unique working directory
    JOBTMP=$(mktemp -d "$SLURM_TMPDIR/guidance_${PREFIX}_XXXX")

    cp "$FILE" "$JOBTMP/input.fasta"

    # PRANK sometimes writes tmp files in cwd
    export PRANK_TMPDIR="$JOBTMP"
    export TMPDIR="$JOBTMP"

    perl "$WORKDIR/guidance/www/Guidance/guidance.pl" \
        --seqFile "$JOBTMP/input.fasta" \
        --msaProgram PRANK \
        --seqType aa \
        --outDir "$JOBTMP/guidance_out" \
        --bootstraps 10 \
        --seqCutoff 0.9 \
        --colCutoff 0.9

    if [[ -f "$JOBTMP/guidance_out/Seqs.Orig.fas.FIXED.Without_low_SP_Seq.With_Names" ]]; then

        cp "$JOBTMP/guidance_out/Seqs.Orig.fas.FIXED.Without_low_SP_Seq.With_Names" \
           "$WORKDIR/$FINAL"

        echo "Finished $PREFIX"

    else

        echo "$PREFIX failed" >> "$WORKDIR/guidance_failures.txt"

    fi

    rm -rf "$JOBTMP"

}

export -f run_guidance
export WORKDIR
export SLURM_TMPDIR

# ----------------------------
# Run parallel GUIDANCE jobs
# ----------------------------

sed -n "${START},${END}p" $INPUT_LIST | \
xargs -n 1 -P $NCPU -I {} bash -c 'run_guidance "$@"' _ {}

echo "Array task $SLURM_ARRAY_TASK_ID finished."

EOF

chmod +x prank_guidance_array.sh
dos2unix prank_guidance_array.sh

sbatch --array=1-3 prank_guidance_array.sh

# Check SLURM_TMP dir
sq
srun --jobid=57781817_1 --pty bash
echo $SLURM_TMPDIR
ls $SLURM_TMPDIR

more guidance_57781817_1.out

sbatch --array=4-724%50 prank_guidance_array.sh
# 57782132

# update to 100 allowed
scontrol update JobId=57782132 ArrayTaskThrottle=100


#####################################
# Check
ls *_MSA.guidanceFiltered.fasta | wc -l 
# starting 1476 at 4:15pm March 13
# 1478 at 4:45pm 
# 1523 at 11:20pm
# 1547 at 1am 
# 1550 at 1pm March 14th
# 1745 at 3pm
# 1810 at 3:24pm
# 1902 at 4pm
# 1994 at 4:40pm
# 2762 at 6:31pm
# 4234 at 11:25pm
# 4900 at 1am March 15th
# 9547 at 12:45pm 
# 10950 at 3pm
# 10961 at 9pm - done



# 81 failed because of time
# relaunch these with faster script for 7 hours

wc -l guidance_failures.txt
# 4809 guidance_failures.txt

# OG0000054 failed
# OG0000353 failed
# OG0000186 failed
# OG0000380 failed
# OG0000425 failed
# OG0000584 failed
# OG0000604 failed
# OG0000717 failed
# OG0000907 failed
# OG0000945 failed
# OG0000997 failed
# OG0000998 failed
# OG0001037 failed
# OG0001179 failed
# OG0001229 failed
# OG0001259 failed
# OG0001187 failed
# OG0001304 failed
# OG0001060 failed
# OG0001036 failed
# OG0001164 failed
# OG0001365 failed
# OG0001182 failed
# OG0001338 failed
# OG0000604 failed
# OG0000054 failed
# OG0000353 failed
# OG0000584 failed
# OG0000380 failed
# OG0000425 failed
# OG0001304 failed
# OG0001259 failed


cd logs/
ls guidance* 

more guidance_57751254_3.out
# OG0000054
# ERROR: Seq: 'Polavi_Chr200001774' contained the character: '*' which is not a standard Amino Acid<br>

# ERROR: Seq: 'maker-lg8-snap-gene-2.214-mRNA-1' contained the character: '*' which is
 # not a standard Amino Acid<br>
# ERROR: Seq: 'snap_masked-lg6-processed-gene-116.35-mRNA-1' contained the character:
# '*' which is not a standard Amino Acid<br>
# ERROR: Seq: 'maker-lg1-snap-gene-49.68-mRNA-1' contained the character: '*' which is
 # not a standard Amino Acid<br>

grep TIME guidance_57782132_*.out

grep TIME guidance_57782132_*.err

#############################
# Get a new list remaining 

# Get completed OG IDs - files larger than 1k
find . -name "*_MSA.guidanceFiltered.fasta" -size +1k \
| sed 's|./||' \
| sed 's/_MSA.guidanceFiltered.fasta//' \
> completed_OGs.txt

cut -d' ' -f1 guidance_failures.txt \
| sort | uniq > failed_ogs.txt

sed 's#.*/##' prank_list.txt \
| sed 's/_tree.txt_Proteins.fasta//' \
| sort > input_OGs.txt

wc -l input_OGs.txt
wc -l completed_OGs.txt
wc -l guidance_failures.txt
# 16328 prank_list.txt
# 9846 completed_OGs.txt
# 4809 guidance_failures.txt

ls *_MSA.guidanceFiltered.fasta | wc -l 
# 10961

#-------------------
# figure out missing

# Combine completed + failed
cat completed_OGs.txt failed_ogs.txt | sort | uniq > done_ogs.txt

wc -l done_ogs.txt
# 14633 done_ogs.txt

# Get missing
comm -23 input_OGs.txt done_ogs.txt > missing_ogs.txt
wc -l missing_ogs.txt
# 1695 missing_ogs.txt


while read OG; do
    echo "${OG}_tree.txt_Proteins.fasta"
done < missing_ogs.txt > missing_protein_files.txt

#------------------------
# rerun with more time

cat << 'EOF' > prank_guidance_missingOG.sh
#!/bin/bash
#SBATCH --account=def-henryg
#SBATCH --time=0-12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=80G
#SBATCH --output=logs/guidance_%A_%a.out
#SBATCH --error=logs/guidance_%A_%a.err

set -euo pipefail

# ----------------------------
# Load modules
# ----------------------------
module load StdEnv/2020
module load gcc/9.3.0
module load seqkit/2.3.1
module load prank/170427
module load perl/5.30.2
module load bioperl/1.7.8

# ----------------------------
# Working directory
# ----------------------------
WORKDIR=/home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees/Arctic_trees
cd $WORKDIR

INPUT_LIST="missing_protein_files.txt"

# ----------------------------
# Parallel settings
# ----------------------------
CHUNK=20
NCPU=$SLURM_CPUS_PER_TASK

TOTAL=$(wc -l < $INPUT_LIST)

START=$(( (SLURM_ARRAY_TASK_ID - 1) * CHUNK + 1 ))
END=$(( START + CHUNK - 1 ))

if [ $START -gt $TOTAL ]; then
    exit 0
fi

if [ $END -gt $TOTAL ]; then
    END=$TOTAL
fi

# Ensure tmp is set
export TMPDIR=$SLURM_TMPDIR

# ----------------------------
# GUIDANCE function
# ----------------------------

run_guidance() {

    FILE="$1"

    PREFIX=$(basename "$FILE" _tree.txt_Proteins.fasta)
    FINAL="${PREFIX}_MSA.guidanceFiltered.fasta"

    # unique working directory
    JOBTMP=$(mktemp -d "$SLURM_TMPDIR/guidance_${PREFIX}_XXXX")

    cp "$FILE" "$JOBTMP/input.fasta"

    # PRANK sometimes writes tmp files in cwd
    export PRANK_TMPDIR="$JOBTMP"
    export TMPDIR="$JOBTMP"

    perl "$WORKDIR/guidance/www/Guidance/guidance.pl" \
        --seqFile "$JOBTMP/input.fasta" \
        --msaProgram PRANK \
        --seqType aa \
        --outDir "$JOBTMP/guidance_out" \
        --bootstraps 10 \
        --seqCutoff 0.9 \
        --colCutoff 0.9

    if [[ -f "$JOBTMP/guidance_out/Seqs.Orig.fas.FIXED.Without_low_SP_Seq.With_Names" ]]; then

        cp "$JOBTMP/guidance_out/Seqs.Orig.fas.FIXED.Without_low_SP_Seq.With_Names" \
           "$WORKDIR/$FINAL"

        echo "Finished $PREFIX"

    else

        echo "$PREFIX failed" >> "$WORKDIR/guidance_failures.txt"

    fi

    rm -rf "$JOBTMP"

}

export -f run_guidance
export WORKDIR
export SLURM_TMPDIR

# ----------------------------
# Run parallel GUIDANCE jobs
# ----------------------------

sed -n "${START},${END}p" $INPUT_LIST | \
xargs -n 1 -P $NCPU -I {} bash -c 'run_guidance "$@"' _ {}

echo "Array task $SLURM_ARRAY_TASK_ID finished."

EOF

chmod +x prank_guidance_missingOG.sh
dos2unix prank_guidance_missingOG.sh

sbatch --array=1-85 prank_guidance_missingOG.sh
# Submitted batch job 57824180

ls *_MSA.guidanceFiltered.fasta | wc -l
# 10961 at 11:15pm March 15th
# 11023 at 1:11am March 16th
# 11244 8:20am
# 11291 at 10:20am 
# 11314 - done/ran out of time

#353 more ran

more logs/guidance_57824180_3.out

grep TIME logs/guidance_57824180_*.err | wc -l 
# 44

#------------------------
# Check counts

# Starting
# 1695 missing_ogs.txt
# Original 4809 guidance_failures.txt

ls *_MSA.guidanceFiltered.fasta | wc -l
# 10961 at 11:15pm March 15th
# 11314 - done/ran out of time

# 353 more ran

wc -l guidance_failures.txt
# 5276 guidance_failures.txt
# 467 more failed

# Not enough time:
#1695-353-467 = 875

# Skip these for now

########################################
# Guidance palnal test

cp ./OG0022767_guidance/Seqs.Orig.fas.FIXED.Without_low_SP_Seq.With_Names OG0022767.MSA.guidanceFiltered.fasta

nano pal2nal_test.sh
#!/bin/bash

# Set the file
PROT="OG0022767.MSA.guidanceFiltered.fasta"
PREFIX="${PROT%.MSA.guidanceFiltered.fasta}"
CDS="${PREFIX}_tree.txt_Transcripts.fasta"
PAL2NAL_OUT="${PREFIX}_guidance_pal2nal.fasta"

# Check files
if [[ ! -s "$PROT" || ! -s "$CDS" ]]; then
    echo "Missing protein or CDS — abort"
    exit 1
fi

# Filter CDS to match protein sequences
CDS_FILTERED="${CDS%.fasta}.filtered.fasta"
seqkit grep -f <(seqkit seq -n "$PROT") "$CDS" > "$CDS_FILTERED"

# Run PAL2NAL
/home/celphin/scratch/Oxyria_Positive_Selection_Test/pal2nal.v14/pal2nal.pl \
    "$PROT" "$CDS_FILTERED" -output fasta -nogap > "$PAL2NAL_OUT"

# Check
if [[ -s "$PAL2NAL_OUT" ]]; then
    echo "PAL2NAL completed: $PAL2NAL_OUT"
else
    echo "PAL2NAL failed"
fi

chmod +x pal2nal_test.sh
./pal2nal_test.sh

# [INFO] 4 patterns loaded from file
# PAL2NAL completed: OG0022767_guidance_pal2nal.fasta


#------------------
# Run pal2nal on all guidance output
# Very fast - maybe don't need slurm submission

# Go to the working directory 
cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees/Arctic_trees

# Check how many there are
ls *_MSA.guidanceFiltered.fasta > guidance_list.txt

wc -l guidance_list.txt

cat << 'EOF' > guidance_pal2nal.sh
#!/bin/bash
#SBATCH --account=def-henryg
#SBATCH --time=0-01:00:00      # 1 hour per task
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1      # PAL2NAL is single-threaded but 20 chunks
#SBATCH --mem=4G                 # enough for most alignments
#SBATCH --output=logs/pal2nal_guidance_%A_%a.out
#SBATCH --error=logs/pal2nal_guidance_%A_%a.err

set -euo pipefail

module load StdEnv/2020 seqkit/2.3.1 perl/5.30.2

cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees/Arctic_trees

# List all filtered protein MSAs
TOTAL=$(wc -l < guidance_list.txt)

# Number of OGs per SLURM array task
CHUNK=20

START=$(( (SLURM_ARRAY_TASK_ID - 1) * CHUNK + 1 ))
END=$(( START + CHUNK - 1 ))

if [ $START -gt $TOTAL ]; then
    echo "No files to process for this task."
    exit 0
fi

if [ $END -gt $TOTAL ]; then
    END=$TOTAL
fi

echo "Processing lines $START to $END"

# Process assigned files
sed -n "${START},${END}p" guidance_list.txt | while read PROT; do
    PREFIX="${PROT%_MSA.guidanceFiltered.fasta}"
    CDS="${PREFIX}_tree.txt_Transcripts.fasta"
    PAL2NAL_OUT="${PREFIX}_guidance_pal2nal.fasta"

    # Skip if output exists
    if [[ -f "$PAL2NAL_OUT" ]]; then
        echo "PAL2NAL output already exists for $PREFIX — skipping"
        continue
    fi

    # Check files
    if [[ ! -s "$PROT" || ! -s "$CDS" ]]; then
        echo "Missing protein or CDS for $PREFIX — skipping"
        continue
    fi

    # Filter CDS to match filtered protein sequences
    CDS_FILTERED="${CDS%.fasta}.filtered.fasta"
    seqkit grep -f <(seqkit seq -n "$PROT") "$CDS" > "$CDS_FILTERED"
    if [[ ! -s "$CDS_FILTERED" ]]; then
        echo "CDS filtering failed for $PREFIX"
        continue
    fi

    # Run PAL2NAL
    /home/celphin/scratch/Oxyria_Positive_Selection_Test/pal2nal.v14/pal2nal.pl \
        "$PROT" "$CDS_FILTERED" -output fasta -nogap > "$PAL2NAL_OUT"

    if [[ -s "$PAL2NAL_OUT" ]]; then
        echo "PAL2NAL completed for $PREFIX"
    else
        echo "PAL2NAL failed for $PREFIX"
    fi

done

echo "All files in this array task finished."

EOF

chmod +x guidance_pal2nal.sh
dos2unix guidance_pal2nal.sh


# Calculate number of jobs
TOTAL=$(wc -l < guidance_list.txt)
JOBS=$(( (TOTAL + 19) / 20 ))  # 20 OGs per job
# JOBS=566

sbatch --array=1-$JOBS%100 guidance_pal2nal.sh
# Submitted batch job 57851846

more logs/pal2nal_guidance_57851846_1.out

# PAL2NAL completed for OG0000203
# PAL2NAL completed for OG0000204
# PAL2NAL completed for OG0000206
# All files in this array task finished.

ls *_guidance_pal2nal.fasta | wc -l 
# 11314

# 11314 guidance_list.txt

###################################

# Try getting IQTREE2 for better gene trees?


iqtree2 -s aligned_cds.fasta \
  -spp codon \
  -p split \
  -m GTR+G \
  -B 1000 \
  -alrt 1000 \
  --prefix cds_codon_tree