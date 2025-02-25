##########################
# Running GeneSpace
# Feb-Aug 2024
############################

# Paper: https://elifesciences.org/articles/78526 
# https://github.com/jtlovell/GENESPACE

################################
# Installation already done

################################
# Collecting and renameing annotation files
# Need both the fasta files and the gff3 files
# Consistient naming 

##############################
# Polygonaceae genome downloads/copies

cd /home/celphin/scratch/Oxyria/GeneSpace/OrthoFinder_source/Protein_files
cp /home/celphin/scratch/Oxyria/Polygonaceae_Genomes_Annotations/MS_CE_Fagopyrum_tataricum/Fagopyrum_proteins_AED0.6.fasta .
cp /home/celphin/scratch/Oxyria/Polygonaceae_Genomes_Annotations/MS_CE_Oxyria_digyna/Oxyria_proteins_AED0.6.fasta .
cp /home/celphin/scratch/Oxyria/Polygonaceae_Genomes_Annotations/Fagopyrum_escelentum/F_escelentum_H1.proteins.fasta .
cp /home/celphin/scratch/Oxyria/Polygonaceae_Genomes_Annotations/Fagopyrum_escelentum/F_escelentum_H2.proteins.fasta .
cp /home/celphin/scratch/Oxyria/Polygonaceae_Genomes_Annotations/NCBI_Oxyria_digyna/Oxyria_NCBI_proteins_AED0.6.fasta .
cp /home/celphin/scratch/Oxyria/Polygonaceae_Genomes_Annotations/OTHER_Fagopyrum_tataricum/F_tataricum_H1.proteins.fasta .
cp /home/celphin/scratch/Oxyria/Polygonaceae_Genomes_Annotations/OTHER_Fagopyrum_tataricum/F_tataricum_H2.proteins.fasta .
cp /home/celphin/scratch/Oxyria/Polygonaceae_Genomes_Annotations/Polygonum_aviculare/Polavi_proteins_AED0.6.fasta .
cp /home/celphin/scratch/Oxyria/Polygonaceae_Genomes_Annotations/Rheum_nobile/R_nobile_proteins.fasta .
cp /home/celphin/scratch/Oxyria/Polygonaceae_Genomes_Annotations/Rheum_tangaticum/R_tangaticum.proteins.fasta .
#cp /home/celphin/scratch/Oxyria/Polygonaceae_Genomes_Annotations/Rumex_hastatulus/

#---------
# rename all
# proteins
mv Fagopyrum_proteins_AED0.6.fasta Fagopyrum_tataricum_H0.fa
mv Oxyria_proteins_AED0.6.fasta Oxyria_digyna_H0.fa
mv F_escelentum_H1.proteins.fasta Fagopyrum_escelentum_H1.fa
mv F_escelentum_H2.proteins.fasta Fagopyrum_escelentum_H2.fa
mv Oxyria_NCBI_proteins_AED0.6.fasta Oxyria_digyna_H1.fa
mv F_tataricum_H1.proteins.fasta Fagopyrum_tataricum_H1.fa
mv F_tataricum_H2.proteins.fasta Fagopyrum_tataricum_H2.fa
mv Polavi_proteins_AED0.6.fasta Polygunum_aviculare_H0.fa
mv R_nobile_proteins.fasta Rheum_nobile_H0.fa
mv R_tangaticum.proteins.fasta Rheum_tangaticum_H0.fa 
# mv XX.fasta Rumex_hastalus_H0.fa

#---------
# copy over gff3 files
cd /home/celphin/scratch/Oxyria/Polygonaceae_Genomes_Annotations
tree |grep .gff3

cd /home/celphin/scratch/Oxyria/GeneSpace/OrthoFinder_source/gff3_files

cd /home/celphin/scratch/Oxyria/GeneSpace/OrthoFinder_source/Protein_files
cp /home/celphin/scratch/Oxyria/Polygonaceae_Genomes_Annotations/MS_CE_Fagopyrum_tataricum/FINAL_Fagopyrum.AED_0.6.sorted.gff3 .
cp /home/celphin/scratch/Oxyria/Polygonaceae_Genomes_Annotations/MS_CE_Oxyria_digyna/FINAL_Oxyria.AED_0.6.sorted.gff3 .
cp /home/celphin/scratch/Oxyria/Polygonaceae_Genomes_Annotations/Fagopyrum_escelentum/F_escelentum_H1.gff3 .
cp /home/celphin/scratch/Oxyria/Polygonaceae_Genomes_Annotations/Fagopyrum_escelentum/F_escelentum_H2.gff3 .
cp /home/celphin/scratch/Oxyria/Polygonaceae_Genomes_Annotations/NCBI_Oxyria_digyna/FINAL_Oxyria_NCBI.AED_0.6.sorted.gff3 .
cp /home/celphin/scratch/Oxyria/Polygonaceae_Genomes_Annotations/OTHER_Fagopyrum_tataricum/F_tataricum_H1.gff3 .
cp /home/celphin/scratch/Oxyria/Polygonaceae_Genomes_Annotations/OTHER_Fagopyrum_tataricum/F_tataricum_H2.gff3 .
cp /home/celphin/scratch/Oxyria/Polygonaceae_Genomes_Annotations/Polygonum_aviculare/FINAL_Polavi.AED_0.6.sorted.gff3 .
cp /home/celphin/scratch/Oxyria/Polygonaceae_Genomes_Annotations/Rheum_nobile/R_nobile.gff3 .
cp /home/celphin/scratch/Oxyria/Polygonaceae_Genomes_Annotations/Rheum_tangaticum/R_tangaticum.gff3 .
#cp /home/celphin/scratch/Oxyria/Polygonaceae_Genomes_Annotations/Rumex_hastatulus/


#-----
# gff3

mv FINAL_Fagopyrum.AED_0.6.sorted.gff3 Fagopyrum_tataricum_H0.gff3
mv FINAL_Oxyria.AED_0.6.sorted.gff3 Oxyria_digyna_H0.gff3
mv F_escelentum_H1.gff3 Fagopyrum_escelentum_H1.gff3
mv F_escelentum_H2.gff3 Fagopyrum_escelentum_H2.gff3
mv FINAL_Oxyria_NCBI.AED_0.6.sorted.gff3 Oxyria_digyna_H1.gff3
mv F_tataricum_H1.gff3 Fagopyrum_tataricum_H1.gff3
mv F_tataricum_H2.gff3 Fagopyrum_tataricum_H2.gff3
mv FINAL_Polavi.AED_0.6.sorted.gff3 Polygunum_aviculare_H0.gff3
mv R_nobile.gff3 Rheum_nobile_H0.gff3
mv R_tangaticum.gff3 Rheum_tangaticum_H0.gff3

# mv XX.fasta Rumex_hastalus_H0.fa

#-----------------------------------
# make and move to unique directories
mkdir Polygonaceae_genomes; cd Polygonaceae_genomes

mv Protein_files/ ../Polygonaceae_genomes
mv gff3_files/ ../Polygonaceae_genomes

mkdir Fagopyrum_tataricum_H0
mkdir Oxyria_digyna_H0
mkdir Fagopyrum_escelentum_H1
mkdir Fagopyrum_escelentum_H2
mkdir Oxyria_digyna_H1
mkdir Fagopyrum_tataricum_H1
mkdir Fagopyrum_tataricum_H2
mkdir  Polygunum_aviculare_H0
mkdir  Rheum_nobile_H0
mkdir  Rheum_tangaticum_H0

SPP_Hap=Fagopyrum_tataricum_H0
mv ./Protein_files/${SPP_Hap}.fa ${SPP_Hap}/${SPP_Hap}.fa
mv ./gff3_files/${SPP_Hap}.gff3 ${SPP_Hap}/${SPP_Hap}.gff3

SPP_Hap=Oxyria_digyna_H0
mv ./Protein_files/${SPP_Hap}.fa ${SPP_Hap}/${SPP_Hap}.fa
mv ./gff3_files/${SPP_Hap}.gff3 ${SPP_Hap}/${SPP_Hap}.gff3

SPP_Hap=Fagopyrum_escelentum_H1
mv ./Protein_files/${SPP_Hap}.fa ${SPP_Hap}/${SPP_Hap}.fa
mv ./gff3_files/${SPP_Hap}.gff3 ${SPP_Hap}/${SPP_Hap}.gff3

SPP_Hap=Fagopyrum_escelentum_H2
mv ./Protein_files/${SPP_Hap}.fa ${SPP_Hap}/${SPP_Hap}.fa
mv ./gff3_files/${SPP_Hap}.gff3 ${SPP_Hap}/${SPP_Hap}.gff3

SPP_Hap=Oxyria_digyna_H1
mv ./Protein_files/${SPP_Hap}.fa ${SPP_Hap}/${SPP_Hap}.fa
mv ./gff3_files/${SPP_Hap}.gff3 ${SPP_Hap}/${SPP_Hap}.gff3

SPP_Hap=Fagopyrum_tataricum_H1
mv ./Protein_files/${SPP_Hap}.fa ${SPP_Hap}/${SPP_Hap}.fa
mv ./gff3_files/${SPP_Hap}.gff3 ${SPP_Hap}/${SPP_Hap}.gff3

SPP_Hap=Fagopyrum_tataricum_H2
mv ./Protein_files/${SPP_Hap}.fa ${SPP_Hap}/${SPP_Hap}.fa
mv ./gff3_files/${SPP_Hap}.gff3 ${SPP_Hap}/${SPP_Hap}.gff3

SPP_Hap=Polygunum_aviculare_H0
mv ./Protein_files/${SPP_Hap}.fa ${SPP_Hap}/${SPP_Hap}.fa
mv ./gff3_files/${SPP_Hap}.gff3 ${SPP_Hap}/${SPP_Hap}.gff3

SPP_Hap=Rheum_nobile_H0
mv ./Protein_files/${SPP_Hap}.fa ${SPP_Hap}/${SPP_Hap}.fa
mv ./gff3_files/${SPP_Hap}.gff3 ${SPP_Hap}/${SPP_Hap}.gff3

SPP_Hap=Rheum_tangaticum_H0
mv ./Protein_files/${SPP_Hap}.fa ${SPP_Hap}/${SPP_Hap}.fa
mv ./gff3_files/${SPP_Hap}.gff3 ${SPP_Hap}/${SPP_Hap}.gff3

tree
.
├── Fagopyrum_escelentum_H1
│   ├── Fagopyrum_escelentum_H1.fa
│   └── Fagopyrum_escelentum_H1.gff3
├── Fagopyrum_escelentum_H2
│   ├── Fagopyrum_escelentum_H2.fa
│   └── Fagopyrum_escelentum_H2.gff3
├── Fagopyrum_tataricum_H0
│   ├── Fagopyrum_tataricum_H0.fa
│   └── Fagopyrum_tataricum_H0.gff3
├── Fagopyrum_tataricum_H1
│   ├── Fagopyrum_tataricum_H1.fa
│   └── Fagopyrum_tataricum_H1.gff3
├── Fagopyrum_tataricum_H2
│   ├── Fagopyrum_tataricum_H2.fa
│   └── Fagopyrum_tataricum_H2.gff3
├── Oxyria_digyna_H0
│   ├── Oxyria_digyna_H0.fa
│   └── Oxyria_digyna_H0.gff3
├── Oxyria_digyna_H1
│   ├── Oxyria_digyna_H1.fa
│   └── Oxyria_digyna_H1.gff3
├── Polygunum_aviculare_H0
│   ├── Polygunum_aviculare_H0.fa
│   └── Polygunum_aviculare_H0.gff3
├── Rheum_nobile_H0
│   ├── Rheum_nobile_H0.fa
│   └── Rheum_nobile_H0.gff3
└── Rheum_tangaticum_H0
    ├── Rheum_tangaticum_H0.fa
    └── Rheum_tangaticum_H0.gff3

#---------------------------------------------
cd /home/celphin/scratch/Oxyria/GeneSpace/Polygonaceae_genomes
# move the genomes that are duplicated haplotypes
# keep
Fagopyrum_tataricum_H1  
Fagopyrum_escelentum_H2  
Rheum_tangaticum_H0
Oxyria_digyna_H1
Rheum_nobile_H0
Polygunum_aviculare_H0

# move
mv Fagopyrum_tataricum_H0 ../other_haplotypes
mv Fagopyrum_escelentum_H1  ../other_haplotypes
mv Oxyria_digyna_H0  ../other_haplotypes
mv Fagopyrum_tataricum_H2  ../other_haplotypes

#-------------------------------------------
# make and move to unique directories
 tree
 
 
├── Fagopyrum_escelentum_H2
│   ├── Fagopyrum_escelentum_H2.fa
│   └── Fagopyrum_escelentum_H2.gff3
├── Fagopyrum_tataricum_H1
│   ├── Fagopyrum_tataricum_H1.fa
│   └── Fagopyrum_tataricum_H1.gff3
├── Oxyria_digyna_H1
│   ├── Oxyria_digyna_H1.fa
│   └── Oxyria_digyna_H1.gff3
├── Polygunum_aviculare_H0
│   ├── Polygunum_aviculare_H0.fa
│   └── Polygunum_aviculare_H0.gff3
├── Rheum_nobile_H0
│   ├── Rheum_nobile_H0.fa
│   └── Rheum_nobile_H0.gff3
├── Rheum_tangaticum_H0
│   ├── Rheum_tangaticum_H0.fa
│   └── Rheum_tangaticum_H0.gff3


#--------------------------
# check gff file formats

SPP_Hap="Fagopyrum_escelentum_H2"
head ./${SPP_Hap}/${SPP_Hap}.gff3
head ./${SPP_Hap}/${SPP_Hap}.fa


SPP_Hap="Oxyria_digyna_H1"
head ./${SPP_Hap}/${SPP_Hap}.gff3
head ./${SPP_Hap}/${SPP_Hap}.fa


SPP_Hap="Fagopyrum_tataricum_H1"
head ./${SPP_Hap}/${SPP_Hap}.gff3
head ./${SPP_Hap}/${SPP_Hap}.fa


SPP_Hap="Polygunum_aviculare_H0"
head ./${SPP_Hap}/${SPP_Hap}.gff3
head ./${SPP_Hap}/${SPP_Hap}.fa


SPP_Hap="Rheum_nobile_H0"
head ./${SPP_Hap}/${SPP_Hap}.gff3
head ./${SPP_Hap}/${SPP_Hap}.fa

SPP_Hap="Rheum_tangaticum_H0"
head ./${SPP_Hap}/${SPP_Hap}.gff3
head ./${SPP_Hap}/${SPP_Hap}.fa


#######################################
# Rosaceae genome files

cd /home/celphin/scratch/Oxyria/GeneSpace/Rosaceae_genomes

# download other genomes

mkdir Pyrus_bretschneideri; cd Pyrus_bretschneideri
# download protein files from other genomes
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_019419815.1/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/019/419/815/GCF_019419815.1_Pyrus_bretschneideri_v1/GCF_019419815.1_Pyrus_bretschneideri_v1_translated_cds.faa.gz
gunzip GCF_019419815.1_Pyrus_bretschneideri_v1_translated_cds.faa.gz
mv GCF_019419815.1_Pyrus_bretschneideri_v1_translated_cds.faa  Pyrus_bretschneideri.fa

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/019/419/815/GCF_019419815.1_Pyrus_bretschneideri_v1/GCF_019419815.1_Pyrus_bretschneideri_v1_genomic.gff.gz
gunzip GCF_019419815.1_Pyrus_bretschneideri_v1_genomic.gff.gz
mv GCF_019419815.1_Pyrus_bretschneideri_v1_genomic.gff  Pyrus_bretschneideri.gff

cd ..

mkdir Malus_sylvestris; cd Malus_sylvestris
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_916048215.2/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/916/048/215/GCF_916048215.2_drMalSylv7.2/GCF_916048215.2_drMalSylv7.2_translated_cds.faa.gz
gunzip GCF_916048215.2_drMalSylv7.2_translated_cds.faa.gz
mv   GCF_916048215.2_drMalSylv7.2_translated_cds.faa  Malus_sylvestris.fa

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/916/048/215/GCF_916048215.2_drMalSylv7.2/GCF_916048215.2_drMalSylv7.2_genomic.gff.gz
gunzip GCF_916048215.2_drMalSylv7.2_genomic.gff.gz
mv   GCF_916048215.2_drMalSylv7.2_genomic.gff  Malus_sylvestris.gff

cd ..

mkdir Prunus_persica; cd Prunus_persica
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000346465.2/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/346/465/GCF_000346465.2_Prunus_persica_NCBIv2/GCF_000346465.2_Prunus_persica_NCBIv2_translated_cds.faa.gz
gunzip GCF_000346465.2_Prunus_persica_NCBIv2_translated_cds.faa.gz
mv   GCF_000346465.2_Prunus_persica_NCBIv2_translated_cds.faa  Prunus_persica.fa

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/346/465/GCF_000346465.2_Prunus_persica_NCBIv2/GCF_000346465.2_Prunus_persica_NCBIv2_genomic.gff.gz
gunzip GCF_000346465.2_Prunus_persica_NCBIv2_genomic.gff.gz
mv   GCF_000346465.2_Prunus_persica_NCBIv2_genomic.gff  Prunus_persica.gff

cd ..

mkdir Fragaria_vesca; cd Fragaria_vesca
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000184155.1/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/184/155/GCF_000184155.1_FraVesHawaii_1.0/GCF_000184155.1_FraVesHawaii_1.0_translated_cds.faa.gz
gunzip GCF_000184155.1_FraVesHawaii_1.0_translated_cds.faa.gz
mv GCF_000184155.1_FraVesHawaii_1.0_translated_cds.faa Fragaria_vesca.fa

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/184/155/GCF_000184155.1_FraVesHawaii_1.0/GCF_000184155.1_FraVesHawaii_1.0_genomic.gff.gz
gunzip GCF_000184155.1_FraVesHawaii_1.0_genomic.gff.gz
mv  GCF_000184155.1_FraVesHawaii_1.0_genomic.gff   Fragaria_vesca.gff

cd ..

mkdir Rosa_rugosa; cd Rosa_rugosa
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_958449725.1/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/958/449/725/GCF_958449725.1_drRosRugo1.1/GCF_958449725.1_drRosRugo1.1_translated_cds.faa.gz
gunzip GCF_958449725.1_drRosRugo1.1_translated_cds.faa.gz
mv  GCF_958449725.1_drRosRugo1.1_translated_cds.faa   Rosa_rugosa.fa

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/958/449/725/GCF_958449725.1_drRosRugo1.1/GCF_958449725.1_drRosRugo1.1_genomic.gff.gz
gunzip GCF_958449725.1_drRosRugo1.1_genomic.gff.gz
mv  GCF_958449725.1_drRosRugo1.1_genomic.gff   Rosa_rugosa.gff

cd ..

mkdir Argentina_anserina; cd Argentina_anserina
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_933775445.1/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/933/775/445/GCF_933775445.1_drPotAnse1.1/GCF_933775445.1_drPotAnse1.1_translated_cds.faa.gz
gunzip GCF_933775445.1_drPotAnse1.1_translated_cds.faa.gz
mv   GCF_933775445.1_drPotAnse1.1_translated_cds.faa  Argentina_anserina.fa

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/933/775/445/GCF_933775445.1_drPotAnse1.1/GCF_933775445.1_drPotAnse1.1_genomic.gff.gz
gunzip GCF_933775445.1_drPotAnse1.1_genomic.gff.gz
mv   GCF_933775445.1_drPotAnse1.1_genomic.gff   Argentina_anserina.gff
cd ..

mkdir Rubus_argutus; cd Rubus_argutus
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_040183295.1/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/040/183/295/GCA_040183295.1_ASM4018329v1/GCA_040183295.1_ASM4018329v1_genomic.gff.gz
gunzip GCA_040183295.1_ASM4018329v1_genomic.gff.gz
mv GCA_040183295.1_ASM4018329v1_genomic.gff     Rubus_argutus.gff

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/040/183/295/GCA_040183295.1_ASM4018329v1/GCA_040183295.1_ASM4018329v1_translated_cds.faa.gz
gunzip GCA_040183295.1_ASM4018329v1_translated_cds.faa.gz
mv GCA_040183295.1_ASM4018329v1_translated_cds.faa     Rubus_argutus.fa
cd ..

#----------------------------
# copy over Dryas annotation we made
mkdir Dryas_octopetala; cd Dryas_octopetala
cp /home/celphin/scratch/Annotation/CAP_Snakemake/DoctH0/FINAL_ANNOTATION/FINAL_DoctH0.AED_0.6.sorted.gff3 \
/home/celphin/scratch/Oxyria/GeneSpace/Rosaceae_genomes/Dryas_octopetala/Dryas_octopetala.gff3
cp /home/celphin/scratch/Annotation/CAP_Snakemake/DoctH0/Summary_data/DoctH0.AED_0.6_protein.fasta \
/home/celphin/scratch/Oxyria/GeneSpace/Rosaceae_genomes/Dryas_octopetala/Dryas_octopetala.fa
cd ..

#------------------------------
# check gff and fasta file formats
SPP_Hap="Dryas_octopetala"
head -n 5 ./${SPP_Hap}/${SPP_Hap}.gff3
grep ">" ./${SPP_Hap}/${SPP_Hap}.fa | head 

###################################
# make and move to unique directories
 tree

├── Argentina_anserina
│   ├── Argentina_anserina.fa
│   └── Argentina_anserina.gff
├── Dryas_octopetala
│   ├── Dryas_octopetala.fa
│   └── Dryas_octopetala.gff3
├── Fragaria_vesca
│   ├── Fragaria_vesca.fa
│   └── Fragaria_vesca.gff
├── Malus_sylvestris
│   ├── Malus_sylvestris.fa
│   └── Malus_sylvestris.gff
├── Prunus_persica
│   ├── Prunus_persica.fa
│   └── Prunus_persica.gff
├── Pyrus_bretschneideri
│   ├── Pyrus_bretschneideri.fa
│   └── Pyrus_bretschneideri.gff
└── Rosa_rugosa
    ├── Rosa_rugosa.fa
    └── Rosa_rugosa.gff


#####################################################
# Brassicaceae
cd /home/celphin/scratch/Oxyria/GeneSpace/Brassicaceae_genomes/genomes


mkdir Arabis_alpina; cd Arabis_alpina
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_000733195.1/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/733/195/GCA_000733195.1_A_alpina_V4/GCA_000733195.1_A_alpina_V4_genomic.gff.gz
gunzip GCA_000733195.1_A_alpina_V4_genomic.gff.gz
mv GCA_000733195.1_A_alpina_V4_genomic.gff  Arabis_alpina.gff

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/733/195/GCA_000733195.1_A_alpina_V4/GCA_000733195.1_A_alpina_V4_translated_cds.faa.gz
gunzip GCA_000733195.1_A_alpina_V4_translated_cds.faa.gz
mv  GCA_000733195.1_A_alpina_V4_translated_cds.faa Arabis_alpina.fa
cd ..

mkdir Arabidopsis_lyrata; cd Arabidopsis_lyrata # arctic species too
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000004255.2/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/004/255/GCF_000004255.2_v.1.0/GCF_000004255.2_v.1.0_genomic.gff.gz
gunzip GCF_000004255.2_v.1.0_genomic.gff.gz
mv GCF_000004255.2_v.1.0_genomic.gff Arabidopsis_lyrata.gff

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/004/255/GCF_000004255.2_v.1.0/GCF_000004255.2_v.1.0_translated_cds.faa.gz
gunzip GCF_000004255.2_v.1.0_translated_cds.faa.gz
mv GCF_000004255.2_v.1.0_translated_cds.faa Arabidopsis_lyrata.fa
cd ..

mkdir Arabidopsis_thaliana; cd Arabidopsis_thaliana
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001735.4/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.gff.gz
gunzip GCF_000001735.4_TAIR10.1_genomic.gff.gz
mv GCF_000001735.4_TAIR10.1_genomic.gff Arabidopsis_thaliana.gff

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_translated_cds.faa.gz
gunzip GCF_000001735.4_TAIR10.1_translated_cds.faa.gz
mv GCF_000001735.4_TAIR10.1_translated_cds.faa  Arabidopsis_thaliana.fa
cd ..

mkdir Capsella_rubella; cd Capsella_rubella
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000375325.1/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/375/325/GCF_000375325.1_Caprub1_0/GCF_000375325.1_Caprub1_0_genomic.gff.gz
gunzip GCF_000375325.1_Caprub1_0_genomic.gff.gz
mv GCF_000375325.1_Caprub1_0_genomic.gff Capsella_rubella.gff

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/375/325/GCF_000375325.1_Caprub1_0/GCF_000375325.1_Caprub1_0_translated_cds.faa.gz
gunzip GCF_000375325.1_Caprub1_0_translated_cds.faa.gz
mv GCF_000375325.1_Caprub1_0_translated_cds.faa Capsella_rubella.fa
cd ..

mkdir Brassica_oleracea; cd Brassica_oleracea
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000695525.1/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/695/525/GCF_000695525.1_BOL/GCF_000695525.1_BOL_genomic.gff.gz
gunzip GCF_000695525.1_BOL_genomic.gff.gz
mv GCF_000695525.1_BOL_genomic.gff Brassica_oleracea.gff

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/695/525/GCF_000695525.1_BOL/GCF_000695525.1_BOL_translated_cds.faa.gz
gunzip GCF_000695525.1_BOL_translated_cds.faa.gz
mv GCF_000695525.1_BOL_translated_cds.faa Brassica_oleracea.fa
cd ..

mkdir Thlaspi_arvense; cd Thlaspi_arvense
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_911865555.2/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/911/865/555/GCA_911865555.2_T_arvense_v2/GCA_911865555.2_T_arvense_v2_genomic.gff.gz
gunzip GCA_911865555.2_T_arvense_v2_genomic.gff.gz
mv GCA_911865555.2_T_arvense_v2_genomic.gff Thlaspi_arvense.gff

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/911/865/555/GCA_911865555.2_T_arvense_v2/GCA_911865555.2_T_arvense_v2_translated_cds.faa.gz
gunzip GCA_911865555.2_T_arvense_v2_translated_cds.faa.gz
mv GCA_911865555.2_T_arvense_v2_translated_cds.faa Thlaspi_arvense.fa
cd ..

mkdir Cochlearia_groenlandica; cd Cochlearia_groenlandica
# https://springernature.figshare.com/collections/Whole-genome_sequencing_of_13_Arctic_plants_and_draft_genomes_of_Oxyria_digyna_and_Cochlearia_groenlandica/6965802/1
# transfer by globus
rename cg Cochlearia_groenlandica *
cd ..

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

#ERROR    Line 1: FORMAT: "##gff-version" missing from the first line-> scaffold74   maker   gene    284456  285556  .       +       .       ID=snap_masked-scaffold74-processed-gene-0.99;Name=snap_masked-scaffold74-processed-gene-0.99
# worked

#-----------------------
# View genome tree

tree
.
├── Arabidopsis_lyrata
│   ├── Arabidopsis_lyrata.fa
│   └── Arabidopsis_lyrata.gff
├── Arabidopsis_thaliana
│   ├── Arabidopsis_thaliana.fa
│   └── Arabidopsis_thaliana.gff
├── Arabis_alpina
│   ├── Arabis_alpina.fa
│   └── Arabis_alpina.gff
├── Brassica_oleracea
│   ├── Brassica_oleracea.fa
│   └── Brassica_oleracea.gff
├── Capsella_rubella
│   ├── Capsella_rubella.fa
│   └── Capsella_rubella.gff
├── Cochlearia_groenlandica
│   ├── Cochlearia_groenlandica.fa
│   ├── Cochlearia_groenlandica.gff3
├── Draba_nivalis
│   ├── Draba_nivalis.fa
│   └── Draba_nivalis.gff
└── Thlaspi_arvense
    ├── Thlaspi_arvense.fa
    └── Thlaspi_arvense.gff


##############################
# For total genomes 

# Copy over all genomes from the three families above

# Rosaceae
cp -rv /home/celphin/scratch/Oxyria/GeneSpace/Rosaceae_genomes/genomes/* /home/celphin/scratch/Oxyria/GeneSpace/Total_genomes
# Polygonaceae
cp  -rv /home/celphin/scratch/Oxyria/GeneSpace/Polygonaceae_genomes/genomes/* /home/celphin/scratch/Oxyria/GeneSpace/Total_genomes
# Brassicaceae
cp  -rv  /home/celphin/scratch/Oxyria/GeneSpace/Brassicaceae_genomes/genomes/* /home/celphin/scratch/Oxyria/GeneSpace/Total_genomes

###################################
# check tree
cd /home/celphin/scratch/Oxyria/GeneSpace/Total_genomes

tree

├── Alnus_glutinosa
│   ├── Alnus_glutinosa.fa
│   └── Alnus_glutinosa.gff
├── Arabidopsis_lyrata
│   ├── Arabidopsis_lyrata.fa
│   └── Arabidopsis_lyrata.gff
├── Arabidopsis_thaliana
│   ├── Arabidopsis_thaliana.fa
│   └── Arabidopsis_thaliana.gff
├── Arabis_alpina
│   ├── Arabis_alpina.fa
│   └── Arabis_alpina.gff
├── Argentina_anserina
│   ├── Argentina_anserina.fa
│   └── Argentina_anserina.gff
├── Betula_nana
│   ├── Betula_nana.fa
│   └── Betula_nana.gff
├── Brassica_oleracea
│   ├── Brassica_oleracea.fa
│   └── Brassica_oleracea.gff
├── Capsella_rubella
│   ├── Capsella_rubella.fa
│   └── Capsella_rubella.gff
├── Carpinus_fangiana
│   ├── Carpinus_fangiana.fa
│   └── Carpinus_fangiana.gff
 Cochlearia_groenlandica
│   ├── Cochlearia_groenlandica.fa
│   └── Cochlearia_groenlandica.gtf
├── Corylus_avellana
│   ├── Corylus_avellana.fa
│   └── Corylus_avellana.gff
├── Draba_nivalis
│   ├── Draba_nivalis.fa
│   └── Draba_nivalis.gff
├── Dryas_octopetala
│   ├── Dryas_octopetala.fa
│   └── Dryas_octopetala.gff3
├── Fagopyrum_escelentum_H2
│   ├── Fagopyrum_escelentum_H2.fa
│   └── Fagopyrum_escelentum_H2.gff3
├── Fagopyrum_tataricum_H1
│   ├── Fagopyrum_tataricum_H1.fa
│   └── Fagopyrum_tataricum_H1.gff3
├── Fragaria_vesca
│   ├── Fragaria_vesca.fa
│   └── Fragaria_vesca.gff
├── Malus_sylvestris
│   ├── Malus_sylvestris.fa
│   └── Malus_sylvestris.gff
├── Oxyria_digyna_H1
│   ├── Oxyria_digyna_H1.fa
│   └── Oxyria_digyna_H1.gff3
├── Polygunum_aviculare_H0
│   ├── Polygunum_aviculare_H0.fa
│   └── Polygunum_aviculare_H0.gff3
├── Prunus_persica
│   ├── Prunus_persica.fa
│   └── Prunus_persica.gff
├── Pyrus_bretschneideri
│   ├── Pyrus_bretschneideri.fa
│   └── Pyrus_bretschneideri.gff
├── Rheum_nobile_H0
│   ├── Rheum_nobile_H0.fa
│   └── Rheum_nobile_H0.gff3
├── Rheum_tangaticum_H0
│   ├── Rheum_tangaticum_H0.fa
│   └── Rheum_tangaticum_H0.gff3
├── Rosa_rugosa
│   ├── Rosa_rugosa.fa
│   └── Rosa_rugosa.gff
└── Thlaspi_arvense
    ├── Thlaspi_arvense.fa
    └── Thlaspi_arvense.gff

####################################





















###################################
# to run

tmux new-session -s GeneSpace1
tmux attach-session -t GeneSpace1

module load StdEnv/2020 python/3.11.5 scipy-stack/2021a
module load StdEnv/2020 java/13.0.2
module load StdEnv/2020 r/4.2.2 glpk/5.0

export PATH=$PATH:/home/celphin/scratch/Oxyria/GeneSpace/OrthoFinder_source/orthofinder.py
alias orthofinder='python /home/celphin/scratch/Oxyria/GeneSpace/OrthoFinder_source/orthofinder.py'

R

library(GENESPACE)
gpar <- init_genespace(
  wd = "/home/celphin/scratch/Oxyria/GeneSpace/Polygonaceae_genomes", 
  path2mcscanx = "/home/celphin/scratch/Oxyria/GeneSpace/MCScanX/")
out <- run_genespace(gsParam = gpar) 

# Could not find a valid path to the orthofinder program from R. To run         
# orthofinder, ensure that the orthofinder program is in the         
# $PATH, then call the following from the shell: 
# orthofinder -f /home/celphin/scratch/Oxyria/GeneSpace/Polygonaceae_genomes/tmp -t 16 -a 1 -X -o /home/celphin/scratch/Oxyria/GeneSpace/orthofinder
# Error in run_orthofinder(gsParam = gsParam, verbose = TRUE) :
# Once OrthoFinder has been run, re-call run_genespace

######################################
# try complete rerun into orthoFinder using more time total

# bash
salloc -c16 --time 5:00:00 --mem 191000M --account def-cronk

module load StdEnv/2020  python/3.10.2 scipy-stack/2021a diamond/2.1.7 fastme/2.1.6.2 gcc/9.3.0 mcl/14.137 java/13.0.2 r/4.2.2 glpk/5.0

# other version of orthofinder
export PATH=$PATH:/home/celphin/scratch/Oxyria/GeneSpace/OrthoFinder/bin
export PATH=$PATH:/home/celphin/scratch/Oxyria/GeneSpace/OrthoFinder_source
alias orthofinder='orthofinder.py'

orthofinder -h

cd /home/celphin/scratch/Oxyria/GeneSpace/
orthofinder -f /home/celphin/scratch/Oxyria/GeneSpace/Polygonaceae_genomes/tmp -t 16 -a 1 -X -o /home/celphin/scratch/Oxyria/GeneSpace/Polygonaceae_genomes/orthofinder

# finished

#######################################
# Exploring results

# https://davidemms.github.io/orthofinder_tutorials/exploring-orthofinders-results.html
# download tree and view in http://etetoolkit.org/treeview/
more /home/celphin/scratch/Oxyria/GeneSpace/Polygonaceae_genomes/orthofinder/Results_Feb24/Species_Tree/SpeciesTree_rooted_node_labels.txt

#################################
# restart in R

tmux new-session -s GeneSpace1
tmux attach-session -t GeneSpace1
tmux kill-session -t GeneSpace1

# bash
salloc -c20 --time 2:40:00 --mem 191000M --account def-cronk

module load StdEnv/2020  python/3.10.2 scipy-stack/2021a diamond/2.1.7 fastme/2.1.6.2 gcc/9.3.0 mcl/14.137 java/13.0.2 r/4.2.2 glpk/5.0

# other version of orthofinder
export PATH=$PATH:/home/celphin/scratch/Oxyria/GeneSpace/OrthoFinder/bin
export PATH=$PATH:/home/celphin/scratch/Oxyria/GeneSpace/OrthoFinder_source
alias orthofinder='orthofinder.py'

#orthofinder -h
# works
cd /home/celphin/scratch/Oxyria/GeneSpace/Polygonaceae_genomes
R

library(GENESPACE)

# https://rdrr.io/github/jtlovell/GENESPACE/man/init_genespace.html
gpar <- init_genespace(
  wd = "/home/celphin/scratch/Oxyria/GeneSpace", 
  path2mcscanx = "/home/celphin/scratch/Oxyria/GeneSpace/MCScanX/",
  rawOrthofinderDir = "/home/celphin/scratch/Oxyria/GeneSpace/Polygonaceae_genomes/orthofinder/Results_Feb24")

# works

####################################
run_genespace(gsParam = gpar)
# runs until first set of dot plots


####################################################
# Try removing  dotplots directly in R
# https://github.com/jtlovell/GENESPACE/tree/master/R

# https://github.com/jtlovell/GENESPACE/blob/master/R/run_genespace.R

# try whole function without dotplots (below)
gpar <- init_genespace(
  wd = "/home/celphin/scratch/Oxyria/GeneSpace", 
  path2mcscanx = "/home/celphin/scratch/Oxyria/GeneSpace/MCScanX/",
  rawOrthofinderDir = "/home/celphin/scratch/Oxyria/GeneSpace/Polygonaceae_genomes/orthofinder/Results_Feb24")

library(data.table)
out <- run_genespace_2(gsParam = gpar)

############################
# GENESPACE run complete!  All results are stored in
# /home/celphin/scratch/Oxyria/GeneSpace in the following subdirectories:
        # syntenic block dotplots: /dotplots (...synHits.pdf)
        # annotated blast files  : /syntenicHits
        # annotated/combined bed : /results/combBed.txt
        # syntenic block coords. : /results/blkCoords.txt
        # syn. blk. by ref genome: /riparian/refPhasedBlkCoords.txt
        # pan-genome annotations : /pangenes (...pangenes.txt.gz)
        # riparian plots         : /riparian
        # genespace param. list  : /results/gsParams.rda
# ############################
# **NOTE** the genespace parameter object is returned or can be loaded
        # into R via
        # `load('/home/celphin/scratch/Oxyria/GeneSpace/results/gsParams.rda',
        # verbose = TRUE)`. Then you can customize your riparian plots by
        # calling `plot_riparian(gsParam = gsParam, ...)`. The source
        # data and ggplot2 objects are also stored in the /riparian
        # directory and can also be accessed by `load(...)`.
# **NOTE** To query genespace results by position or gene, use
        # `query_genespace(...)`. See specifications in ?query_genespace
        # for details.

##############################################
# more in depth use and plots:
# https://htmlpreview.github.io/?https://github.com/jtlovell/tutorials/blob/main/genespaceGuide.html

tmux new-session -s GeneSpace
tmux attach-session -t GeneSpace

cd /home/celphin/scratch/Oxyria/GeneSpace/

module load StdEnv/2020  python/3.10.2 scipy-stack/2021a diamond/2.1.7 fastme/2.1.6.2 gcc/9.3.0 mcl/14.137 java/13.0.2 r/4.2.2 glpk/5.0

# other version of orthofinder
export PATH=$PATH:/home/celphin/scratch/Oxyria/GeneSpace/OrthoFinder/bin
export PATH=$PATH:/home/celphin/scratch/Oxyria/GeneSpace/OrthoFinder_source
alias orthofinder='orthofinder.py'

#orthofinder -h
# works

R

library(GENESPACE)

# need to move to main Genespace folder to use
out <- load('/home/celphin/scratch/Oxyria/GeneSpace/results/gsParams.rda', verbose = TRUE)

# https://htmlpreview.github.io/?https://github.com/jtlovell/tutorials/blob/main/riparianGuide.html
# remake riparian plot

# gsParam$genomeIDs
 # [1] "Oxyria_digyna_H1"        "Oxyria_digyna_H0"
 # [3] "Rheum_nobile_H0"         "Rheum_tangaticum_H0"
 # [5] "Polygunum_aviculare_H0"  "Fagopyrum_tataricum_H1"
 # [7] "Fagopyrum_tataricum_H2"  "Fagopyrum_tataricum_H0"
 # [9] "Fagopyrum_escelentum_H1" "Fagopyrum_escelentum_H2"

# invert some chromosomes
invchr <- data.frame(
  genome = c("Fagopyrum_escelentum_H2", "Fagopyrum_escelentum_H2", "Fagopyrum_escelentum_H2",
			"Fagopyrum_escelentum_H1",
			"Fagopyrum_tataricum_H2", "Fagopyrum_tataricum_H2",
			"Fagopyrum_tataricum_H1", "Fagopyrum_tataricum_H1",
			"Polygunum_aviculare_H0", "Polygunum_aviculare_H0", "Polygunum_aviculare_H0", "Polygunum_aviculare_H0", "Polygunum_aviculare_H0",
			"Rheum_tangaticum_H0", "Rheum_tangaticum_H0","Rheum_tangaticum_H0","Rheum_tangaticum_H0","Rheum_tangaticum_H0","Rheum_tangaticum_H0","Rheum_tangaticum_H0",
			 "Rheum_nobile_H0", "Rheum_nobile_H0"),

 chr = c("2", "7", "5",
		"5",
		"4", "6",
		"1", "4",
		"polavi-7", "polavi-10", "polavi-9", "polavi-3", "polavi-6",
		"rta02", "rta05", "rta06", "rta01", "rta07", "rta08", "rta04",
		"rno02", "rno06"))

ggthemes <- ggplot2::theme(
  panel.background = ggplot2::element_rect(fill = "white"))
customPal <- colorRampPalette(
  c("darkorange", "skyblue", "darkblue", "purple", "darkred", "salmon"))
  
  
grDevices::pdf("/home/celphin/scratch/Oxyria/GeneSpace/Polygonaceae_genomes/OxydigH1_riparian.pdf")
print(plot_riparian(
  gsParam = gsParam, 
  highlightBed = roi, 
  backgroundColor = NULL,
  useOrder = TRUE, useRegions = FALSE, reorderBySynteny = TRUE,
  genomeIDs = c("Polygunum_aviculare_H0","Fagopyrum_tataricum_H1", "Fagopyrum_escelentum_H2",  "Rheum_tangaticum_H0", "Rheum_nobile_H0", "Oxyria_digyna_H1"),
  refGenome = "Oxyria_digyna_H1",   
  palette = customPal,
  braidAlpha = .75,
  chrFill = "lightgrey",
  addThemes = ggthemes))
grDevices::dev.off()

#----------------------
roi <- data.frame(
  genome = c("Oxyria_digyna_H1"), 
  chr = c("Oxyrt-1-86582034", "Oxyrt-2-79714091", "Oxyrt-3-79472951","Oxyrt-4-78410798","Oxyrt-5-76064323","Oxyrt-6-73303751","Oxyrt-7-72361354"))
  #chr = c("Oxyrt-4-78410798","Oxyrt-5-76064323","Oxyrt-6-73303751"))
  
ggthemes <- ggplot2::theme(
  panel.background = ggplot2::element_rect(fill = "white"))
customPal <- colorRampPalette(
  c("darkorange", "skyblue", "darkblue", "purple", "darkred", "salmon"))
  
grDevices::pdf("/home/celphin/scratch/Oxyria/GeneSpace/OxydigH1_riparian_4spp.pdf")
print(plot_riparian(
  gsParam = gsParam, 
  useOrder = FALSE, useRegions = FALSE, reorderBySynteny = TRUE,
  highlightBed = roi, 
  backgroundColor = NULL,
  palette = customPal,
  braidAlpha = .75,
  chrFill = "lightgrey",
  addThemes = ggthemes,
  genomeIDs = c( "Fagopyrum_escelentum_H2",  "Polygunum_aviculare_H0", "Rheum_nobile_H0", "Oxyria_digyna_H1"),
  refGenome = "Oxyria_digyna_H1"))
grDevices::dev.off()

grDevices::pdf("/home/celphin/scratch/Oxyria/GeneSpace/Polygonaceae_genomes/OxydigH1_riparian_Rnob.pdf")
print(plot_riparian(
  gsParam = gsParam, 
  useOrder = FALSE, useRegions = FALSE, reorderBySynteny = TRUE,
  highlightBed = roi, 
  backgroundColor = NULL,
  palette = customPal,
  braidAlpha = .75,
  chrFill = "lightgrey",
  addThemes = ggthemes,
  genomeIDs = c("Rheum_nobile_H0","Oxyria_digyna_H1", "Polygunum_aviculare_H0"),
  refGenome = "Oxyria_digyna_H1"))
grDevices::dev.off()

grDevices::pdf("/home/celphin/scratch/Oxyria/GeneSpace/Polygonaceae_genomes/OxydigH1_riparian_Polavi.pdf")
print(plot_riparian(
  gsParam = gsParam, 
  useOrder = FALSE, useRegions = FALSE, reorderBySynteny = TRUE,
  highlightBed = roi, 
  backgroundColor = NULL,
  palette = customPal,
  braidAlpha = .75,
  chrFill = "lightgrey",
  addThemes = ggthemes,
  genomeIDs = c("Oxyria_digyna_H1", "Polygunum_aviculare_H0"),
  refGenome = "Oxyria_digyna_H1"))
grDevices::dev.off()

grDevices::pdf("/home/celphin/scratch/Oxyria/GeneSpace/Polygonaceae_genomes/OxydigH1_H0_riparian.pdf")
print(plot_riparian(
  gsParam = gsParam, 
  useOrder = FALSE, useRegions = FALSE, reorderBySynteny = TRUE,
  highlightBed = roi, 
  backgroundColor = NULL,
  palette = customPal,
  braidAlpha = .75,
  chrFill = "lightgrey",
  addThemes = ggthemes,
  genomeIDs = c("Oxyria_digyna_H1", "Oxyria_digyna_H0"),
  refGenome = "Oxyria_digyna_H1"))
grDevices::dev.off()

#########################################
# plot specific chromosomes
R

library(GENESPACE)

out <- load('/home/celphin/scratch/Oxyria/GeneSpace/results/gsParams.rda', verbose = TRUE)

roi <- data.frame(
  genome = c("Oxyria_digyna_H1"), 
  chr = c("Oxyrt-6-73303751"), 
  color = c("#17B5C5"))

# invert some chromosomes
invchr <- data.frame(
  genome = c("Rheum_nobile_H0", "Rheum_nobile_H0", "Rheum_nobile_H0"),
  chr = c("rno05", "rno08", "rno09"))

grDevices::pdf("/home/celphin/scratch/Oxyria/GeneSpace/OxydigH1_riparian_subsetchr6only.pdf")
print(
ripDat <- plot_riparian(
  gsParam = gsParam, 
  highlightBed = roi, 
  useOrder = FALSE, useRegions = FALSE, #reorderBySynteny = TRUE,
  backgroundColor = NULL,
  invertTheseChrs = invchr,  
  genomeIDs = c("Oxyria_digyna_H1","Rheum_tangaticum_H0", "Rheum_nobile_H0", "Polygunum_aviculare_H0", "Fagopyrum_escelentum_H2", "Fagopyrum_tataricum_H1"),
  refGenome = "Oxyria_digyna_H1"))
grDevices::dev.off()


grDevices::pdf("/home/celphin/scratch/Oxyria/GeneSpace/OxydigH1_riparian_subsetchr6only_geneorder.pdf")
print(
ripDat <- plot_riparian(
  gsParam = gsParam, 
  highlightBed = roi, 
  useOrder = TRUE, useRegions = TRUE, reorderBySynteny = TRUE,
  backgroundColor = NULL, 
  invertTheseChrs = invchr,
  genomeIDs = c("Oxyria_digyna_H1","Rheum_tangaticum_H0", "Rheum_nobile_H0", "Polygunum_aviculare_H0", "Fagopyrum_escelentum_H2", "Fagopyrum_tataricum_H1"),
  refGenome = "Oxyria_digyna_H1"))
grDevices::dev.off()

#---------------------------------
# plot region of interest
roibed <- roi[,c("genome", "chr")]
roibed$color <- c("pink", "cyan", "green")
ripd <- plot_riparian(
  gsParam = gsParam, 
  useRegions = FALSE, 
  highlightBed = roibed)

# find pangenes
# genes in same orthogroup and are found together
test <- query_pangenes(
  gsParam = gsParam, bed = roi)


##############################


  
  

          
                 chrLab
 1:                    7
 2:                    2
 3:                    1
 4:                    8
 5:                    4
 6:                    5
 7:                    6
 8:                    3
 9:                    7
10:                    2
11:                    8
12:                    1
13:                    4
14:                    5
15:                    6
16:                    3
17: fagopyrum-6-52287906
18: fagopyrum-2-61235386
19: fagopyrum-1-68031765
20: fagopyrum-4-56655744
21: fagopyrum-3-57706077
22: fagopyrum-8-49982843
23: fagopyrum-5-53883329
24: fagopyrum-7-51545819
25:                    7
26:                    3
27:                    2
28:                    1
29:                    6
30:                    4
31:                    8
32:                    5
33:                    7
34:                    3
35:                    2
36:                    1
37:                    4
38:                    8
39:                    6
40:                    5
41:       oxy-1-88146246
42:       oxy-2-80787722
43:       oxy-3-77175000
44:       oxy-4-76036264
45:        oxy-9-3540000
46:       oxy-7-70136240
47:       oxy-11-2000000
48:       oxy-12-1350245
49:       oxy-5-73842794
50:       oxy-6-73289246
51:     oxyrt-1-86582034
52:     oxyrt-2-79714091
53:     oxyrt-3-79472951
54:     oxyrt-4-78410798
55:     oxyrt-5-76064323
56:     oxyrt-6-73303751
57:     oxyrt-7-72361354
58:             polavi-8
59:             polavi-5
60:             polavi-7
61:            polavi-10
62:             polavi-2
63:             polavi-1
64:             polavi-9
65:             polavi-4
66:             polavi-3
67:             polavi-6
68:                rno02
69:                rno07
70:                rno03
71:                rno05
72:                rno10
73:                rno01
74:                rno06
75:                rno09
76:                rno08
77:                rno11
78:                rno04
79:                rta03
80:                rta09
81:                rta02
82:                rta05
83:                rta06
84:                rta01
85:                rta07
86:                rta10
87:                rta08
88:                rta11
89:                rta04



# https://htmlpreview.github.io/?https://github.com/jtlovell/tutorials/blob/main/genomeVizGuide.html

#########################################
run_genespace_2 <- function(gsParam,
                          overwrite = FALSE,
                          overwriteBed = overwrite,
                          overwriteSynHits = overwrite,
                          overwriteInBlkOF = TRUE,
                          makePairwiseFiles = FALSE){

   gsParam$paths$rawOrthofinder <- gsParam$paths$orthofinder
  ##############################################################################
  # 1. Run orthofinder ...
  cat("\n############################", strwrap(
    "1. Running orthofinder (or parsing existing results)",
    indent = 0, exdent = 8), sep = "\n")

  ##############################################################################
  # -- 1.1 Check for existing parsed orthofinder results
  cat("\tChecking for existing orthofinder results ...\n")
  gsParam <- set_syntenyParams(gsParam)

  if("synteny" %in% names(gsParam)){
    noResults <- is.na(gsParam$synteny$SpeciesIDs)
  }else{
    noResults <- TRUE
  }

  ##############################################################################
  # -- 1.2 If no results exist, check for raw orthofinder run
  if(noResults){
    if(dir.exists(gsParam$paths$rawOrthofinder)){
      chkOf <- find_ofFiles(gsParam$paths$rawOrthofinder)
      noOrthofinder <- is.na(chkOf[[1]])
    }else{
      noOrthofinder <- TRUE
    }
  }else{
    noOrthofinder <- FALSE
  }

  ##############################################################################
  # -- 1.3 If raw results exist, copy them over
  if(!noOrthofinder && noResults){
    with(gsParam, copy_of2results(
      orthofinderDir = paths$rawOrthofinder, resultsDir = paths$results))

    if(dir.exists(gsParam$paths$rawOrthofinder)){
      chkOf <- find_ofFiles(gsParam$paths$rawOrthofinder)
      noOrthofinder <- is.na(chkOf[[1]])
      noResults <- noOrthofinder
    }else{
      noOrthofinder <- TRUE
    }
  }

  print(noResults)

  if(!noResults){
    spids <- names(read_orthofinderSpeciesIDs(
      file.path(gsParam$paths$results, "SpeciesIDs.txt")))
    gid <- unique(c(gsParam$genomeIDs, gsParam$outgroup))
    gid <- gid[!is.na(gid)]
    ps <- all(gid %in% spids) && all(spids %in% gid)
    if(ps){
      cat("\t... found existing run, not re-running orthofinder\n")
    }else{
      stop("genomes in the existing orthofinder run do not exactly match specified genomeIDs\n")
    }
  }


  ##############################################################################
  # -- 1.4 if no orthofinder run, make one
  if(noResults)
    tmp <- run_orthofinder(gsParam = gsParam, verbose = TRUE)

  gsParam <- set_syntenyParams(gsParam)
  gsParam <<- gsParam
  ##############################################################################
  # -- 1.5 get the files in order if the run is complete
  if(noResults){
    chkOf <- find_ofFiles(gsParam$paths$orthofinder)
    noOrthofinder <- is.na(chkOf[[1]])
    if(noOrthofinder)
      stop("could not find orthofinder files!")
    with(gsParam, copy_of2results(
      orthofinderDir = paths$orthofinder,
      resultsDir = paths$results))
  }
  gsParam <- set_syntenyParams(gsParam)

  ##############################################################################
  # -- 1.6 if the species tree exists, re-order the genomeIDs
  tmp <- gsParam$synteny$speciesTree

  if(requireNamespace("ape", quietly = T)){
    if(!is.na(tmp) && !is.null(tmp)){
      if(file.exists(tmp) && length(gsParam$genomeIDs) > 2){
        treLabs <- get_orderedTips(
          treFile = gsParam$synteny$speciesTree,
          ladderize = TRUE,
          genomeIDs = gsParam$genomeIDs)
        cat(strwrap(sprintf(
          "re-ordering genomeIDs by the species tree: %s",
          paste(treLabs, collapse = ", ")), indent = 8, exdent = 16),
          sep = "\n")
        gsParam$genomeIDs <- treLabs
      }
    }
  }

  # -- 1.7 if useHOGs, check if the N0.tsv file exists
  useHOGs <- gsParam$params$useHOGs
  if(useHOGs){
    if(is.na(gsParam$synteny$hogs)){
      useHOGs <- FALSE
    }else{
      if(!file.exists(gsParam$synteny$hogs)){
        useHOGs <- FALSE
      }
    }
  }
  gsParam$params$useHOGs <- useHOGs
  useHOGs <- NULL

  ##############################################################################
  # 2. Get the data ready for synteny
  hasBed <- FALSE
  bedf <- gsParam$synteny$combBed
  if(file.exists(bedf))
    hasBed <- is.data.table(read_combBed(bedf))
  if(overwriteBed)
    hasBed <- FALSE

  if(!hasBed){
    cat("\n############################", strwrap(
      "2. Combining and annotating bed files w/ OGs and tandem array info ... ",
      indent = 0, exdent = 8), sep = "\n")
    bed <- annotate_bed(gsParam = gsParam)
  }else{
    cat("\n############################", strwrap(
      "2. Annotated/concatenated bed file exists", indent = 0, exdent = 8),
      sep = "\n")
  }

  ##############################################################################
  # 3. Annotate the blast files ...
  # -- First make sure that the blast files are all there, then go through
  # and annotate them with the combined bed file
  # -- This also makes the first round of dotplots
  # -- 3.1 check if all the synHits files exist. If so, and !overwriteSynHits
  # don't re-annotate
  hasHits <- FALSE
  if(all(file.exists(gsParam$synteny$blast$synHits)))
    if(!overwriteSynHits)
      hasHits <- TRUE

  # -- 3.2 iterate through and annotate all synHits files
  if(!hasHits){
    cat("\n############################", strwrap(
      "3. Combining and annotating the blast files with orthogroup info ...",
      indent = 0, exdent = 8), sep = "\n")
    gsParam <- annotate_blast(gsParam = gsParam)
  }else{
    cat("\n############################", strwrap(
      "3. Annotated/blast files exists", indent = 0, exdent = 8),
      sep = "\n")
  }

  dpFiles <- with(gsParam$synteny$blast, file.path(
    file.path(gsParam$paths$wd, "dotplots",
    sprintf("%s_vs_%s.rawHits.pdf",
            query, target))))
# if(!all(file.exists(dpFiles)) || overwrite){
    # cat("\t##############\n\tGenerating dotplots for all hits ... ")
    # nu <- plot_hits(gsParam = gsParam, type = "raw")
    # cat("Done!\n")
  # }

  ##############################################################################
  # 4. Run synteny
  # -- goes through each pair of genomes and pulls syntenic anchors and the hits
  # nearby. This is the main engine of genespace
  bed <- read_combBed(bedf)
  hasSynFiles <- all(file.exists(gsParam$synteny$blast$synHits))
  hasOg <- all(!is.na(bed$og))
  if(!hasOg || !hasSynFiles || overwrite){
    cat("\n############################", strwrap(
      "4. Flagging synteny for each pair of genomes ...",
      indent = 0, exdent = 8), sep = "\n")
    gsParam <- synteny(gsParam = gsParam)
  }

  ##############################################################################
  # 5. Build syntenic orthogroups
  cat("\n############################", strwrap(
    "5. Building synteny-constrained orthogroups ... ",
    indent = 0, exdent = 8), sep = "\n")
  # -- in the case of polyploid genomes, this also runs orthofinder in blocks,
  # then re-runs synteny and re-aggregates blocks,.
  if(gsParam$params$orthofinderInBlk & overwriteInBlkOF){

    # -- returns the gsparam obj and overwrites the bed file with a new og col
    cat("\t##############\n\tRunning Orthofinder within syntenic regions\n")
    tmp <- run_orthofinderInBlk(
      gsParam = gsParam, overwrite = overwriteSynHits)

    # -- adds a new column to the bed file
    cat("\t##############\n\tPulling syntenic orthogroups\n")
  }

  gsParam <- syntenic_orthogroups(
    gsParam, updateArrays = gsParam$params$orthofinderInBlk)
  cat("\tDone!\n")
  gsParam <<- gsParam

  ##############################################################################
  # 6. Make dotplots
  cat("\n############################", strwrap(
    "6. Integrating syntenic positions across genomes ... ",
    indent = 0, exdent = 8), sep = "\n")

  dpFiles <- with(gsParam$synteny$blast, file.path(
    file.path(gsParam$paths$wd, "dotplots",
              sprintf("%s_vs_%s.synHits.pdf",
                      query, target))))
  # if(!all(file.exists(dpFiles)) || overwrite){
    # cat("\t##############\n\tGenerating syntenic dotplots ... ")
    # nu <- plot_hits(gsParam = gsParam, type = "syntenic")
    # cat("Done!\n")
  # }

  ##############################################################################
  # 7. Interpolate syntenic positions
  cat("\t##############\n\tInterpolating syntenic positions of genes ... \n")
  nsynPos <- interp_synPos(gsParam)
  cat("\tDone!\n")
  gsParam <<- gsParam

  ##############################################################################
  # 8. Phase syntenic blocks against reference chromosomes
  cat("\n############################", strwrap(
    "7. Final block coordinate calculation and riparian plotting ... ",
    indent = 0, exdent = 8), sep = "\n")
  cat("\t##############\n\tCalculating syntenic blocks by reference chromosomes ... \n")
  reg <- nophase_blks(gsParam = gsParam, useRegions = T)
  cat(sprintf("\t\tn regions (aggregated by %s gene radius): %s\n",
              gsParam$params$blkRadius, nrow(reg)))
  fwrite(reg, file = file.path(gsParam$paths$results, "syntenicRegion_coordinates.csv"))
  blk <- nophase_blks(gsParam = gsParam, useRegions = F)
  cat(sprintf("\t\tn blocks (collinear sets of > %s genes): %s\n",
              gsParam$params$blkSize, nrow(blk)))
  fwrite(blk, file = file.path(gsParam$paths$results, "syntenicBlock_coordinates.csv"))
  hapGenomes <- names(gsParam$ploidy)[gsParam$ploidy == 1]

  bed <- read_combBed(file.path(gsParam$paths$results, "combBed.txt"))
  minChrSize <- gsParam$params$blkSize * 2
  isOK <- og <- chr <- genome <- NULL
  bed[,isOK := uniqueN(og) >= minChrSize, by = c("genome", "chr")]
  ok4pg <- bed[,list(propOK = (sum(isOK) / .N) > .75), by = "genome"]
  ok4rip <- bed[,list(nGood = (uniqueN(chr[isOK]) < 100)), by = "genome"]

  if(length(hapGenomes) == 0){
    cat(strwrap("NOTE!!! No genomes provided with ploidy < 2. Phasing of polyploid references is not currently supported internally. You will need to make custom riparian plots",
                indent = 8, exdent = 8), sep = "\n")
  }else{
    okg <- subset(ok4rip, genome %in% hapGenomes)
    if(any(!okg$nGood)){
      cat(strwrap(sprintf(
        "**WARNING**: genomes %s have > 100 chromosomes in the synteny map. This is too complex to make riparian plots.\n",
        paste(okg$genome[!okg$nGood], collapse = ", ")), indent = 8, exdent = 16),
        sep = "\n")
      hapGenomes <- hapGenomes[!hapGenomes %in% okg$genome[!okg$nGood]]
    }
    if(length(hapGenomes) > 0){
      cat("\t##############\n\tBuilding ref.-phased blks and riparian plots for haploid genomes:\n")
      labs <- align_charLeft(hapGenomes)
      names(labs) <- hapGenomes
      for(i in hapGenomes){
        plotf <- file.path(gsParam$paths$riparian,
                           sprintf("%s_geneOrder.rip.pdf", i))

        srcf <- file.path(gsParam$paths$riparian,
                          sprintf("%s_geneOrder_rSourceData.rda", i))
        blkf <- file.path(gsParam$paths$riparian,
                          sprintf("%s_phasedBlks.csv", i))

        rip <- plot_riparian(
          gsParam = gsParam, useRegions = TRUE, refGenome = i, pdfFile = plotf)
        cat(sprintf("\t\t%s: %s phased blocks\n", labs[i], nrow(rip$blks)))

        srcd <- rip$plotData
        save(srcd, file = srcf)
        fwrite(rip$blks, file = blkf)

        plotf <- file.path(gsParam$paths$riparian,
                           sprintf("%s_bp.rip.pdf", i))
        srcf <- file.path(gsParam$paths$riparian,
                          sprintf("%s_bp_rSourceData.rda", i))
        rip <- plot_riparian(
          gsParam = gsParam, useOrder = FALSE, useRegions = TRUE,
          refGenome = i, pdfFile = plotf)
        srcd <- rip$plotData
        save(srcd, file = srcf)
      }
      cat("\tDone!\n")
    }
  }

  gsParam <<- gsParam

  ##############################################################################
  # 9. Build pan-genes (aka pan-genome annotations)
  cat("\n############################", strwrap(
    "8. Constructing syntenic pan-gene sets ... ",
    indent = 0, exdent = 8), sep = "\n")
  gids <- gsParam$genomeIDs
  tp <- paste(ok4pg$genome[!ok4pg$propOK], collapse = ", ")
  if(any(!ok4pg$propOK))
    cat(strwrap(sprintf(
      "**WARNING**: genomes %s have < 75%% of genes on chromosomes that contain > %s genes. Synteny is not a useful metric for these genomes. Be very careful with your pan-gene sets.\n",
      tp, minChrSize), indent = 8, exdent = 16),
      sep = "\n")
  labs <- align_charLeft(gids)
  names(labs) <- gids
  for(i in gids){
    cat(sprintf("\t%s: ", labs[i]))
    pgref <- syntenic_pangenes(gsParam = gsParam, refGenome = i)
    with(pgref, cat(sprintf(
      "n pos. = %s, synOgs = %s, array mem. = %s, NS orthos %s\n",
      uniqueN(pgID), sum(flag == "PASS"), sum(flag == "array"), sum(flag == "NSOrtho"))))
  }



  ##############################################################################
  # 8. Print summaries and return

  # --- make pairwise files
  if(makePairwiseFiles){
      cat("Building pairiwse hits files ...\n")
    pull_pairwise(gsParam, verbose = TRUE)
      cat("Done!")
  }

  gpFile <- file.path(gsParam$paths$results, "gsParams.rda")
  cat("\n############################", strwrap(sprintf(
    "GENESPACE run complete!\n All results are stored in %s in the following subdirectories:",
    gsParam$paths$wd), indent = 0, exdent = 0),
    "\tsyntenic block dotplots: /dotplots (...synHits.pdf)",
    "\tannotated blast files  : /syntenicHits",
    "\tannotated/combined bed : /results/combBed.txt",
    "\tsyntenic block coords. : /results/blkCoords.txt",
    "\tsyn. blk. by ref genome: /riparian/refPhasedBlkCoords.txt",
    "\tpan-genome annotations : /pangenes (...pangenes.txt.gz)",
    "\triparian plots         : /riparian",
    "\tgenespace param. list  : /results/gsParams.rda",
    "############################",
    strwrap(sprintf(
      "**NOTE** the genespace parameter object is returned or can be loaded
      into R via `load('%s', verbose = TRUE)`. Then you can customize your
      riparian plots by calling `plot_riparian(gsParam = gsParam, ...)`. The
      source data and ggplot2 objects are also stored in the /riparian
      directory and can also be accessed by `load(...)`. ",
      gpFile), indent = 0, exdent = 8),
    strwrap(
      "**NOTE** To query genespace results by position or gene,
      use `query_genespace(...)`. See specifications in ?query_genespace for
      details.",  indent = 0, exdent = 8),
    "############################",
    sep = "\n")
  save(gsParam, file = gpFile)
  return(gsParam)
}


###############################################
# ExampleData

# -- change paths to those valid on your system
genomeRepo <- "~/path/to/store/rawGenomes"
wd <- "~/path/to/genespace/workingDirectory"
path2mcscanx <- "~/path/to/MCScanX/"
###############################################

# -- download raw data from NCBI for human and chicken genomes
dir.create(genomeRepo)
rawFiles <- download_exampleData(filepath = genomeRepo)

# -- parse the annotations to fastas with headers that match a gene bed file
parsedPaths <- parse_annotations(
  rawGenomeRepo = genomeRepo,
  genomeDirs = c("human", "chicken"),
  genomeIDs = c("human", "chicken"),
  presets = "ncbi",
  genespaceWd = wd)

# -- initalize the run and QC the inputs
gpar <- init_genespace(
  wd = wd, 
  path2mcscanx = path2mcscanx)

# -- accomplish the run
out <- run_genespace(gpar)
