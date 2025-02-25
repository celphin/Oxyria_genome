#############################################################
#January 2024
#Cleaned version of all genome fasta downloads and links
#See Notes_May2023_Oxyria_Orthofinder.bash
    #For original downloads, as well as throughout Annotation Notes

#############################################################
# For our annotation: Fagopyrum tataricum 2n = 2x = 16
# ~/projects/def-rieseber/Dryas_shared_data/Oxyria/Polygonaceae_Genomes_Annotations/MS_CE_Fagopyrum_tataricum/Fagopyrum_tataricum_Main.fasta
# https://www.ncbi.nlm.nih.gov/assembly/GCA_002319775.1

#Download genome file
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/319/775/GCA_002319775.1_Ft1.0/GCA_002319775.1_Ft1.0_genomic.fna.gz
gunzip GCA_002319775.1_Ft1.0_genomic.fna.gz
mv GCA_002319775.1_Ft1.0_genomic.fna Fagopyrum_tataricum_genome.fasta

###########################################################
# Polygonum aviculare n = x = 10
# ~/projects/def-rieseber/Dryas_shared_data/Oxyria/Polygonaceae_Genomes_Annotations/Polygonum_aviculare
# https://www.ncbi.nlm.nih.gov/assembly/GCA_934048045.1

#Download fna file
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/934/048/045/GCA_934048045.1_dcPolAvic1.1/GCA_934048045.1_dcPolAvic1.1_genomic.fna.gz
gunzip GCA_934048045.1_dcPolAvic1.1_genomic.fna.gz
mv GCA_934048045.1_dcPolAvic1.1_genomic.fna PolAvi_genome.fasta
#############################################################
# Available annotation: Fagopyrum tataricum 2n = 2x = 16
# ~/projects/def-rieseber/Dryas_shared_data/Oxyria/Polygonaceae_Genomes_Annotations/Other_Fagopyrum_tataricum/Fagopyrum_tataricum.fasta
# https://figshare.com/articles/dataset/Fagopyrum_genome_data/21617562/2
wget https://figshare.com/ndownloader/files/39865297
mv 39865297 F_tataricum_H1.fasta
wget https://figshare.com/ndownloader/files/39867931
mv 39867931 F_tataricum_H1.proteins.fasta
wget https://figshare.com/ndownloader/files/39867925
mv 39867925 F_tataricum_H1.genes.cds.fasta
wget https://figshare.com/ndownloader/files/39867928
mv 39867928 F_tataricum_H1.gff3

#Haploid 2:
wget https://figshare.com/ndownloader/files/39865312
mv 39865312 F_tataricum_H2.fasta
wget https://figshare.com/ndownloader/files/39867940
mv 39867940 F_tataricum_H2.proteins.fasta
wget https://figshare.com/ndownloader/files/39867934
mv 39867934 F_tataricum_H2.genes.cds.fasta
wget https://figshare.com/ndownloader/files/39867937
mv 39867937 F_tataricum_H2.gff3


#############################################################
#Fagopyrum esculentum : 
#Haploid 1,2
#~/projects/def-rieseber/Dryas_shared_data/Oxyria/Polygonaceae_Genomes_Annotations/Fagopyrum_escelentum/
# https://figshare.com/articles/dataset/Fagopyrum_genome_data/21617562/2

#Haploid 1:
wget https://figshare.com/ndownloader/files/39865477
mv 39865477 F_escelentum_H1.fasta
wget https://figshare.com/ndownloader/files/39867913
mv 39867913 F_escelentum_H1.proteins.fasta
wget https://figshare.com/ndownloader/files/39867907
mv 39867907 F_escelentum_H1.genes.cds.fasta
wget https://figshare.com/ndownloader/files/39867910
mv 39867910 F_escelentum_H1.gff3

#Haploid 2:
wget https://figshare.com/ndownloader/files/39866572
mv 39866572 F_escelentum_H2.fasta
wget https://figshare.com/ndownloader/files/39867922
mv 39867922 F_escelentum_H2.proteins.fasta
wget https://figshare.com/ndownloader/files/39867916
mv 39867916 F_escelentum_H2.genes.cds.fasta
wget https://figshare.com/ndownloader/files/39867919
mv 39867919 F_escelentum_H2.gff3


#################################################
#Rheum nobile: (not chromosome level)
#~/projects/def-rieseber/Dryas_shared_data/Oxyria/Polygonaceae_Genomes_Annotations/Rheum_nobile
#https://figshare.com/articles/dataset/Rheum_nobile_genome/19662933

#Rno.genome.Chr.fa.gz
wget https://figshare.com/ndownloader/files/38250711
mv 38250711 Rno.genome.Chr.fa.gz
gunzip Rno.genome.Chr.fa.gz
mv Rno.genome.Chr.fa R_nobile.fasta

#Rno.genomic.Chr.gff.gz
wget https://figshare.com/ndownloader/files/38246958
mv 38246958 Rno.genomic.Chr.gff.gz
gunzip Rno.genomic.Chr.gff.gz
mv Rno.genomic.Chr.gff R_nobile.gff3

###############################################
#Oxyria digyna (not chromosome level)
#https://www.ncbi.nlm.nih.gov/assembly/GCA_029168935.1#/st
#
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/029/168/935/GCA_029168935.1_ASM2916893v1/GCA_029168935.1_ASM2916893v1_genomic.fna.gz
gunzip GCA_029168935.1_ASM2916893v1_genomic.fna.gz
mv GCA_029168935.1_ASM2916893v1_genomic.fna Oxyria_digyna.fasta
###############################################
#Rheum tangaticum:
#https://figshare.com/articles/dataset/Rheum_tanguticum_Genome/19663062
#Download of proteins also in Notes5_Oxyria_Annotation_Rerun.bash
wget https://figshare.com/ndownloader/files/41810904
mv 41810904 R_tangaticum.fasta.gz
gunzip 41910904
wget https://figshare.com/ndownloader/files/39014225
mv 39014225 R_tangaticum.gff3
wget https://figshare.com/ndownloader/files/39014219
mv 39014219 R_tangaticum.genes.cds.fasta
################################################
#Rumex Hastalus:
#https://datadryad.org/stash/dataset/doi:10.5061/dryad.s7h44j14h

wget https://datadryad.org/stash/downloads/file_stream/630777
mv 630777 R_hastalus_linkage.fasta
wget https://datadryad.org/stash/downloads/file_stream/630774
mv 630774 R_hastalus.gff3