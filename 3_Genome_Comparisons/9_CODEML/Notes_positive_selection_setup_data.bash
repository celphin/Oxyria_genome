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
cp /home/celphin/scratch/Oxyria/Polygonaceae_Genomes_Annotations/NCBI_Oxyria_digyna/Oxyria_NCBI.AED_0.6.genes.fasta .
cp /home/celphin/scratch/Oxyria/Polygonaceae_Genomes_Annotations/OTHER_Fagopyrum_tataricum/F_tataricum_H1.genes.fasta .
cp /home/celphin/scratch/Oxyria/Polygonaceae_Genomes_Annotations/OTHER_Fagopyrum_tataricum/F_tataricum_H2.genes.fasta .
cp /home/celphin/scratch/Oxyria/Polygonaceae_Genomes_Annotations/Polygonum_aviculare/Polavi.AED_0.6.genes.fasta .
cp /home/celphin/scratch/Oxyria/Polygonaceae_Genomes_Annotations/Rheum_nobile/R_nobile.genes.fasta .
cp /home/celphin/scratch/Oxyria/Polygonaceae_Genomes_Annotations/Rheum_tangaticum/R_tangaticum.genes.fasta .


#-------------------------
# rename files

mkdir Oxyria_digyna_H1; cd Oxyria_digyna_H1
mv ../Oxyria_NCBI.AED_0.6.genes.gasta  Oxyria_digyna_H1.fna
cd ..

mkdir Polygunum_aviculare_H0; cd Polygunum_aviculare_H0
mv ../Pol-avi.AED_0.6.genes.fasta  Polygunum_aviculare_H0.fna
cd ..

mkdir Rheum_nobile_H0; cd Rheum_nobile_H0
mv ../R_nobile.genes.cds.fasta  Rheum_nobile_H0.fna
cd ..

mkdir Rheum_tangaticum_H0; cd Rheum_tangaticum_H0
mv ../R_tangaticum.genes.cds.fasta  Rheum_tangaticum_H0.fna
cd ..

mkdir Fagopyrum_escelentum_H2; cd Fagopyrum_escelentum_H2
mv ../F_escelentum_H2.genes.cds.fasta  Fagopyrum_escelentum_H2.fna
cd ..

mkdir Fagopyrum_tataricum_H1; cd Fagopyrum_tataricum_H1
mv ../F_tataricum_H1.genes.cds.fasta  Fagopyrum_tataricum_H1.fna
cd ..


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
cat /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/proteins/*fasta >AllCleanProteins.fasta

# Make a giant fasta file containing all transcripts from every species
cat /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/peptide/*fasta >AllCleanTranscripts.fasta

#---
# index files
cdbfasta AllCleanTranscripts.fasta
cdbfasta AllCleanProteins.fasta

#------------------
# Extract orthologue protein and transcript names by looping through 
# the tree files and formatting to get gene names.

for f in /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees/*tree.txt; do sed 's/(//g' $f |sed 's/,/\n/g' |sed 's/:.*//g' |sort|uniq|sed 's/Proteins_/\t/1'|cut -f 2  |cdbyank AllCleanProteins.fasta.cidx >${f}_Proteins.fasta;done
for f in /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees/*; do sed 's/(//g' $f |sed 's/,/\n/g' |sed 's/:.*//g' |sort|uniq|sed 's/Proteins_/\t/1'|cut -f 2  |cdbyank AllCleanTranscripts.fasta.cidx >${f}_Transcripts.fasta;done

#------------------------
# Generate protein alignments with MAFFT

cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/

for f in *_Proteins.fasta; do
  echo "mafft \"$f\" > \"${f%.*}_alignment.fasta\""
done > generateMAFFTAlignment.sh

#------------------------
# Can we do protein alignments with PRANK?
# More phylogenetically correct but more variable too
# https://gensoft.pasteur.fr/docs/prank/170427/
# https://ariloytynoja.github.io/prank-msa/

module load StdEnv/2020  gcc/9.3.0  prank/170427

for f in *_Proteins.fasta; do
  echo "prank -protein -d=\"$f\" -o=\"${f%.*}_alignment.fasta\""
done > generatePRANKAlignment.sh

#---
# Advanced usage: 'prank [optional parameters] -d=sequence_file [optional parameters]'

 # input/output parameters:
  # -d=sequence_file (in FASTA format)
  # -t=tree_file [default: no tree, generate approximate NJ tree]
  # -o=output_file [default: 'output']
  # -f=output_format ['fasta' (default), 'phylipi', 'phylips', 'paml', 'nexus']
  # -showxml [output xml-files]
  # -showtree [output dnd-files]
  # -showanc [output ancestral sequences]
  # -showevents [output evolutioanry events]
  # -showall [output all of these]
  # -support [compute posterior support]
  # -njtree [estimate tree from input alignment (and realign)]
  # -treeonly [estimate tree only]
  # -quiet

 # model parameters:
  # +F or -F [force insertions to be always skipped]
  # -gaprate=# [gap opening rate; default: dna 0.025 / prot 0.005]
  # -gapext=# [gap extension probability; default: dna 0.75 / prot 0.5]
  # -codon [for coding DNA: use empirical codon model]
  # -DNA / -protein [no autodetection: use dna or protein model]
  # -termgap [penalise terminal gaps normally]
  # -nomissing [no missing data, use -F for terminal gaps ]

 # other parameters:
  # -keep [keep alignment "as is" (e.g. for ancestor inference)]
  # -iterate=# [rounds of re-alignment iteration]
  # -once [run only once; same as -iterate=1]
  # -prunetree [prune guide tree branches with no sequence data]
  # -prunedata [prune sequence data with no guide tree leaves]
  # -uselogs [slower but should work for a greater number of sequences]
  # -translate [translate to protein]
  # -mttranslate [translate to protein using mt table]

 # other:
  # -convert [no alignment, just convert to another format]
  # -version [check for updates]
  # -verbose [print progress etc. during runtime]

  # -help [show more options]


#---------------------------
# Run pal2nal – convert alignments into paml format
# The suggestion is to use pal2nal with “-nogap”, though codeml has a way to deal with these.
# However, if you do not have sequence in your pal2nal output, likely you had gaps
# along the entire sequence across all of the species analyzed. 

cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/

wget http://www.bork.embl.de/pal2nal/distribution/pal2nal.v14.tar.gz
tar -zxvf pal2nal.v14.tar.gz

cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees

for f in *alignment; \
do echo "/home/celphin/scratch/Oxyria_Positive_Selection_Test/pal2nal.v14/pal2nal.pl \
 "$f" "${f%_*}"_Transcripts.fasta \
 -output paml -nogap >"${f%_*}"_pal2nal" ; done >pal2nal.sh

#------------------------












