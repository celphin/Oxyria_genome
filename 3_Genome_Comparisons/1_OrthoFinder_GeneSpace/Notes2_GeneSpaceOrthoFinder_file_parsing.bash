##########################
##########################
# Running GeneSpace file parsing and formatting 
# June-August 2024
############################

# Paper: https://elifesciences.org/articles/78526 
# https://github.com/jtlovell/GENESPACE

####################################
# run parse annotations
# https://rdrr.io/github/jtlovell/GENESPACE/man/parse_annotations.html

tmux new-session -s GeneSpace
tmux attach-session -t GeneSpace

module load StdEnv/2020 r/4.2.2 glpk/5.0

R

library(GENESPACE)

#---------------------------
# Rosaceae

parsing_files_other <- function(SPP_Hap, GeneID, Wd){
  parsedPaths_other <- parse_annotations(
    rawGenomeRepo = "/home/celphin/scratch/Oxyria/GeneSpace/Rosaceae_genomes/genomes", 
    genomeDirs = SPP_Hap,
    genomeIDs = SPP_Hap,
    gffString = "gff",
    faString = "fa",
    headerEntryIndex = 1, 
    gffIdColumn = GeneID,
    genespaceWd = "/home/celphin/scratch/Oxyria/GeneSpace/Rosaceae_genomes")
}

parse_ncbi_other <- function(SPP_Hap, gffID, gffStrip){
  parse_ncbi(
    rawGenomeRepo="/home/celphin/scratch/Oxyria/GeneSpace/Rosaceae_genomes/genomes",
    genomeDirs = SPP_Hap,
    genomeIDs = SPP_Hap,
    gffString = "gff",
    faString = "fa",
    genespaceWd="/home/celphin/scratch/Oxyria/GeneSpace/Rosaceae_genomes",
    troubleShoot = FALSE
  )
}

# Dryas_octopetala
parsing_files_other("Dryas_octopetala", "ID")
# Dryas_octopetala: n unique sequences = 47946, n matched to gff = 39696

#Pyrus_bretschneideri
parse_ncbi_other("Pyrus_bretschneideri")
# Pyrus_bretschneideri: n unique sequences = 35293, n matched to gff = 35293

parse_ncbi_other("Prunus_persica")
# Prunus_persica: n unique sequences = 23128, n matched to gff = 23128

parse_ncbi_other("Malus_sylvestris")
# n unique sequences = 37467, n matched to gff = 37467

parse_ncbi_other("Argentina_anserina")
# Argentina_anserina: n unique sequences = 19620, n matched to gff = 19620

parse_ncbi_other("Rosa_rugosa")
# Rosa_rugosa: n unique sequences = 29146, n matched to gff = 29146

#-----------------------
# Troubleshooting
# in bash
SPP_Hap="Fragaria_vesca"

grep "ID=gene-" ./${SPP_Hap}/${SPP_Hap}.gff|head -n 20 | tail
grep ">" ./${SPP_Hap}/${SPP_Hap}.fa | head -n 5

sed -i 's/db_xref\=GeneID\:/gene\=/' ./${SPP_Hap}/${SPP_Hap}.fa
sed -i 's/Name\=LOC/Name\=/' ./${SPP_Hap}/${SPP_Hap}.gff
sed -i 's/ID\=gene\-LOC/ID\=gene\-/' ./${SPP_Hap}/${SPP_Hap}.gff
sed -i 's/gene\=LOC/gene\=/' ./${SPP_Hap}/${SPP_Hap}.gff

# back in R
parse_ncbi_other <- function(SPP_Hap, gffID){
  parse_ncbi(
    rawGenomeRepo="/home/celphin/scratch/Oxyria/GeneSpace/Rosaceae_genomes/genomes",
    genomeDirs = SPP_Hap,
    genomeIDs = SPP_Hap,
    gffString = "gff",
    faString = "fa",
    genespaceWd="/home/celphin/scratch/Oxyria/GeneSpace/Rosaceae_genomes",
    troubleShoot = FALSE,
    gffIdColumn = gffID
  )
}

parse_ncbi_other("Fragaria_vesca", "ID")
# Fragaria_vesca: n unique sequences = 22376, n matched to gff = 1057
# Fragaria_vesca: n unique sequences = 22376, n matched to gff = 21397

####################################
# Polygonaceae

parsing_files_other <- function(SPP_Hap, GeneID){
  parsedPaths_other <- parse_annotations(
    rawGenomeRepo = "/home/celphin/scratch/Oxyria/GeneSpace/Polygonaceae_genomes/genomes", 
    genomeDirs = SPP_Hap,
    genomeIDs = SPP_Hap,
    gffString = "gff3",
    faString = "fa",
    headerEntryIndex = 1, 
    gffIdColumn = GeneID,
    genespaceWd = "/home/celphin/scratch/Oxyria/GeneSpace/Polygonaceae_genomes")
}


#parsing_files_other("Fagopyrum_tataricum_H0", "ID")
#Fagopyrum_tataricum_H0: n unique sequences = 30903, n matched to gff = 28502

#parsing_files_other("Oxyria_digyna_H0", "ID")
#Oxyria_digyna_H0: n unique sequences = 42052, n matched to gff = 37871

#parsing_files_other("Fagopyrum_escelentum_H1", "ID")
#Fagopyrum_escelentum_H1: n unique sequences = 52725, n matched to gff = 52725

parsing_files_other("Fagopyrum_escelentum_H2", "ID")
# Fagopyrum_escelentum_H2: n unique sequences = 52149, n matched to gff = 52149

parsing_files_other("Oxyria_digyna_H1", "ID")
# Oxyria_digyna_H1: n unique sequences = 37263, n matched to gff = 33799

parsing_files_other("Fagopyrum_tataricum_H1", "ID")
# Fagopyrum_tataricum_H1: n unique sequences = 39560, n matched to gff = 39560

#parsing_files_other("Fagopyrum_tataricum_H2", "ID")
#Fagopyrum_tataricum_H2: n unique sequences = 37200, n matched to gff = 37200

parsing_files_other("Polygunum_aviculare_H0", "ID")
# Polygunum_aviculare_H0: n unique sequences = 27714, n matched to gff = 26201

parsing_files_other("Rheum_nobile_H0", "ID")
# Rheum_nobile_H0: n unique sequences = 34698, n matched to gff = 34698

parsing_files_other("Rheum_tangaticum_H0", "ID")
# Rheum_tangaticum_H0: n unique sequences = 31898, n matched to gff = 30938

# it worked!

#################################
# Brassicaceae


parsing_files_other <- function(SPP_Hap, GeneID, Wd){
  parsedPaths_other <- parse_annotations(
    rawGenomeRepo = "/home/celphin/scratch/Oxyria/GeneSpace/Brassicaceae_genomes/genomes", 
    genomeDirs = SPP_Hap,
    genomeIDs = SPP_Hap,
    gffString = "gff",
    faString = "fa",
    headerEntryIndex = 1, 
    gffIdColumn = GeneID,
    genespaceWd = "/home/celphin/scratch/Oxyria/GeneSpace/Brassicaceae_genomes")
}

parse_ncbi_other <- function(SPP_Hap, gffID, gffStrip){
  parse_ncbi(
    rawGenomeRepo="/home/celphin/scratch/Oxyria/GeneSpace/Brassicaceae_genomes/genomes",
    genomeDirs = SPP_Hap,
    genomeIDs = SPP_Hap,
    gffString = "gff",
    faString = "fa",
    genespaceWd="/home/celphin/scratch/Oxyria/GeneSpace/Brassicaceae_genomes",
    troubleShoot = FALSE
  )
}

# Draba_nivalis
parsing_files_other("Draba_nivalis", "ID")
# Draba_nivalis: n unique sequences = 33557, n matched to gff = 33557

parse_ncbi_other("Arabidopsis_lyrata")
# Arabidopsis_lyrata: n unique sequences = 29812, n matched to gff = 29811 

#----------------
# Did not work troubleshoot FALSE, try TRUE - no better
# https://rdrr.io/github/jtlovell/GENESPACE/man/parse_annotations.html
# https://rdrr.io/github/jtlovell/GENESPACE/src/R/parse_annotations.R

  # gff <- data.table(rtracklayer::readGFF(path2gff, tags = gffIdColumn))
  # setnames(gff, gffIdColumn, "id")
  # gff[,id := gsub(gffStripText, "", id)]
  # if(troubleShoot){
    # cat("\n### first 6 gff lines after parsing ... \n")
    # print(head(gff))
  # }

parse_ncbi_other <- function(SPP_Hap, gffID){
  parse_ncbi(
    rawGenomeRepo="/home/celphin/scratch/Oxyria/GeneSpace/Brassicaceae_genomes/genomes",
    genomeDirs = SPP_Hap,
    genomeIDs = SPP_Hap,
    gffString = "gff",
    faString = "fa",
    genespaceWd="/home/celphin/scratch/Oxyria/GeneSpace/Brassicaceae_genomes",
    troubleShoot = FALSE,
    gffIdColumn = gffID
  )
}

#-----------------------------
# Troubleshooting
# does not work when locus tag not gene in gff matches ID=gene- in fasta
# works with edits (sed) below

SPP_Hap="Arabis_alpina"

sed -i 's/locus\_tag\=/gene\=/' ./${SPP_Hap}/${SPP_Hap}.gff
sed -i 's/locus\_tag\=/gene\=/' ./${SPP_Hap}/${SPP_Hap}.fa

parse_ncbi_other("Arabis_alpina", "ID")
# Arabis_alpina: n unique sequences = 21609, n matched to gff = 0
# Arabis_alpina: n unique sequences = 21609, n matched to gff = 21609
#----------------------------
SPP_Hap="Arabidopsis_thaliana"

sed -i 's/gene\=/Genes\=/' ./${SPP_Hap}/${SPP_Hap}.gff
sed -i 's/gene\=/Genes\=/' ./${SPP_Hap}/${SPP_Hap}.fa

sed -i 's/locus\_tag\=/gene\=/' ./${SPP_Hap}/${SPP_Hap}.gff
sed -i 's/locus\_tag\=/gene\=/' ./${SPP_Hap}/${SPP_Hap}.fa

parse_ncbi_other("Arabidopsis_thaliana", "ID")
# Arabidopsis_thaliana: n unique sequences = 27354, n matched to gff = 10250
# Arabidopsis_thaliana: n unique sequences = 27354, n matched to gff = 17070
#---------------------
SPP_Hap="Brassica_oleracea"

grep "ID=gene-" ./${SPP_Hap}/${SPP_Hap}.gff|head -n 20 | tail
grep ">" ./${SPP_Hap}/${SPP_Hap}.fa | head 
# matching Dbxref=GeneID:
sed -i 's/gene\=/Genes\=/' ./${SPP_Hap}/${SPP_Hap}.gff
sed -i 's/gene\=/Genes\=/' ./${SPP_Hap}/${SPP_Hap}.fa

sed -i 's/Dbxref\=GeneID\:/gene\=/' ./${SPP_Hap}/${SPP_Hap}.gff
sed -i 's/db_xref\=GeneID\:/gene\=/' ./${SPP_Hap}/${SPP_Hap}.fa

sed -i 's/gene\=LOC/gene\=/' ./${SPP_Hap}/${SPP_Hap}.gff


parse_ncbi_other("Brassica_oleracea", "ID")
# Brassica_oleracea: n unique sequences = 44382, n matched to gff = 1968
# Brassica_oleracea: n unique sequences = 44382, n matched to gff = 42491

#--------------------------
SPP_Hap="Thlaspi_arvense"

grep "ID=gene-" ./${SPP_Hap}/${SPP_Hap}.gff|head -n 20 | tail
grep ">" ./${SPP_Hap}/${SPP_Hap}.fa | head -n 5

sed -i 's/locus\_tag\=/gene\=/' ./${SPP_Hap}/${SPP_Hap}.gff
sed -i 's/locus\_tag\=/gene\=/' ./${SPP_Hap}/${SPP_Hap}.fa

parse_ncbi_other("Thlaspi_arvense", "ID")
# Thlaspi_arvense: n unique sequences = 26392, n matched to gff = 0
# Thlaspi_arvense: n unique sequences = 26392, n matched to gff = 26392

#--------------------------
SPP_Hap="Cochlearia_groenlandica"

grep "ID=gene-" ./${SPP_Hap}/${SPP_Hap}.gtf|head -n 20 | tail
grep ">" ./${SPP_Hap}/${SPP_Hap}.fa | head -n 5

grep "ID=g" ./${SPP_Hap}/${SPP_Hap}.gff3

#-----------------------------
# Convert gtf to gff
# https://agat.readthedocs.io/en/latest/tools/agat_convert_sp_gxf2gxf.html
# AGAT

cd /home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/genomes/Cochlearia_groenlandica

# Load modules on graham or cedar (or use instructions for conda at https://github.com/NBISweden/AGAT if using your own computer)
module load StdEnv/2023
module load apptainer/1.2.4

# Install `AGAT` via `Apptainer` into a temp directory
apptainer pull docker://quay.io/biocontainers/agat:1.4.0--pl5321hdfd78af_0
apptainer run agat_1.4.0--pl5321hdfd78af_0.sif

# Kick the tires, was it successfully installed?
agat_convert_sp_gxf2gxf.pl --help | less # It should show tool specific information. Use 'q' to return the prompt
agat_<tab>  # tab-completion should result in a list of all the possible scripts; hit space to advance to next page, or 'q' to return the prompt
agat_convert_sp_gxf2gxf.pl --help

agat_convert_sp_gxf2gxf.pl -gtf Cochlearia_groenlandica.gtf -o Cochlearia_groenlandica.gff3

# => Number of line in file: 602809
# => Number of comment lines: 0
# => Fasta included: No
# => Number of features lines: 602809
# => Number of feature type (3rd column): 7
        # * Level1: 1 => gene
        # * level2: 1 => transcript
        # * level3: 5 => CDS stop_codon exon start_codon intron
        # * unknown: 0 =>
# There is a problem we found several formats in this file: 1,2
# Let's see what we can do...


# WARNING gff3 reader: Hmmm, be aware that your feature doesn't contain any Parent and locus tag. No worries, we will handle it by considering it as strictly sequential. If you disagree, please provide an ID or a comon tag by locus. @ the feature is:
# h2tg000015l     AUGUSTUS        transcript      36789   37112   0.96    +       .       ID "g10.t1"

#-------------------------------
##gff-version 3
h2tg000001l     AUGUSTUS        gene    1       900     0.53    -       .       ID=g6166
h2tg000001l     AUGUSTUS        transcript      1       900     0.53    -       .       ID=g6166.t1;
Parent=g6166
h2tg000001l     AUGUSTUS        exon    1       165     .       -       .       ID=agat-exon-34710;P
arent=g6166.t1;gene_id=g6166;transcript_id=g6166.t1
h2tg000001l     AUGUSTUS        exon    342     418     .       -       .       ID=agat-exon-34711;P
arent=g6166.t1;gene_id=g6166;transcript_id=g6166.t1
h2tg000001l     AUGUSTUS        exon    500     570     .       -       .       ID=agat-exon-34712;P
arent=g6166.t1;gene_id=g6166;transcript_id=g6166.t1
h2tg000001l     AUGUSTUS        exon    662     748     .       -       .       ID=agat-exon-34713;P
arent=g6166.t1;gene_id=g6166;transcript_id=g6166.t1

#-------------------------
# remove all but transcripts

grep "ID=g" Cochlearia_groenlandica2.gff3 | grep "transcript" > Cochlearia_groenlandica.gff3
sed -i 's/g/ g/' Cochlearia_groenlandica.fa

#---------------------------
parse_ncbi_other("Cochlearia_groenlandica", "ID")
# Cochlearia_groenlandica: n unique sequences = 31127, n matched to gff = 31127

#-----------------------
# Check BUSCO for protein file?

tmux new-session -s BUSCO
tmux attach-session -t BUSCO

cd /home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/genomes/Cochlearia_groenlandica

salloc -c10 --time 2:55:00 --mem 120000m --account def-rieseber

module load StdEnv/2020
module load gcc/9.3.0 python augustus hmmer blast+ metaeuk prodigal r bbmap
source ~/busco_env/bin/activate

busco --offline --in Cochlearia_groenlandica.fa \
--out  BUSCO_Cochlearia_groenlandica_pep  --lineage_dataset eudicots_odb10 --mode protein --cpu 10 \
--download_path ~/BUSCO_downloads/

# remove spaces
sed -i 's/ //g' Cochlearia_groenlandica.fa

# C:91.0%[S:79.7%,D:11.3%],F:0.6%,M:8.4%,n:2326


#######################################
# try to parse Dryas spp. genome annotations

cd /home/celphin/scratch/Oxyria/GeneSpace/Dryas_genomes/genomes/

cp /home/celphin/scratch/Dryas/Dryas_genomes/Liftoff/Dry-*_proteins.fa .
cp /home/celphin/scratch/Dryas/Dryas_genomes/Liftoff/Dry-*_liftoffpolishslurm.gff3 .
cp /home/celphin/scratch/Dryas/Dryas_genomes/Liftoff/FINAL* . 
mv FINAL_DoctH0.AED_0.6.sorted.gff3 Dry-octo-H0.gff3


# for all Dryas files
 rename _liftoffpolishslurm.gff3 .gff3 *
 rename _proteins.fa .fa *
 rename - _ *

mkdir Dry_alask
mv Dry_alask* Dry_alask/

mkdir Dry_int
mv Dry_int* Dry_int/

mkdir Dry_drumm
mv Dry_drumm* Dry_drumm/

mkdir Dry_octo_H0
mv Dry_octo_H0* Dry_octo_H0/

mkdir Dry_octo_H1
mv Dry_octo_H1* Dry_octo_H1/

mkdir Dry_ajan
mv Dry_ajan* Dry_ajan/


#_--------------------------
tmux new-session -s GeneSpace2
tmux attach-session -t GeneSpace2

module load StdEnv/2020 r/4.2.2 glpk/5.0

R

library(GENESPACE)

#---------------------------
# Dryas

parsing_files_other <- function(SPP_Hap, GeneID, Wd){
  parsedPaths_other <- parse_annotations(
    rawGenomeRepo = "/home/celphin/scratch/Oxyria/GeneSpace/Dryas_genomes/genomes", 
    genomeDirs = SPP_Hap,
    genomeIDs = SPP_Hap,
    gffString = "gff3",
    faString = "fa",
    headerEntryIndex = 1, 
    gffIdColumn = GeneID,
    genespaceWd = "/home/celphin/scratch/Oxyria/GeneSpace/Dryas_genomes")
}

# Dryas_octopetala
parsing_files_other("Dry_octo_H0", "ID")
# Dry_octo_H0: n unique sequences = 39696, n matched to gff = 39696

parsing_files_other("Dry_octo_H1", "ID")
# Dry_octo_H1: n unique sequences = 37395, n matched to gff = 37395 (dropped 94791 ./- chars)

parsing_files_other("Dry_ajan", "ID")
# Dry_ajan: n unique sequences = 34179, n matched to gff = 34179 (dropped 112295 ./- chars)

parsing_files_other("Dry_int", "ID")
# Dry_int: n unique sequences = 34097, n matched to gff = 34097 (dropped 112851 ./_ chars)

parsing_files_other("Dry_drumm", "ID")
# Dry_drumm: n unique sequences = 34456, n matched to gff = 34456 (dropped 152512 ./_ chars)

parsing_files_other("Dry_alask", "ID")
# Dry_alask: n unique sequences = 34226, n matched to gff = 34226 (dropped 114241 ./- chars)

#---------------------------
# copy over annotation files for Rosa, pyrus and Malus_sylvestris

cd /home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/peptide
cp Malus_sylvestris.fa /home/celphin/scratch/Oxyria/GeneSpace/Dryas_genomes/peptide/
cp Pyrus_bretschneideri.fa /home/celphin/scratch/Oxyria/GeneSpace/Dryas_genomes/peptide/
cp Rosa_rugosa.fa /home/celphin/scratch/Oxyria/GeneSpace/Dryas_genomes/peptide/

cd /home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/bed
cp Malus_sylvestris.bed /home/celphin/scratch/Oxyria/GeneSpace/Dryas_genomes/bed/
cp Pyrus_bretschneideri.bed /home/celphin/scratch/Oxyria/GeneSpace/Dryas_genomes/bed/
cp Rosa_rugosa.bed /home/celphin/scratch/Oxyria/GeneSpace/Dryas_genomes/bed/

cd /home/celphin/scratch/Oxyria/GeneSpace/Dryas_genomes/bed/

#######################################
# try to parse Oxyria spp. genome annotations

cd /home/celphin/scratch/Oxyria/GeneSpace/Oxyria_genomes/genomes

rename _pep.fasta .fa *

mkdir DToL_h1
mv DToL_h1* DToL_h1/

mkdir DToL_h2
mv DToL_h2* DToL_h2/

mkdir Oxy_Elles_Hap1
mv Oxy_Elles_Hap1* Oxy_Elles_Hap1/

mkdir Oxy_Elles_Hap2
mv Oxy_Elles_Hap2* Oxy_Elles_Hap2/

mkdir Oxy_Sval_h1
mv Oxy_Sval_h1* Oxy_Sval_h1/

mkdir Oxy_Sval_h2
mv Oxy_Sval_h2* Oxy_Sval_h2/

#---------------------------
tmux new-session -s GeneSpace
tmux attach-session -t GeneSpace

module load StdEnv/2020 r/4.2.2 glpk/5.0

R

library(GENESPACE)

#---------------------------
# Oxyria

parsing_files_other <- function(SPP_Hap, GeneID, Wd){
  parsedPaths_other <- parse_annotations(
    rawGenomeRepo = "/home/celphin/scratch/Oxyria/GeneSpace/Oxyria_genomes/genomes", 
    genomeDirs = SPP_Hap,
    genomeIDs = SPP_Hap,
    gffString = "gff3",
    faString = "fa",
    headerEntryIndex = 1, 
    gffIdColumn = GeneID,
    genespaceWd = "/home/celphin/scratch/Oxyria/GeneSpace/Oxyria_genomes")
}

parsing_files_other("DToL_h1", "ID")
# DToL_h1: n unique sequences = 31523, n matched to gff = 31523 (dropped 125791 ./- chars)
parsing_files_other("DToL_h2", "ID")
# DToL_h2: n unique sequences = 35477, n matched to gff = 35477 (dropped 127566 ./- chars)

parsing_files_other("Oxy_Elles_Hap1", "ID")
# Oxy_Elles_Hap1: n unique sequences = 38198, n matched to gff = 38198 (dropped 43824 ./- chars)

parsing_files_other("Oxy_Elles_Hap2", "ID")
# Oxy_Elles_Hap2: n unique sequences = 35196, n matched to gff = 35196 (dropped 47302 ./- chars)

parsing_files_other("Oxy_Sval_h1", "ID")
# Oxy_Sval_h1: n unique sequences = 33402, n matched to gff = 33402 (dropped 135951 ./- chars)

parsing_files_other("Oxy_Sval_h2", "ID")
# Oxy_Sval_h2: n unique sequences = 31714, n matched to gff = 31714 (dropped 129380 ./- chars)

#---------------------------
# copy over annotation files for Rosa, pyrus and Malus_sylvestris

cd /home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/peptide
cp Rheum_nobile_H0.fa /home/celphin/scratch/Oxyria/GeneSpace/Oxyria_genomes/peptide/
cp Polygunum_aviculare_H0.fa /home/celphin/scratch/Oxyria/GeneSpace/Oxyria_genomes/peptide/
cp Fagopyrum_tataricum_H1.fa /home/celphin/scratch/Oxyria/GeneSpace/Oxyria_genomes/peptide/

cd /home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/bed
cp Rheum_nobile_H0.bed /home/celphin/scratch/Oxyria/GeneSpace/Oxyria_genomes/bed/
cp Polygunum_aviculare_H0.bed /home/celphin/scratch/Oxyria/GeneSpace/Oxyria_genomes/bed/
cp Fagopyrum_tataricum_H1.bed /home/celphin/scratch/Oxyria/GeneSpace/Oxyria_genomes/bed/

cd /home/celphin/scratch/Oxyria/GeneSpace/Oxyria_genomes/bed/