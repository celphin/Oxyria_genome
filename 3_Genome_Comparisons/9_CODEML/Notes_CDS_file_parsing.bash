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

#narval2
tmux new-session -s GeneSpace
tmux attach-session -t GeneSpace

# move peptide folder to proteins

cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes
mv peptide proteins
mkdir peptide 

#-----------------
cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/genomes

module load StdEnv/2020 r/4.2.2 glpk/5.0

R

# Install
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("jtlovell/GENESPACE")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("Biostrings", "rtracklayer"))

#----------------
# Load
library(GENESPACE)

#---------------------------
# Rosaceae

parsing_files_other <- function(SPP_Hap, GeneID, Wd){
  parsedPaths_other <- parse_annotations(
    rawGenomeRepo = "/home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/genomes", 
    genomeDirs = SPP_Hap,
    genomeIDs = SPP_Hap,
    gffString = "gff",
    faString = "fna",
    headerEntryIndex = 1, 
    gffIdColumn = GeneID,
    genespaceWd = "/home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes")
}

parse_ncbi_other <- function(SPP_Hap, gffID, gffStrip){
  parse_ncbi(
    rawGenomeRepo="/home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/genomes",
    genomeDirs = SPP_Hap,
    genomeIDs = SPP_Hap,
    gffString = "gff",
    faString = "fna",
    genespaceWd="/home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes",
    troubleShoot = FALSE
  )
}

# Dryas_octopetala
parsing_files_other("Dryas_octopetala", "ID")
# Dryas_octopetala: n unique sequences = 47946, n matched to gff = 39696
# Dryas_octopetala: n unique sequences = 39696, n matched to gff = 39696

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

parse_ncbi_other("Capsella_rubella")
# Capsella_rubella: n unique sequences = 26410, n matched to gff = 26410

#-----------------------
# Troubleshooting
# in bash
SPP_Hap="Fragaria_vesca"

grep "ID=gene-" ./${SPP_Hap}/${SPP_Hap}.gff|head -n 20 | tail
grep ">" ./${SPP_Hap}/${SPP_Hap}.fna | head -n 5

sed -i 's/db_xref\=GeneID\:/gene\=/' ./${SPP_Hap}/${SPP_Hap}.fna
sed -i 's/Name\=LOC/Name\=/' ./${SPP_Hap}/${SPP_Hap}.gff
sed -i 's/ID\=gene\-LOC/ID\=gene\-/' ./${SPP_Hap}/${SPP_Hap}.gff
sed -i 's/gene\=LOC/gene\=/' ./${SPP_Hap}/${SPP_Hap}.gff

# back in R
parse_ncbi_other <- function(SPP_Hap, gffID){
  parse_ncbi(
    rawGenomeRepo="/home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/genomes",
    genomeDirs = SPP_Hap,
    genomeIDs = SPP_Hap,
    gffString = "gff",
    faString = "fna",
    genespaceWd="/home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes",
    troubleShoot = FALSE,
    gffIdColumn = gffID
  )
}

parse_ncbi_other("Fragaria_vesca", "ID")
# Fragaria_vesca: n unique sequences = 22376, n matched to gff = 78
# Fragaria_vesca: n unique sequences = 22376, n matched to gff = 21397

####################################
# Polygonaceae

parsing_files_other <- function(SPP_Hap, GeneID){
  parsedPaths_other <- parse_annotations(
    rawGenomeRepo = "/home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/genomes", 
    genomeDirs = SPP_Hap,
    genomeIDs = SPP_Hap,
    gffString = "gff3",
    faString = "fna",
    headerEntryIndex = 1, 
    gffIdColumn = GeneID,
    genespaceWd = "/home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes")
}


#parsing_files_other("Fagopyrum_tataricum_H0", "ID")
#Fagopyrum_tataricum_H0: n unique sequences = 30903, n matched to gff = 28502

#parsing_files_other("Oxyria_digyna_H0", "ID")
#Oxyria_digyna_H0: n unique sequences = 42052, n matched to gff = 37871

#parsing_files_other("Fagopyrum_escelentum_H1", "ID")
#Fagopyrum_escelentum_H1: n unique sequences = 52725, n matched to gff = 52725

parsing_files_other("Fagopyrum_escelentum_H2", "ID")
# Fagopyrum_escelentum_H2: n unique sequences = 52149, n matched to gff = 52149

parsing_files_other("Fagopyrum_tataricum_H1", "ID")
# Fagopyrum_tataricum_H1: n unique sequences = 39560, n matched to gff = 39560

#parsing_files_other("Fagopyrum_tataricum_H2", "ID")
#Fagopyrum_tataricum_H2: n unique sequences = 37200, n matched to gff = 37200

parsing_files_other("Rheum_tangaticum_H0", "ID")
# Rheum_tangaticum_H0: n unique sequences = 31898, n matched to gff = 30938

# it worked!


#-------------
parsing_files_other("Rheum_nobile_H0", "ID")
# Rheum_nobile_H0: n unique sequences = 34698, n matched to gff = 34698

parsing_files_other("Polygunum_aviculare_H0", "ID")
# Polygunum_aviculare_H0: n unique sequences = 27714, n matched to gff = 26201

parsing_files_other("Oxyria_digyna_H1", "ID")
# Oxyria_digyna_H1: n unique sequences = 33799, n matched to gff = 33799

#################################
# Brassicaceae


parsing_files_other <- function(SPP_Hap, GeneID, Wd){
  parsedPaths_other <- parse_annotations(
    rawGenomeRepo = "/home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/genomes", 
    genomeDirs = SPP_Hap,
    genomeIDs = SPP_Hap,
    gffString = "gff",
    faString = "fna",
    headerEntryIndex = 1, 
    gffIdColumn = GeneID,
    genespaceWd = "/home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes")
}

parse_ncbi_other <- function(SPP_Hap, gffID, gffStrip){
  parse_ncbi(
    rawGenomeRepo="/home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/genomes",
    genomeDirs = SPP_Hap,
    genomeIDs = SPP_Hap,
    gffString = "gff",
    faString = "fna",
    genespaceWd="/home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes",
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
    rawGenomeRepo="/home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/genomes",
    genomeDirs = SPP_Hap,
    genomeIDs = SPP_Hap,
    gffString = "gff",
    faString = "fna",
    genespaceWd="/home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes",
    troubleShoot = FALSE,
    gffIdColumn = gffID
  )
}

#-----------------------------
# Troubleshooting
# does not work when locus tag not gene in gff matches ID=gene- in fasta
# works with edits (sed) below

SPP_Hap="Arabis_alpina"

sed -i 's/locus\_tag\=/gene\=/' ./${SPP_Hap}/${SPP_Hap}.fna

parse_ncbi_other("Arabis_alpina", "ID")
# Arabis_alpina: n unique sequences = 21609, n matched to gff = 0
# Arabis_alpina: n unique sequences = 21609, n matched to gff = 21609

#----------------------------
SPP_Hap="Arabidopsis_thaliana"
sed -i 's/gene\=/Genes\=/' ./${SPP_Hap}/${SPP_Hap}.fna
sed -i 's/locus\_tag\=/gene\=/' ./${SPP_Hap}/${SPP_Hap}.fna

parse_ncbi_other("Arabidopsis_thaliana", "ID")
# Arabidopsis_thaliana: n unique sequences = 27354, n matched to gff = 10250
# Arabidopsis_thaliana: n unique sequences = 27354, n matched to gff = 17070
#---------------------
SPP_Hap="Brassica_oleracea"
grep ">" ./${SPP_Hap}/${SPP_Hap}.fna | head 
# matching Dbxref=GeneID:
sed -i 's/gene\=/Genes\=/' ./${SPP_Hap}/${SPP_Hap}.fna
sed -i 's/db_xref\=GeneID\:/gene\=/' ./${SPP_Hap}/${SPP_Hap}.fna

parse_ncbi_other("Brassica_oleracea", "ID")
# Brassica_oleracea: n unique sequences = 44382, n matched to gff = 1968
# Brassica_oleracea: n unique sequences = 44382, n matched to gff = 42414

#--------------------------
SPP_Hap="Thlaspi_arvense"
grep ">" ./${SPP_Hap}/${SPP_Hap}.fna | head -n 5
sed -i 's/locus\_tag\=/gene\=/' ./${SPP_Hap}/${SPP_Hap}.fna

parse_ncbi_other("Thlaspi_arvense", "ID")
# Thlaspi_arvense: n unique sequences = 26392, n matched to gff = 0
# Thlaspi_arvense: n unique sequences = 26392, n matched to gff = 26392

#--------------------------
SPP_Hap="Cochlearia_groenlandica"
grep ">" ./${SPP_Hap}/${SPP_Hap}.fna | head -n 5
sed -i 's/g/ g/' ./${SPP_Hap}/${SPP_Hap}.fna

parse_ncbi_other("Cochlearia_groenlandica", "ID")
# Cochlearia_groenlandica: n unique sequences = 31127, n matched to gff = 31127


#######################################
# Missing
# Capsella_rubella _ now added
# Rubus_argutus - no protein file in the total genomes...


