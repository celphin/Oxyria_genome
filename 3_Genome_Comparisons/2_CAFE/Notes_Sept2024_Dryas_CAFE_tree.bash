##########################
# Running CAFE - gene enrichment analyses on all genomes on beluga 
# Sept 2024
############################

# Paper: https://academic.oup.com/bioinformatics/article/22/10/1269/237347
#https://github.com/hahnlab/CAFE5 

################################
# Installation

cd /home/celphin/scratch/Oxyria/CAFE

# from https://github.com/hahnlab/CAFE5/releases
wget https://github.com/hahnlab/CAFE5/releases/download/v5.1/CAFE5-5.1.0.tar.gz

tar -xvf CAFE5-5.1.0.tar.gz

module load StdEnv/2020 flexiblas/3.0.4 
module load StdEnv/2020  intel/2021.2.0 openblas/0.3.20

# To get the library search paths for the GNU c++ compiler (g++) one can type:
g++ -print-search-dirs | sed '/^lib/b 1;d;:1;s,/[^/.][^/]*/\.\./,/,;t 1;s,:[^=]*=,:;,;s,;,;  ,g' | tr \; \\012

#To get the Include search paths type:
g++ -E -Wp,-v -xc++ /dev/null

# Install
cd /home/celphin/scratch/Oxyria/CAFE/CAFE5

./configure
make

# test
cafe5 -h

#####################################
# Input files needed

# The count data file is a tab-delimited family "counts" file that contains a column for a 
# description of the gene family, the unique ID for each family, and a column for each taxon 
# that has count data for each family. This file is acquired by first peforming a clustering 
# analysis, often using software such as OrthoMCL, SwiftOrtho, FastOrtho, OrthAgogue, or OrthoFinder 
# and then parsing the output into a table like the one below 
# (Note: if a functional description is not desired, include this column anyway 
# with a place holder as below (null)).

# The file is tab-separated with a header line giving the order of the species. 
# Each line thereafter consists of a description, a family ID, and counts for each species in the tree.

# Desc	Family ID	human	chimp	orang	baboon	gibbon	macaque	marmoset rat	mouse	cat	horse	cow
# ATPase	ORTHOMCL1	 52	 55	 54	 57	 54	  56	  56	 53	 52	57	55	 54
# (null)	ORTHOMCL2	 76	 51	 41	 39	 45	  36	  37	 67	 79	37	41	 49
# HMG box	ORTHOMCL3	 50	 49	 48	 48	 46	  49	  48	 55	 52	51	47	 55
# (null)	ORTHOMCL4	 43	 43	 47	 53	 44	  47	  46	 59	 58	51	50	 55
# Dynamin	ORTHOMCL5	 43	 40	 43	 44	 31	  46	  33	 79	 70	43	49	 50
# ......
# ....
# ..
# DnaJ	ORTHOMCL10016	 45	 46	 50	 46	 46 	  47	  46	 48	 49	45	44	 48


#---------------------------
# parsing file from OrthoFinder

#https://www.biostars.org/p/9553290/

tmux new-session -s CAFE3
tmux attach-session -t CAFE3

# make directories
cd /home/celphin/scratch/Oxyria/CAFE
mkdir Dryas_genomes; cd Dryas_genomes
mkdir input_data
mkdir output
mkdir analysis
mkdir expanded
mkdir contracted

# copy over important files
cd /home/celphin/scratch/Oxyria/GeneSpace/Dryas_genomes/orthofinder/Results_Oct30
cp Species_Tree/SpeciesTree_rooted.txt /lustre04/scratch/celphin/Oxyria/CAFE/Dryas_genomes/input_data
cp Phylogenetic_Hierarchical_Orthogroups/N0.tsv /lustre04/scratch/celphin/Oxyria/CAFE/Dryas_genomes/input_data

#-------------------------------
# in R
cd /lustre04/scratch/celphin/Oxyria/CAFE/Dryas_genomes/input_data

module load StdEnv/2020 r/4.2.2 
R

#install.packages("ape")
# https://rdrr.io/cran/ape/man/chronos.html
# https://phylobotanist.blogspot.com/2016/12/how-to-set-multiple-calibration-points.html

library(ape)
library(data.table)

# make orthofinder's species tree ultrametric. It should be already binary and rooted as required by cafe.
setwd("/home/celphin/scratch/Oxyria/CAFE/Dryas_genomes/input_data")
tre <- read.tree('./SpeciesTree_rooted.txt')

stopifnot(is.binary(tre))
stopifnot(is.rooted(tre))

# Rosa diverged from Prunus 80mya, total tree 104 mya

# set calibrations
cal <- makeChronosCalib(tre, node="root", age.min=90, age.max=110) 

if(is.ultrametric(tre)) {
    utre <- tre
} else{
    utre <- chronos(tre, lambda = 1, model = "discrete", calibration = cal, control = chronos.control())
}

# Setting initial dates...
# Fitting in progress... get a first set of estimates
         # (Penalised) log-lik = -2.385786
# Optimising rates... frequencies... dates... -2.385786
# Optimising rates... frequencies... dates... -2.254234
# Optimising rates... frequencies... dates... -2.253594
# Optimising rates... frequencies... dates... -2.253559
# Optimising rates... frequencies... dates... -2.253558

# log-Lik = -2.253558
# PHIIC = 58.51

write.tree(utre, './SpeciesTree_rooted_ultra.txt')

#######################

#----------------------------------------
# Then prepare the table of counts. I use the N0 output since the original orthogroups are 
# deprecated in favour of the phylogenetic hierarchical orthogroups:
library(ape)
library(data.table)

hog <- fread('./N0.tsv')
hog[, OG := NULL]
hog[, `Gene Tree Parent Clade` := NULL]
hog <- melt(hog, id.vars='HOG', variable.name='species', value.name='pid')
hog <- hog[pid != '']
hog$n <- sapply(hog$pid, function(x) length(strsplit(x, ', ')[[1]]))

# Exclude HOGs with lots of genes in a one or more species. 
# See also cafe tutorial about filtering gene families
keep <- hog[, list(n_max=max(n)), HOG][n_max < 100]$HOG
hog <- hog[HOG %in% keep]

# remove Arabdiopsis thaliana because low gene counts due to poor parsing of annotation file
#hog <- hog[-which(hog$species=="Arabidopsis_thaliana"),]

# Exclude HOGs not present in at least 9 species
keep <- hog[, .N, HOG][N > 8]$HOG
hog <- hog[HOG %in% keep]

counts <- dcast(hog, HOG ~ species, value.var='n', fill=0)
counts[, Desc := 'n/a']
setcolorder(counts, 'Desc')
fwrite(counts, 'hog_gene_counts.tsv', sep='\t')


#-----------------------------
# Keep HOGs with lots of genes in another file. 

hog <- fread('./N0.tsv')
hog[, OG := NULL]
hog[, `Gene Tree Parent Clade` := NULL]
hog <- melt(hog, id.vars='HOG', variable.name='species', value.name='pid')
hog <- hog[pid != '']
hog$n <- sapply(hog$pid, function(x) length(strsplit(x, ', ')[[1]]))

# Exclude HOGs with lots of genes in a one or more species. 
# See also cafe tutorial about filtering gene families
keep <- hog[, list(n_max=max(n)), HOG][n_max < 100]$HOG
hog <- hog[HOG %in% keep]

# remove Arabdiopsis thaliana because low gene counts due to poor parsing of annotation file
hog <- hog[-which(hog$species=="Arabidopsis_thaliana"),]


# See also cafe tutorial about filtering gene families
keep <- hog[, list(n_max=max(n)), HOG][n_max > 100]$HOG
hog <- hog[HOG %in% keep]

# Exclude HOGs not present in all 20 species
keep <- hog[, .N, HOG][N > 19]$HOG
hog <- hog[HOG %in% keep]

counts <- dcast(hog, HOG ~ species, value.var='n', fill=0)
counts[, Desc := 'n/a']
setcolorder(counts, 'Desc')
fwrite(counts, 'large_hog_gene_counts.tsv', sep='\t')

#--------------------------
wc -l hog_gene_counts.tsv
11 183 hog_gene_counts.tsv

#---------------------------------
# The tree file should contain a binary, rooted, ultrametric, tree in Newick format. 
# Typically one obtains this tree using one of several molecular dating methods. 
# If you are unsure if your tree is binary, rooted, or ultrametric CAFE will 
# report this when you try to use it for an analysis. Alternatively, you can use the R package,
# Ape with its included functions: is.ultrametric, is.rooted, and is.binary.
# view here: http://etetoolkit.org/treeview/

more SpeciesTree_rooted_ultra.txt

((Rosa_rugosa:80.35408113,(Pyrus_bretschneideri:13.4018725,Malus_sylvestris:13.4018725)0.9
83537:66.95220863)0.980262:29.64591887,(Dry_drumm:32.0847978,((Dry_octo_H1:5.37717023,Dry_
octo_H0:5.37717023)0.420959:18.92554367,(Dry_int:18.02740367,(Dry_ajan:11.85396054,Dry_ala
sk:11.85396054)0.315012:6.17344313)0.281908:6.275310226)0.464241:7.782083899)0.980262:77.9
152022);


########################################
# Running CAFE
# https://github.com/hahnlab/CAFE5?tab=readme-ov-file

# Tutorial https://github.com/elsemikk/Willisornis_Genome_Assembly/blob/master/4.2_Gene_Family_Expansions.md

# Compare scenarios in which the whole phylogeny shares the same (global) lambda vs. scenarios in which different parts of the phylogeny share different (local) lambdas.

# Classify specific gene families as "fastly evolving” in at least two different ways (see documentation below).

# Infer gene family counts at all internal nodes (ancestral populations) of the phylogeny provided as input. By comparing different nodes in the phylogeny, the user will be able to tell along which branches gene families have contracted or expanded.

# Account for non-biological factors (e.g., genome sequencing and coverage differences, gene family clustering errors, etc.) leading to incorrect gene family counts in input files. This is done with an error model.

#-----------------------------------
# Make an error model

# all models with 5 spp stored in: /lustre04/scratch/celphin/Oxyria/CAFE/Dryas_genomes/5spp_output/
# run main model with all spp here: /lustre04/scratch/celphin/Oxyria/CAFE/Dryas_genomes/output/

tmux new-session -s CAFE
tmux attach-session -t CAFE

salloc -c30 --time 3:00:00 --mem 191000M --account def-cronk

module load StdEnv/2020 intel/2020.1.217 cafe5/5.1.0

cd /lustre04/scratch/celphin/Oxyria/CAFE
cafe5 -i /lustre04/scratch/celphin/Oxyria/CAFE/Dryas_genomes/input_data/hog_gene_counts.tsv \
-t /lustre04/scratch/celphin/Oxyria/CAFE/Dryas_genomes/input_data/SpeciesTree_rooted_ultra.txt \
-o /lustre04/scratch/celphin/Oxyria/CAFE/Dryas_genomes/output/error_model -e --cores 30

#-------------------------------
# Completed 63 iterations
# Time: 0H 0M 17S
# Best matches are: 0.0015245772530848,-0.0024665595472228
# Final -lnL: 106251.9058827

# 125 values were attempted (0% rejected)

# Inferring processes for Base model
# Score (-lnL):  106251.9058827
# Maximum possible lambda for this topology: 0.012444918614428
# Computing pvalues...
# done!

#----------------------------------
cd /lustre04/scratch/celphin/Oxyria/CAFE/Dryas_genomes/output
more error_model/Base_error_model.txt #view results


##########################################
# Explore data
cd /lustre04/scratch/celphin/Oxyria/CAFE/Dryas_genomes/output/error_model/

more Base_asr.tre 
#gives reconstructed numbers of each gene family at each node

#--------------------------------
more Base_branch_probabilities.tab 
#gives calculated likelihood of the gene family size at each node for each gene family

#------------------------------
more Base_clade_results.txt

# #Taxon_ID       Increase        Decrease
# Rosa_rugosa<1>  589     357
# Pyrus_bretschneideri<2> 1034    194
# Malus_sylvestris<3>     1001    201
# Dry_drumm<4>    67      128
# Dry_octo_H1<5>  212     129
# Dry_octo_H0<6>  232     8
# Dry_int<7>      68      153
# Dry_ajan<8>     76      123
# Dry_alask<9>    88      111
# <11>    27      24
# <12>    5541    110
# <13>    169     962
# <14>    32      4
# <15>    205     11
# <16>    1       107
# <17>    9       48



#--------------------------------------------------
cd /lustre04/scratch/celphin/Oxyria/CAFE/Dryas_genomes/output/error_model/

grep "Lambda" ./*/*_results.txt


grep "Epsilon" ./*/*_results.txt


#------------------------------------
# run many times to check lambdas 

printf '%s\n' {1..50} | parallel cafe5 -i /lustre04/scratch/celphin/Oxyria/CAFE/Dryas_genomes/input_data/hog_gene_counts.tsv \
-t /lustre04/scratch/celphin/Oxyria/CAFE/Dryas_genomes/input_data/SpeciesTree_rooted_ultra.txt -e \
-o /lustre04/scratch/celphin/Oxyria/CAFE/Dryas_genomes/output/error_model_{1} >> \
error_model_100runs.log 2>> error_model_100runs.errorlog

#-----------------------------------------
#Summarise the results of all runs
cd /lustre04/scratch/celphin/Oxyria/CAFE/Dryas_genomes/output/

grep "Lambda" ./error_model_*/Base_results.txt | sed 's/^.* //g' > lambdas.txt


grep "Epsilon" ./error_model_*/Base_results.txt | sed 's/^.* //g' > epsilons.txt

#-------------------------q---
module load StdEnv/2020 r/4.2.2 

R

epsilons <- read.table("epsilons.txt")
lambdas <- read.table("lambdas.txt")

mean(epsilons[,1])
# 0.001346833
sd(epsilons[,1])
# 0.007842287

mean(lambdas[,1])
#  0.004702567
sd(lambdas[,1])
# 0.0001320786

q()
n

#####################
# Plotting results
# https://github.com/moshi4/CafePlotter

cd /lustre04/scratch/celphin/Oxyria/CAFE/

module load StdEnv/2023 python/3.11.5 scipy-stack/2023b
pip install cafeplotter

cafeplotter -i /lustre04/scratch/celphin/Oxyria/CAFE/Dryas_genomes/output/error_model \
-o /lustre04/scratch/celphin/Oxyria/CAFE/Dryas_genomes/Cafe_plots_error_model_pdf  --format pdf


######################################
# Now we can investigate which orthogroups are expanded/contracted
# See notes file specifc to this!