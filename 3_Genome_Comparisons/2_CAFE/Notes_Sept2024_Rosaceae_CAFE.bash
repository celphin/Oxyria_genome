##########################
# Running CAFE - gene enrichment analyses on Rosaceae genomes on beluga 
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

tmux new-session -s CAFE2
tmux attach-session -t CAFE2

# make directories
cd /home/celphin/scratch/Oxyria/CAFE
mkdir Rosaceae_data; cd Rosaceae_data
mkdir input_data
mkdir output
mkdir analysis
mkdir expanded
mkdir contracted

# copy over important files
cd /home/celphin/scratch/Oxyria/GeneSpace/Rosaceae_genomes/orthofinder/Results_Aug16
cp Species_Tree/SpeciesTree_rooted.txt /lustre04/scratch/celphin/Oxyria/CAFE/Rosaceae_data/input_data
cp Phylogenetic_Hierarchical_Orthogroups/N0.tsv /lustre04/scratch/celphin/Oxyria/CAFE/Rosaceae_data/input_data

#----------------------------------
# in R
cd /lustre04/scratch/celphin/Oxyria/CAFE/Rosaceae_data/input_data

module load StdEnv/2020 r/4.2.2 
R

#install.packages("ape")
# https://rdrr.io/cran/ape/man/chronos.html
# https://phylobotanist.blogspot.com/2016/12/how-to-set-multiple-calibration-points.html

library(ape)
library(data.table)

# make orthofinder's species tree ultrametric. It should be already binary and rooted as required by cafe.
setwd("/home/celphin/scratch/Oxyria/CAFE/Rosaceae_data/input_data")
tre <- read.tree('./SpeciesTree_rooted.txt')
stopifnot(is.binary(tre))
stopifnot(is.rooted(tre))

# Rosa divergence 74 Mya
# https://www.researchgate.net/figure/The-estimation-of-divergence-time-and-expansion-or-contraction-of-gene-families-in-Prunus_fig2_357148509
# https://bmcplantbiol.biomedcentral.com/articles/10.1186/1471-2229-8-67
# Prunus (about 13my old) and diverged from fragaria 29mya
# https://bmcplantbiol.biomedcentral.com/articles/10.1186/1471-2229-8-67

# set calibrations
cal <- makeChronosCalib(tre, node="root", age.min=40, age.max=70) 

if(is.ultrametric(tre)) {
    utre <- tre
} else{
    utre <- chronos(tre, lambda = 1, model = "clock", calibration = cal, control = chronos.control())
}

# Setting initial dates...
# Fitting in progress... get a first set of estimates
         # (Penalised) log-lik = -2.454595
# Optimising rates... dates... -2.454595

# log-Lik = -2.454595
# PHIIC = 18.91

write.tree(utre, './SpeciesTree_rooted_ultra.txt')

#-----------------------------
# Then prepare the table of counts. I use the N0 output since the original orthogroups are 
# deprecated in favour of the phylogenetic hierarchical orthogroups:

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

# Exclude HOGs present in in all 7 species
keep <- hog[, .N, HOG][N > 6]$HOG
hog <- hog[HOG %in% keep]

counts <- dcast(hog, HOG ~ species, value.var='n', fill=0)
counts[, Desc := 'n/a']
setcolorder(counts, 'Desc')
fwrite(counts, 'hog_gene_counts.tsv', sep='\t')

#-----------------------------
# Keep HOGs with lots of genes in another file. 
# See also cafe tutorial about filtering gene families
keep <- hog[, list(n_max=max(n)), HOG][n_max > 100]$HOG
hog <- hog[HOG %in% keep]

# Exclude HOGs present in all 7 species
keep <- hog[, .N, HOG][N > 6]$HOG
hog <- hog[HOG %in% keep]

counts <- dcast(hog, HOG ~ species, value.var='n', fill=0)
counts[, Desc := 'n/a']
setcolorder(counts, 'Desc')
fwrite(counts, 'large_hog_gene_counts.tsv', sep='\t')
#----------------------------
wc -l hog_gene_counts.tsv
12917
wc -l large_hog_gene_counts.tsv
1

#---------------------------------
# The tree file should contain a binary, rooted, ultrametric, tree in Newick format. 
# Typically one obtains this tree using one of several molecular dating methods. 
# If you are unsure if your tree is binary, rooted, or ultrametric CAFE will 
# report this when you try to use it for an analysis. Alternatively, you can use the R package,
# Ape with its included functions: is.ultrametric, is.rooted, and is.binary.
# download tree and view in http://etetoolkit.org/treeview/

(Dryas_octopetala:62.04146729,((Rosa_rugosa:22.40944509,(Argentina_anserina:18.71836684,Fragaria_ve
sca:18.71836684)0.615002:3.691078256)0.97422:23.19883189,((Pyrus_bretschneideri:6.774856001,Malus_s
ylvestris:6.774856001)0.9861:26.86055074,Prunus_persica:33.63540674)0.817829:11.97287024)1:16.43319
03);

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

tmux new-session -s CAFE2
tmux attach-session -t CAFE2

salloc -c30 --time 3:00:00 --mem 191000M --account def-cronk

module load StdEnv/2020 intel/2020.1.217 cafe5/5.1.0

cd /lustre04/scratch/celphin/Oxyria/CAFE
cafe5 -i /lustre04/scratch/celphin/Oxyria/CAFE/Rosaceae_data/input_data/hog_gene_counts.tsv \
-t /lustre04/scratch/celphin/Oxyria/CAFE/Rosaceae_data/input_data/SpeciesTree_rooted_ultra.txt \
-o /lustre04/scratch/celphin/Oxyria/CAFE/Rosaceae_data/output/error_model -e --cores 30

# Best matches are: 0.0027606901630401,0.00033345133146483
# Final -lnL: 129700.06255048
# 128 values were attempted (1% rejected)
# Score (-lnL): 129700.06255048
# Maximum possible lambda for this topology: 0.01611825193182

#----------------------------------
cd /lustre04/scratch/celphin/Oxyria/CAFE/Rosaceae_data/output
more error_model/Base_error_model.txt #view results

# maxcnt: 170
# cntdiff: -1 0 1
# 0 0 0.999667 0.000333451
# 1 0.000333451 0.999333 0.000333451

##########################################
# Run again with error model ( note no space after -e)

cd /lustre04/scratch/celphin/Oxyria/CAFE
#no error model
cafe5 -i /lustre04/scratch/celphin/Oxyria/CAFE/Rosaceae_data/input_data/hog_gene_counts.tsv \
-t /lustre04/scratch/celphin/Oxyria/CAFE/Rosaceae_data/input_data/SpeciesTree_rooted_ultra.txt \
-o /lustre04/scratch/celphin/Oxyria/CAFE/Rosaceae_data/output/single_lambda_noe --cores 30 

# Best match is: 0.0027688392610675
# Final -lnL: 129700.26838166
# Score (-lnL): 129700.26838166
# Maximum possible lambda for this topology: 0.01611825193182

#-------------
#No among family rate variation, try with and without error model
cafe5 -i /lustre04/scratch/celphin/Oxyria/CAFE/Rosaceae_data/input_data/hog_gene_counts.tsv \
-t /lustre04/scratch/celphin/Oxyria/CAFE/Rosaceae_data/input_data/SpeciesTree_rooted_ultra.txt \
-o /lustre04/scratch/celphin/Oxyria/CAFE/Rosaceae_data/output/single_lambda \
-e/lustre04/scratch/celphin/Oxyria/CAFE/Rosaceae_data/output/error_model/Base_error_model.txt --cores 30 

# Best match is: 0.0027606603304254
# Final -lnL: 129700.07140915
# Score (-lnL): 129700.07140915
# Maximum possible lambda for this topology: 0.01611825193182

#--------------------
# #Yes among family rate variation, estimate lambda and alpha and three discrete gamma rate categories
cafe5 -i /lustre04/scratch/celphin/Oxyria/CAFE/Rosaceae_data/input_data/hog_gene_counts.tsv \
-t /lustre04/scratch/celphin/Oxyria/CAFE/Rosaceae_data/input_data/SpeciesTree_rooted_ultra.txt \
-k 3 -o /lustre04/scratch/celphin/Oxyria/CAFE/Rosaceae_data/output/three_lambda --cores 30 

# Best matches are: 0.0027948489079349,4.4995382788564
# Final -lnL: 129383.23520309
# Score (-lnL): 129383.23520309
# Maximum possible lambda for this topology: 0.01611825193182

#--------------------------
# #with error model
cafe5 -i /lustre04/scratch/celphin/Oxyria/CAFE/Rosaceae_data/input_data/hog_gene_counts.tsv \
-t /lustre04/scratch/celphin/Oxyria/CAFE/Rosaceae_data/input_data/SpeciesTree_rooted_ultra.txt \
-k 3 -o /lustre04/scratch/celphin/Oxyria/CAFE/Rosaceae_data/output/three_lambda_error \
-e/lustre04/scratch/celphin/Oxyria/CAFE/Rosaceae_data/output/error_model/Base_error_model.txt --cores 30

# Best matches are: 0.0027865957370893,4.4630399317598
# Final -lnL: 129381.49737525
# Score (-lnL): 129381.49737525
# Maximum possible lambda for this topology: 0.01611825193182

#-----------------------
#Now analyze the large families using the previous values of lambda (skipped)
# cafe5 -i large_/lustre04/scratch/celphin/Oxyria/CAFE/Rosaceae_data/input_data/hog_gene_counts.tsv \
# -t /lustre04/scratch/celphin/Oxyria/CAFE/Rosaceae_data/input_data/SpeciesTree_rooted_ultra.txt \
# -l 0.001319516 -o /lustre04/scratch/celphin/Oxyria/CAFE/Rosaceae_data/output/large_results --cores 30 >


##########################################
# Explore data
cd /lustre04/scratch/celphin/Oxyria/CAFE/Rosaceae_data/output/error_model/

more Base_asr.tre 
#gives reconstructed numbers of each gene family at each node

#--------------------------------
more Base_branch_probabilities.tab 
#gives calculated likelihood of the gene family size at each node for each gene family

#------------------------------
more Base_clade_results.txt 
#how many families are expanded or contracted at each node
# #Taxon_ID                    Increase  Decrease
                       # <12>    284     44
       # Dryas_octopetala<1>     378     795
   # Pyrus_bretschneideri<5>     1315    106
       # Malus_sylvestris<6>     1226    153
                       # <9>     15      44
                       # <10>    225     294
            # Rosa_rugosa<2>     1063    161
                       # <11>    26      176
     # Argentina_anserina<3>     333     460
         # Fragaria_vesca<4>     561     218
                       # <13>    6747    102
         # Prunus_persica<7>     345     679

#--------------------------------------------------
cd /lustre04/scratch/celphin/Oxyria/CAFE/Rosaceae_data/output/

grep "Lambda" ./*/*_results.txt
# ./error_model/Base_results.txt:Lambda: 0.0027606901630401
# ./single_lambda/Base_results.txt:Lambda: 0.0027606603304254
# ./single_lambda_noe/Base_results.txt:Lambda: 0.0027688392610675
# ./three_lambda_error/Gamma_results.txt:Lambda: 0.0027865957370893
# ./three_lambda/Gamma_results.txt:Lambda: 0.0027948489079349


grep "Epsilon" ./*/*_results.txt
# ./error_model/Base_results.txt:Epsilon: 0.000333451
# ./single_lambda/Base_results.txt:Epsilon: 0.000333451
# ./three_lambda_error/Gamma_results.txt:Epsilon: 0.000333451


#------------------------------------
# run many times to check lambdas 

printf '%s\n' {1..50} | parallel cafe5 -i /lustre04/scratch/celphin/Oxyria/CAFE/Rosaceae_data/input_data/hog_gene_counts.tsv \
-t /lustre04/scratch/celphin/Oxyria/CAFE/Rosaceae_data/input_data/SpeciesTree_rooted_ultra.txt \
-o /lustre04/scratch/celphin/Oxyria/CAFE/Rosaceae_data/output/error_model_{1} \
-e --cores 30 >> \
error_model_100runs.log 2>> error_model_100runs.errorlog


#Now analyze the large families using the previous values of lambda
# printf '%s\n' {1..50} | parallel cafe5 -i large_/lustre04/scratch/celphin/Oxyria/CAFE/Rosaceae_data/input_data/hog_gene_counts.tsv \
# -t /lustre04/scratch/celphin/Oxyria/CAFE/Rosaceae_data/input_data/SpeciesTree_rooted_ultra.txt -l 0.001488216 -o large_results_{1} >> \
# large_results_100runs.log


#-----------------------------------------
#Summarise the results of all runs

grep "Lambda" ./error_model_*/Base_results.txt | sed 's/^.* //g' > lambdas.txt
grep "Epsilon" ./error_model_*/Base_results.txt | sed 's/^.* //g' > epsilons.txt

#----------------------------
module load StdEnv/2020 r/4.2.2 

R

epsilons <- read.table("epsilons.txt")
lambdas <- read.table("lambdas.txt")

mean(epsilons[,1])
sd(epsilons[,1])

mean(lambdas[,1])
sd(lambdas[,1])
q()
n

# > mean(epsilons[,1])
# [1] 0.0003320915
# > sd(epsilons[,1])
# [1] 1.477355e-05
# >
# > mean(lambdas[,1])
# [1] 0.002760916
# > sd(lambdas[,1])
# [1] 1.367513e-06

#####################
# Plotting results
# https://github.com/moshi4/CafePlotter

cd /lustre04/scratch/celphin/Oxyria/CAFE/

module load StdEnv/2023 python/3.11.5 scipy-stack/2023b
pip install cafeplotter

cafeplotter -i /lustre04/scratch/celphin/Oxyria/CAFE/Rosaceae_data/output/error_model \
-o /lustre04/scratch/celphin/Oxyria/CAFE/Rosaceae_data/output/error_model/Cafe_plots


#####################################################
# Now we can investigate which orthogroups are expanded/contracted
# See next notes file