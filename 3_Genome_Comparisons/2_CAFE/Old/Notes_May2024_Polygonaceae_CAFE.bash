##########################
# Running CAFE - gene enrichment analyses on Polygonaceae genomes on beluga 
# Mar 2024
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

tmux new-session -s CAFE
tmux attach-session -t CAFE

# copy over important files

cd /home/celphin/scratch/Oxyria/GeneSpace/Polygonaceae_genomes/orthofinder/Results_Jul03

cp Species_Tree/SpeciesTree_rooted.txt /lustre04/scratch/celphin/Oxyria/CAFE/Polygonaceae_data/input_data
cp Phylogenetic_Hierarchical_Orthogroups/N0.tsv /lustre04/scratch/celphin/Oxyria/CAFE/Polygonaceae_data/input_data

# in R
cd /lustre04/scratch/celphin/Oxyria/CAFE/Polygonaceae_data/input_data

module load StdEnv/2020 r/4.2.2 
R

#install.packages("ape")

library(ape)
library(data.table)

# make orthofinder's species tree ultrametric. It should be already binary and rooted as required by cafe.
setwd("/home/celphin/scratch/Oxyria/CAFE/Polygonaceae_data/input_data")
tre <- read.tree('./SpeciesTree_rooted.txt')
stopifnot(is.binary(tre))
stopifnot(is.rooted(tre))

if(is.ultrametric(tre)) {
    utre <- tre
} else{
    utre <- chronos(tre)
}

# Setting initial dates...
# Fitting in progress... get a first set of estimates
         # (Penalised) log-lik = -2.614781
# Optimising rates... dates... -2.614781
# Optimising rates... dates... -2.614727

# log-Lik = -2.613039
# PHIIC = 57.24

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

# Exclude HOGs present in only 1 species
keep <- hog[, .N, HOG][N > 1]$HOG
hog <- hog[HOG %in% keep]

counts <- dcast(hog, HOG ~ species, value.var='n', fill=0)
counts[, Desc := 'n/a']
setcolorder(counts, 'Desc')
fwrite(counts, 'hog_gene_counts.tsv', sep='\t')

#---------------------------------
# The tree file should contain a binary, rooted, ultrametric, tree in Newick format. 
# Typically one obtains this tree using one of several molecular dating methods. 
# If you are unsure if your tree is binary, rooted, or ultrametric CAFE will 
# report this when you try to use it for an analysis. Alternatively, you can use the R package,
# Ape with its included functions: is.ultrametric, is.rooted, and is.binary.

# Feb23
# ((Polygunum_aviculare_H0:0.170533,((Rheum_nobile_H0:0.0423251,Rheum_tangaticum_H0:0.050222)N6:0.0642639
# ,(Oxyria_digyna_H0:0.0166409,Oxyria_digyna_H1:0.0169317)N7:0.136957)N3:0.0364266)N1:0.0634775,((Fagopyr
# um_escelentum_H1:0.00875405,Fagopyrum_escelentum_H2:0.00914042)N4:0.038296,(Fagopyrum_tataricum_H0:0.04
# 94432,(Fagopyrum_tataricum_H1:0.00106153,Fagopyrum_tataricum_H2:0.000814587)N8:0.0163544)N5:0.0381724)N
# 2:0.0634775)N0;

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

tmux new-session -s CAFE
tmux attach-session -t CAFE

salloc -c30 --time 3:00:00 --mem 191000M --account def-cronk

module load StdEnv/2020 intel/2020.1.217 cafe5/5.1.0

cd /lustre04/scratch/celphin/Oxyria/CAFE
cafe5 -i /lustre04/scratch/celphin/Oxyria/CAFE/Polygonaceae_data/input_data/hog_gene_counts.tsv \
-t /lustre04/scratch/celphin/Oxyria/CAFE/Polygonaceae_data/input_data/SpeciesTree_rooted_ultra.txt \
-o /lustre04/scratch/celphin/Oxyria/CAFE/Polygonaceae_data/output/error_model -e --cores 30

# Filtering families not present at the root from: 42683 to 22790
# No root family size distribution specified, using uniform distribution
# Optimizer strategy: Nelder-Mead with similarity cutoff
# Iterations: 300
# Expansion: 2
# Reflection: 1
# Starting Search for Initial Parameter Values
# Calculating probability: epsilon=0.1, lambda=0.054274828693924
# Score (-lnL): 357905.70862791
# Calculating probability: epsilon=0.1, lambda=0.054274828693924
# Completed 67 iterations
# Time: 0H 0M 33S
# Best matches are: 0.31132230351936,0.038010282712594
# Final -lnL: 325352.13580164

# 131 values were attempted (0% rejected)
# Inferring processes for Base model
# Score (-lnL): 325352.13580164
# Maximum possible lambda for this topology: 1.3735308813446
# Computing pvalues...
# done!
# Starting reconstruction processes for Base model
# Done!

#----------------------------------
cd /lustre04/scratch/celphin/Oxyria/CAFE/Polygonaceae_data/output
more error_model/Base_error_model.txt #view results

# maxcnt: 170
# cntdiff: -1 0 1
# 0 0 0.96199 0.0380103
# 1 0.0380103 0.923979 0.0380103

##########################################
# Run again with error model ( note no space after -e)

cd /lustre04/scratch/celphin/Oxyria/CAFE
#no error model
cafe5 -i /lustre04/scratch/celphin/Oxyria/CAFE/Polygonaceae_data/input_data/hog_gene_counts.tsv \
-t /lustre04/scratch/celphin/Oxyria/CAFE/Polygonaceae_data/input_data/SpeciesTree_rooted_ultra.txt \
-o /lustre04/scratch/celphin/Oxyria/CAFE/Polygonaceae_data/output/single_lambda_noe --cores 30 

#No among family rate variation, try with an without error model
cafe5 -i /lustre04/scratch/celphin/Oxyria/CAFE/Polygonaceae_data/input_data/hog_gene_counts.tsv \
-t /lustre04/scratch/celphin/Oxyria/CAFE/Polygonaceae_data/input_data/SpeciesTree_rooted_ultra.txt \
-o /lustre04/scratch/celphin/Oxyria/CAFE/Polygonaceae_data/output/single_lambda \
-e/lustre04/scratch/celphin/Oxyria/CAFE/Polygonaceae_data/output/error_model/Base_error_model.txt --cores 30 

# #Yes among family rate variation, estimate lambda and alpha and three discrete gamma rate categories
cafe5 -i /lustre04/scratch/celphin/Oxyria/CAFE/Polygonaceae_data/input_data/hog_gene_counts.tsv \
-t /lustre04/scratch/celphin/Oxyria/CAFE/Polygonaceae_data/input_data/SpeciesTree_rooted_ultra.txt \
-k 3 -o /lustre04/scratch/celphin/Oxyria/CAFE/Polygonaceae_data/output/three_lambda --cores 30 

# #with error model
cafe5 -i /lustre04/scratch/celphin/Oxyria/CAFE/Polygonaceae_data/input_data/hog_gene_counts.tsv \
-t /lustre04/scratch/celphin/Oxyria/CAFE/Polygonaceae_data/input_data/SpeciesTree_rooted_ultra.txt \
-k 3 -o /lustre04/scratch/celphin/Oxyria/CAFE/Polygonaceae_data/output/three_lambda_error \
-e/lustre04/scratch/celphin/Oxyria/CAFE/Polygonaceae_data/output/error_model/Base_error_model.txt --cores 30

#Now analyze the large families using the previous values of lambda (skipped)
# cafe5 -i large_/lustre04/scratch/celphin/Oxyria/CAFE/Polygonaceae_data/input_data/hog_gene_counts.tsv \
# -t /lustre04/scratch/celphin/Oxyria/CAFE/Polygonaceae_data/input_data/SpeciesTree_rooted_ultra.txt \
# -l 0.001319516 -o /lustre04/scratch/celphin/Oxyria/CAFE/Polygonaceae_data/output/large_results --cores 30 >


##########################################
# Explore data
cd /lustre04/scratch/celphin/Oxyria/CAFE/Polygonaceae_data/output/

more Base_asr.tre 
#gives reconstructed numbers of each gene family at each node

#nexus
# BEGIN TREES;
  # TREE N0.HOG0000000 = ((Polygunum_aviculare_H0<1>_0:0.728051,((Oxyria_digyna_H1<2>*_2:0.0649694,Oxyria_digyna_H0<3>_1:0.0649694)<14>_1:0.507763,(Rheum_nobile_H0<4>_0:0.238872,Rheum_tangaticum_H0<5>_0:0.238872)<15>_0:0.333861)<13>_1:0.155318)<12>_1:0.271949,((Fagopyrum_tataricum_H0<6>*_8:0.256333,(Fagopyrum_tataricum_H1<7>_0:0.01234,Fagopyrum_tataricum_H2<8>_0:0.01234)<18>_0:0.243993)<17>_1:0.267586,(Fagopyrum_escelentum_H1<9>_0:0.0966194,Fagopyrum_escelentum_H2<10>_0:0.0966194)<19>_0:0.427299)<16>_1:0.476081)<11>_1;
  # TREE N0.HOG0000010 = ((Polygunum_aviculare_H0<1>_0:0.728051,((Oxyria_digyna_H1<2>*_0:0.0649694,Oxyria_digyna_H0<3>_1:0.0649694)<14>_1:0.507763,(Rheum_nobile_H0<4>_0:0.238872,Rheum_tangaticum_H0<5>_0:0.238872)<15>_0:0.333861)<13>_1:0.155318)<12>_1:0.271949,((Fagopyrum_tataricum_H0<6>*_42:0.256333,(Fagopyrum_tataricum_H1<7>_1:0.01234,Fagopyrum_tataricum_H2<8>_1:0.01234)<18>*_1:0.243993)<17>*_6:0.267586,(Fagopyrum_escelentum_H1<9>_2:0.0966194,Fagopyrum_escelentum_H2<10>*_0:0.0966194)<19>_2:0.427299)<16>_3:0.476081)<11>_2;
  # TREE N0.HOG0000011 = ((Polygunum_aviculare_H0<1>_0:0.728051,((Oxyria_digyna_H1<2>*_0:0.0649694,Oxyria_digyna_H0<3>_1:0.0649694)<14>_1:0.507763,(Rheum_nobile_H0<4>_0:0.238872,Rheum_tangaticum_H0<5>_0:0.238872)<15>_0:0.333861)<13>_1:0.155318)<12>_1:0.271949,((Fagopyrum_tataricum_H0<6>*_40:0.256333,(Fagopyrum_tataricum_H1<7>_0:0.01234,Fagopyrum_tataricum_H2<8>_0:0.01234)<18>*_0:0.243993)<17>*_3:0.267586,(Fagopyrum_escelentum_H1<9>_0:0.0966194,Fagopyrum_escelentum_H2<10>_0:0.0966194)<19>_0:0.427299)<16>_1:0.476081)<11>_1;
  # TREE N0.HOG0000014 = ((Polygunum_aviculare_H0<1>_0:0.728051,((Oxyria_digyna_H1<2>_1:0.0649694,Oxyria_digyna_H0<3>*_0:0.0649694)<14>_1:0.507763,(Rheum_nobile_H0<4>_0:0.238872,Rheum_tangaticum_H0<5>_0:0.238872)<15>_0:0.333861)<13>_1:0.155318)<12>_1:0.271949,((Fagopyrum_tataricum_H0<6>*_29:0.256333,(Fagopyrum_tataricum_H1<7>*_0:0.01234,Fagopyrum_tataricum_H2<8>_1:0.01234)<18>*_1:0.243993)<17>*_3:0.267586,(Fagopyrum_escelentum_H1<9>_0:0.0966194,Fagopyrum_escelentum_H2<10>_0:0.0966194)<19>_0:0.427299)<16>_1:0.476081)<11>_1;
  # TREE N0.HOG0000017 = ((Polygunum_aviculare_H0<1>_0:0.728051,((Oxyria_digyna_H1<2>_1:0.0649694,Oxyria_digyna_H0<3>*_0:0.0649694)<14>_1:0.507763,(Rheum_nobile_H0<4>_0:0.238872,Rheum_tangaticum_H0<5>_0:0.238872)<15>_0:0.333861)<13>_1:0.155318)<12>*_1:0.271949,((Fagopyrum_tataricum_H0<6>*_88:0.256333,(Fagopyrum_tataricum_H1<7>_1:0.01234,Fagopyrum_tataricum_H2<8>_1:0.01234)<18>*_1:0.243993)<17>*_14:0.267586,(Fagopyrum_escelentum_H1<9>*_12:0.0966194,Fagopyrum_escelentum_H2<10>*_0:0.0966194)<19>_8:0.427299)<16>*_9:0.476081)<11>_4;
  # TREE N0.HOG0000018 = ((Polygunum_aviculare_H0<1>_0:0.728051,((Oxyria_digyna_H1<2>*_0:0.0649694,Oxyria_digyna_H0<3>_1:0.0649694)<14>_1:0.507763,(Rheum_nobile_H0<4>_0:0.238872,Rheum_tangaticum_H0<5>_0:0.238872)<15>_0:0.333861)<13>_1:0.155318)<12>_1:0.271949,((Fagopyrum_tataricum_H0<6>*_6:0.256333,(Fagopyrum_tataricum_H1<7>_0:0.01234,Fagopyrum_tataricum_H2<8>_0:0.01234)<18>_0:0.243993)<17>_1:0.267586,(Fagopyrum_escelentum_H1<9>_0:0.0966194,Fagopyrum_escelentum_H2<10>_0:0.0966194)<19>_0:0.427299)<16>_1:0.476081)<11>_1;


#--------------------------------
more Base_branch_probabilities.tab 
#gives calculated likelihood of the gene family size at each node for each gene family

# FamilyID        Polygunum_aviculare_H0<1>       Oxyria_digyna_H1<2>     Oxyria_digyna_H0<3>     Rheum_nobile_H0<4>    Rheum_tangaticum_H0<5>  Fagopyrum_tataricum_H0<6>       Fagopyrum_tataricum_H1<7>    Fagopyrum_tataricum_H2<8>        Fagopyrum_escelentum_H1<9>      Fagopyrum_escelentum_H2<10>     <11> <12>     <13>    <14>    <15>    <16>    <17>    <18>    <19>
# N0.HOG0000000   0.287782        0.0124962       0.524695        0.5     0.5     2.99402e-08     0.5  0.5      0.5     0.5     N/A     0.59317 0.556777        0.654618        0.162752        0.647373     0.591991 0.125242        0.197834
# N0.HOG0000010   0.287782        0.0368868       0.524695        0.5     0.5     1.08253e-32     0.504775      0.504775        0.568808        0.000869829     N/A     0.242626        0.556777        0.654618      0.162752        0.15786 0.00423242      2.78041e-05     0.423383
# N0.HOG0000011   0.287782        0.0368868       0.524695        0.5     0.5     1.75022e-36     0.5  0.5      0.5     0.5     N/A     0.59317 0.556777        0.654618        0.162752        0.647373     0.00462824       0.000964007     0.197834
# N0.HOG0000014   0.287782        0.524695        0.0368868       0.5     0.5     1.9808e-25      0.00715651    0.504775        0.5     0.5     N/A     0.59317 0.556777        0.654618        0.162752     0.647373 0.00462824      0.0151824       0.197834
# N0.HOG0000017   0.287782        0.524695        0.0368868       0.5     0.5     2.05592e-63     0.504775      0.504775        0.000307081     5.13341e-12     N/A     0.00334109      0.556777        0.654618      0.162752        0.00192581      0.00207076      1.74354e-13     0.656342

#------------------------------
more Base_clade_results.txt 
#how many families are expanded or contracted at each node

# #Taxon_ID       Increase        Decrease
# Oxyria_digyna_H0<3>     925     760
                # <15>    661     1972
                # <12>    56      269
# Fagopyrum_escelentum_H1<9>      1926    1035
# Rheum_tangaticum_H0<5>  2402    3644
# Polygunum_aviculare_H0<1>       1613    6095
                     # <16>    1580    972
# Fagopyrum_tataricum_H2<8>       832     684
                       # <17>    528     1279
# Fagopyrum_tataricum_H0<6>       914     3300
                     # <13>    66      914
                     # <14>    1044    4394
# Oxyria_digyna_H1<2>     665     926
# Rheum_nobile_H0<4>      1797    1587
              # <18>    1366    1404
# Fagopyrum_tataricum_H1<7>       760     563
                      # <19>    1713    2066
# Fagopyrum_escelentum_H2<10>     1574    1365


#--------------------------------------------------
cd /lustre04/scratch/celphin/Oxyria/CAFE/Polygonaceae_data/output/

grep "Lambda" ./*/*_results.txt
# ./error_model/Base_results.txt:Lambda: 0.24217528611495
# ./single_lambda/Base_results.txt:Lambda: 0.24217638329046
# ./single_lambda_noe/Base_results.txt:Lambda: 0.32892330721424
# ./three_lambda_error/Gamma_results.txt:Lambda: 0.49042723523936
# ./three_lambda/Gamma_results.txt:Lambda: 0.16734859913189


grep "Epsilon" ./*/*_results.txt
# ./error_model/Base_results.txt:Epsilon: 0.0577946
# ./single_lambda/Base_results.txt:Epsilon: 0.0577946
# ./three_lambda_error/Gamma_results.txt:Epsilon: 0.0577946


#------------------------------------
# run many times to check lambdas 

printf '%s\n' {1..10} | parallel cafe5 -i /lustre04/scratch/celphin/Oxyria/CAFE/Polygonaceae_data/input_data/hog_gene_counts.tsv \
-t /lustre04/scratch/celphin/Oxyria/CAFE/Polygonaceae_data/input_data/SpeciesTree_rooted_ultra.txt -e \
-o /lustre04/scratch/celphin/Oxyria/CAFE/Polygonaceae_data/output/error_model_{1} >> \
error_model_100runs.log 2>> error_model_100runs.errorlog

#MORE!
printf '%s\n' {11..20} | parallel cafe5 -i /lustre04/scratch/celphin/Oxyria/CAFE/Polygonaceae_data/input_data/hog_gene_counts.tsv \
-t /lustre04/scratch/celphin/Oxyria/CAFE/Polygonaceae_data/input_data/SpeciesTree_rooted_ultra.txt -e \
-o /lustre04/scratch/celphin/Oxyria/CAFE/Polygonaceae_data/output/error_model_{1} >> \
error_model_100runs.log 2>> error_model_100runs.errorlog

#MORE
printf '%s\n' {21..50} | parallel cafe5 -i /lustre04/scratch/celphin/Oxyria/CAFE/Polygonaceae_data/input_data/hog_gene_counts.tsv \
-t /lustre04/scratch/celphin/Oxyria/CAFE/Polygonaceae_data/input_data/SpeciesTree_rooted_ultra.txt -e \
-o /lustre04/scratch/celphin/Oxyria/CAFE/Polygonaceae_data/output/error_model_{1} >> \
error_model_100runs.log 2>> error_model_100runs.errorlog

#Now analyze the large families using the previous values of lambda
# printf '%s\n' {1..50} | parallel cafe5 -i large_/lustre04/scratch/celphin/Oxyria/CAFE/Polygonaceae_data/input_data/hog_gene_counts.tsv \
# -t /lustre04/scratch/celphin/Oxyria/CAFE/Polygonaceae_data/input_data/SpeciesTree_rooted_ultra.txt -l 0.001488216 -o large_results_{1} >> \
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

# mean(epsilons[,1])
# 0.05741533
# sd(epsilons[,1])
# 0.001128798

# mean(lambdas[,1])
# 0.2424507
# sd(lambdas[,1])
# 0.001008741


#####################################################
# Now we can investigate which orthogroups are expanded/contracted
# See next notes