##########################
# Running CAFE - gene enrichment analyses on Brassicaceae genomes on beluga 
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

tmux new-session -s CAFE1
tmux attach-session -t CAFE1

# make directories
cd /home/celphin/scratch/Oxyria/CAFE
mkdir Brassicaceae_data; cd Brassicaceae_data
mkdir input_data
mkdir output
mkdir analysis
mkdir expanded
mkdir contracted

# copy over important files
cd /home/celphin/scratch/Oxyria/GeneSpace/Brassicaceae_genomes/orthofinder/Results_Oct09
cp Species_Tree/SpeciesTree_rooted.txt /lustre04/scratch/celphin/Oxyria/CAFE/Brassicaceae_data/input_data
cp Phylogenetic_Hierarchical_Orthogroups/N0.tsv /lustre04/scratch/celphin/Oxyria/CAFE/Brassicaceae_data/input_data


#--------------------------
# in R
cd /lustre04/scratch/celphin/Oxyria/CAFE/Brassicaceae_data/input_data

# remove Arabidopsis thaliana since poor parsing of annotation

module load StdEnv/2020 r/4.2.2 
R

#install.packages("ape")
# https://rdrr.io/cran/ape/man/chronos.html
# https://phylobotanist.blogspot.com/2016/12/how-to-set-multiple-calibration-points.html

library(ape)
library(data.table)

# make orthofinder's species tree ultrametric. It should be already binary and rooted as required by cafe.
setwd("/home/celphin/scratch/Oxyria/CAFE/Brassicaceae_data/input_data")
tre <- read.tree('./SpeciesTree_rooted.txt')

tre <- drop.tip(tre, "Arabidopsis_thaliana") ## taxon can be a single taxon or a vector of taxon names

stopifnot(is.binary(tre))
stopifnot(is.rooted(tre))

#  Brassica and Arabidopsis diverged about 43MYA
#  A. thaliana and Arabidopsis lyrata occurred about 13 Mya
# https://pubmed.ncbi.nlm.nih.gov/20921408/
# https://www.sciencedirect.com/science/article/abs/pii/S1055790316000658

# set calibrations
cal <- makeChronosCalib(tre, node="root", age.min=38, age.max=54) 
chronos.control() 
# $tol= 1e-08
# $iter.max=10000
# $eval.max=10000
# $nb.rate.cat= 10
# $dual.iter.max= 20
# $epsilon=1e-06

if(is.ultrametric(tre)) {
    utre <- tre
} else{
    utre <- chronos(tre, lambda = 1, model = "discrete", calibration = cal, control = chronos.control())
}

# Setting initial dates...
# Fitting in progress... get a first set of estimates
         # (Penalised) log-lik = -2.133754
# Optimising rates... frequencies... dates... -2.133754
# Optimising rates... frequencies... dates... -1.919909
# Optimising rates... frequencies... dates... -1.894968
# Optimising rates... frequencies... dates... -1.894298
# Optimising rates... frequencies... dates... -1.894289

# log-Lik = -1.894289
# PHIIC = 51.79


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

# remove Arabdiopsis thaliana because low gene counts due to poor parsing of annotation file
hog <- hog[-which(hog$species=="Arabidopsis_thaliana"),]

# Exclude HOGs present in all 7 species
keep <- hog[, .N, HOG][N > 5]$HOG
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

#--------------------------
wc -l hog_gene_counts.tsv
8471

#---------------------------------
# The tree file should contain a binary, rooted, ultrametric, tree in Newick format. 
# Typically one obtains this tree using one of several molecular dating methods. 
# If you are unsure if your tree is binary, rooted, or ultrametric CAFE will 
# report this when you try to use it for an analysis. Alternatively, you can use the R package,
# Ape with its included functions: is.ultrametric, is.rooted, and is.binary.
# download tree and view in http://etetoolkit.org/treeview/

more SpeciesTree_rooted_ultra.txt

((Brassica_oleracea:34.45426048,Thlaspi_arvense:34.45426048)0.363985:3.545739519,(Cochlearia_groenl
andica:35.7725719,(Arabidopsis_lyrata:31.37649352,(Draba_nivalis:23.16672676,Arabis_alpina:23.16672
676)0.767655:8.209766754)0.330489:4.396078379)0.363985:2.227428105);


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

tmux new-session -s CAFE1
tmux attach-session -t CAFE1

salloc -c30 --time 3:00:00 --mem 191000M --account def-cronk

module load StdEnv/2020 intel/2020.1.217 cafe5/5.1.0

cd /lustre04/scratch/celphin/Oxyria/CAFE
cafe5 -i /lustre04/scratch/celphin/Oxyria/CAFE/Brassicaceae_data/input_data/hog_gene_counts.tsv \
-t /lustre04/scratch/celphin/Oxyria/CAFE/Brassicaceae_data/input_data/SpeciesTree_rooted_ultra.txt \
-o /lustre04/scratch/celphin/Oxyria/CAFE/Brassicaceae_data/output/error_model -e --cores 30

# Best matches are: 0.0047988816744054,-0.087780928308963
# Final -lnL: 75214.762709139
# 144 values were attempted (2% rejected)
# Score (-lnL): 75214.762709139
# Maximum possible lambda for this topology: 0.027954378085966

#----------------------------------
cd /lustre04/scratch/celphin/Oxyria/CAFE/Brassicaceae_data/output
more error_model/Base_error_model.txt #view results

# maxcnt: 170
# cntdiff: -1 0 1
# 0 0 1.08778 -0.0877809
# 1 -0.0877809 1.17556 -0.0877809

##########################################
# Run again with error model ( note no space after -e)

cd /lustre04/scratch/celphin/Oxyria/CAFE
#no error model
cafe5 -i /lustre04/scratch/celphin/Oxyria/CAFE/Brassicaceae_data/input_data/hog_gene_counts.tsv \
-t /lustre04/scratch/celphin/Oxyria/CAFE/Brassicaceae_data/input_data/SpeciesTree_rooted_ultra.txt \
-o /lustre04/scratch/celphin/Oxyria/CAFE/Brassicaceae_data/output/single_lambda_noe --cores 30 

# Best match is: 0.0028764690589508
# Final -lnL: 75710.562329231
# Score (-lnL): 75710.562329231
# Maximum possible lambda for this topology: 0.027954378085966

#-----------------
#No among family rate variation, try with an without error model
cafe5 -i /lustre04/scratch/celphin/Oxyria/CAFE/Brassicaceae_data/input_data/hog_gene_counts.tsv \
-t /lustre04/scratch/celphin/Oxyria/CAFE/Brassicaceae_data/input_data/SpeciesTree_rooted_ultra.txt \
-o /lustre04/scratch/celphin/Oxyria/CAFE/Brassicaceae_data/output/single_lambda \
-e/lustre04/scratch/celphin/Oxyria/CAFE/Brassicaceae_data/output/error_model/Base_error_model.txt --cores 30 

# Best match is: 0.004798922266262
# Final -lnL: 75214.854185751
# Score (-lnL): 75214.854185751
# Maximum possible lambda for this topology: 0.027954378085966

#----------------------
# #Yes among family rate variation, estimate lambda and alpha and three discrete gamma rate categories
cafe5 -i /lustre04/scratch/celphin/Oxyria/CAFE/Brassicaceae_data/input_data/hog_gene_counts.tsv \
-t /lustre04/scratch/celphin/Oxyria/CAFE/Brassicaceae_data/input_data/SpeciesTree_rooted_ultra.txt \
-k 2 -o /lustre04/scratch/celphin/Oxyria/CAFE/Brassicaceae_data/output/three_lambda --cores 30 

# Best matches are: 0.0028825354771294,3.2041748501872
# Final -lnL: 75626.160901245
# Score (-lnL): 75626.160901245
# Maximum possible lambda for this topology: 0.027954378085966

#--------------------------------
# #with error model
cafe5 -i /lustre04/scratch/celphin/Oxyria/CAFE/Brassicaceae_data/input_data/hog_gene_counts.tsv \
-t /lustre04/scratch/celphin/Oxyria/CAFE/Brassicaceae_data/input_data/SpeciesTree_rooted_ultra.txt \
-k 3 -o /lustre04/scratch/celphin/Oxyria/CAFE/Brassicaceae_data/output/three_lambda_error \
-e/lustre04/scratch/celphin/Oxyria/CAFE/Brassicaceae_data/output/error_model/Base_error_model.txt --cores 30

# Best matches are: 0.0038969423002203,0.95876428593321
# Final -lnL: 75158.37890783
# Score (-lnL):  75158.37890783
# Maximum possible lambda for this topology: 0.027954378085966

#-------------------------------
#Now analyze the large families using the previous values of lambda (skipped)
 # cafe5 -i /lustre04/scratch/celphin/Oxyria/CAFE/Brassicaceae_data/input_data/large_hog_gene_counts.tsv \
 # -t /lustre04/scratch/celphin/Oxyria/CAFE/Brassicaceae_data/input_data/SpeciesTree_rooted_ultra.txt \
 # -l 0.01487 -o /lustre04/scratch/celphin/Oxyria/CAFE/Brassicaceae_data/output/large_results --cores 30 

# large_hog_gene_counts.tsv is empty

##########################################
# Explore data
cd /lustre04/scratch/celphin/Oxyria/CAFE/Brassicaceae_data/output/error_model

more Base_asr.tre 
#gives reconstructed numbers of each gene family at each node

#--------------------------------
more Base_branch_probabilities.tab 
#gives calculated likelihood of the gene family size at each node for each gene family

#------------------------------
more Base_clade_results.txt 
#how many families are expanded or contracted at each node


# #Taxon_ID                 Increase      Decrease
                        # <8>     8       5
         # Thlaspi_arvense<2>     630     412
       # Brassica_oleracea<1>     3182    151
                        # <9>      0       69
 # Cochlearia_groenlandica<3>      1182    260
                        # <10>      8       88
      # Arabidopsis_lyrata<4>       581     202
                        # <11>       70      75
           # Draba_nivalis<5>        809     92
           # Arabis_alpina<6>        362     189

#--------------------------------------------------
cd /lustre04/scratch/celphin/Oxyria/CAFE/Brassicaceae_data/output/

grep "Lambda" ./*/*_results.txt
# ./error_model/Base_results.txt:Lambda: 0.0047988816744054
# ./single_lambda/Base_results.txt:Lambda: 0.004798922266262
# ./single_lambda_noe/Base_results.txt:Lambda: 0.0028764690589508
# ./three_lambda_error/Gamma_results.txt:Lambda: 0.0038969423002203
# ./three_lambda/Gamma_results.txt:Lambda: 0.0028825354771294

grep "Epsilon" ./*/*_results.txt
# ./error_model/Base_results.txt:Epsilon: -0.0877809
# ./single_lambda/Base_results.txt:Epsilon: -0.0877809
# ./three_lambda_error/Gamma_results.txt:Epsilon: -0.0877809

#------------------------------------
# run many times to check lambdas 
cd /lustre04/scratch/celphin/Oxyria/CAFE/Brassicaceae_data/output
module load StdEnv/2020 intel/2020.1.217 cafe5/5.1.0

printf '%s\n' {1..10} | parallel cafe5 -i /lustre04/scratch/celphin/Oxyria/CAFE/Brassicaceae_data/input_data/hog_gene_counts.tsv \
-t /lustre04/scratch/celphin/Oxyria/CAFE/Brassicaceae_data/input_data/SpeciesTree_rooted_ultra.txt -e \
-o /lustre04/scratch/celphin/Oxyria/CAFE/Brassicaceae_data/output/error_model_{1} >> \
error_model_100runs.log 2>> error_model_100runs.errorlog

#MORE!
printf '%s\n' {11..20} | parallel cafe5 -i /lustre04/scratch/celphin/Oxyria/CAFE/Brassicaceae_data/input_data/hog_gene_counts.tsv \
-t /lustre04/scratch/celphin/Oxyria/CAFE/Brassicaceae_data/input_data/SpeciesTree_rooted_ultra.txt -e \
-o /lustre04/scratch/celphin/Oxyria/CAFE/Brassicaceae_data/output/error_model_{1} >> \
error_model_100runs.log 2>> error_model_100runs.errorlog

#MORE
printf '%s\n' {21..50} | parallel cafe5 -i /lustre04/scratch/celphin/Oxyria/CAFE/Brassicaceae_data/input_data/hog_gene_counts.tsv \
-t /lustre04/scratch/celphin/Oxyria/CAFE/Brassicaceae_data/input_data/SpeciesTree_rooted_ultra.txt -e \
-o /lustre04/scratch/celphin/Oxyria/CAFE/Brassicaceae_data/output/error_model_{1} >> \
error_model_100runs.log 2>> error_model_100runs.errorlog

#Now analyze the large families using the previous values of lambda
# printf '%s\n' {1..50} | parallel cafe5 -i large_/lustre04/scratch/celphin/Oxyria/CAFE/Brassicaceae_data/input_data/hog_gene_counts.tsv \
# -t /lustre04/scratch/celphin/Oxyria/CAFE/Brassicaceae_data/input_data/SpeciesTree_rooted_ultra.txt -l 0.001488216 -o large_results_{1} >> \
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
# [1] -0.08675126
# > sd(epsilons[,1])
# [1] 0.002076295

# > mean(lambdas[,1])
# [1] 0.00477913
# > sd(lambdas[,1])
# [1] 4.010276e-05

#####################
# Plotting results
# https://github.com/moshi4/CafePlotter

cd /lustre04/scratch/celphin/Oxyria/CAFE/

module load StdEnv/2023 python/3.11.5 scipy-stack/2023b
pip install cafeplotter

cafeplotter -i /lustre04/scratch/celphin/Oxyria/CAFE/Brassicaceae_data/output/error_model \
-o /lustre04/scratch/celphin/Oxyria/CAFE/Brassicaceae_data/output/error_model/Cafe_plots

#####################################################
# Now we can investigate which orthogroups are expanded/contracted
# See other notes files

