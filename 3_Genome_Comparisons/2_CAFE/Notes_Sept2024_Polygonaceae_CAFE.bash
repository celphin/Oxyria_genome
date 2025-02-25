##########################
# Running CAFE - gene enrichment analyses on Polygonaceae genomes on beluga 
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
# parsing file from OrthoFinder for CAFE
#https://www.biostars.org/p/9553290/

tmux new-session -s CAFE
tmux attach-session -t CAFE

# make directories
cd /home/celphin/scratch/Oxyria/CAFE
mkdir Polygonaceae_data; cd Polygonaceae_data
mkdir input_data
mkdir output
mkdir analysis
mkdir expanded
mkdir contracted

# copy over important files
cd /home/celphin/scratch/Oxyria/GeneSpace/Polygonaceae_genomes/orthofinder/Results_Oct07
cp Species_Tree/SpeciesTree_rooted.txt /lustre04/scratch/celphin/Oxyria/CAFE/Polygonaceae_data/input_data
cp Phylogenetic_Hierarchical_Orthogroups/N0.tsv /lustre04/scratch/celphin/Oxyria/CAFE/Polygonaceae_data/input_data


#--------------------------
# in R
cd /lustre04/scratch/celphin/Oxyria/CAFE/Polygonaceae_data/input_data

module load StdEnv/2020 r/4.2.2 
R

#install.packages("ape")
# https://rdrr.io/cran/ape/man/chronos.html
# https://phylobotanist.blogspot.com/2016/12/how-to-set-multiple-calibration-points.html

library(ape)
library(data.table)

# make orthofinder's species tree ultrametric. It should be already binary and rooted as required by cafe.
setwd("/home/celphin/scratch/Oxyria/CAFE/Polygonaceae_data/input_data")
tre <- read.tree('./SpeciesTree_rooted.txt')
stopifnot(is.binary(tre))
stopifnot(is.rooted(tre))

# Fagopyrum and Rheum (rhubarb) shared a common ancestor about 48 million years ago
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10123988/ # 48 mya - Fagopyrum and Rheum (rhubarb)
# https://www.nature.com/articles/s42003-023-05248-5 # 28mya - Fagopyrum and Rheum (rhubarb)
# https://www.frontiersin.org/journals/plant-science/articles/10.3389/fpls.2022.893201/full
# https://pubmed.ncbi.nlm.nih.gov/23585884/ # 90-150mya - whole Polygonaceae
# https://www.researchgate.net/figure/BEAST-derived-chronograms-of-Polygonaceae-based-on-74-shared-protein-coding-genes-PCGs_fig4_355973888
# Oxyria and Rheum - ~27mya

# set calibrations
cal <- makeChronosCalib(tre, node="root", age.min=39, age.max=68) 
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
    utre <- chronos(tre, lambda = 1, model = "discrete", calibration = cal, control = chronos.control() )
    attr(utre, "rates")
}

# Setting initial dates...
# Fitting in progress... get a first set of estimates
         # (Penalised) log-lik = -2.388797
# Optimising rates... frequencies... dates... -2.388797
# Optimising rates... frequencies... dates... -2.250021

# log-Lik = -2.250021
# PHIIC = 52.5
 # [1] 0.003769139 

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

# Exclude HOGs present in all 6 species
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

# Exclude HOGs present in only 1 species
keep <- hog[, .N, HOG][N > 5]$HOG
hog <- hog[HOG %in% keep]

counts <- dcast(hog, HOG ~ species, value.var='n', fill=0)
counts[, Desc := 'n/a']
setcolorder(counts, 'Desc')
fwrite(counts, 'large_hog_gene_counts.tsv', sep='\t')

#---------------------------------
wc -l hog_gene_counts.tsv
11401

#---------------------------------
# The tree file should contain a binary, rooted, ultrametric, tree in Newick format. 
# Typically one obtains this tree using one of several molecular dating methods. 
# If you are unsure if your tree is binary, rooted, or ultrametric CAFE will 
# report this when you try to use it for an analysis. Alternatively, you can use the R package,
# Ape with its included functions: is.ultrametric, is.rooted, and is.binary.
# download tree and view in http://etetoolkit.org/treeview/

more SpeciesTree_rooted_ultra.txt
((Fagopyrum_tataricum_H1:16.35614806,Fagopyrum_escelentum_H2:16.35614806)0.98348:34.30301222,(Polyg
unum_aviculare_H0:38.77384923,(Oxyria_digyna_H1:30.99240788,(Rheum_nobile_H0:12.94302104,Rheum_tang
aticum_H0:12.94302104)0.878113:18.04938684)0.658735:7.781441346)0.98348:11.88531105);

########################################
# Running CAFE
# https://github.com/hahnlab/CAFE5?tab=readme-ov-file
# https://github.com/hahnlab/CAFE5/blob/master/docs/tutorial/tutorial.md
# Tutorial https://github.com/elsemikk/Willisornis_Genome_Assembly/blob/master/4.2_Gene_Family_Expansions.md

# Compare scenarios in which the whole phylogeny shares the same (global) lambda vs. scenarios in 
# which different parts of the phylogeny share different (local) lambdas.

# Classify specific gene families as "fastly evolving” in at least two different ways (see documentation below).

# Infer gene family counts at all internal nodes (ancestral populations) of the phylogeny provided as input.
# By comparing different nodes in the phylogeny, the user will be able to tell along which branches 
# gene families have contracted or expanded.

# Account for non-biological factors (e.g., genome sequencing and coverage differences, gene family
# clustering errors, etc.) leading to incorrect gene family counts in input files. This is done with 
# an error model.


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

# Best matches are: 0.002383227827052,0.011967534422874
# Final -lnL: 102704.51782527
# Score (-lnL): 102704.51782527
# Maximum possible lambda for this topology: 0.025790578440334


#----------------------------------
cd /lustre04/scratch/celphin/Oxyria/CAFE/Polygonaceae_data/output
more error_model/Base_error_model.txt #view results

# maxcnt: 170
# cntdiff: -1 0 1
# 0 0 0.988032 0.0119675
# 1 0.0119675 0.976065 0.0119675


##########################################
# Run again with error model ( note no space after -e)

cd /lustre04/scratch/celphin/Oxyria/CAFE
#no error model
cafe5 -i /lustre04/scratch/celphin/Oxyria/CAFE/Polygonaceae_data/input_data/hog_gene_counts.tsv \
-t /lustre04/scratch/celphin/Oxyria/CAFE/Polygonaceae_data/input_data/SpeciesTree_rooted_ultra.txt \
-o /lustre04/scratch/celphin/Oxyria/CAFE/Polygonaceae_data/output/single_lambda_noe --cores 30 

# Best match is: 0.0027165417373796
# Final -lnL: 102786.71895722
# Score (-lnL): 102786.71895722
# Maximum possible lambda for this topology: 0.025790578440334

#------------------------
#No among family rate variation, try with an without error model
cafe5 -i /lustre04/scratch/celphin/Oxyria/CAFE/Polygonaceae_data/input_data/hog_gene_counts.tsv \
-t /lustre04/scratch/celphin/Oxyria/CAFE/Polygonaceae_data/input_data/SpeciesTree_rooted_ultra.txt \
-o /lustre04/scratch/celphin/Oxyria/CAFE/Polygonaceae_data/output/single_lambda \
-e/lustre04/scratch/celphin/Oxyria/CAFE/Polygonaceae_data/output/error_model/Base_error_model.txt --cores 30 

# Best match is: 0.0023831902736841
# Final -lnL: 102704.5178237
# Score (-lnL):  102704.5178237
# Maximum possible lambda for this topology: 0.025790578440334

#--------------------------
# #Yes among family rate variation, estimate lambda and alpha and three discrete gamma rate categories
cafe5 -i /lustre04/scratch/celphin/Oxyria/CAFE/Polygonaceae_data/input_data/hog_gene_counts.tsv \
-t /lustre04/scratch/celphin/Oxyria/CAFE/Polygonaceae_data/input_data/SpeciesTree_rooted_ultra.txt \
-k 3 -o /lustre04/scratch/celphin/Oxyria/CAFE/Polygonaceae_data/output/three_lambda --cores 30 

# Best matches are: 0.002807174772633,1.0699069462293
# Final -lnL: 101945.31223627
# 76 values were attempted (0% rejected)
# Score (-lnL): 101945.31223627
# Maximum possible lambda for this topology: 0.025790578440334

#---------------------------
# #with error model
cafe5 -i /lustre04/scratch/celphin/Oxyria/CAFE/Polygonaceae_data/input_data/hog_gene_counts.tsv \
-t /lustre04/scratch/celphin/Oxyria/CAFE/Polygonaceae_data/input_data/SpeciesTree_rooted_ultra.txt \
-k 3 -o /lustre04/scratch/celphin/Oxyria/CAFE/Polygonaceae_data/output/three_lambda_error \
-e/lustre04/scratch/celphin/Oxyria/CAFE/Polygonaceae_data/output/error_model/Base_error_model.txt --cores 30

# Best matches are: 0.0024360829102418,0.62049499002915
# Final -lnL: 101603.04113835
# 135 values were attempted (1% rejected)
# Score (-lnL): 101603.04113835
# Maximum possible lambda for this topology: 0.025790578440334

#------------------------------
#Now analyze the large families using the previous values of lambda (skipped)
# cafe5 -i /lustre04/scratch/celphin/Oxyria/CAFE/Polygonaceae_data/input_data/large_hog_gene_counts.tsv \
 # -t /lustre04/scratch/celphin/Oxyria/CAFE/Polygonaceae_data/input_data/SpeciesTree_rooted_ultra.txt \
 # -l 0.0154 -o /lustre04/scratch/celphin/Oxyria/CAFE/Polygonaceae_data/output/large_results --cores 30 

# large gene count file empty

##########################################
# Explore data
cd /lustre04/scratch/celphin/Oxyria/CAFE/Polygonaceae_data/output/error_model/

more Base_asr.tre 
#gives reconstructed numbers of each gene family at each node

#--------------------------------
more Base_branch_probabilities.tab 
#gives calculated likelihood of the gene family size at each node for each gene family

#------------------------------
cd /lustre04/scratch/celphin/Oxyria/CAFE/Polygonaceae_data/output/
more ./*/Base_clade_results.txt 


#how many families are expanded or contracted at each node

# #Taxon_ID                     Increase  Decrease
                       # <8>       812     664
# Fagopyrum_tataricum_H1 <1>       1176    208
# Fagopyrum_escelentum_H2<2>       1195    182
                       # <9>       5       76
# Polygunum_aviculare_H0 <3>       607     521
                       # <10>      22      80
       # Oxyria_digyna_H1<4>      596     367
                       # <11>     454     174
        # Rheum_nobile_H0<5>      1036    192
    # Rheum_tangaticum_H0<6>      1557    252

#--------------------------------------------------
cd /lustre04/scratch/celphin/Oxyria/CAFE/Polygonaceae_data/output/

grep "Lambda" ./*/*_results.txt
# ./error_model/Base_results.txt:Lambda: 0.002383227827052
# ./single_lambda/Base_results.txt:Lambda: 0.0023831902736841
# ./single_lambda_noe/Base_results.txt:Lambda: 0.0027165417373796
# ./three_lambda_error/Gamma_results.txt:Lambda: 0.0024360829102418
# ./three_lambda/Gamma_results.txt:Lambda: 0.002807174772633

grep "Epsilon" ./*/*_results.txt
# ./error_model/Base_results.txt:Epsilon: 0.0119675
# ./single_lambda/Base_results.txt:Epsilon: 0.0119675
# ./three_lambda_error/Gamma_results.txt:Epsilon: 0.0119675

#------------------------------------
# run many times to check lambdas 
cd /lustre04/scratch/celphin/Oxyria/CAFE/Polygonaceae_data/output
module load StdEnv/2020 intel/2020.1.217 cafe5/5.1.0

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

# > mean(epsilons[,1])
# [1] 0.01196818
# > sd(epsilons[,1])
# [1] 1.084367e-05

# > mean(lambdas[,1])
# [1] 0.002383249
# > sd(lambdas[,1])
# [1] 2.492556e-07

#########################################
# Plotting CAFE results - does not work yet
# https://github.com/LKremer/CAFE_fig
# https://github.com/etetoolkit/ete/issues/354

cd /lustre04/scratch/celphin/Oxyria/CAFE/

module load StdEnv/2023 python/3.11.5 scipy-stack/2023b qt/6.5.3 
pip3 install 'ete3==3.0.0b35'

wget https://raw.githubusercontent.com/LKremer/CAFE_fig/refs/heads/master/CAFE_fig.py

python3 /lustre04/scratch/celphin/Oxyria/CAFE/CAFE_fig.py \
/lustre04/scratch/celphin/Oxyria/CAFE/Polygonaceae_data/output/error_model/Base_report.cafe \
-pb 0.05 -pf 0.05 --dump test/ -g .pdf --count_all_expansions

#------------------------------
# Try another plot type
# https://github.com/moshi4/CafePlotter

cd /lustre04/scratch/celphin/Oxyria/CAFE/

module load StdEnv/2023 python/3.11.5 scipy-stack/2023b
pip install cafeplotter

cafeplotter -i /lustre04/scratch/celphin/Oxyria/CAFE/Polygonaceae_data/output/error_model \
-o /lustre04/scratch/celphin/Oxyria/CAFE/Polygonaceae_data/output/error_model/Cafe_plots

#####################################################
# Now we can investigate which orthogroups are expanded/contracted
# See next notes