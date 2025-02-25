##########################
# Running CAFE - gene enrichment analyses on all genomes on beluga 
# July 2024
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

cd /home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/orthofinder/Results_Jul03

cp Species_Tree/SpeciesTree_rooted.txt /lustre04/scratch/celphin/Oxyria/CAFE/Total_data/input_data
cp Phylogenetic_Hierarchical_Orthogroups/N0.tsv /lustre04/scratch/celphin/Oxyria/CAFE/Total_data/input_data

# in R
cd /lustre04/scratch/celphin/Oxyria/CAFE/Total_data/input_data

module load StdEnv/2020 r/4.2.2 
R

#install.packages("ape")

library(ape)
library(data.table)

# make orthofinder's species tree ultrametric. It should be already binary and rooted as required by cafe.
setwd("/home/celphin/scratch/Oxyria/CAFE/Total_data/input_data")
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
cafe5 -i /lustre04/scratch/celphin/Oxyria/CAFE/Total_data/input_data/hog_gene_counts.tsv \
-t /lustre04/scratch/celphin/Oxyria/CAFE/Total_data/input_data/SpeciesTree_rooted_ultra.txt \
-o /lustre04/scratch/celphin/Oxyria/CAFE/Total_data/output/error_model -e --cores 30

# Polygonaceae
# Filtering families not present at the root from: 42683 to 22790

# Total
# Filtering families not present at the root from: 29835 to 17839

# Completed 54 iterations
# Time: 0H 0M 43S
# Best matches are: 0.55721745244227,0.047621493185943
# Final -lnL: 478490.60580415

# 112 values were attempted (0% rejected)

# Inferring processes for Base model
# Score (-lnL): 478490.60580415
# Maximum possible lambda for this topology: 1.8875541699733
# Computing pvalues...
# done!

# Starting reconstruction processes for Base model
# Done!


#----------------------------------
cd /lustre04/scratch/celphin/Oxyria/CAFE/Total_data/output
more error_model/Base_error_model.txt #view results

# maxcnt: 170
# cntdiff: -1 0 1
# 0 0 0.96199 0.0380103
# 1 0.0380103 0.923979 0.0380103

# maxcnt: 170
# cntdiff: -1 0 1
# 0 0 0.952379 0.0476215
# 1 0.0476215 0.904757 0.0476215

##########################################
# Run again with error model ( note no space after -e)

cd /lustre04/scratch/celphin/Oxyria/CAFE
#no error model
cafe5 -i /lustre04/scratch/celphin/Oxyria/CAFE/Total_data/input_data/hog_gene_counts.tsv \
-t /lustre04/scratch/celphin/Oxyria/CAFE/Total_data/input_data/SpeciesTree_rooted_ultra.txt \
-o /lustre04/scratch/celphin/Oxyria/CAFE/Total_data/output/single_lambda_noe --cores 30 

#No among family rate variation, try with an without error model
cafe5 -i /lustre04/scratch/celphin/Oxyria/CAFE/Total_data/input_data/hog_gene_counts.tsv \
-t /lustre04/scratch/celphin/Oxyria/CAFE/Total_data/input_data/SpeciesTree_rooted_ultra.txt \
-o /lustre04/scratch/celphin/Oxyria/CAFE/Total_data/output/single_lambda \
-e/lustre04/scratch/celphin/Oxyria/CAFE/Total_data/output/error_model/Base_error_model.txt --cores 30 

# #Yes among family rate variation, estimate lambda and alpha and three discrete gamma rate categories
cafe5 -i /lustre04/scratch/celphin/Oxyria/CAFE/Total_data/input_data/hog_gene_counts.tsv \
-t /lustre04/scratch/celphin/Oxyria/CAFE/Total_data/input_data/SpeciesTree_rooted_ultra.txt \
-k 3 -o /lustre04/scratch/celphin/Oxyria/CAFE/Total_data/output/three_lambda --cores 30 

# #with error model
cafe5 -i /lustre04/scratch/celphin/Oxyria/CAFE/Total_data/input_data/hog_gene_counts.tsv \
-t /lustre04/scratch/celphin/Oxyria/CAFE/Total_data/input_data/SpeciesTree_rooted_ultra.txt \
-k 3 -o /lustre04/scratch/celphin/Oxyria/CAFE/Total_data/output/three_lambda_error \
-e/lustre04/scratch/celphin/Oxyria/CAFE/Total_data/output/error_model/Base_error_model.txt --cores 30

#Now analyze the large families using the previous values of lambda (skipped)
# cafe5 -i large_/lustre04/scratch/celphin/Oxyria/CAFE/Total_data/input_data/hog_gene_counts.tsv \
# -t /lustre04/scratch/celphin/Oxyria/CAFE/Total_data/input_data/SpeciesTree_rooted_ultra.txt \
# -l 0.001319516 -o /lustre04/scratch/celphin/Oxyria/CAFE/Total_data/output/large_results --cores 30 >


##########################################
# Explore data
cd /lustre04/scratch/celphin/Oxyria/CAFE/Total_data/output/

more Base_asr.tre 
#gives reconstructed numbers of each gene family at each node

# BEGIN TREES;
# TREE N0.HOG0000000 = (((Polygunum_aviculare_H0<1>_0:0.529786,(Oxyria_digyna_H1<2>_0:0.412583,(Rheum_t
# angaticum_H0<3>_0:0.156531,Rheum_nobile_H0<4>_0:0.156531)<25>_0:0.256052)<24>_0:0.117203)<23>_0:0.116
# 818,(Fagopyrum_tataricum_H1<5>_1:0.167024,Fagopyrum_escelentum_H2<6>_0:0.167024)<26>_1:0.479579)<22>_
# 1:0.353396,((((Capsella_rubella<7>_0:0.170537,(Arabidopsis_lyrata<8>_0:0.113768,Arabidopsis_thaliana<9>
# _0:0.113768)<31>_0:0.0567696)<30>_0:0.0815547,(Draba_nivalis<10>_0:0.184379,Arabis_alpina<11>_0:0.184379)<3
# 2>_0:0.0677131)<29>_0:0.0348375,(Brassica_oleracea<12>_0:0.236399,Thlaspi_arvense<13>_0:0.236399)<33>_0
# :0.0505304)<28>_0:0.483808,((Rosa_rugosa<14>_0:0.207543,(Argentina_anserina<15>_0:0.16543,Fragaria_vesc
# a<16>_0:0.16543)<36>_0:0.0421135)<35>_0:0.238568,(Dryas_octopetala<17>_4:0.366067,(Prunus_persica<18>_0
# :0.270153,(Malus_sylvestris<19>_0:0.050394,Pyrus_bretschneideri<20>_0:0.050394)<39>_0:0.219759)<38>_0:0
# .0959135)<37>_1:0.0800447)<34>_1:0.324626)<27>_1:0.229262)<21>_1;

  # TREE N0.HOG0000001 = (((Polygunum_aviculare_H0<1>_1:0.529786,(Oxyria_digyna_H1<2>_0:0.412583,(Rheum_t
# angaticum_H0<3>_0:0.156531,Rheum_nobile_H0<4>_0:0.156531)<25>_0:0.256052)<24>_0:0.117203)<23>_1:0.11681
# 8,(Fagopyrum_tataricum_H1<5>_1:0.167024,Fagopyrum_escelentum_H2<6>_0:0.167024)<26>_1:0.479579)<22>_1:0.
# 353396,((((Capsella_rubella<7>_0:0.170537,(Arabidopsis_lyrata<8>_0:0.113768,Arabidopsis_thaliana<9>_0:0
# .113768)<31>_0:0.0567696)<30>_0:0.0815547,(Draba_nivalis<10>_0:0.184379,Arabis_alpina<11>_0:0.184379)<3
# 2>_0:0.0677131)<29>_0:0.0348375,(Brassica_oleracea<12>_0:0.236399,Thlaspi_arvense<13>_0:0.236399)<33>_0
# :0.0505304)<28>_0:0.483808,((Rosa_rugosa<14>_0:0.207543,(Argentina_anserina<15>_0:0.16543,Fragaria_vesc
# a<16>_0:0.16543)<36>_0:0.0421135)<35>_0:0.238568,(Dryas_octopetala<17>_2:0.366067,(Prunus_persica<18>_0
# :0.270153,(Malus_sylvestris<19>_0:0.050394,Pyrus_bretschneideri<20>_0:0.050394)<39>_0:0.219759)<38>_0:0
# .0959135)<37>_1:0.0800447)<34>_1:0.324626)<27>_1:0.229262)<21>_1;

  # TREE N0.HOG0000002 = (((Polygunum_aviculare_H0<1>*_70:0.529786,(Oxyria_digyna_H1<2>*_94:0.412583,(Rhe
# um_tangaticum_H0<3>_13:0.156531,Rheum_nobile_H0<4>*_4:0.156531)<25>*_14:0.256052)<24>_48:0.117203)<23>_
# 50:0.116818,(Fagopyrum_tataricum_H1<5>_90:0.167024,Fagopyrum_escelentum_H2<6>*_97:0.167024)<26>*_87:0.4
# 79579)<22>*_49:0.353396,((((Capsella_rubella<7>_0:0.170537,(Arabidopsis_lyrata<8>_0:0.113768,Arabidopsi
# s_thaliana<9>_0:0.113768)<31>_0:0.0567696)<30>_0:0.0815547,(Draba_nivalis<10>_0:0.184379,Arabis_alpina<
# 11>_0:0.184379)<32>_0:0.0677131)<29>_0:0.0348375,(Brassica_oleracea<12>_0:0.236399,Thlaspi_arvense<13>_
# 0:0.236399)<33>_0:0.0505304)<28>*_0:0.483808,((Rosa_rugosa<14>_0:0.207543,(Argentina_anserina<15>_0:0.1
# 6543,Fragaria_vesca<16>_1:0.16543)<36>_1:0.0421135)<35>_1:0.238568,(Dryas_octopetala<17>_1:0.366067,(Pr
# unus_persica<18>_1:0.270153,(Malus_sylvestris<19>_0:0.050394,Pyrus_bretschneideri<20>_0:0.050394)<39>_0
# :0.219759)<38>_1:0.0959135)<37>_1:0.0800447)<34>*_1:0.324626)<27>*_4:0.229262)<21>_20;
  

#--------------------------------
more Base_branch_probabilities.tab 
#gives calculated likelihood of the gene family size at each node for each gene family

# FamilyID        Polygunum_aviculare_H0<1>       Oxyria_digyna_H1<2>     Rheum_tangaticum_H0<3>  Rheum_nobile_H0<4>     Fagopyrum_tataricum_H1<5>       Fagopyrum_escelentum_H2<6>      Capsella_rubella<7>   Arabidopsis_lyrata<8>    Arabidopsis_thaliana<9> Draba_nivalis<10>       Arabis_alpina<11>       Brassica_oleracea<12>  Thlaspi_arvense<13>     Rosa_rugosa<14> Argentina_anserina<15>  Fragaria_vesca<16>    Dryas_octopetala<17>     Prunus_persica<18>      Malus_sylvestris<19>    Pyrus_bretschneideri<20>      <21>     <22>    <23>    <24>    <25>    <26>    <27>    <28>    <29>    <30>    <31>    <32>    <33>  <34>     <35>    <36>    <37>    <38>    <39>
# N0.HOG0000002   0.00062883      9.56077e-15     0.617332        4.06033e-09     0.417097        0.0162343      0.5     0.5     0.5     0.5     0.5     0.5     0.5     0.144429        0.11921 0.580655      0.655048 0.622223        0.5     0.5     N/A     2.37248e-13     0.614107        0.488347        2.59671e-21    5.16218e-10     1.7386e-12      0.0028256       0.5     0.5     0.5     0.5     0.5     0.0128645      0.610235        0.522606        0.541764        0.549011        0.15131
# N0.HOG0000020   0.701746        0.66928 0.5     0.5     0.581509        0.120452        0.5     0.5   0.5      0.0460947       0.130843        4.79951e-07     0.736814        0.5     0.5     0.5     0.0822671      0.179059        0.526738        0.0399238       N/A     0.650859        0.55887 0.559331      0.171674 0.688484        0.60676 0.689577        0.51842 0.0629132       0.5     0.535342        0.000367089    0.641235        0.161924        0.5     0.541764        0.549011        0.602844
# N0.HOG0000024   1.26285e-06     0.245204        0.5     0.5     0.781494        0.235152        0.5   0.5      0.5     0.58867 0.130843        0.0573251       0.160823        0.5     0.5     0.5     0.5   0.5      0.5     0.5     N/A     0.161876        0.648006        0.00691933      0.171674        0.0327642      0.272306        0.689577        0.51842 0.0629132       0.5     0.535342        0.526738      0.206006 0.5     0.5     0.5     0.5     0.5
# N0.HOG0000025   0.791173        0.0737203       0.5     0.5     0.581509        0.581509        0.127783       0.557481        0.0853448       0.000409933     0.00112825      2.93337e-20     0.00237108    0.0511548        0.580655        0.11921 0.5     0.5     0.5     0.5     N/A     0.737336        0.607556       0.608321        0.0132481       0.425267        0.680653        0.207735        0.552288      0.611655 0.00156104      0.595351        0.574026        0.341858        0.610235        0.522606       0.0621914       0.5     0.5
# N0.HOG0000031   8.11584e-07     1.02529e-06     0.5     0.5     0.5     0.5     0.122308        0.0853448      0.557481        0.5     0.5     0.000370428     0.011268        5.24759e-10     3.12709e-09   0.814369 1.35439e-19     0.0147136       0.526738        0.0399238       N/A     0.861062        0.0633784      0.772701        5.13436e-08     9.57687e-05     0.81906 0.021287        0.0537498       0.542253       0.529802        0.0526896       0.551351        0.128139        0.373868        0.242216       0.722888        3.18438e-07     0.263891
# N0.HOG0000035   0.222697        3.55304e-15     0.5     0.5     0.581509        0.120452        0.5   0.5      0.5     0.5     0.5     0.5     0.5     0.5     0.5     0.5     0.0139356       0.5     0.5   0.5      N/A     0.737336        0.607556        0.608321        0.0132481       0.425267        0.272306       0.273123        0.5     0.5     0.5     0.5     0.5     0.641235        0.161924        0.5   0.541764 0.072884        0.5

#------------------------------
more Base_clade_results.txt 
#how many families are expanded or contracted at each node

# #Taxon_ID       Increase        Decrease
# <22>    704     714
# <25>    1108    972
# Rheum_tangaticum_H0<3>  2125    2687
# <23>    81      1584
# <26>    2840    2452
# Polygunum_aviculare_H0<1>       1597    4230
# <34>    677     3003
# <35>    480     2072
# Rosa_rugosa<14> 1612    727
# <27>    121     1630
# <28>    2011    3274
# <29>    41      1044
# <24>    169     467
# Oxyria_digyna_H1<2>     1480    3523
# <31>    141     234
# Fagopyrum_tataricum_H1<5>       2207    1659
# Fagopyrum_escelentum_H2<6>      2726    1609
# Rheum_nobile_H0<4>      1923    1280

# Arabidopsis_lyrata<8>   1337    196
# Arabidopsis_thaliana<9> 180     6341
# <32>    161     1042
# Argentina_anserina<15>  358     1436
# Draba_nivalis<10>       1626    1017
# <33>    497     954
# <36>    12      748
# <30>    200     976
# Capsella_rubella<7>     727     667
# Arabis_alpina<11>       522     4437
# Thlaspi_arvense<13>     987     2938
# Brassica_oleracea<12>   6276    1068

# Prunus_persica<18>      473     1536
# <37>    205     528
# Dryas_octopetala<17>    1257    1961
# <38>    365     955
# Fragaria_vesca<16>      632     1510
# <39>    7574    377
# Malus_sylvestris<19>    1548    498
# Pyrus_bretschneideri<20>        1573    702

#--------------------------------------------------
cd /lustre04/scratch/celphin/Oxyria/CAFE/Total_data/output/

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

printf '%s\n' {1..10} | parallel cafe5 -i /lustre04/scratch/celphin/Oxyria/CAFE/Total_data/input_data/hog_gene_counts.tsv \
-t /lustre04/scratch/celphin/Oxyria/CAFE/Total_data/input_data/SpeciesTree_rooted_ultra.txt -e \
-o /lustre04/scratch/celphin/Oxyria/CAFE/Total_data/output/error_model_{1} >> \
error_model_100runs.log 2>> error_model_100runs.errorlog

#MORE!
printf '%s\n' {11..20} | parallel cafe5 -i /lustre04/scratch/celphin/Oxyria/CAFE/Total_data/input_data/hog_gene_counts.tsv \
-t /lustre04/scratch/celphin/Oxyria/CAFE/Total_data/input_data/SpeciesTree_rooted_ultra.txt -e \
-o /lustre04/scratch/celphin/Oxyria/CAFE/Total_data/output/error_model_{1} >> \
error_model_100runs.log 2>> error_model_100runs.errorlog

#MORE
printf '%s\n' {21..50} | parallel cafe5 -i /lustre04/scratch/celphin/Oxyria/CAFE/Total_data/input_data/hog_gene_counts.tsv \
-t /lustre04/scratch/celphin/Oxyria/CAFE/Total_data/input_data/SpeciesTree_rooted_ultra.txt -e \
-o /lustre04/scratch/celphin/Oxyria/CAFE/Total_data/output/error_model_{1} >> \
error_model_100runs.log 2>> error_model_100runs.errorlog

#Now analyze the large families using the previous values of lambda
# printf '%s\n' {1..50} | parallel cafe5 -i large_/lustre04/scratch/celphin/Oxyria/CAFE/Total_data/input_data/hog_gene_counts.tsv \
# -t /lustre04/scratch/celphin/Oxyria/CAFE/Total_data/input_data/SpeciesTree_rooted_ultra.txt -l 0.001488216 -o large_results_{1} >> \
# large_results_100runs.log


#-----------------------------------------
#Summarise the results of all runs

grep "Lambda" ./error_model_*/Base_results.txt | sed 's/^.* //g' > lambdas.txt

# 0.55722341575828
# 0.5572225405751
# 0.55722276889862
# 0.55722504362198
# 0.55722375130075

grep "Epsilon" ./error_model_*/Base_results.txt | sed 's/^.* //g' > epsilons.txt

# 0.0476209
# 0.0476235
# 0.0476211

#-------------------------q---
module load StdEnv/2020 r/4.2.2 

R

epsilons <- read.table("epsilons.txt")
lambdas <- read.table("lambdas.txt")

mean(epsilons[,1])
# 0.04762016
sd(epsilons[,1])
# 8.623014e-06

mean(lambdas[,1])
# 0.5572264

sd(lambdas[,1])
# 2.285778e-05

q()
n

#-----------------
# Polygonaceae
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
# See notes file specifc to this!