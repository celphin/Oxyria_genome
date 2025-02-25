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
mkdir Total_genomes; cd Total_genomes
mkdir input_data
mkdir output
mkdir analysis
mkdir expanded
mkdir contracted

# copy over important files
cd /home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/orthofinder/Results_Aug18
cp Species_Tree/SpeciesTree_rooted.txt /lustre04/scratch/celphin/Oxyria/CAFE/Total_genomes/input_data
cp Phylogenetic_Hierarchical_Orthogroups/N0.tsv /lustre04/scratch/celphin/Oxyria/CAFE/Total_genomes/input_data

#-------------------------------
# in R
cd /lustre04/scratch/celphin/Oxyria/CAFE/Total_genomes/input_data

module load StdEnv/2020 r/4.2.2 
R

#install.packages("ape")
# https://rdrr.io/cran/ape/man/chronos.html
# https://phylobotanist.blogspot.com/2016/12/how-to-set-multiple-calibration-points.html

library(ape)
library(data.table)

# make orthofinder's species tree ultrametric. It should be already binary and rooted as required by cafe.
setwd("/home/celphin/scratch/Oxyria/CAFE/Total_genomes/input_data")
tre <- read.tree('./SpeciesTree_rooted.txt')

tre <- drop.tip(tre, "Arabidopsis_thaliana") ## taxon can be a single taxon or a vector of taxon names

stopifnot(is.binary(tre))
stopifnot(is.rooted(tre))

# Brassica and Arabidopsis diverged about 43MYA
# https://pubmed.ncbi.nlm.nih.gov/20921408/
# Rosa diverged from Prunus 74mya - but crown is 100-110mya
# Fagopyrum and Rheum diverged 48mya

# set calibrations
cal <- makeChronosCalib(tre, node="root", age.min=90, age.max=130) 

if(is.ultrametric(tre)) {
    utre <- tre
} else{
    utre <- chronos(tre, lambda = 1, model = "discrete", calibration = cal, control = chronos.control())
}

# Setting initial dates...
# Fitting in progress... get a first set of estimates
         # (Penalised) log-lik = -8.095333
# Optimising rates... frequencies... dates... -8.095333
# Optimising rates... frequencies... dates... -7.741562
# Optimising rates... frequencies... dates... -7.739322
# Optimising rates... frequencies... dates... -7.738777
# Optimising rates... frequencies... dates... -7.738646
# Optimising rates... frequencies... dates... -7.738615
# Optimising rates... frequencies... dates... -7.73861

# log-Lik = -7.73861
# PHIIC = 91.48

write.tree(utre, './SpeciesTree_rooted_ultra2.txt')

#----------------------------------------
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

# Exclude HOGs not present in at least 6 species
keep <- hog[, .N, HOG][N > 5]$HOG
hog <- hog[HOG %in% keep]

counts <- dcast(hog, HOG ~ species, value.var='n', fill=0)
counts[, Desc := 'n/a']
setcolorder(counts, 'Desc')
fwrite(counts, 'hog_gene_counts_5spp.tsv', sep='\t')

#----------------------------------
# Run again for HOGs present in all spp

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

# Exclude HOGs not present in all 20 species
keep <- hog[, .N, HOG][N > 19]$HOG
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
wc -l hog_gene_counts_5spp.tsv
# 18663

wc -l hog_gene_counts.tsv
# 5234 

#---------------------------------
# The tree file should contain a binary, rooted, ultrametric, tree in Newick format. 
# Typically one obtains this tree using one of several molecular dating methods. 
# If you are unsure if your tree is binary, rooted, or ultrametric CAFE will 
# report this when you try to use it for an analysis. Alternatively, you can use the R package,
# Ape with its included functions: is.ultrametric, is.rooted, and is.binary.
# view here: http://etetoolkit.org/treeview/

more SpeciesTree_rooted_ultra.txt

(((Fagopyrum_escelentum_H2:9.5333339,Fagopyrum_tataricum_H1:9.5333339)0.973441:32.93660324,(Polygu
num_aviculare_H0:34.3137747,((Rheum_nobile_H0:9.537871213,Rheum_tangaticum_H0:9.537871213)0.880931
:16.94451109,Oxyria_digyna_H1:26.48238231)0.593853:7.831392398)0.56789:8.156162436)0.966876:32.530
06286,(((Rosa_rugosa:16.062299,(Fragaria_vesca:12.96039051,Argentina_anserina:12.96039051)0.504625
:3.101908489)0.958221:18.18780382,(Dryas_octopetala:27.97813247,(Prunus_persica:20.39123076,(Malus
_sylvestris:3.793886609,Pyrus_bretschneideri:3.793886609)0.982393:16.59734415)0.714115:7.586901711
)0.626082:6.271970344)0.961504:25.72918221,(Cochlearia_groenlandica:26.85314704,(((Arabidopsis_lyr
ata:12.0916722,Capsella_rubella:12.0916722)0.797672:8.199098375,(Arabis_alpina:15.30631759,Draba_n
ivalis:15.30631759)0.719188:4.984452987)0.183229:2.367998821,(Thlaspi_arvense:19.04005471,Brassica
_oleracea:19.04005471)0.281409:3.61871468)0.289466:4.194377652)0.965085:33.12613798)0.966876:15.02
071498);

 more SpeciesTree_rooted_ultra2.txt
(((Fagopyrum_escelentum_H2:14.21232417,Fagopyrum_tataricum_H1:14.21232417)0.973441:49.0218
6106,(Polygunum_aviculare_H0:51.10984779,((Rheum_nobile_H0:14.22049366,Rheum_tangaticum_H0
:14.22049366)0.880931:25.24358696,Oxyria_digyna_H1:39.46408062)0.593853:11.64576717)0.5678
9:12.12433743)0.966876:48.21786211,(((Rosa_rugosa:23.93803053,(Fragaria_vesca:19.31643358,
Argentina_anserina:19.31643358)0.504625:4.621596953)0.958221:27.09223951,(Dryas_octopetala
:41.6975285,(Prunus_persica:30.40795111,(Malus_sylvestris:5.659707504,Pyrus_bretschneideri
:5.659707504)0.982393:24.74824361)0.714115:11.28957739)0.626082:9.332741539)0.961504:38.19
880119,(Cochlearia_groenlandica:40.0145861,(((Arabidopsis_lyrata:18.03179591,Capsella_rube
lla:18.03179591)0.797672:12.22233548,(Arabis_alpina:22.82923984,Draba_nivalis:22.82923984)
0.719188:7.424891548)0.183229:3.519275738,(Thlaspi_arvense:28.3699873,Brassica_oleracea:28
.3699873)0.281409:5.403419823)0.289466:6.241178972)0.965085:49.21448513)0.966876:22.222976
1);

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

# all models with 5 spp stored in: /lustre04/scratch/celphin/Oxyria/CAFE/Total_genomes/5spp_output/
# run main model with all spp here: /lustre04/scratch/celphin/Oxyria/CAFE/Total_genomes/output/

tmux new-session -s CAFE
tmux attach-session -t CAFE

salloc -c30 --time 3:00:00 --mem 191000M --account def-cronk

module load StdEnv/2020 intel/2020.1.217 cafe5/5.1.0

cd /lustre04/scratch/celphin/Oxyria/CAFE
cafe5 -i /lustre04/scratch/celphin/Oxyria/CAFE/Total_genomes/input_data/hog_gene_counts.tsv \
-t /lustre04/scratch/celphin/Oxyria/CAFE/Total_genomes/input_data/SpeciesTree_rooted_ultra2.txt \
-o /lustre04/scratch/celphin/Oxyria/CAFE/Total_genomes/output/error_model -e --cores 30

#-------------------------------
# Total with 5 spp
# Filtering families not present at the root from: 18662 to 15157

# Completed 52 iterations
# Best matches are: 0.0064360341120366,0.026301913603283
# Final -lnL: 405112.40881178
# 106 values were attempted (0% rejected)
# Score (-lnL): 405112.40881178
# Maximum possible lambda for this topology: 0.029142815348729

#------------------------------
# total with all spp
# all families present at root

# Completed 52 iterations
# Best matches are: 0.0047212248287059,0.00023411162372211
# Final -lnL: 125973.64463733
# 117 values were attempted (0% rejected)
# Score (-lnL): 125973.64463733
# Maximum possible lambda for this topology: 0.029142815348729


#----------------------------------
cd /lustre04/scratch/celphin/Oxyria/CAFE/Total_genomes/output
more error_model/Base_error_model.txt #view results

# 5 spp
# maxcnt: 170
# cntdiff: -1 0 1
# 0 0 0.973698 0.0263019
# 1 0.0263019 0.947396 0.0263019

# all spp
# maxcnt: 170
# cntdiff: -1 0 1
# 0 0 0.999766 0.000234112
# 1 0.000234112 0.999532 0.000234112


##########################################
# Run again with error model ( note no space after -e)

cd /lustre04/scratch/celphin/Oxyria/CAFE
#no error model
cafe5 -i /lustre04/scratch/celphin/Oxyria/CAFE/Total_genomes/input_data/hog_gene_counts.tsv \
-t /lustre04/scratch/celphin/Oxyria/CAFE/Total_genomes/input_data/SpeciesTree_rooted_ultra.txt \
-o /lustre04/scratch/celphin/Oxyria/CAFE/Total_genomes/output/single_lambda_noe --cores 30 

# 5 spp
# Filtering families not present at the root from: 18662 to 15157
# Best match is: 0.0072746135262669
# Final -lnL: 408404.73068398
# Score (-lnL): 408404.73068398
# Maximum possible lambda for this topology: 0.029142815348729

# all spp
# Best match is: 0.0047278163131363
# Final -lnL: 125973.69625303
# Score (-lnL): 125973.69625303
# Maximum possible lambda for this topology: 0.029142815348729

#-----------------------------
#No among family rate variation, try with an without error model
cafe5 -i /lustre04/scratch/celphin/Oxyria/CAFE/Total_genomes/input_data/hog_gene_counts.tsv \
-t /lustre04/scratch/celphin/Oxyria/CAFE/Total_genomes/input_data/SpeciesTree_rooted_ultra.txt \
-o /lustre04/scratch/celphin/Oxyria/CAFE/Total_genomes/output/single_lambda \
-e/lustre04/scratch/celphin/Oxyria/CAFE/Total_genomes/output/error_model/Base_error_model.txt --cores 30 

#------------------------------
# #Yes among family rate variation, estimate lambda and alpha and three discrete gamma rate categories
cafe5 -i /lustre04/scratch/celphin/Oxyria/CAFE/Total_genomes/input_data/hog_gene_counts.tsv \
-t /lustre04/scratch/celphin/Oxyria/CAFE/Total_genomes/input_data/SpeciesTree_rooted_ultra.txt \
-k 3 -o /lustre04/scratch/celphin/Oxyria/CAFE/Total_genomes/output/three_lambda --cores 30 

#------------------------------------
# #with error model
cafe5 -i /lustre04/scratch/celphin/Oxyria/CAFE/Total_genomes/input_data/hog_gene_counts.tsv \
-t /lustre04/scratch/celphin/Oxyria/CAFE/Total_genomes/input_data/SpeciesTree_rooted_ultra.txt \
-k 3 -o /lustre04/scratch/celphin/Oxyria/CAFE/Total_genomes/output/three_lambda_error \
-e/lustre04/scratch/celphin/Oxyria/CAFE/Total_genomes/output/error_model/Base_error_model.txt --cores 30


##########################################
# Explore data
cd /lustre04/scratch/celphin/Oxyria/CAFE/Total_genomes/output/error_model/

more Base_asr.tre 
#gives reconstructed numbers of each gene family at each node

#--------------------------------
more Base_branch_probabilities.tab 
#gives calculated likelihood of the gene family size at each node for each gene family

#------------------------------
more Base_clade_results.txt

# Fagopyrum_escelentum_H2    <1>       737     191
# Fagopyrum_tataricum_H1     <2>       725     176
# Polygunum_aviculare_H0     <3>       487     564
# Rheum_nobile_H0            <4>       530     184
# Rheum_tangaticum_H0        <5>       756     396
# Oxyria_digyna_H1           <6>       397     446

# Rosa_rugosa                <7>       467     76
# Fragaria_vesca             <8>       215     148
# Argentina_anserina         <9>       108     283
# Dryas_octopetala           <10>      280     344
# Prunus_persica             <11>      123     294
# Malus_sylvestris           <12>      486     57
# Pyrus_bretschneideri       <13>      586     39

# Cochlearia_groenlandica    <14>      1124    589
# Arabidopsis_lyrata         <15>      402     86
# Capsella_rubella           <16>      161     164
# Arabis_alpina              <17>      231     674
# Draba_nivalis              <18>      566     227
# Thlaspi_arvense            <19>      320     640
# Brassica_oleracea          <20>      2624    85

# <22>    693     261
# <23>    778     192
# <32>    109     13
# <28>    161     834
# <27>    11      98
# <24>    59      61
# <25>    45      85
# <26>    345     170
# <30>    7       118
# <29>    121     210
# <33>    3199    35
# <35>    31      7
# <31>    46      13
# <36>    13      76
# <37>    137     125
# <34>    889     207
# <38>    60      246
# <39>    198     18

#----------------------------
more /lustre04/scratch/celphin/Oxyria/CAFE/Total_genomes/5spp_output/error_model/Base_clade_results.txt 

#how many families are expanded or contracted at each node
# List with HOGs in at least 6 spp 
#                        Taxon_ID       Increase  Decrease
# Fagopyrum_escelentum_H2        <1>      2298    1091
# Fagopyrum_tataricum_H1         <2>      1978    925
# Polygunum_aviculare_H0         <3>      1374    2917
# Rheum_nobile_H0                <4>      1598    818
# Rheum_tangaticum_H0            <5>      1852    2240
# Oxyria_digyna_H1               <6>      1230    2557

# Rosa_rugosa                    <7>      1532    498
# Fragaria_vesca                 <8>       613     1390
# Argentina_anserina             <9>       352     1313
# Dryas_octopetala               <10>      1059    1497
# Prunus_persica                 <11>      440     1105
# Malus_sylvestris               <12>      1499    322
# Pyrus_bretschneideri           <13>      1582    461

# Cochlearia_groenlandica        <14>     2662    2820
# Arabidopsis_lyrata             <15>     1069    233
# Capsella_rubella               <16>     430     529
# Arabis_alpina                  <17>     461     4075
# Draba_nivalis                  <18>     1487    855
# Thlaspi_arvense                <19>     883     2793
# Brassica_oleracea              <20>     6152    627

# <22>    2064    802
# <23>    2385    1169
# <24>    70      683
# <25>    168     338
# <26>    933     657
# <27>    54      392
# <28>    679     1968
# <29>    422     1048
# <30>    34      480
# <31>    204     165
# <32>    345     307
# <33>    7542    276
# <34>    2023    2086
# <35>    145     58
# <36>    39      424
# <37>    405     479
# <38>    136     897
# <39>    420     264

#--------------------------------------------------
cd /lustre04/scratch/celphin/Oxyria/CAFE/Total_genomes/output/error_model/

grep "Lambda" ./*/*_results.txt
# ./error_model/Base_results.txt:Lambda: 0.0047212248287059
# ./single_lambda/Base_results.txt:Lambda: 0.0047212679379012
# ./single_lambda_noe/Base_results.txt:Lambda: 0.0047278163131363
# ./three_lambda_error/Gamma_results.txt:Lambda: 0.004950998138347
# ./three_lambda/Gamma_results.txt:Lambda: 0.0049572400620943

grep "Epsilon" ./*/*_results.txt
# ./error_model/Base_results.txt:Epsilon: 0.000234112
# ./single_lambda/Base_results.txt:Epsilon: 0.000234112
# ./three_lambda_error/Gamma_results.txt:Epsilon: 0.000234112

#------------------------------------
# run many times to check lambdas 

printf '%s\n' {1..50} | parallel cafe5 -i /lustre04/scratch/celphin/Oxyria/CAFE/Total_genomes/input_data/hog_gene_counts.tsv \
-t /lustre04/scratch/celphin/Oxyria/CAFE/Total_genomes/input_data/SpeciesTree_rooted_ultra.txt -e \
-o /lustre04/scratch/celphin/Oxyria/CAFE/Total_genomes/output/error_model_{1} >> \
error_model_100runs.log 2>> error_model_100runs.errorlog

#-----------------------------------------
#Summarise the results of all runs
cd /lustre04/scratch/celphin/Oxyria/CAFE/Total_genomes/output/

grep "Lambda" ./error_model_*/Base_results.txt | sed 's/^.* //g' > lambdas.txt

# 0.004719143550386
# 0.0047212040192077
# 0.0047213948768954
# 0.0047212996652848
# 0.0047212430903214
# 0.0047213965327005

# more than 5 spp
# 0.006438551061351
# 0.0064385791880088
# 0.0064385810414676
# 0.0064385298168018

grep "Epsilon" ./error_model_*/Base_results.txt | sed 's/^.* //g' > epsilons.txt

# 0.000254736
# 0.000234536
# 0.000234125
# 0.000233262
# 0.000232297
# 0.000235579
# 0.000232998
# 0.000234114


# more than 5 spp
# 0.0264686
# 0.0264735
# 0.0264791
# 0.0264764
# 0.0264767


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

cafeplotter -i /lustre04/scratch/celphin/Oxyria/CAFE/Total_genomes/output/error_model \
-o /lustre04/scratch/celphin/Oxyria/CAFE/Total_genomes/Cafe_plots_error_model0_pdf  --format pdf


######################################
# Now we can investigate which orthogroups are expanded/contracted
# See notes file specifc to this!

