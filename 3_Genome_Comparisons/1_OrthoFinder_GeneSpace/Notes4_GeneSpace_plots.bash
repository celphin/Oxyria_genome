##########################
# Making GeneSpace plots of synteny using orthoFinder results in R
# June 2024
############################

# Polygonaceae
# restart in R

tmux new-session -s GeneSpace
tmux attach-session -t GeneSpace
tmux kill-session -t GeneSpace

# bash
salloc -c20 --time 2:40:00 --mem 191000M --account def-cronk

module load StdEnv/2020  python/3.10.2 scipy-stack/2021a diamond/2.1.7 fastme/2.1.6.2 gcc/9.3.0 mcl/14.137 java/13.0.2 r/4.2.2 glpk/5.0

# other version of orthofinder
export PATH=$PATH:/home/celphin/scratch/Oxyria/GeneSpace/OrthoFinder/bin
export PATH=$PATH:/home/celphin/scratch/Oxyria/GeneSpace/OrthoFinder_source
alias orthofinder='orthofinder.py'

#orthofinder -h
# works
cd /home/celphin/scratch/Oxyria/GeneSpace/

R

library(GENESPACE)

# https://rdrr.io/github/jtlovell/GENESPACE/man/init_genespace.html

###############################
# Total
gpar <- init_genespace(
  wd = "/home/celphin/scratch/Oxyria/GeneSpace/Total_genomes", 
  path2mcscanx = "/home/celphin/scratch/Oxyria/GeneSpace/MCScanX/",
  rawOrthofinderDir = "/home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/orthofinder/Results_Aug18")

     # Flagging chrs. w/ < 10 unique orthogroups
        # ...Arabidopsis_lyrata     :  694 genes on  376 small chrs.
        # ...Arabidopsis_thaliana   :    0 genes on    0 small chrs.
        # ...Arabis_alpina          : 2415 genes on 1734 small chrs. ***
        # ...Argentina_anserina     :   13 genes on    7 small chrs.
        # ...Brassica_oleracea      :    0 genes on    0 small chrs.
        # ...Capsella_rubella       :  217 genes on  120 small chrs.
        # ...Cochlearia_groenlandica:   72 genes on   24 small chrs.
        # ...Draba_nivalis          :  138 genes on   60 small chrs.
        # ...Dryas_octopetala       :   70 genes on   15 small chrs.
        # ...Fagopyrum_escelentum_H2:  162 genes on   63 small chrs.
        # ...Fagopyrum_tataricum_H1 : 2791 genes on  881 small chrs. ***
        # ...Fragaria_vesca         :    0 genes on    0 small chrs.
        # ...Malus_sylvestris       :    5 genes on    2 small chrs.
        # ...Oxyria_digyna_H1       : 1608 genes on  273 small chrs.
        # ...Polygunum_aviculare_H0 :  590 genes on   75 small chrs.
        # ...Prunus_persica         :   62 genes on   38 small chrs.
        # ...Pyrus_bretschneideri   :  520 genes on  175 small chrs.
        # ...Rheum_nobile_H0        :    0 genes on    0 small chrs.
        # ...Rheum_tangaticum_H0    :    0 genes on    0 small chrs.
        # ...Rosa_rugosa            :   70 genes on   42 small chrs.
        # ...Thlaspi_arvense        : 2146 genes on  652 small chrs. ***
        # NOTE! Genomes flagged *** have > 5% of genes on small chrs.
                # These are likely not great assemblies and should be
                # examined carefully

       ##############
        # Flagging over-dispered OGs
        # ...Arabidopsis_lyrata     :  2250 genes in  87 OGs hit > 8 unique places ***
        # ...Arabidopsis_thaliana   :   704 genes in  21 OGs hit > 8 unique places
        # ...Arabis_alpina          :  2143 genes in  70 OGs hit > 8 unique places ***
        # ...Argentina_anserina     :   390 genes in  23 OGs hit > 8 unique places
        # ...Brassica_oleracea      :  5989 genes in 235 OGs hit > 8 unique places ***
        # ...Capsella_rubella       :  1231 genes in  50 OGs hit > 8 unique places
        # ...Cochlearia_groenlandica:  2458 genes in  98 OGs hit > 8 unique places ***
        # ...Draba_nivalis          :  3906 genes in 122 OGs hit > 8 unique places ***
        # ...Dryas_octopetala       : 11208 genes in 352 OGs hit > 8 unique places ***
        # ...Fagopyrum_escelentum_H2: 13011 genes in 338 OGs hit > 8 unique places ***
        # ...Fagopyrum_tataricum_H1 :  4438 genes in 150 OGs hit > 8 unique places ***
        # ...Fragaria_vesca         :   911 genes in  43 OGs hit > 8 unique places
        # ...Malus_sylvestris       :  2739 genes in 100 OGs hit > 8 unique places ***
        # ...Oxyria_digyna_H1       :  7171 genes in 243 OGs hit > 8 unique places ***
        # ...Polygunum_aviculare_H0 :  2633 genes in  78 OGs hit > 8 unique places ***
        # ...Prunus_persica         :  1075 genes in  40 OGs hit > 8 unique places
        # ...Pyrus_bretschneideri   :  1339 genes in  58 OGs hit > 8 unique places
        # ...Rheum_nobile_H0        :  2769 genes in 161 OGs hit > 8 unique places ***
        # ...Rheum_tangaticum_H0    :  1801 genes in 109 OGs hit > 8 unique places ***
        # ...Rosa_rugosa            :  3301 genes in 139 OGs hit > 8 unique places ***
        # ...Thlaspi_arvense        :  1073 genes in  53 OGs hit > 8 unique places
        # NOTE! Genomes flagged *** have > 5% of genes in over-dispersed
                # orthogroups. These are likely not great annotations, or
                # the synteny run contains un-specified WGDs. Regardless,
                # these should be examined carefully

#----------------------
# Polygonaceae
gpar <- init_genespace(
  wd = "/home/celphin/scratch/Oxyria/GeneSpace/Polygonaceae_genomes", 
  path2mcscanx = "/home/celphin/scratch/Oxyria/GeneSpace/MCScanX/",
  rawOrthofinderDir = "/home/celphin/scratch/Oxyria/GeneSpace/Polygonaceae_genomes/orthofinder/Results_Oct07")

# Flagging over-dispered OGs
        # ...Fagopyrum_escelentum_H2: 12819 genes in 285 OGs hit > 8 unique places ***
        # ...Fagopyrum_tataricum_H1 :  4214 genes in 140 OGs hit > 8 unique places ***
        # ...Oxyria_digyna_H1       :  6187 genes in 218 OGs hit > 8 unique places ***
        # ...Polygunum_aviculare_H0 :  2025 genes in  70 OGs hit > 8 unique places ***
        # ...Rheum_nobile_H0        :  2384 genes in 144 OGs hit > 8 unique places ***
        # ...Rheum_tangaticum_H0    :  1647 genes in 101 OGs hit > 8 unique places ***
        # NOTE! Genomes flagged *** have > 5% of genes in over-dispersed
                # orthogroups. These are likely not great annotations, or
                # the synteny run contains un-specified WGDs. Regardless,
                # these should be examined carefully
#------------------------------
# Brassicaceae
gpar <- init_genespace(
  wd = "/home/celphin/scratch/Oxyria/GeneSpace/Brassicaceae_genomes", 
  path2mcscanx = "/home/celphin/scratch/Oxyria/GeneSpace/MCScanX/",
  rawOrthofinderDir = "/home/celphin/scratch/Oxyria/GeneSpace/Brassicaceae_genomes/orthofinder/Results_Oct09")

# Flagging chrs. w/ < 10 unique orthogroups
        # ...Arabidopsis_lyrata     :  694 genes on  376 small chrs.
        # ...Arabidopsis_thaliana   :    0 genes on    0 small chrs.
        # ...Arabis_alpina          : 2414 genes on 1734 small chrs. ***
        # ...Brassica_oleracea      :    0 genes on    0 small chrs.
        # ...Cochlearia_groenlandica:   72 genes on   24 small chrs.
        # ...Draba_nivalis          :  138 genes on   60 small chrs.
        # ...Thlaspi_arvense        : 2163 genes on  654 small chrs. ***
        # NOTE! Genomes flagged *** have > 5% of genes on small chrs.
                # These are likely not great assemblies and should be
                # examined carefully

#-------------------------------
# Rosaceae

gpar <- init_genespace(
  wd = "/home/celphin/scratch/Oxyria/GeneSpace/Rosaceae_genomes", 
  path2mcscanx = "/home/celphin/scratch/Oxyria/GeneSpace/MCScanX/",
  rawOrthofinderDir = "/home/celphin/scratch/Oxyria/GeneSpace/Rosaceae_genomes/orthofinder/Results_Aug16")

# Checking annotation files (.bed and peptide .fa):
        # Argentina_anserina  : 19620 / 19620 geneIDs exactly match (PASS)
        # Dryas_octopetala    : 39696 / 39696 geneIDs exactly match (PASS)
        # Fragaria_vesca      : 21397 / 21397 geneIDs exactly match (PASS)
        # Malus_sylvestris    : 37467 / 37467 geneIDs exactly match (PASS)
        # Prunus_persica      : 23128 / 23128 geneIDs exactly match (PASS)
        # Pyrus_bretschneideri: 35293 / 35293 geneIDs exactly match (PASS)
        # Rosa_rugosa         : 29146 / 29146 geneIDs exactly match (PASS)

#-----------------------------
# Dryas
gpar <- init_genespace(
  wd = "/home/celphin/scratch/Oxyria/GeneSpace/Dryas_genomes", 
  path2mcscanx = "/home/celphin/scratch/Oxyria/GeneSpace/MCScanX/",
  rawOrthofinderDir = "/home/celphin/scratch/Oxyria/GeneSpace/Dryas_genomes/orthofinder/Results_Oct30")

        # Flagging over-dispered OGs
        # ...Dry_ajan            :  936 genes in  55 OGs hit > 8 unique places
        # ...Dry_alask           : 1011 genes in  63 OGs hit > 8 unique places
        # ...Dry_drumm           : 2160 genes in 124 OGs hit > 8 unique places ***
        # ...Dry_int             :  869 genes in  52 OGs hit > 8 unique places
        # ...Dry_octo_H0         : 4428 genes in 187 OGs hit > 8 unique places ***
        # ...Dry_octo_H1         : 3485 genes in 158 OGs hit > 8 unique places ***
        # ...Malus_sylvestris    : 1962 genes in  71 OGs hit > 8 unique places ***
        # ...Pyrus_bretschneideri:  479 genes in  30 OGs hit > 8 unique places
        # ...Rosa_rugosa         : 2300 genes in 100 OGs hit > 8 unique places ***

      # Annotation summaries (after exclusions):
        # ...Dry_ajan            : 27779 genes in 25524 OGs || 1093 genes in  513 arrays
        # ...Dry_alask           : 27527 genes in 25272 OGs || 1086 genes in  508 arrays
        # ...Dry_drumm           : 31834 genes in 27031 OGs || 2143 genes in  949 arrays
        # ...Dry_int             : 27882 genes in 25611 OGs || 1096 genes in  512 arrays
        # ...Dry_octo_H0         : 35214 genes in 28218 OGs || 4447 genes in 1744 arrays
        # ...Dry_octo_H1         : 32431 genes in 26570 OGs || 4027 genes in 1728 arrays
        # ...Malus_sylvestris    : 35504 genes in 20904 OGs || 4236 genes in 1702 arrays
        # ...Pyrus_bretschneideri: 34366 genes in 20137 OGs || 4791 genes in 1837 arrays
        # ...Rosa_rugosa         : 26784 genes in 19903 OGs || 5490 genes in 2027 arrays

#-------------------------------
run_genespace(gsParam = gpar)
# runs until first set of dot plots then freezes

####################################################
# Try removing  dotplots directly in R
# https://github.com/jtlovell/GENESPACE/tree/master/R

# https://github.com/jtlovell/GENESPACE/blob/master/R/run_genespace.R

# try whole function without dotplots (below)

# Total
gpar <- init_genespace(
  wd = "/home/celphin/scratch/Oxyria/GeneSpace/Total_genomes", 
  path2mcscanx = "/home/celphin/scratch/Oxyria/GeneSpace/MCScanX/",
  rawOrthofinderDir = "/home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/orthofinder/Results_Aug18")

#-------------------------
# Polygonaceae
gpar <- init_genespace(
  wd = "/home/celphin/scratch/Oxyria/GeneSpace/Polygonaceae_genomes", 
  path2mcscanx = "/home/celphin/scratch/Oxyria/GeneSpace/MCScanX/",
  rawOrthofinderDir = "/home/celphin/scratch/Oxyria/GeneSpace/Polygonaceae_genomes/orthofinder/Results_Oct07")

#-------------------------------
# Brassicaceae

gpar <- init_genespace(
  wd = "/home/celphin/scratch/Oxyria/GeneSpace/Brassicaceae_genomes", 
  path2mcscanx = "/home/celphin/scratch/Oxyria/GeneSpace/MCScanX/",
  rawOrthofinderDir = "/home/celphin/scratch/Oxyria/GeneSpace/Brassicaceae_genomes/orthofinder/Results_Oct09")

#-------------------------------
# Rosaceae

gpar <- init_genespace(
  wd = "/home/celphin/scratch/Oxyria/GeneSpace/Rosaceae_genomes", 
  path2mcscanx = "/home/celphin/scratch/Oxyria/GeneSpace/MCScanX/",
  rawOrthofinderDir = "/home/celphin/scratch/Oxyria/GeneSpace/Rosaceae_genomes/orthofinder/Results_Aug16")

#-------------------------------
# Dryas

gpar <- init_genespace(
  wd = "/home/celphin/scratch/Oxyria/GeneSpace/Dryas_genomes", 
  path2mcscanx = "/home/celphin/scratch/Oxyria/GeneSpace/MCScanX/",
  rawOrthofinderDir = "/home/celphin/scratch/Oxyria/GeneSpace/Dryas_genomes/orthofinder/Results_Oct30")


#-------------------------------
# Oxyria

gpar <- init_genespace(
  wd = "/home/celphin/scratch/Oxyria/GeneSpace/Oxyria_genomes", 
  path2mcscanx = "/home/celphin/scratch/Oxyria/GeneSpace/MCScanX/",
  rawOrthofinderDir = "/home/celphin/scratch/Oxyria/GeneSpace/Oxyria_genomes/orthofinder/Results_Oct31")



#-------------------
library(data.table)
# using code at bottom
out <- run_genespace_2(gsParam = gpar)

############################
# GENESPACE run complete!  All results are stored in
# /home/celphin/scratch/Oxyria/GeneSpace in the following subdirectories:
        # syntenic block dotplots: /dotplots (...synHits.pdf)
        # annotated blast files  : /syntenicHits
        # annotated/combined bed : /results/combBed.txt
        # syntenic block coords. : /results/blkCoords.txt
        # syn. blk. by ref genome: /riparian/refPhasedBlkCoords.txt
        # pan-genome annotations : /pangenes (...pangenes.txt.gz)
        # riparian plots         : /riparian
        # genespace param. list  : /results/gsParams.rda
# ############################
# **NOTE** the genespace parameter object is returned or can be loaded
        # into R via
        # `load('/home/celphin/scratch/Oxyria/GeneSpace/results/gsParams.rda',
        # verbose = TRUE)`. Then you can customize your riparian plots by
        # calling `plot_riparian(gsParam = gsParam, ...)`. The source
        # data and ggplot2 objects are also stored in the /riparian
        # directory and can also be accessed by `load(...)`.
# **NOTE** To query genespace results by position or gene, use
        # `query_genespace(...)`. See specifications in ?query_genespace
        # for details.

##############################################
# more in depth use and plots:
# https://htmlpreview.github.io/?https://github.com/jtlovell/tutorials/blob/main/genespaceGuide.html

tmux new-session -s GeneSpace
tmux attach-session -t GeneSpace

cd /home/celphin/scratch/Oxyria/GeneSpace/

module load StdEnv/2020  python/3.10.2 scipy-stack/2021a diamond/2.1.7 fastme/2.1.6.2 gcc/9.3.0 mcl/14.137 java/13.0.2 r/4.2.2 glpk/5.0

# other version of orthofinder
export PATH=$PATH:/home/celphin/scratch/Oxyria/GeneSpace/OrthoFinder/bin
export PATH=$PATH:/home/celphin/scratch/Oxyria/GeneSpace/OrthoFinder_source
alias orthofinder='orthofinder.py'

#orthofinder -h
# works

R

library(GENESPACE)

out <- load('/home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/results/gsParams.rda', verbose = TRUE)

# https://htmlpreview.github.io/?https://github.com/jtlovell/tutorials/blob/main/riparianGuide.html
# remake riparian plot

# gsParam$genomeIDs
 # [1] "Arabidopsis_lyrata"      "Arabidopsis_thaliana"
 # [3] "Capsella_rubella"        "Arabis_alpina"
 # [5] "Draba_nivalis"           "Thlaspi_arvense"
 # [7] "Brassica_oleracea"       "Cochlearia_groenlandica"
 # [9] "Malus_sylvestris"        "Pyrus_bretschneideri"
# [11] "Prunus_persica"          "Dryas_octopetala"
# [13] "Fragaria_vesca"          "Argentina_anserina"
# [15] "Rosa_rugosa"             "Rheum_nobile_H0"
# [17] "Rheum_tangaticum_H0"     "Oxyria_digyna_H1"
# [19] "Polygunum_aviculare_H0"  "Fagopyrum_escelentum_H2"
# [21] "Fagopyrum_tataricum_H1"

ggthemes <- ggplot2::theme(
  panel.background = ggplot2::element_rect(fill = "white"))
customPal <- colorRampPalette(
  c("darkorange", "skyblue", "darkblue", "purple", "darkred", "salmon"))

#--------------------------------------
grDevices::pdf("/home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/Draba_Coch_riparian.pdf")
print(plot_riparian(
  gsParam = gsParam, 
  backgroundColor = NULL,
  useOrder = TRUE, useRegions = FALSE, reorderBySynteny = TRUE,
  genomeIDs = c("Draba_nivalis", "Cochlearia_groenlandica" ),
  refGenome = "Draba_nivalis",   
  palette = customPal,
  braidAlpha = .75,
  chrFill = "lightgrey",
  addThemes = ggthemes))
grDevices::dev.off()

#--------------------------------------
grDevices::pdf("/home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/Dryas_Fragaria_riparian.pdf")
print(plot_riparian(
  gsParam = gsParam, 
  backgroundColor = NULL,
  useOrder = TRUE, useRegions = FALSE, reorderBySynteny = TRUE,
  genomeIDs = c("Dryas_octopetala", "Fragaria_vesca" ),
  refGenome = "Dryas_octopetala",   
  palette = customPal,
  braidAlpha = .75,
  chrFill = "lightgrey",
  addThemes = ggthemes))
grDevices::dev.off()

#--------------------------------------
grDevices::pdf("/home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/Dryas_Oxyria_riparian.pdf")
print(plot_riparian(
  gsParam = gsParam, 
  backgroundColor = NULL,
  useOrder = TRUE, useRegions = FALSE, reorderBySynteny = TRUE,
  genomeIDs = c("Dryas_octopetala", "Oxyria_digyna_H1", "Draba_nivalis" ),
  refGenome = "Dryas_octopetala",   
  palette = customPal,
  braidAlpha = .75,
  chrFill = "lightgrey",
  addThemes = ggthemes))
grDevices::dev.off()

#--------------------------------------
# invert some chromosomes
# invchr <- data.frame(
  # genome = c("Fagopyrum_escelentum_H2", "Fagopyrum_escelentum_H2", "Fagopyrum_escelentum_H2"),
 # chr = c("2", "7", "5"))

# roi <- data.frame(
  # genome = c("Oxyria_digyna_H1"), 
  # chr = c("Oxyrt-1-86582034", "Oxyrt-2-79714091", "Oxyrt-3-79472951","Oxyrt-4-78410798","Oxyrt-5-76064323","Oxyrt-6-73303751","Oxyrt-7-72361354"))
  # #chr = c("Oxyrt-4-78410798","Oxyrt-5-76064323","Oxyrt-6-73303751"))
 

grDevices::pdf("/home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/OxydigH1_riparian.pdf")
print(plot_riparian(
  gsParam = gsParam, 
  backgroundColor = NULL,
  useOrder = TRUE, useRegions = FALSE, reorderBySynteny = TRUE,
  genomeIDs = c("Oxyria_digyna_H1", "Rheum_nobile_H0" ),
  refGenome = "Oxyria_digyna_H1",   
  palette = customPal,
  braidAlpha = .75,
  chrFill = "lightgrey",
  addThemes = ggthemes))
grDevices::dev.off()

#----------------------
roi <- data.frame(
  genome = c("Oxyria_digyna_H1"), 
  chr = c("Oxyrt-1-86582034", "Oxyrt-2-79714091", "Oxyrt-3-79472951","Oxyrt-4-78410798","Oxyrt-5-76064323","Oxyrt-6-73303751","Oxyrt-7-72361354"))
  #chr = c("Oxyrt-4-78410798","Oxyrt-5-76064323","Oxyrt-6-73303751"))
  
ggthemes <- ggplot2::theme(
  panel.background = ggplot2::element_rect(fill = "white"))
customPal <- colorRampPalette(
  c("darkorange", "skyblue", "darkblue", "purple", "darkred", "salmon"))
  
grDevices::pdf("/home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/OxydigH1_riparian_4spp.pdf")
print(plot_riparian(
  gsParam = gsParam, 
  useOrder = FALSE, useRegions = FALSE, reorderBySynteny = TRUE,
  highlightBed = roi, 
  backgroundColor = NULL,
  palette = customPal,
  braidAlpha = .75,
  chrFill = "lightgrey",
  addThemes = ggthemes,
  genomeIDs = c("Oxyria_digyna_H1","Draba_nivalis"),
  refGenome = "Oxyria_digyna_H1"))
grDevices::dev.off()

grDevices::pdf("/home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/OxydigH1_riparian_Rnob.pdf")
print(plot_riparian(
  gsParam = gsParam, 
  useOrder = FALSE, useRegions = FALSE, reorderBySynteny = TRUE,
  highlightBed = roi, 
  backgroundColor = NULL,
  palette = customPal,
  braidAlpha = .75,
  chrFill = "lightgrey",
  addThemes = ggthemes,
  genomeIDs = c("Oxyria_digyna_H1", "Cochlearia_groenlandica"),
  refGenome = "Oxyria_digyna_H1"))
grDevices::dev.off()

grDevices::pdf("/home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/OxydigH1_riparian_Polavi.pdf")
print(plot_riparian(
  gsParam = gsParam, 
  useOrder = FALSE, useRegions = FALSE, reorderBySynteny = TRUE,
  highlightBed = roi, 
  backgroundColor = NULL,
  palette = customPal,
  braidAlpha = .75,
  chrFill = "lightgrey",
  addThemes = ggthemes,
  genomeIDs = c("Oxyria_digyna_H1", "Dryas_octopetala"),
  refGenome = "Oxyria_digyna_H1"))
grDevices::dev.off()

#########################################
# plot specific chromosomes
R

library(GENESPACE)

out <- load('/home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/results/gsParams.rda', verbose = TRUE)

roi <- data.frame(
  genome = c("Oxyria_digyna_H1"), 
  chr = c("Oxyrt-6-73303751"), 
  color = c("#17B5C5"))

# invert some chromosomes
invchr <- data.frame(
  genome = c("Rheum_nobile_H0", "Rheum_nobile_H0", "Rheum_nobile_H0"),
  chr = c("rno05", "rno08", "rno09"))

grDevices::pdf("/home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/OxydigH1_riparian_subsetchr6only.pdf")
print(
ripDat <- plot_riparian(
  gsParam = gsParam, 
  highlightBed = roi, 
  useOrder = FALSE, useRegions = FALSE, #reorderBySynteny = TRUE,
  backgroundColor = NULL,
  invertTheseChrs = invchr,  
  genomeIDs = c("Oxyria_digyna_H1","Rheum_tangaticum_H0", "Rheum_nobile_H0", "Polygunum_aviculare_H0", "Fagopyrum_escelentum_H2", "Fagopyrum_tataricum_H1"),
  refGenome = "Oxyria_digyna_H1"))
grDevices::dev.off()


grDevices::pdf("/home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/OxydigH1_riparian_subsetchr6only_geneorder.pdf")
print(
ripDat <- plot_riparian(
  gsParam = gsParam, 
  highlightBed = roi, 
  useOrder = TRUE, useRegions = TRUE, reorderBySynteny = TRUE,
  backgroundColor = NULL, 
  invertTheseChrs = invchr,
  genomeIDs = c("Oxyria_digyna_H1","Rheum_tangaticum_H0", "Rheum_nobile_H0", "Polygunum_aviculare_H0", "Fagopyrum_escelentum_H2", "Fagopyrum_tataricum_H1"),
  refGenome = "Oxyria_digyna_H1"))
grDevices::dev.off()

#---------------------------------
# plot region of interest
roibed <- roi[,c("genome", "chr")]
roibed$color <- c("pink", "cyan", "green")
ripd <- plot_riparian(
  gsParam = gsParam, 
  useRegions = FALSE, 
  highlightBed = roibed)

# find pangenes
# genes in same orthogroup and are found together
test <- query_pangenes(
  gsParam = gsParam, bed = roi)

##########################################
# For Dryas Genomes

tmux new-session -s GeneSpace
tmux attach-session -t GeneSpace

cd /home/celphin/scratch/Oxyria/GeneSpace/

module load StdEnv/2020  python/3.10.2 scipy-stack/2021a diamond/2.1.7 fastme/2.1.6.2 gcc/9.3.0 mcl/14.137 java/13.0.2 r/4.2.2 glpk/5.0

# other version of orthofinder
export PATH=$PATH:/home/celphin/scratch/Oxyria/GeneSpace/OrthoFinder/bin
export PATH=$PATH:/home/celphin/scratch/Oxyria/GeneSpace/OrthoFinder_source
alias orthofinder='orthofinder.py'

#orthofinder -h
# works
R

library(GENESPACE)

out <- load('/home/celphin/scratch/Oxyria/GeneSpace/Dryas_genomes/results/gsParams.rda', verbose = TRUE)

# https://htmlpreview.github.io/?https://github.com/jtlovell/tutorials/blob/main/riparianGuide.html
# remake riparian plot

# gsParam$genomeIDs
 # [1] "Arabidopsis_lyrata"      "Arabidopsis_thaliana"
 # [3] "Capsella_rubella"        "Arabis_alpina"
 # [5] "Draba_nivalis"           "Thlaspi_arvense"
 # [7] "Brassica_oleracea"       "Cochlearia_groenlandica"
 # [9] "Malus_sylvestris"        "Pyrus_bretschneideri"
# [11] "Prunus_persica"          "Dryas_octopetala"
# [13] "Fragaria_vesca"          "Argentina_anserina"
# [15] "Rosa_rugosa"             "Rheum_nobile_H0"
# [17] "Rheum_tangaticum_H0"     "Oxyria_digyna_H1"
# [19] "Polygunum_aviculare_H0"  "Fagopyrum_escelentum_H2"
# [21] "Fagopyrum_tataricum_H1"

ggthemes <- ggplot2::theme(
  panel.background = ggplot2::element_rect(fill = "white"))
customPal <- colorRampPalette(
  c("darkorange", "skyblue", "darkblue", "purple", "darkred", "salmon"))

#--------------------------------------
grDevices::pdf("/home/celphin/scratch/Oxyria/GeneSpace/Dryas_genomes/Dryas_riparian_edited.pdf")
print(plot_riparian(
  gsParam = gsParam, 
  backgroundColor = NULL,
  useOrder = TRUE, useRegions = FALSE, reorderBySynteny = TRUE,
  genomeIDs = c("Dry_alask", "Dry_ajan", "Dry_int",  "Dry_octo_H0", "Dry_octo_H1", "Dry_drumm", "Rosa_rugosa", "Malus_sylvestris" ),
  refGenome = "Dry_octo_H0",   
  palette = customPal,
  braidAlpha = .75,
  chrFill = "lightgrey",
  addThemes = ggthemes))
grDevices::dev.off()

##############################
# Chromosome labels for Polygonaceae

                 chrLab
 1:                    7
 2:                    2
 3:                    1
 4:                    8
 5:                    4
 6:                    5
 7:                    6
 8:                    3
 9:                    7
10:                    2
11:                    8
12:                    1
13:                    4
14:                    5
15:                    6
16:                    3
17: fagopyrum-6-52287906
18: fagopyrum-2-61235386
19: fagopyrum-1-68031765
20: fagopyrum-4-56655744
21: fagopyrum-3-57706077
22: fagopyrum-8-49982843
23: fagopyrum-5-53883329
24: fagopyrum-7-51545819
25:                    7
26:                    3
27:                    2
28:                    1
29:                    6
30:                    4
31:                    8
32:                    5
33:                    7
34:                    3
35:                    2
36:                    1
37:                    4
38:                    8
39:                    6
40:                    5
41:       oxy-1-88146246
42:       oxy-2-80787722
43:       oxy-3-77175000
44:       oxy-4-76036264
45:        oxy-9-3540000
46:       oxy-7-70136240
47:       oxy-11-2000000
48:       oxy-12-1350245
49:       oxy-5-73842794
50:       oxy-6-73289246
51:     oxyrt-1-86582034
52:     oxyrt-2-79714091
53:     oxyrt-3-79472951
54:     oxyrt-4-78410798
55:     oxyrt-5-76064323
56:     oxyrt-6-73303751
57:     oxyrt-7-72361354
58:             polavi-8
59:             polavi-5
60:             polavi-7
61:            polavi-10
62:             polavi-2
63:             polavi-1
64:             polavi-9
65:             polavi-4
66:             polavi-3
67:             polavi-6
68:                rno02
69:                rno07
70:                rno03
71:                rno05
72:                rno10
73:                rno01
74:                rno06
75:                rno09
76:                rno08
77:                rno11
78:                rno04
79:                rta03
80:                rta09
81:                rta02
82:                rta05
83:                rta06
84:                rta01
85:                rta07
86:                rta10
87:                rta08
88:                rta11
89:                rta04



# https://htmlpreview.github.io/?https://github.com/jtlovell/tutorials/blob/main/genomeVizGuide.html

















#########################################
# New editable verion of the overall function

run_genespace_2 <- function(gsParam,
                          overwrite = FALSE,
                          overwriteBed = overwrite,
                          overwriteSynHits = overwrite,
                          overwriteInBlkOF = TRUE,
                          makePairwiseFiles = FALSE){

   gsParam$paths$rawOrthofinder <- gsParam$paths$orthofinder
  ##############################################################################
  # 1. Run orthofinder ...
  cat("\n############################", strwrap(
    "1. Running orthofinder (or parsing existing results)",
    indent = 0, exdent = 8), sep = "\n")

  ##############################################################################
  # -- 1.1 Check for existing parsed orthofinder results
  cat("\tChecking for existing orthofinder results ...\n")
  gsParam <- set_syntenyParams(gsParam)

  if("synteny" %in% names(gsParam)){
    noResults <- is.na(gsParam$synteny$SpeciesIDs)
  }else{
    noResults <- TRUE
  }

  ##############################################################################
  # -- 1.2 If no results exist, check for raw orthofinder run
  if(noResults){
    if(dir.exists(gsParam$paths$rawOrthofinder)){
      chkOf <- find_ofFiles(gsParam$paths$rawOrthofinder)
      noOrthofinder <- is.na(chkOf[[1]])
    }else{
      noOrthofinder <- TRUE
    }
  }else{
    noOrthofinder <- FALSE
  }

  ##############################################################################
  # -- 1.3 If raw results exist, copy them over
  if(!noOrthofinder && noResults){
    with(gsParam, copy_of2results(
      orthofinderDir = paths$rawOrthofinder, resultsDir = paths$results))

    if(dir.exists(gsParam$paths$rawOrthofinder)){
      chkOf <- find_ofFiles(gsParam$paths$rawOrthofinder)
      noOrthofinder <- is.na(chkOf[[1]])
      noResults <- noOrthofinder
    }else{
      noOrthofinder <- TRUE
    }
  }

  print(noResults)

  if(!noResults){
    spids <- names(read_orthofinderSpeciesIDs(
      file.path(gsParam$paths$results, "SpeciesIDs.txt")))
    gid <- unique(c(gsParam$genomeIDs, gsParam$outgroup))
    gid <- gid[!is.na(gid)]
    ps <- all(gid %in% spids) && all(spids %in% gid)
    if(ps){
      cat("\t... found existing run, not re-running orthofinder\n")
    }else{
      stop("genomes in the existing orthofinder run do not exactly match specified genomeIDs\n")
    }
  }


  ##############################################################################
  # -- 1.4 if no orthofinder run, make one
  if(noResults)
    tmp <- run_orthofinder(gsParam = gsParam, verbose = TRUE)

  gsParam <- set_syntenyParams(gsParam)
  gsParam <<- gsParam
  ##############################################################################
  # -- 1.5 get the files in order if the run is complete
  if(noResults){
    chkOf <- find_ofFiles(gsParam$paths$orthofinder)
    noOrthofinder <- is.na(chkOf[[1]])
    if(noOrthofinder)
      stop("could not find orthofinder files!")
    with(gsParam, copy_of2results(
      orthofinderDir = paths$orthofinder,
      resultsDir = paths$results))
  }
  gsParam <- set_syntenyParams(gsParam)

  ##############################################################################
  # -- 1.6 if the species tree exists, re-order the genomeIDs
  tmp <- gsParam$synteny$speciesTree

  if(requireNamespace("ape", quietly = T)){
    if(!is.na(tmp) && !is.null(tmp)){
      if(file.exists(tmp) && length(gsParam$genomeIDs) > 2){
        treLabs <- get_orderedTips(
          treFile = gsParam$synteny$speciesTree,
          ladderize = TRUE,
          genomeIDs = gsParam$genomeIDs)
        cat(strwrap(sprintf(
          "re-ordering genomeIDs by the species tree: %s",
          paste(treLabs, collapse = ", ")), indent = 8, exdent = 16),
          sep = "\n")
        gsParam$genomeIDs <- treLabs
      }
    }
  }

  # -- 1.7 if useHOGs, check if the N0.tsv file exists
  useHOGs <- gsParam$params$useHOGs
  if(useHOGs){
    if(is.na(gsParam$synteny$hogs)){
      useHOGs <- FALSE
    }else{
      if(!file.exists(gsParam$synteny$hogs)){
        useHOGs <- FALSE
      }
    }
  }
  gsParam$params$useHOGs <- useHOGs
  useHOGs <- NULL

  ##############################################################################
  # 2. Get the data ready for synteny
  hasBed <- FALSE
  bedf <- gsParam$synteny$combBed
  if(file.exists(bedf))
    hasBed <- is.data.table(read_combBed(bedf))
  if(overwriteBed)
    hasBed <- FALSE

  if(!hasBed){
    cat("\n############################", strwrap(
      "2. Combining and annotating bed files w/ OGs and tandem array info ... ",
      indent = 0, exdent = 8), sep = "\n")
    bed <- annotate_bed(gsParam = gsParam)
  }else{
    cat("\n############################", strwrap(
      "2. Annotated/concatenated bed file exists", indent = 0, exdent = 8),
      sep = "\n")
  }

  ##############################################################################
  # 3. Annotate the blast files ...
  # -- First make sure that the blast files are all there, then go through
  # and annotate them with the combined bed file
  # -- This also makes the first round of dotplots
  # -- 3.1 check if all the synHits files exist. If so, and !overwriteSynHits
  # don't re-annotate
  hasHits <- FALSE
  if(all(file.exists(gsParam$synteny$blast$synHits)))
    if(!overwriteSynHits)
      hasHits <- TRUE

  # -- 3.2 iterate through and annotate all synHits files
  if(!hasHits){
    cat("\n############################", strwrap(
      "3. Combining and annotating the blast files with orthogroup info ...",
      indent = 0, exdent = 8), sep = "\n")
    gsParam <- annotate_blast(gsParam = gsParam)
  }else{
    cat("\n############################", strwrap(
      "3. Annotated/blast files exists", indent = 0, exdent = 8),
      sep = "\n")
  }

  dpFiles <- with(gsParam$synteny$blast, file.path(
    file.path(gsParam$paths$wd, "dotplots",
    sprintf("%s_vs_%s.rawHits.pdf",
            query, target))))
# if(!all(file.exists(dpFiles)) || overwrite){
    # cat("\t##############\n\tGenerating dotplots for all hits ... ")
    # nu <- plot_hits(gsParam = gsParam, type = "raw")
    # cat("Done!\n")
  # }

  ##############################################################################
  # 4. Run synteny
  # -- goes through each pair of genomes and pulls syntenic anchors and the hits
  # nearby. This is the main engine of genespace
  bed <- read_combBed(bedf)
  hasSynFiles <- all(file.exists(gsParam$synteny$blast$synHits))
  hasOg <- all(!is.na(bed$og))
  if(!hasOg || !hasSynFiles || overwrite){
    cat("\n############################", strwrap(
      "4. Flagging synteny for each pair of genomes ...",
      indent = 0, exdent = 8), sep = "\n")
    gsParam <- synteny(gsParam = gsParam)
  }

  ##############################################################################
  # 5. Build syntenic orthogroups
  cat("\n############################", strwrap(
    "5. Building synteny-constrained orthogroups ... ",
    indent = 0, exdent = 8), sep = "\n")
  # -- in the case of polyploid genomes, this also runs orthofinder in blocks,
  # then re-runs synteny and re-aggregates blocks,.
  if(gsParam$params$orthofinderInBlk & overwriteInBlkOF){

    # -- returns the gsparam obj and overwrites the bed file with a new og col
    cat("\t##############\n\tRunning Orthofinder within syntenic regions\n")
    tmp <- run_orthofinderInBlk(
      gsParam = gsParam, overwrite = overwriteSynHits)

    # -- adds a new column to the bed file
    cat("\t##############\n\tPulling syntenic orthogroups\n")
  }

  gsParam <- syntenic_orthogroups(
    gsParam, updateArrays = gsParam$params$orthofinderInBlk)
  cat("\tDone!\n")
  gsParam <<- gsParam

  ##############################################################################
  # 6. Make dotplots
  cat("\n############################", strwrap(
    "6. Integrating syntenic positions across genomes ... ",
    indent = 0, exdent = 8), sep = "\n")

  dpFiles <- with(gsParam$synteny$blast, file.path(
    file.path(gsParam$paths$wd, "dotplots",
              sprintf("%s_vs_%s.synHits.pdf",
                      query, target))))
  # if(!all(file.exists(dpFiles)) || overwrite){
    # cat("\t##############\n\tGenerating syntenic dotplots ... ")
    # nu <- plot_hits(gsParam = gsParam, type = "syntenic")
    # cat("Done!\n")
  # }

  ##############################################################################
  # 7. Interpolate syntenic positions
  cat("\t##############\n\tInterpolating syntenic positions of genes ... \n")
  nsynPos <- interp_synPos(gsParam)
  cat("\tDone!\n")
  gsParam <<- gsParam

  ##############################################################################
  # 8. Phase syntenic blocks against reference chromosomes
  cat("\n############################", strwrap(
    "7. Final block coordinate calculation and riparian plotting ... ",
    indent = 0, exdent = 8), sep = "\n")
  cat("\t##############\n\tCalculating syntenic blocks by reference chromosomes ... \n")
  reg <- nophase_blks(gsParam = gsParam, useRegions = T)
  cat(sprintf("\t\tn regions (aggregated by %s gene radius): %s\n",
              gsParam$params$blkRadius, nrow(reg)))
  fwrite(reg, file = file.path(gsParam$paths$results, "syntenicRegion_coordinates.csv"))
  blk <- nophase_blks(gsParam = gsParam, useRegions = F)
  cat(sprintf("\t\tn blocks (collinear sets of > %s genes): %s\n",
              gsParam$params$blkSize, nrow(blk)))
  fwrite(blk, file = file.path(gsParam$paths$results, "syntenicBlock_coordinates.csv"))
  hapGenomes <- names(gsParam$ploidy)[gsParam$ploidy == 1]

  bed <- read_combBed(file.path(gsParam$paths$results, "combBed.txt"))
  minChrSize <- gsParam$params$blkSize * 2
  isOK <- og <- chr <- genome <- NULL
  bed[,isOK := uniqueN(og) >= minChrSize, by = c("genome", "chr")]
  ok4pg <- bed[,list(propOK = (sum(isOK) / .N) > .75), by = "genome"]
  ok4rip <- bed[,list(nGood = (uniqueN(chr[isOK]) < 100)), by = "genome"]

  if(length(hapGenomes) == 0){
    cat(strwrap("NOTE!!! No genomes provided with ploidy < 2. Phasing of polyploid references is not currently supported internally. You will need to make custom riparian plots",
                indent = 8, exdent = 8), sep = "\n")
  }else{
    okg <- subset(ok4rip, genome %in% hapGenomes)
    if(any(!okg$nGood)){
      cat(strwrap(sprintf(
        "**WARNING**: genomes %s have > 100 chromosomes in the synteny map. This is too complex to make riparian plots.\n",
        paste(okg$genome[!okg$nGood], collapse = ", ")), indent = 8, exdent = 16),
        sep = "\n")
      hapGenomes <- hapGenomes[!hapGenomes %in% okg$genome[!okg$nGood]]
    }
    if(length(hapGenomes) > 0){
      cat("\t##############\n\tBuilding ref.-phased blks and riparian plots for haploid genomes:\n")
      labs <- align_charLeft(hapGenomes)
      names(labs) <- hapGenomes
      for(i in hapGenomes){
        plotf <- file.path(gsParam$paths$riparian,
                           sprintf("%s_geneOrder.rip.pdf", i))

        srcf <- file.path(gsParam$paths$riparian,
                          sprintf("%s_geneOrder_rSourceData.rda", i))
        blkf <- file.path(gsParam$paths$riparian,
                          sprintf("%s_phasedBlks.csv", i))

        rip <- plot_riparian(
          gsParam = gsParam, useRegions = TRUE, refGenome = i, pdfFile = plotf)
        cat(sprintf("\t\t%s: %s phased blocks\n", labs[i], nrow(rip$blks)))

        srcd <- rip$plotData
        save(srcd, file = srcf)
        fwrite(rip$blks, file = blkf)

        plotf <- file.path(gsParam$paths$riparian,
                           sprintf("%s_bp.rip.pdf", i))
        srcf <- file.path(gsParam$paths$riparian,
                          sprintf("%s_bp_rSourceData.rda", i))
        rip <- plot_riparian(
          gsParam = gsParam, useOrder = FALSE, useRegions = TRUE,
          refGenome = i, pdfFile = plotf)
        srcd <- rip$plotData
        save(srcd, file = srcf)
      }
      cat("\tDone!\n")
    }
  }

  gsParam <<- gsParam

  ##############################################################################
  # 9. Build pan-genes (aka pan-genome annotations)
  cat("\n############################", strwrap(
    "8. Constructing syntenic pan-gene sets ... ",
    indent = 0, exdent = 8), sep = "\n")
  gids <- gsParam$genomeIDs
  tp <- paste(ok4pg$genome[!ok4pg$propOK], collapse = ", ")
  if(any(!ok4pg$propOK))
    cat(strwrap(sprintf(
      "**WARNING**: genomes %s have < 75%% of genes on chromosomes that contain > %s genes. Synteny is not a useful metric for these genomes. Be very careful with your pan-gene sets.\n",
      tp, minChrSize), indent = 8, exdent = 16),
      sep = "\n")
  labs <- align_charLeft(gids)
  names(labs) <- gids
  for(i in gids){
    cat(sprintf("\t%s: ", labs[i]))
    pgref <- syntenic_pangenes(gsParam = gsParam, refGenome = i)
    with(pgref, cat(sprintf(
      "n pos. = %s, synOgs = %s, array mem. = %s, NS orthos %s\n",
      uniqueN(pgID), sum(flag == "PASS"), sum(flag == "array"), sum(flag == "NSOrtho"))))
  }



  ##############################################################################
  # 8. Print summaries and return

  # --- make pairwise files
  if(makePairwiseFiles){
      cat("Building pairiwse hits files ...\n")
    pull_pairwise(gsParam, verbose = TRUE)
      cat("Done!")
  }

  gpFile <- file.path(gsParam$paths$results, "gsParams.rda")
  cat("\n############################", strwrap(sprintf(
    "GENESPACE run complete!\n All results are stored in %s in the following subdirectories:",
    gsParam$paths$wd), indent = 0, exdent = 0),
    "\tsyntenic block dotplots: /dotplots (...synHits.pdf)",
    "\tannotated blast files  : /syntenicHits",
    "\tannotated/combined bed : /results/combBed.txt",
    "\tsyntenic block coords. : /results/blkCoords.txt",
    "\tsyn. blk. by ref genome: /riparian/refPhasedBlkCoords.txt",
    "\tpan-genome annotations : /pangenes (...pangenes.txt.gz)",
    "\triparian plots         : /riparian",
    "\tgenespace param. list  : /results/gsParams.rda",
    "############################",
    strwrap(sprintf(
      "**NOTE** the genespace parameter object is returned or can be loaded
      into R via `load('%s', verbose = TRUE)`. Then you can customize your
      riparian plots by calling `plot_riparian(gsParam = gsParam, ...)`. The
      source data and ggplot2 objects are also stored in the /riparian
      directory and can also be accessed by `load(...)`. ",
      gpFile), indent = 0, exdent = 8),
    strwrap(
      "**NOTE** To query genespace results by position or gene,
      use `query_genespace(...)`. See specifications in ?query_genespace for
      details.",  indent = 0, exdent = 8),
    "############################",
    sep = "\n")
  save(gsParam, file = gpFile)
  return(gsParam)
}
