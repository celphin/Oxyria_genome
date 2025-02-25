##########################
# Running GeneSpace and orthoFinder on all genomes on beluga 
# August 2024
############################

# Paper: https://elifesciences.org/articles/78526 
# https://github.com/jtlovell/GENESPACE


###################################
# to run GeneSpace and setup all files

tmux new-session -s GeneSpace1
tmux attach-session -t GeneSpace1

module load StdEnv/2020 python/3.11.5 scipy-stack/2021a
module load StdEnv/2020 java/13.0.2
module load StdEnv/2020 r/4.2.2 glpk/5.0

export PATH=$PATH:/home/celphin/scratch/Oxyria/GeneSpace/OrthoFinder_source/orthofinder.py
alias orthofinder='python /home/celphin/scratch/Oxyria/GeneSpace/OrthoFinder_source/orthofinder.py'

R

library(GENESPACE)

#--------------------------------------
# Total
gpar <- init_genespace(
  wd = "/home/celphin/scratch/Oxyria/GeneSpace/Total_genomes", 
  path2mcscanx = "/home/celphin/scratch/Oxyria/GeneSpace/MCScanX/")
out <- run_genespace(gsParam = gpar) 

        # Arabidopsis_lyrata     : 29811 / 29811 geneIDs exactly match (PASS)
        # Arabidopsis_thaliana   : 17070 / 17070 geneIDs exactly match (PASS)
        # Arabis_alpina          : 21609 / 21609 geneIDs exactly match (PASS)
        # Argentina_anserina     : 19620 / 19620 geneIDs exactly match (PASS)
        # Brassica_oleracea      : 42491 / 42491 geneIDs exactly match (PASS)
        # Capsella_rubella       : 26410 / 26410 geneIDs exactly match (PASS)
		# Cochlearia_groenlandica: 31127 / 31127 geneIDs exactly match (PASS)
        # Draba_nivalis          : 33557 / 33557 geneIDs exactly match (PASS)
        # Dryas_octopetala       : 39696 / 39696 geneIDs exactly match (PASS)
        # Fagopyrum_escelentum_H2: 52149 / 52149 geneIDs exactly match (PASS)
        # Fagopyrum_tataricum_H1 : 39560 / 39560 geneIDs exactly match (PASS)
        # Fragaria_vesca         : 21397 / 21397 geneIDs exactly match (PASS)
        # Malus_sylvestris       : 37467 / 37467 geneIDs exactly match (PASS)
        # Oxyria_digyna_H1       : 33799 / 33799 geneIDs exactly match (PASS)
        # Polygunum_aviculare_H0 : 26201 / 26201 geneIDs exactly match (PASS)
        # Prunus_persica         : 23128 / 23128 geneIDs exactly match (PASS)
        # Pyrus_bretschneideri   : 35293 / 35293 geneIDs exactly match (PASS)
        # Rheum_nobile_H0        : 34698 / 34698 geneIDs exactly match (PASS)
        # Rheum_tangaticum_H0    : 30938 / 30938 geneIDs exactly match (PASS)
        # Rosa_rugosa            : 29146 / 29146 geneIDs exactly match (PASS)
        # Thlaspi_arvense        : 26392 / 26392 geneIDs exactly match (PASS)


# Could not find a valid path to the orthofinder program from R. To run         
# orthofinder, ensure that the orthofinder program is in the $PATH, 
# then call the following from the shell: 
# orthofinder -f /home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/tmp -t 16 -a 1 -X -o /home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/orthofinder

# Error in run_orthofinder(gsParam = gsParam, verbose = TRUE) :
  # Once OrthoFinder has been run, re-call run_genespace

#-------------------------------------
# Rosaceae

gpar <- init_genespace(
  wd = "/home/celphin/scratch/Oxyria/GeneSpace/Rosaceae_genomes", 
  path2mcscanx = "/home/celphin/scratch/Oxyria/GeneSpace/MCScanX/")
out <- run_genespace(gsParam = gpar)

# Checking annotation files (.bed and peptide .fa):
        # Argentina_anserina  : 19620 / 19620 geneIDs exactly match (PASS)
        # Dryas_octopetala    : 39696 / 39696 geneIDs exactly match (PASS)
        # Fragaria_vesca      : 21397 / 21397 geneIDs exactly match (PASS)
        # Malus_sylvestris    : 37467 / 37467 geneIDs exactly match (PASS)
        # Prunus_persica      : 23128 / 23128 geneIDs exactly match (PASS)
        # Pyrus_bretschneideri: 35293 / 35293 geneIDs exactly match (PASS)
        # Rosa_rugosa         : 29146 / 29146 geneIDs exactly match (PASS)


# Could not find a valid path to the orthofinder program from R. To run         orthofinder, ensure th
# at the orthofinder program is in the         $PATH, then call the following from the shell: 
orthofinder -f /home/celphin/scratch/Oxyria/GeneSpace/Rosaceae_genomes/tmp -t 16 -a 1 -X -o /home/celphin/scratch/Oxyria/GeneSpace/Rosaceae_genomes/orthofinderError in run_orthofinder(gsParam = gsParam, verbose = TRUE)
  # Once OrthoFinder has been run, re-call run_genespace

#----------------------------------
# Brassicaceae

gpar <- init_genespace(
  wd = "/home/celphin/scratch/Oxyria/GeneSpace/Brassicaceae_genomes", 
  path2mcscanx = "/home/celphin/scratch/Oxyria/GeneSpace/MCScanX/")
out <- run_genespace(gsParam = gpar)

# Checking annotation files (.bed and peptide .fa):
        # Arabidopsis_lyrata     : 29811 / 29811 geneIDs exactly match (PASS)
        # Arabidopsis_thaliana   : 10250 / 10250 geneIDs exactly match (PASS)
        # Arabis_alpina          : 21609 / 21609 geneIDs exactly match (PASS)
        # Brassica_oleracea      : 42491 / 42491 geneIDs exactly match (PASS)
        # Cochlearia_groenlandica: 31127 / 31127 geneIDs exactly match (PASS)
        # Draba_nivalis          : 33557 / 33557 geneIDs exactly match (PASS)
        # Thlaspi_arvense        : 26392 / 26392 geneIDs exactly match (PASS)


# Could not find a valid path to the orthofinder program from R. To run         orthofinder, ensure th
# at the orthofinder program is in the         $PATH, then call the following from the shell: 
orthofinder -f /home/celphin/scratch/Oxyria/GeneSpace/Brassicaceae_genomes/tmp -t 16 -a 1 -X -o /home/celphin/scratch/Oxyria/GeneSpace/Brassicaceae_genomes/orthofinder
# Error in run_orthofinder(gsParam = gsParam, verbose = TRUE)
  # Once OrthoFinder has been run, re-call run_genespace

#-----------------------------------
# Polygonaceae

gpar <- init_genespace(
  wd = "/home/celphin/scratch/Oxyria/GeneSpace/Polygonaceae_genomes", 
  path2mcscanx = "/home/celphin/scratch/Oxyria/GeneSpace/MCScanX/")
out <- run_genespace(gsParam = gpar)

# Checking annotation files (.bed and peptide .fa):
        # Fagopyrum_escelentum_H2: 52149 / 52149 geneIDs exactly match (PASS)
        # Fagopyrum_tataricum_H1 : 39560 / 39560 geneIDs exactly match (PASS)
        # Oxyria_digyna_H1       : 33799 / 33799 geneIDs exactly match (PASS)
        # Polygunum_aviculare_H0 : 26201 / 26201 geneIDs exactly match (PASS)
        # Rheum_nobile_H0        : 34698 / 34698 geneIDs exactly match (PASS)
        # Rheum_tangaticum_H0    : 30938 / 30938 geneIDs exactly match (PASS)

# Could not find a valid path to the orthofinder program from R. To run         orthofinder, ensure th
# at the orthofinder program is in the         $PATH, then call the following from the shell: 
orthofinder -f /home/celphin/scratch/Oxyria/GeneSpace/Polygonaceae_genomes/tmp -t 16 -a 1 -X -o /home/celphin/scratch/Oxyria/GeneSpace/Polygonaceae_genomes/orthofinder
# Error in run_orthofinder(gsParam = gsParam, verbose = TRUE)
  # Once OrthoFinder has been run, re-call run_genespace

#-----------------
# Dryas
gpar <- init_genespace(
  wd = "/home/celphin/scratch/Oxyria/GeneSpace/Dryas_genomes", 
  path2mcscanx = "/home/celphin/scratch/Oxyria/GeneSpace/MCScanX/")
out <- run_genespace(gsParam = gpar)

# Checking annotation files (.bed and peptide .fa):
        # Dry_ajan            : 34179 / 34179 geneIDs exactly match (PASS)
        # Dry_alask           : 34226 / 34226 geneIDs exactly match (PASS)
        # Dry_drumm           : 34456 / 34456 geneIDs exactly match (PASS)
        # Dry_int             : 34097 / 34097 geneIDs exactly match (PASS)
        # Dry_octo_H0         : 39696 / 39696 geneIDs exactly match (PASS)
        # Dry_octo_H1         : 37395 / 37395 geneIDs exactly match (PASS)
        # Malus_sylvestris    : 37467 / 37467 geneIDs exactly match (PASS)
        # Pyrus_bretschneideri: 35293 / 35293 geneIDs exactly match (PASS)
        # Rosa_rugosa         : 29146 / 29146 geneIDs exactly match (PASS)


# Could not find a valid path to the orthofinder program from R. To run         
# orthofinder, ensure that the orthofinder program is in the         
# $PATH, then call the following from the shell: 
orthofinder -f /home/celphin/scratch/Oxyria/GeneSpace/Dryas_genomes/tmp -t 16 -a 1 -X -o /home/celphin/scratch/Oxyria/GeneSpace/Dryas_genomes/orthofinder
# Error in run_orthofinder(gsParam = gsParam, verbose = TRUE) :

#-----------------
# Oxyria
gpar <- init_genespace(
  wd = "/home/celphin/scratch/Oxyria/GeneSpace/Oxyria_genomes", 
  path2mcscanx = "/home/celphin/scratch/Oxyria/GeneSpace/MCScanX/")
out <- run_genespace(gsParam = gpar)

# Checking annotation files (.bed and peptide .fa):
        # DToL_h1               : 31523 / 31523 geneIDs exactly match (PASS)
        # DToL_h2               : 35477 / 35477 geneIDs exactly match (PASS)
        # Fagopyrum_tataricum_H1: 39560 / 39560 geneIDs exactly match (PASS)
        # Oxy_Elles_Hap1        : 38198 / 38198 geneIDs exactly match (PASS)
        # Oxy_Elles_Hap2        : 35196 / 35196 geneIDs exactly match (PASS)
        # Oxy_Sval_h1           : 33402 / 33402 geneIDs exactly match (PASS)
        # Oxy_Sval_h2           : 31714 / 31714 geneIDs exactly match (PASS)
        # Polygunum_aviculare_H0: 26201 / 26201 geneIDs exactly match (PASS)
        # Rheum_nobile_H0       : 34698 / 34698 geneIDs exactly match (PASS

# Could not find a valid path to the orthofinder program from R. To run         
# orthofinder, ensure that the orthofinder program is in the         
# $PATH, then call the following from the shell: 
orthofinder -f /home/celphin/scratch/Oxyria/GeneSpace/Oxyria_genomes/tmp -t 16 -a 1 -X -o /home/celphin/scratch/Oxyria/GeneSpace/Oxyria_genomes/orthofinder
#Error in run_orthofinder(gsParam = gsParam, verbose = TRUE) :


######################################
# Run orthoFinder outside GeneSpace

# bash
salloc -c40 --time 5:00:00 --mem 191000M --account def-rieseber

module load StdEnv/2020  python/3.10.2 scipy-stack/2021a diamond/2.1.7 fastme/2.1.6.2 gcc/9.3.0 mcl/14.137 java/13.0.2 r/4.2.2 glpk/5.0

# other version of orthofinder
export PATH=$PATH:/home/celphin/scratch/Oxyria/GeneSpace/OrthoFinder/bin
export PATH=$PATH:/home/celphin/scratch/Oxyria/GeneSpace/OrthoFinder_source
alias orthofinder='orthofinder.py'

orthofinder -h

#-------------------------------
# Total genomes
cd /home/celphin/scratch/Oxyria/GeneSpace/Total_genomes

orthofinder -f /home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/tmp \
-t 40 -a 1 -X -o /home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/orthofinder

# START: /scratch (user celphin)            5461G/20T           634k/1000k
# END: /scratch (user celphin)            5470G/20T           775k/1000k

#--------------------------------------------
# Rosaceae
cd /home/celphin/scratch/Oxyria/GeneSpace/Rosaceae_genomes/

orthofinder -f /home/celphin/scratch/Oxyria/GeneSpace/Rosaceae_genomes/tmp \
-t 40 -a 1 -X -o /home/celphin/scratch/Oxyria/GeneSpace/Rosaceae_genomes/orthofinder

# Start: /scratch (user celphin)            5456G/20T           504k/1000k
# End: /scratch (user celphin)            5461G/20T           634k/1000k

#-----------------------------------
# Brassicaceae
cd /home/celphin/scratch/Oxyria/GeneSpace/Brassicaceae_genomes

orthofinder -f /home/celphin/scratch/Oxyria/GeneSpace/Brassicaceae_genomes/tmp \
-t 40 -a 1 -X -o /home/celphin/scratch/Oxyria/GeneSpace/Brassicaceae_genomes/orthofinder


#--------------------------------
# Polygonaceae

cd /home/celphin/scratch/Oxyria/GeneSpace/Polygonaceae_genomes

orthofinder -f /home/celphin/scratch/Oxyria/GeneSpace/Polygonaceae_genomes/tmp \
-t 40 -a 1 -X -o /home/celphin/scratch/Oxyria/GeneSpace/Polygonaceae_genomes/orthofinder

# Done: /home/celphin/scratch/Oxyria/GeneSpace/Polygonaceae_genomes/orthofinder/Results_Oct07/

#-------------------------
# Dryas

cd /home/celphin/scratch/Oxyria/GeneSpace/Dryas_genomes

orthofinder -f /home/celphin/scratch/Oxyria/GeneSpace/Dryas_genomes/tmp \
-t 40 -a 1 -X -o /home/celphin/scratch/Oxyria/GeneSpace/Dryas_genomes/orthofinder

#-------------------------
# Oxyria

cd /home/celphin/scratch/Oxyria/GeneSpace/Oxyria_genomes

orthofinder -f /home/celphin/scratch/Oxyria/GeneSpace/Oxyria_genomes/tmp \
-t 40 -a 1 -X -o /home/celphin/scratch/Oxyria/GeneSpace/Oxyria_genomes/orthofinder

# takes too much space - did not finish

#######################################
# Exploring results

# https://davidemms.github.io/orthofinder_tutorials/exploring-orthofinders-results.html
# download tree and view in http://etetoolkit.org/treeview/
more /home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/orthofinder/Results_Aug18/Species_Tree/SpeciesTree_rooted_node_labels.txt

# (((Polygunum_aviculare_H0:0.156996,(Oxyria_digyna_H1:0.136637,(Rheum_tangaticum_H0:0.0440109,Rheum_
# nobile_H0:0.0354046)N12:0.0630934)N7:0.0341106)N3:0.0347272,(Fagopyrum_tataricum_H1:0.0360394,Fagop
# yrum_escelentum_H2:0.037885)N4:0.103084)N1:0.093549,((((Capsella_rubella:0.0516886,(Arabidopsis_lyr
# ata:0.0229793,Arabidopsis_thaliana:0.077621)N17:0.02683)N13:0.0314861,(Draba_nivalis:0.0862001,Arab
# is_alpina:0.0831463)N14:0.0326396)N8:0.0158126,(Brassica_oleracea:0.0837281,Thlaspi_arvense:0.10468
# 4)N9:0.0201726)N5:0.215788,((Rosa_rugosa:0.0480735,(Argentina_anserina:0.0780935,Fragaria_vesca:0.0
# 594828)N15:0.0187218)N10:0.0798625,(Dryas_octopetala:0.143586,(Prunus_persica:0.0789428,(Malus_sylv
# estris:0.0172918,Pyrus_bretschneideri:0.0191209)N18:0.0817642)N16:0.0315978)N11:0.029253)N6:0.11457
# 8)N2:0.093549)N0;

more /home/celphin/scratch/Oxyria/GeneSpace/Rosaceae_genomes/orthofinder/Results_Aug16/Species_Tree/SpeciesTree_rooted_node_labels.txt

more /home/celphin/scratch/Oxyria/GeneSpace/Polygonaceae_genomes/orthofinder/Results_Oct07/Species_Tree/SpeciesTree_rooted_node_labels.txt
# ((Fagopyrum_tataricum_H1:0.0467894,Fagopyrum_escelentum_H2:0.0465864)N1:0.06654,(Polygunum_avicular
# e_H0:0.17436,(Oxyria_digyna_H1:0.156198,(Rheum_nobile_H0:0.0430757,Rheum_tangaticum_H0:0.052304)N4:
# 0.0649788)N3:0.0379018)N2:0.06654)N0;

more /home/celphin/scratch/Oxyria/GeneSpace/Brassicaceae_genomes/orthofinder/Results_Aug18/Species_Tree/SpeciesTree_rooted_node_labels.txt
# (Brassica_oleracea:0.135322,(Thlaspi_arvense:0.107401,(Cochlearia_groenlandica:0.109195,((Arabis_al
# pina:0.0810288,Draba_nivalis:0.0660331)N4:0.0424301,(Arabidopsis_thaliana:0.0527601,Arabidopsis_lyr
# ata:0.0115863)N5:0.0546759)N3:0.0278312)N2:0.051759)N1:0.135322)N0;

more /home/celphin/scratch/Oxyria/GeneSpace/Dryas_genomes/orthofinder/Results_Oct30/Species_Tree/SpeciesTree_rooted_node_labels.txt
((Rosa_rugosa:0.172356,(Pyrus_bretschneideri:0.0265023,Malus_sylvestris:0.0235352)N3:0.135
978)N1:0.073315,(Dry_drumm:0.0311242,((Dry_octo_H1:0.0116492,Dry_octo_H0:6.90812e-07)N5:0.
00886338,(Dry_int:0.027865,(Dry_ajan:0.0271139,Dry_alask:0.0293075)N7:0.0188818)N6:0.01821
96)N4:0.0130102)N2:0.073315)N0;


#################################
# Explore orthofinder results
# https://davidemms.github.io/orthofinder_tutorials/exploring-orthofinders-results.html

more  /lustre04/scratch/celphin/Oxyria/GeneSpace/Brassicaceae_genomes/orthofinder/Results_Oct09/Comparative_Genomics_Statistics/Statistics_Overall.tsv
# Number of species       7
# Number of genes 195237
# Number of genes in orthogroups  183039
# Number of unassigned genes      12198
# Percentage of genes in orthogroups      93.8
# Percentage of unassigned genes  6.2
# Number of orthogroups   23464
# Number of species-specific orthogroups  2856
# Number of genes in species-specific orthogroups 13643
# Percentage of genes in species-specific orthogroups     7.0
# Mean orthogroup size    7.8
# Median orthogroup size  7.0
# G50 (assigned genes)    8
# G50 (all genes) 8
# O50 (assigned genes)    6536
# O50 (all genes) 7299
# Number of orthogroups with all species present  4687
# Number of single-copy orthogroups       1061
# Date    2024-10-09
# Orthogroups file        Orthogroups.tsv
# Unassigned genes file   Orthogroups_UnassignedGenes.tsv
# Per-species statistics  Statistics_PerSpecies.tsv
# Overall statistics      Statistics_Overall.tsv
# Orthogroups shared between species      Orthogroups_SpeciesOverlaps.tsv
# Number of species in orthogroup Number of orthogroups
# 1       2856
# 2       1531
# 3       1006
# 4       1705
# 5       4198
# 6       7481
# 7       4687


more  /lustre04/scratch/celphin/Oxyria/GeneSpace/Brassicaceae_genomes/orthofinder/Results_Oct09/Comparative_Genomics_Statistics/Statistics_PerSpecies.tsv
        # Arabidopsis_lyrata      Arabidopsis_thaliana    Arabis_alpina   Brassica_oleracea       Cochlearia_groenlandica    Draba_nivalis   Thlaspi_arvense
# Number of genes 29811   10250   21609   42491   31127   33557   26392
# Number of genes in orthogroups  29068   10213   20323   39752   28383   30406   24894
# Number of unassigned genes      743     37      1286    2739    2744    3151    1498
# Percentage of genes in orthogroups      97.5    99.6    94.0    93.6    91.2    90.6    94.3
# Percentage of unassigned genes  2.5     0.4     6.0     6.4     8.8     9.4     5.7
# Number of orthogroups containing species        19516   8286    14229   19240   17534   18855   16781
# Percentage of orthogroups containing species    83.2    35.3    60.6    82.0    74.7    80.4    71.5
# Number of species-specific orthogroups  262     5       212     492     901     610     374
# Number of genes in species-specific orthogroups 1238    14      1622    2669    3993    2845    1262
# Percentage of genes in species-specific orthogroups     4.2     0.1     7.5     6.3     12.8    8.5        4.8

more  /lustre04/scratch/celphin/Oxyria/GeneSpace/Polygonaceae_genomes/orthofinder/Results_Oct07/Comparative_Genomics_Statistics/Statistics_Overall.tsv
# Number of species       6
# Number of genes 217345
# Number of genes in orthogroups  199612
# Number of unassigned genes      17733
# Percentage of genes in orthogroups      91.8
# Percentage of unassigned genes  8.2
# Number of orthogroups   28237
# Number of species-specific orthogroups  5383
# Number of genes in species-specific orthogroups 29056
# Percentage of genes in species-specific orthogroups     13.4
# Mean orthogroup size    7.1
# Median orthogroup size  6.0
# G50 (assigned genes)    8
# G50 (all genes) 7
# O50 (assigned genes)    6534
# O50 (all genes) 7664
# Number of orthogroups with all species present  11683
# Number of single-copy orthogroups       4168
# Date    2024-10-07
# Orthogroups file        Orthogroups.tsv
# Unassigned genes file   Orthogroups_UnassignedGenes.tsv
# Per-species statistics  Statistics_PerSpecies.tsv
# Overall statistics      Statistics_Overall.tsv
# Orthogroups shared between species      Orthogroups_SpeciesOverlaps.tsv
# Number of species in orthogroup Number of orthogroups
# 1       5383
# 2       4810
# 3       1910
# 4       1554
# 5       2897
# 6       11683

more  /lustre04/scratch/celphin/Oxyria/GeneSpace/Polygonaceae_genomes/orthofinder/Results_Oct07/Comparative_Genomics_Statistics/Statistics_PerSpecies.tsv

                    # Fagopyrum_escelentum_H2           Fagopyrum_tataricum_H1          Oxyria_digyna_H1        Polygunum_aviculare_H0     Rheum_nobile_H0     Rheum_tangaticum_H0
# Number of genes                      52149             39560                          33799                    26201                        34698                30938
# Number of genes in orthogroups        47768                35547                      30486                    24021                        32321                29469
# Number of unassigned genes             4381                4013                        3313                    2180                           2377               1469
# Percentage of genes in orthogroups      91.6               89.9                         90.2                   91.7                           93.1                95.3
# Percentage of unassigned genes            8.4               10.1                        9.8                    8.3                             6.9                4.7
# Number of orthogroups containing species  21007             20312                       16838                  16167                           19337               17871
# Percentage of orthogroups containing species    74.4        71.9                        59.6                     57.3                           68.5                 63.3
# Number of species-specific orthogroups        1780            949                       722                      256                             921                  755
# Number of genes in species-specific orthogroups 13334         3460                      4514                      1299                           3699                  2750
# Percentage of genes in species-specific orthogroups  25.6      8.7                      13.4                    5.0                                 10.7                8.9

more  /lustre04/scratch/celphin/Oxyria/GeneSpace/Rosaceae_genomes/orthofinder/Results_Aug16/Comparative_Genomics_Statistics/Statistics_Overall.tsv
# Number of species       7
# Number of genes 205747
# Number of genes in orthogroups  197391
# Number of unassigned genes      8356
# Percentage of genes in orthogroups      95.9
# Percentage of unassigned genes  4.1
# Number of orthogroups   22986
# Number of species-specific orthogroups  2513
# Number of genes in species-specific orthogroups 15660
# Percentage of genes in species-specific orthogroups     7.6
# Mean orthogroup size    8.6
# Median orthogroup size  8.0
# G50 (assigned genes)    9
# G50 (all genes) 9
# O50 (assigned genes)    6883
# O50 (all genes) 7347
# Number of orthogroups with all species present  12878
# Number of single-copy orthogroups       2493
# Date    2024-08-16
# Orthogroups file        Orthogroups.tsv
# Unassigned genes file   Orthogroups_UnassignedGenes.tsv
# Per-species statistics  Statistics_PerSpecies.tsv
# Overall statistics      Statistics_Overall.tsv
# Orthogroups shared between species      Orthogroups_SpeciesOverlaps.tsv
# Number of species in orthogroup Number of orthogroups
# 1       2513
# 2       2279
# 3       1098
# 4       827
# 5       782
# 6       2609
# 7       12878

more  /lustre04/scratch/celphin/Oxyria/GeneSpace/Rosaceae_genomes/orthofinder/Results_Aug16/Comparative_Genomics_Statistics/Statistics_PerSpecies.tsv
#                                            Argentina_anserina      Dryas_octopetala        Fragaria_vesca  Malus_sylvestris        Prunus_persica     Pyrus_bretschneideri    Rosa_rugosa
# Number of genes                                        19620   39696   21397   37467   23128   35293   29146
# Number of genes in orthogroups                         19405   36345   21045   35867   22664   33962   28103
# Number of unassigned genes                              215     3351    352     1600    464     1331    1043
# Percentage of genes in orthogroups                      98.9    91.6    98.4    95.7    98.0    96.2    96.4
# Percentage of unassigned genes                          1.1     8.4     1.6     4.3     2.0     3.8     3.6
# Number of orthogroups containing species               16313   18261   16299   18595   17489   18163   18263
# Percentage of orthogroups containing species            71.0    79.4    70.9    80.9    76.1    79.0    79.5
# Number of species-specific orthogroups                   46      1243    96      313     174     155     486
# Number of genes in species-specific orthogroups         213     10810   320     1391    768     461     1697
# Percentage of genes in species-specific orthogroups     1.1     27.2    1.5     3.7     3.3     1.3        5.8

more /lustre04/scratch/celphin/Oxyria/GeneSpace/Total_genomes/orthofinder/Results_Aug18/Comparative_Genomics_Statistics/Statistics_PerSpecies.tsv





more  /lustre04/scratch/celphin/Oxyria/GeneSpace/Total_genomes/orthofinder/Results_Aug18/Comparative_Genomics_Statistics/Statistics_Overall.tsv
# Number of species       21
# Number of genes 651559
# Number of genes in orthogroups  624240
# Number of unassigned genes      27319
# Percentage of genes in orthogroups      95.8
# Percentage of unassigned genes  4.2
# Number of orthogroups   33157
# Number of species-specific orthogroups  8168
# Number of genes in species-specific orthogroups 43928
# Percentage of genes in species-specific orthogroups     6.7
# Mean orthogroup size    18.8
# Median orthogroup size  8.0
# G50 (assigned genes)    35
# G50 (all genes) 33
# O50 (assigned genes)    4642
# O50 (all genes) 5046
# Number of orthogroups with all species present  3351
# Number of single-copy orthogroups       84
# Date    2024-08-18
# Orthogroups file        Orthogroups.tsv
# Unassigned genes file   Orthogroups_UnassignedGenes.tsv
# Per-species statistics  Statistics_PerSpecies.tsv
# Overall statistics      Statistics_Overall.tsv
# Orthogroups shared between species      Orthogroups_SpeciesOverlaps.tsv

# Average number of genes per-species in orthogroup       Number of orthogroups   Percentage
 # of orthogroups Number of genes Percentage of genes
# <1      21359   64.4    125605  20.1
# '1      8461    25.5    235694  37.8
# '2      1909    5.8     95984   15.4
# '3      684     2.1     48542   7.8
# '4      274     0.8     25362   4.1
# '5      116     0.3     13114   2.1
# '6      84      0.3     11409   1.8
# '7      69      0.2     10825   1.7
# '8      35      0.1     6208    1.0
# '9      31      0.1     6136    1.0
# '10     27      0.1     5901    0.9
# 11-15   68      0.2     18417   3.0
# 16-20   21      0.1     8014    1.3
# 21-50   17      0.1     9868    1.6
# 51-100  2       0.0     3161    0.5
# 101-150 0       0.0     0       0.0
# 151-200 0       0.0     0       0.0
# 201-500 0       0.0     0       0.0
# 501-1000        0       0.0     0       0.0
# '1001+  0       0.0     0       0.0

# Number of species in orthogroup Number of orthogroups
# 1       8168
# 2       5308
# 3       2034
# 4       1097
# 5       949
# 6       1227
# 7       1241
# 8       568
# 9       135
# 10      138
# 11      153
# 12      211
# 13      358
# 14      305
# 15      275
# 16      224
# 17      431
# 18      791
# 19      1956
# 20      4237
# 21      3351



####################################
# Tar orthofinder directories

tmux new-session -s tar_files
tmux attach-session -t tar_files

cd /home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/orthofinder

tar -zcvf Total_Results_Aug18.tar.gz ./Results_Aug18
rm -r Results_Aug18/

#-------------------
cd /home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/orthofinder_old

tar -zcvf Results_Jul03.tar.gz ./Results_Jul03
rm -r Results_Jul03

tar -zcvf Results_Jul04.tar.gz ./Results_Jul04
rm -r Results_Jul04

#-----------------------------
tmux new-session -s tar_files1
tmux attach-session -t tar_files1

cd /home/celphin/scratch/Oxyria/GeneSpace/Rosaceae_genomes/orthofinder

tar -zcvf Rosaceae_Results_Aug16.tar.gz ./Results_Aug16
rm -r Results_Aug16

#-----------------------------
tmux new-session -s tar_files2
tmux attach-session -t tar_files2

cd /home/celphin/scratch/Oxyria/GeneSpace/Brassicaceae_genomes/orthofinder

tar -zcvf Brassicaceae_Results_Aug18.tar.gz ./Results_Aug18
rm -r Results_Aug18

tar -zcvf Brassicaceae_Results_Oct09.tar.gz ./Results_Oct09
rm -r Results_Oct09

#-----------------------------
tmux new-session -s tar_files
tmux attach-session -t tar_files

cd /home/celphin/scratch/Oxyria/GeneSpace/Polygonaceae_genomes/orthofinder

tar -zcvf Polygonaceae_Results_Oct07.tar.gz ./Results_Oct07
rm -r Results_Oct07

#-----------------------------
tmux new-session -s tar_files
tmux attach-session -t tar_files

cd /home/celphin/scratch/Oxyria/GeneSpace/Dryas_genomes/orthofinder

tar -zcvf Dryas_Results_Oct29.tar.gz ./Results_Oct30
rm -r Results_Oct30

#-----------------------------
# to untar
cd /home/celphin/scratch/Oxyria/GeneSpace/Rosaceae_genomes/orthofinder
tar -xvzf Rosaceae_Results_Aug16.tar.gz

cd /home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/orthofinder
tar -xvzf Total_Results_Aug18.tar.gz

cd /home/celphin/scratch/Oxyria/GeneSpace/Polygonaceae_genomes/orthofinder
tar -xvzf Polygonaceae_Results_Oct07.tar.gz

cd /home/celphin/scratch/Oxyria/GeneSpace/Brassicaceae_genomes/orthofinder
tar -xvzf Brassicaceae_Results_Oct09.tar.gz

###################################

      Dry_ajan        Dry_alask       Dry_drumm       Dry_int Dry_octo_H0     Dry_octo_H1     Malus_sylvestris        Pyrus_bretschneideri    Rosa_rugosa
Number of genes 34179   34226   34456   34097   39696   37395   37467   35293   29146
Number of genes in orthogroups  32973   32986   32981   32978   39412   36482   35993   34093   27682
Number of unassigned genes      1206    1240    1475    1119    284     913     1474    1200    1464
Percentage of genes in orthogroups      96.5    96.4    95.7    96.7    99.3    97.6    96.1    96.6    95.0
Percentage of unassigned genes  3.5     3.6     4.3     3.3     0.7     2.4     3.9     3.4     5.0
Number of orthogroups containing species        26535   26470   24922   26614   26800   25491   18442   18016   17401
Percentage of orthogroups containing species    76.0    75.8    71.3    76.2    76.7    73.0    52.8    51.6    49.8
Number of species-specific orthogroups  44      50      57      43      39      22      310     144     795
Number of genes in species-specific orthogroups 102     102     143     99      506     97      1597    399     3739
Percentage of genes in species-specific orthogroups     0.3     0.3     0.4     0.3     1.3     0.3     4.3     1.1     12.8

