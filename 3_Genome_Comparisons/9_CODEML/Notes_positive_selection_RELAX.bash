#########################
# Running HyPhy programs
# Jan 2026
############################

# follow instructions here
# https://github.com/siribi/RELAX-workflow

# other tutorial: https://bioinformaticsworkbook.org/dataAnalysis/ComparativeGenomics/Finding_Positive_Selection_With_Codeml.html#gsc.tab=0

###############################
# Steps
# Run RELAX and absrel
# Explore results

######################

# Load modules
module load cufflinks/2.2.1 
module load cdbfasta/2017-03-16
module load orthofinder/2.5.2 
module load diamond/2.0.4 
module load clustalo/1.2.4
module load pal2nal.v14 
module load paml/4.9h

##########################
# label Arctic branches as foreground in species tree

((A,B{foreground}),C,(D,E{foreground}));

cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Species_Tree

# Node labels
(((Fagopyrum_escelentum_H2:0.0377569,Fagopyrum_tataricum_H1:0.0354292)N3:0.10263
9,(Polygunum_aviculare_H0:0.156221,((Rheum_nobile_H0:0.035189,Rheum_tangaticum_H
0:0.0440871)N11:0.0634952,Oxyria_digyna_H1:0.136143)N7:0.0338668)N4:0.0351737)N1
:0.093076,(((Rosa_rugosa:0.0478462,(Fragaria_vesca:0.0593422,Argentina_anserina:
0.076546)N12:0.0183635)N8:0.0789374,(Dryas_octopetala:0.142376,(Prunus_persica:0
.0787626,(Malus_sylvestris:0.0169615,Pyrus_bretschneideri:0.0188564)N16:0.080969
9)N13:0.0316866)N9:0.0295128)N5:0.115343,(Cochlearia_groenlandica:0.121011,((((A
rabidopsis_lyrata:0.022788,Arabidopsis_thaliana:0.0793045)N19:0.0265357,Capsella
_rubella:0.0514998)N17:0.0309554,(Arabis_alpina:0.0831942,Draba_nivalis:0.085671
7)N18:0.0322681)N14:0.0134594,(Thlaspi_arvense:0.102948,Brassica_oleracea:0.0829
899)N15:0.0188365)N10:0.0265485)N6:0.207906)N2:0.093076)N0;

# No node labels
(((Fagopyrum_escelentum_H2:0.0377569,Fagopyrum_tataricum_H1:0.0354292)0.973441:0
.102639,(Polygunum_aviculare_H0:0.156221,((Rheum_nobile_H0:0.035189,Rheum_tangat
icum_H0:0.0440871)0.880931:0.0634952,Oxyria_digyna_H1:0.136143)0.593853:0.033866
8)0.56789:0.0351737)0.966876:0.093076,(((Rosa_rugosa:0.0478462,(Fragaria_vesca:0
.0593422,Argentina_anserina:0.076546)0.504625:0.0183635)0.958221:0.0789374,(Drya
s_octopetala:0.142376,(Prunus_persica:0.0787626,(Malus_sylvestris:0.0169615,Pyru
s_bretschneideri:0.0188564)0.982393:0.0809699)0.714115:0.0316866)0.626082:0.0295
128)0.961504:0.115343,(Cochlearia_groenlandica:0.121011,((((Arabidopsis_lyrata:0
.022788,Arabidopsis_thaliana:0.0793045)0.80185:0.0265357,Capsella_rubella:0.0514
998)0.797672:0.0309554,(Arabis_alpina:0.0831942,Draba_nivalis:0.0856717)0.719188
:0.0322681)0.183229:0.0134594,(Thlaspi_arvense:0.102948,Brassica_oleracea:0.0829
899)0.281409:0.0188365)0.289466:0.0265485)0.965085:0.207906)0.966876:0.093076);


#####################
# Run absrel

hyphy absrel \
  --alignment geneX.codon.aln \
  --tree geneX.treefile
  
  
# should I use the gene tree or the species tree here?








#######################
# RELAX info

###############################
# Prepare files and run RELAX



###############################
# Checking results

# Which orthogroups are demonstrating positive selection







