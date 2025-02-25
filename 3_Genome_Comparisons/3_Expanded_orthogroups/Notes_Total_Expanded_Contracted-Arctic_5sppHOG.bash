#################################################################################################################
# Sept 2024 - Arctic Species expanded gene families
# Tutorial https://github.com/elsemikk/Willisornis_Genome_Assembly/blob/master/4.2_Gene_Family_Expansions.md
#################################################################################################################
tmux new-session -s CAFE
tmux attach-session -t CAFE

cd /lustre04/scratch/celphin/Oxyria/CAFE/Total_genomes/5spp_analysis

head -n 3 ../5spp_output/error_model/Base_clade_results.txt

#get a list of the significantly expanded/contracted/rapid genes in Oxy/Rhuem/Dryas/Draba

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

head -n 3 ../5spp_output/error_model/Base_asr.tre

#nexus
# BEGIN TREES;
  # TREE N0.HOG0000110 = (((Fagopyrum_escelentum_H2<1>_13:9.53333,Fagopyrum_tataricum_H1<2>*_18:9.53333)<23>_14
# :32.9366,(Polygunum_aviculare_H0<3>_11:34.3138,((Rheum_nobile_H0<4>_10:9.53787,Rheum_tangaticum_H0<5>_10:9.53
# 787)<26>_10:16.9445,Oxyria_digyna_H1<6>_11:26.4824)<25>_11:7.83139)<24>_11:8.15616)<22>_11:32.5301,(((Rosa_ru
# gosa<7>*_14:16.0623,(Fragaria_vesca<8>*_4:12.9604,Argentina_anserina<9>_5:12.9604)<30>*_6:3.10191)<29>_9:18.1
# 878,(Dryas_octopetala<10>*_23:27.9781,(Prunus_persica<11>*_17:20.3912,(Malus_sylvestris<12>*_6:3.79389,Pyrus_
# bretschneideri<13>_11:3.79389)<33>*_10:16.5973)<32>_13:7.5869)<31>*_13:6.27197)<28>_11:25.7292,(Cochlearia_gr
# oenlandica<14>_3:26.8531,(((Arabidopsis_lyrata<15>_4:12.0917,Capsella_rubella<16>_3:12.0917)<37>_4:8.1991,(Ar
# abis_alpina<17>_4:15.3063,Draba_nivalis<18>_7:15.3063)<38>_5:4.98445)<36>_5:2.368,(Thlaspi_arvense<19>*_2:19.
# 0401,Brassica_oleracea<20>*_8:19.0401)<39>_5:3.61871)<35>_5:4.19438)<34>*_5:33.1261)<27>_9:15.0207)<21>_10;

##################################
# Repeat with 5 spp HOGs

# High Arctic
#get a list of the significantly expanded/contracted/rapid genes in Oxyria_digyna_H1<6>
printf '%s\n' {1..50} | while read number ; do grep "Oxyria_digyna_H1<6>\*" ../5spp_output/error_model_$number/Base_asr.tre | \
sed 's/^.*OG/OG/g' | sed 's/ = .*$//g' > Oxydig_sig_changes_error_model_"$number" ; done 

#get a list of the significantly expanded/contracted/rapid genes in Draba_nivalis<18>
printf '%s\n' {1..50} | while read number ; do grep "Draba_nivalis<18>\*" ../5spp_output/error_model_$number/Base_asr.tre | \
sed 's/^.*OG/OG/g' | sed 's/ = .*$//g' > Drabaniv_sig_changes_error_model_"$number" ; done 

#get a list of the significantly expanded/contracted/rapid genes in Dryas_octopetala<10>
printf '%s\n' {1..50} | while read number ; do grep "Dryas_octopetala<10>\*" ../5spp_output/error_model_$number/Base_asr.tre | \
sed 's/^.*OG/OG/g' | sed 's/ = .*$//g' > Dryasoct_sig_changes_error_model_"$number" ; done 

#get a list of the significantly expanded/contracted/rapid genes in Cochlearia_groenlandica<14>
printf '%s\n' {1..50} | while read number ; do grep "Cochlearia_groenlandica<14>\*" ../5spp_output/error_model_$number/Base_asr.tre | \
sed 's/^.*OG/OG/g' | sed 's/ = .*$//g' > Cochgro_sig_changes_error_model_"$number" ; done 

# Alpine
#get a list of the significantly expanded/contracted/rapid genes in Rheum_nobile_H0<4>
printf '%s\n' {1..50} | while read number ; do grep "Rheum_nobile_H0<4>\*" ../5spp_output/error_model_$number/Base_asr.tre | \
sed 's/^.*OG/OG/g' | sed 's/ = .*$//g' > Rheumnob_sig_changes_error_model_"$number" ; done 

#get a list of the significantly expanded/contracted/rapid genes in Arabis_alpina<17>
printf '%s\n' {1..50} | while read number ; do grep "Arabis_alpina<17>\*" ../5spp_output/error_model_$number/Base_asr.tre | \
sed 's/^.*OG/OG/g' | sed 's/ = .*$//g' > Arabalp_sig_changes_error_model_"$number" ; done 

# Clade
#get a list of the significantly expanded/contracted/rapid genes in Arabis/Draba clase <38>
printf '%s\n' {1..50} | while read number ; do grep "<38>\*" ../5spp_output/error_model_$number/Base_asr.tre | \
sed 's/^.*OG/OG/g' | sed 's/ = .*$//g' > Dra_Ara_sig_changes_error_model_"$number" ; done 

#--------------------------
#Data exploration:
wc -l Oxydig_sig_changes_error_model_* | head # 132
wc -l Rheumnob_sig_changes_error_model_* | head # 190
wc -l Dra_Ara_sig_changes_error_model_* | head # 89
wc -l Drabaniv_sig_changes_error_model_* | head # 171
wc -l Dryasoct_sig_changes_error_model_* | head # 175
wc -l Arabalp_sig_changes_error_model_* | head # 166
wc -l Cochgro_sig_changes_error_model_* | head # 145

#---------------------
#see what is uniquely found in run 6 (this was the time I had no Hypoch/Rhemel)
printf '%s\n' {1..10} | while read number ; do grep -v -f Oxydig_sig_changes_error_model_"$number" Oxydig_sig_changes_error_model_10 ; done
printf '%s\n' {1..10} | while read number ; do grep -v -f Oxydig_sig_changes_error_model_"$number" Oxydig_sig_changes_error_model_1 ; done
# many differences between runs in Oxyria

##################################################################################################################
#Get list of families that have expanded in Oxyria, whether or not significant
# value should be one greater than ID to get right column

printf '%s\n' {1..50} | while read number ; do cut -f 1,7 ../5spp_output/error_model_$number/Base_change.tab | \
awk -F'\t' '$2 != 0' | grep -v "-" > Oxydig_expandedfams_error_model_"$number" ; done 

printf '%s\n' {1..50} | while read number ; do cut -f 1,19 ../5spp_output/error_model_$number/Base_change.tab | \
awk -F'\t' '$2 != 0' | grep -v "-" > Drabaniv_expandedfams_error_model_"$number" ; done 

printf '%s\n' {1..50} | while read number ; do cut -f 1,11 ../5spp_output/error_model_$number/Base_change.tab | \
awk -F'\t' '$2 != 0' | grep -v "-" > Dryasoct_expandedfams_error_model_"$number" ; done 

printf '%s\n' {1..50} | while read number ; do cut -f 1,15 ../5spp_output/error_model_$number/Base_change.tab | \
awk -F'\t' '$2 != 0' | grep -v "-" > Cochgro_expandedfams_error_model_"$number" ; done 


printf '%s\n' {1..50} | while read number ; do cut -f 1,18 ../5spp_output/error_model_$number/Base_change.tab | \
awk -F'\t' '$2 != 0' | grep -v "-" > Arabalp_expandedfams_error_model_"$number" ; done 

printf '%s\n' {1..50} | while read number ; do cut -f 1,5 ../5spp_output/error_model_$number/Base_change.tab | \
awk -F'\t' '$2 != 0' | grep -v "-" > Rheumnob_expandedfams_error_model_"$number" ; done 


printf '%s\n' {1..50} | while read number ; do cut -f 1,39 ../5spp_output/error_model_$number/Base_change.tab | \
awk -F'\t' '$2 != 0' | grep -v "-" > Dra_Ara_expandedfams_error_model_"$number" ; done 


wc -l *expandedfams_error_model_1 #remember the header is one line!

   # 462 Arabalp_expandedfams_error_model_1
  # 2663 Cochgro_expandedfams_error_model_1
   # 137 Dra_Ara_expandedfams_error_model_1
  # 1488 Drabaniv_expandedfams_error_model_1
  # 1060 Dryasoct_expandedfams_error_model_1
  # 1231 Oxydig_expandedfams_error_model_1
  # 1599 Rheumnob_expandedfams_error_model_1

##################################################################################################################
#Get list of families that have contracted in Oxyria, whether or not significant
# value should be one greater than ID to get right column

printf '%s\n' {1..50} | while read number ; do cut -f 1,7 ../5spp_output/error_model_$number/Base_change.tab | \
awk -F'\t' '$2 != 0' | grep "-" > Oxydig_contractedfams_error_model_"$number" ; done 

printf '%s\n' {1..50} | while read number ; do cut -f 1,19 ../5spp_output/error_model_$number/Base_change.tab | \
awk -F'\t' '$2 != 0' | grep "-" > Drabaniv_contractedfams_error_model_"$number" ; done 

printf '%s\n' {1..50} | while read number ; do cut -f 1,11 ../5spp_output/error_model_$number/Base_change.tab | \
awk -F'\t' '$2 != 0' | grep "-" > Dryasoct_contractedfams_error_model_"$number" ; done 

printf '%s\n' {1..50} | while read number ; do cut -f 1,15 ../5spp_output/error_model_$number/Base_change.tab | \
awk -F'\t' '$2 != 0' | grep "-" > Cochgro_contractedfams_error_model_"$number" ; done 


printf '%s\n' {1..50} | while read number ; do cut -f 1,18 ../5spp_output/error_model_$number/Base_change.tab | \
awk -F'\t' '$2 != 0' | grep "-" > Arabalp_contractedfams_error_model_"$number" ; done 

printf '%s\n' {1..50} | while read number ; do cut -f 1,5 ../5spp_output/error_model_$number/Base_change.tab | \
awk -F'\t' '$2 != 0' | grep "-" > Rheumnob_contractedfams_error_model_"$number" ; done 


printf '%s\n' {1..50} | while read number ; do cut -f 1,39 ../5spp_output/error_model_$number/Base_change.tab | \
awk -F'\t' '$2 != 0' | grep "-" > Dra_Ara_contractedfams_error_model_"$number" ; done 


wc -l *contractedfams_error_model_1 #remember the header is one line!
  # 4075 Arabalp_contractedfams_error_model_1
  # 2820 Cochgro_contractedfams_error_model_1
   # 897 Dra_Ara_contractedfams_error_model_1
   # 855 Drabaniv_contractedfams_error_model_1
  # 1497 Dryasoct_contractedfams_error_model_1
  # 2557 Oxydig_contractedfams_error_model_1
   # 818 Rheumnob_contractedfams_error_model_1

##################################################################################################################
#Running python script:
module load StdEnv/2020 python/3.9.6 scipy-stack/2021a
python

import os
from os.path import join as jn
import numpy as np

folderpath= os.getcwd()

#############################
nameroot='Dra_Ara_sig_changes_error_model_'

files = [x for x in os.listdir(folderpath) if nameroot in x]
orthologs=''
superstring=''
for file in files:
    with open(jn(folderpath,file),'r') as fl:
        fl=fl.readlines()
        fl=fl[1:]
        for f in fl:
            superstring+=f.strip().split('\t')[0]+'\n'

print(superstring)
export= open(jn(folderpath,'output_Dra_Ara'),'w')
export.write(superstring)
export.close()
#----------------
nameroot='Oxydig_sig_changes_error_model_'
files = [x for x in os.listdir(folderpath) if nameroot in x]
orthologs=''
superstring=''
for file in files:
    with open(jn(folderpath,file),'r') as fl:
        fl=fl.readlines()
        fl=fl[1:]
        for f in fl:
            superstring+=f.strip().split('\t')[0]+'\n'

print(superstring)
export= open(jn(folderpath,'output_Oxydig'),'w')
export.write(superstring)
export.close()
#---------------------
nameroot='Rheumnob_sig_changes_error_model_'
files = [x for x in os.listdir(folderpath) if nameroot in x]
orthologs=''
superstring=''
for file in files:
    with open(jn(folderpath,file),'r') as fl:
        fl=fl.readlines()
        fl=fl[1:]
        for f in fl:
            superstring+=f.strip().split('\t')[0]+'\n'

print(superstring)
export= open(jn(folderpath,'output_Rheumnob'),'w')
export.write(superstring)
export.close()

#----------------
nameroot='Drabaniv_sig_changes_error_model_'
files = [x for x in os.listdir(folderpath) if nameroot in x]
orthologs=''
superstring=''
for file in files:
    with open(jn(folderpath,file),'r') as fl:
        fl=fl.readlines()
        fl=fl[1:]
        for f in fl:
            superstring+=f.strip().split('\t')[0]+'\n'

print(superstring)
export= open(jn(folderpath,'output_Drabaniv'),'w')
export.write(superstring)
export.close()
#---------------------
nameroot='Dryasoct_sig_changes_error_model_'
files = [x for x in os.listdir(folderpath) if nameroot in x]
orthologs=''
superstring=''
for file in files:
    with open(jn(folderpath,file),'r') as fl:
        fl=fl.readlines()
        fl=fl[1:]
        for f in fl:
            superstring+=f.strip().split('\t')[0]+'\n'

print(superstring)
export= open(jn(folderpath,'output_Dryasoct'),'w')
export.write(superstring)
export.close()

#---------------------
nameroot='Arabalp_sig_changes_error_model_'
files = [x for x in os.listdir(folderpath) if nameroot in x]
orthologs=''
superstring=''
for file in files:
    with open(jn(folderpath,file),'r') as fl:
        fl=fl.readlines()
        fl=fl[1:]
        for f in fl:
            superstring+=f.strip().split('\t')[0]+'\n'

print(superstring)
export= open(jn(folderpath,'output_Arabalp'),'w')
export.write(superstring)
export.close()

#---------------------
nameroot='Cochgro_sig_changes_error_model_'
files = [x for x in os.listdir(folderpath) if nameroot in x]
orthologs=''
superstring=''
for file in files:
    with open(jn(folderpath,file),'r') as fl:
        fl=fl.readlines()
        fl=fl[1:]
        for f in fl:
            superstring+=f.strip().split('\t')[0]+'\n'

print(superstring)
export= open(jn(folderpath,'output_Arabalp'),'w')
export.write(superstring)
export.close()

####################################
nameroot='Dra_Ara_expandedfams_error_model_'

files = [x for x in os.listdir(folderpath) if nameroot in x]
orthologs=''
superstring=''
for file in files:
    with open(jn(folderpath,file),'r') as fl:
        fl=fl.readlines()
        fl=fl[1:]
        for f in fl:
            superstring+=f.strip().split('\t')[0]+'\n'

print(superstring)
export= open(jn(folderpath,'output_Dra_Ara_expandedfams'),'w')
export.write(superstring)
export.close()

#----------------
nameroot='Oxydig_expandedfams_error_model_'
files = [x for x in os.listdir(folderpath) if nameroot in x]
orthologs=''
superstring=''
for file in files:
    with open(jn(folderpath,file),'r') as fl:
        fl=fl.readlines()
        fl=fl[1:]
        for f in fl:
            superstring+=f.strip().split('\t')[0]+'\n'

print(superstring)
export= open(jn(folderpath,'output_Oxydig_expandedfams'),'w')
export.write(superstring)
export.close()

#---------------------
nameroot='Rheumnob_expandedfams_error_model_'
files = [x for x in os.listdir(folderpath) if nameroot in x]
orthologs=''
superstring=''
for file in files:
    with open(jn(folderpath,file),'r') as fl:
        fl=fl.readlines()
        fl=fl[1:]
        for f in fl:
            superstring+=f.strip().split('\t')[0]+'\n'

print(superstring)
export= open(jn(folderpath,'output_Rheumnob_expandedfams'),'w')
export.write(superstring)
export.close()

#----------------
nameroot='Drabaniv_expandedfams_error_model_'
files = [x for x in os.listdir(folderpath) if nameroot in x]
orthologs=''
superstring=''
for file in files:
    with open(jn(folderpath,file),'r') as fl:
        fl=fl.readlines()
        fl=fl[1:]
        for f in fl:
            superstring+=f.strip().split('\t')[0]+'\n'

print(superstring)
export= open(jn(folderpath,'output_Drabaniv_expandedfams'),'w')
export.write(superstring)
export.close()

#---------------------
nameroot='Dryasoct_expandedfams_error_model_'
files = [x for x in os.listdir(folderpath) if nameroot in x]
orthologs=''
superstring=''
for file in files:
    with open(jn(folderpath,file),'r') as fl:
        fl=fl.readlines()
        fl=fl[1:]
        for f in fl:
            superstring+=f.strip().split('\t')[0]+'\n'

print(superstring)
export= open(jn(folderpath,'output_Dryasoct_expandedfams'),'w')
export.write(superstring)
export.close()

#---------------------
nameroot='Arabalp_expandedfams_error_model_'
files = [x for x in os.listdir(folderpath) if nameroot in x]
orthologs=''
superstring=''
for file in files:
    with open(jn(folderpath,file),'r') as fl:
        fl=fl.readlines()
        fl=fl[1:]
        for f in fl:
            superstring+=f.strip().split('\t')[0]+'\n'

print(superstring)
export= open(jn(folderpath,'output_Arabalp_expandedfams'),'w')
export.write(superstring)
export.close()

#---------------------
nameroot='Cochgro_expandedfams_error_model_'
files = [x for x in os.listdir(folderpath) if nameroot in x]
orthologs=''
superstring=''
for file in files:
    with open(jn(folderpath,file),'r') as fl:
        fl=fl.readlines()
        fl=fl[1:]
        for f in fl:
            superstring+=f.strip().split('\t')[0]+'\n'

print(superstring)
export= open(jn(folderpath,'output_Cochgro_expandedfams'),'w')
export.write(superstring)
export.close()

################################
nameroot='Dra_Ara_contractedfams_error_model_'
files = [x for x in os.listdir(folderpath) if nameroot in x]
orthologs=''
superstring=''
for file in files:
    with open(jn(folderpath,file),'r') as fl:
        fl=fl.readlines()
        fl=fl[1:]
        for f in fl:
            superstring+=f.strip().split('\t')[0]+'\n'

print(superstring)
export= open(jn(folderpath,'output_Dra_Ara_contractedfams'),'w')
export.write(superstring)
export.close()

#----------------
nameroot='Oxydig_contractedfams_error_model_'
files = [x for x in os.listdir(folderpath) if nameroot in x]
orthologs=''
superstring=''
for file in files:
    with open(jn(folderpath,file),'r') as fl:
        fl=fl.readlines()
        fl=fl[1:]
        for f in fl:
            superstring+=f.strip().split('\t')[0]+'\n'

print(superstring)
export= open(jn(folderpath,'output_Oxydig_contractedfams'),'w')
export.write(superstring)
export.close()

#---------------------
nameroot='Rheumnob_contractedfams_error_model_'
files = [x for x in os.listdir(folderpath) if nameroot in x]
orthologs=''
superstring=''
for file in files:
    with open(jn(folderpath,file),'r') as fl:
        fl=fl.readlines()
        fl=fl[1:]
        for f in fl:
            superstring+=f.strip().split('\t')[0]+'\n'

print(superstring)
export= open(jn(folderpath,'output_Rheumnob_contractedfams'),'w')
export.write(superstring)
export.close()

#----------------
nameroot='Drabaniv_contractedfams_error_model_'
files = [x for x in os.listdir(folderpath) if nameroot in x]
orthologs=''
superstring=''
for file in files:
    with open(jn(folderpath,file),'r') as fl:
        fl=fl.readlines()
        fl=fl[1:]
        for f in fl:
            superstring+=f.strip().split('\t')[0]+'\n'

print(superstring)
export= open(jn(folderpath,'output_Drabaniv_contractedfams'),'w')
export.write(superstring)
export.close()

#---------------------
nameroot='Dryasoct_contractedfams_error_model_'
files = [x for x in os.listdir(folderpath) if nameroot in x]
orthologs=''
superstring=''
for file in files:
    with open(jn(folderpath,file),'r') as fl:
        fl=fl.readlines()
        fl=fl[1:]
        for f in fl:
            superstring+=f.strip().split('\t')[0]+'\n'

print(superstring)
export= open(jn(folderpath,'output_Dryasoct_contractedfams'),'w')
export.write(superstring)
export.close()

#---------------------
nameroot='Arabalp_contractedfams_error_model_'
files = [x for x in os.listdir(folderpath) if nameroot in x]
orthologs=''
superstring=''
for file in files:
    with open(jn(folderpath,file),'r') as fl:
        fl=fl.readlines()
        fl=fl[1:]
        for f in fl:
            superstring+=f.strip().split('\t')[0]+'\n'

print(superstring)
export= open(jn(folderpath,'output_Arabalp_contractedfams'),'w')
export.write(superstring)
export.close()

#---------------------
nameroot='Cochgro_contractedfams_error_model_'
files = [x for x in os.listdir(folderpath) if nameroot in x]
orthologs=''
superstring=''
for file in files:
    with open(jn(folderpath,file),'r') as fl:
        fl=fl.readlines()
        fl=fl[1:]
        for f in fl:
            superstring+=f.strip().split('\t')[0]+'\n'

print(superstring)
export= open(jn(folderpath,'output_Cochgro_contractedfams'),'w')
export.write(superstring)
export.close()

exit()

##################################################################################################################
#Clean up:
mkdir sig_genes; mv *sig_changes* sig_genes
mkdir expandedfams_error_models; mv *expandedfams_error_model* expandedfams_error_models
mkdir contractedfams_error_models; mv *contractedfams_error_model* contractedfams_error_models

##################################################################################################################
#Oxyria
cat output_Oxydig_expandedfams | sort | uniq > OGtotaltemp
cat OGtotaltemp | while read line ; do grep -c "$line" output_Oxydig_expandedfams >> tempcounts ; done
paste OGtotaltemp tempcounts >> OGcountstemp
cat OGcountstemp | while read orthogroup counts ; do if [ "$counts" -gt 47 ] ; then echo "$orthogroup" >> Oxydig_expandedfams ; fi ; done
rm -f OGcountstemp OGtotaltemp tempcounts output


cat output_Oxydig_contractedfams | sort | uniq > OGtotaltemp
cat OGtotaltemp | while read line ; do grep -c "$line" output_Oxydig_contractedfams >> tempcounts ; done
paste OGtotaltemp tempcounts >> OGcountstemp
cat OGcountstemp | while read orthogroup counts ; do if [ "$counts" -gt 47 ] ; then echo "$orthogroup" >> Oxydig_contractedfams ; fi ; done
rm -f OGcountstemp OGtotaltemp tempcounts output


cat output_Oxydig | sort | uniq > OGtotaltemp
cat OGtotaltemp | while read line ; do grep -c "$line" output_Oxydig >> tempcounts ; done
paste OGtotaltemp tempcounts >> OGcountstemp
cat OGcountstemp | while read orthogroup counts ; do if [ "$counts" -gt 47 ] ; then echo "$orthogroup" >> Oxydig_sig_changes ; fi ; done
rm -f OGcountstemp OGtotaltemp tempcounts output

#-------------------------------------------
#Rheum nobile
cat output_Rheumnob_expandedfams | sort | uniq > OGtotaltemp
cat OGtotaltemp | while read line ; do grep -c "$line" output_Rheumnob_expandedfams >> tempcounts ; done
paste OGtotaltemp tempcounts >> OGcountstemp
cat OGcountstemp | while read orthogroup counts ; do if [ "$counts" -gt 47 ] ; then echo "$orthogroup" >> Rheumnob_expandedfams ; fi ; done
rm -f OGcountstemp OGtotaltemp tempcounts output 


cat output_Rheumnob_contractedfams | sort | uniq > OGtotaltemp
cat OGtotaltemp | while read line ; do grep -c "$line" output_Rheumnob_contractedfams >> tempcounts ; done
paste OGtotaltemp tempcounts >> OGcountstemp
cat OGcountstemp | while read orthogroup counts ; do if [ "$counts" -gt 47 ] ; then echo "$orthogroup" >> Rheumnob_contractedfams ; fi ; done
rm -f OGcountstemp OGtotaltemp tempcounts output 


cat output_Rheumnob | sort | uniq > OGtotaltemp
cat OGtotaltemp | while read line ; do grep -c "$line" output_Rheumnob >> tempcounts ; done
paste OGtotaltemp tempcounts >> OGcountstemp
cat OGcountstemp | while read orthogroup counts ; do if [ "$counts" -gt 47 ] ; then echo "$orthogroup" >>Rheumnob_sig_changes ; fi ; done
rm -f OGcountstemp OGtotaltemp tempcounts output

#-----------------------------------------------
#Draba/Arabis node
cat output_Dra_Ara_expandedfams | sort | uniq > OGtotaltemp
cat OGtotaltemp | while read line ; do grep -c "$line" output_Dra_Ara_expandedfams >> tempcounts ; done
paste OGtotaltemp tempcounts >> OGcountstemp
cat OGcountstemp | while read orthogroup counts ; do if [ "$counts" -gt 47 ] ; then echo "$orthogroup" >> Dra_Ara_expandedfams ; fi ; done
rm -f OGcountstemp OGtotaltemp tempcounts output #here curently


cat output_Dra_Ara_contractedfams | sort | uniq > OGtotaltemp
cat OGtotaltemp | while read line ; do grep -c "$line" output_Dra_Ara_contractedfams >> tempcounts ; done
paste OGtotaltemp tempcounts >> OGcountstemp
cat OGcountstemp | while read orthogroup counts ; do if [ "$counts" -gt 47 ] ; then echo "$orthogroup" >> Dra_Ara_contractedfams ; fi ; done
rm -f OGcountstemp OGtotaltemp tempcounts output


cat output_Dra_Ara | sort | uniq > OGtotaltemp
cat OGtotaltemp | while read line ; do grep -c "$line" output_Dra_Ara >> tempcounts ; done
paste OGtotaltemp tempcounts >> OGcountstemp
cat OGcountstemp | while read orthogroup counts ; do if [ "$counts" -gt 47 ] ; then echo "$orthogroup" >> Dra_Ara_sig_changes ; fi ; done
rm -f OGcountstemp OGtotaltemp tempcounts output

#--------------------------------------------------
# Draba nivalis
cat output_Drabaniv_expandedfams | sort | uniq > OGtotaltemp
cat OGtotaltemp | while read line ; do grep -c "$line" output_Drabaniv_expandedfams >> tempcounts ; done
paste OGtotaltemp tempcounts >> OGcountstemp
cat OGcountstemp | while read orthogroup counts ; do if [ "$counts" -gt 47 ] ; then echo "$orthogroup" >> Drabaniv_expandedfams ; fi ; done
rm -f OGcountstemp OGtotaltemp tempcounts output


cat output_Drabaniv_contractedfams | sort | uniq > OGtotaltemp
cat OGtotaltemp | while read line ; do grep -c "$line" output_Drabaniv_contractedfams >> tempcounts ; done
paste OGtotaltemp tempcounts >> OGcountstemp
cat OGcountstemp | while read orthogroup counts ; do if [ "$counts" -gt 47 ] ; then echo "$orthogroup" >> Drabaniv_contractedfams ; fi ; done
rm -f OGcountstemp OGtotaltemp tempcounts output


cat output_Drabaniv | sort | uniq > OGtotaltemp
cat OGtotaltemp | while read line ; do grep -c "$line" output_Drabaniv >> tempcounts ; done
paste OGtotaltemp tempcounts >> OGcountstemp
cat OGcountstemp | while read orthogroup counts ; do if [ "$counts" -gt 47 ] ; then echo "$orthogroup" >> Drabaniv_sig_changes ; fi ; done
rm -f OGcountstemp OGtotaltemp tempcounts output

#-------------------------------------------
# Dryas octopetala
cat output_Dryasoct_expandedfams | sort | uniq > OGtotaltemp
cat OGtotaltemp | while read line ; do grep -c "$line" output_Dryasoct_expandedfams >> tempcounts ; done
paste OGtotaltemp tempcounts >> OGcountstemp
cat OGcountstemp | while read orthogroup counts ; do if [ "$counts" -gt 47 ] ; then echo "$orthogroup" >> Dryasoct_expandedfams ; fi ; done
rm -f OGcountstemp OGtotaltemp tempcounts output 


cat output_Dryasoct_contractedfams | sort | uniq > OGtotaltemp
cat OGtotaltemp | while read line ; do grep -c "$line" output_Dryasoct_contractedfams >> tempcounts ; done
paste OGtotaltemp tempcounts >> OGcountstemp
cat OGcountstemp | while read orthogroup counts ; do if [ "$counts" -gt 47 ] ; then echo "$orthogroup" >> Dryasoct_contractedfams ; fi ; done
rm -f OGcountstemp OGtotaltemp tempcounts output 


cat output_Dryasoct | sort | uniq > OGtotaltemp
cat OGtotaltemp | while read line ; do grep -c "$line" output_Dryasoct >> tempcounts ; done
paste OGtotaltemp tempcounts >> OGcountstemp
cat OGcountstemp | while read orthogroup counts ; do if [ "$counts" -gt 47 ] ; then echo "$orthogroup" >>Dryasoct_sig_changes ; fi ; done
rm -f OGcountstemp OGtotaltemp tempcounts output

#-------------------------------------------
# Arabis alpina
cat output_Arabalp_expandedfams | sort | uniq > OGtotaltemp
cat OGtotaltemp | while read line ; do grep -c "$line" output_Arabalp_expandedfams >> tempcounts ; done
paste OGtotaltemp tempcounts >> OGcountstemp
cat OGcountstemp | while read orthogroup counts ; do if [ "$counts" -gt 47 ] ; then echo "$orthogroup" >> Arabalp_expandedfams ; fi ; done
rm -f OGcountstemp OGtotaltemp tempcounts output 


cat output_Arabalp_contractedfams | sort | uniq > OGtotaltemp
cat OGtotaltemp | while read line ; do grep -c "$line" output_Arabalp_contractedfams >> tempcounts ; done
paste OGtotaltemp tempcounts >> OGcountstemp
cat OGcountstemp | while read orthogroup counts ; do if [ "$counts" -gt 47 ] ; then echo "$orthogroup" >> Arabalp_contractedfams ; fi ; done
rm -f OGcountstemp OGtotaltemp tempcounts output 


cat output_Arabalp | sort | uniq > OGtotaltemp
cat OGtotaltemp | while read line ; do grep -c "$line" output_Arabalp >> tempcounts ; done
paste OGtotaltemp tempcounts >> OGcountstemp
cat OGcountstemp | while read orthogroup counts ; do if [ "$counts" -gt 47 ] ; then echo "$orthogroup" >>Arabalp_sig_changes ; fi ; done
rm -f OGcountstemp OGtotaltemp tempcounts output

#-------------------------------------------
# Cochgro
cat output_Cochgro_expandedfams | sort | uniq > OGtotaltemp
cat OGtotaltemp | while read line ; do grep -c "$line" output_Cochgro_expandedfams >> tempcounts ; done
paste OGtotaltemp tempcounts >> OGcountstemp
cat OGcountstemp | while read orthogroup counts ; do if [ "$counts" -gt 47 ] ; then echo "$orthogroup" >> Cochgro_expandedfams ; fi ; done
rm -f OGcountstemp OGtotaltemp tempcounts output 


cat output_Cochgro_contractedfams | sort | uniq > OGtotaltemp
cat OGtotaltemp | while read line ; do grep -c "$line" output_Cochgro_contractedfams >> tempcounts ; done
paste OGtotaltemp tempcounts >> OGcountstemp
cat OGcountstemp | while read orthogroup counts ; do if [ "$counts" -gt 47 ] ; then echo "$orthogroup" >> Cochgro_contractedfams ; fi ; done
rm -f OGcountstemp OGtotaltemp tempcounts output 


cat output_Cochgro | sort | uniq > OGtotaltemp
cat OGtotaltemp | while read line ; do grep -c "$line" output_Cochgro >> tempcounts ; done
paste OGtotaltemp tempcounts >> OGcountstemp
cat OGcountstemp | while read orthogroup counts ; do if [ "$counts" -gt 47 ] ; then echo "$orthogroup" >>Cochgro_sig_changes ; fi ; done
rm -f OGcountstemp OGtotaltemp tempcounts output


#######################################################################################################################################
#Make gene list in orthogroups total expanded
mkdir ~/scratch/Oxyria/CAFE/Total_genomes/5spp_expanded
cd ~/scratch/Oxyria/CAFE/Total_genomes/5spp_expanded

cp ../5spp_analysis/Cochgro_expandedfams .
sed -i 's/N0.H//g' Cochgro_expandedfams 
wc -l Cochgro_expandedfams
# 2662 Cochgro_expandedfams

cp ../5spp_analysis/Drabaniv_expandedfams .
sed -i 's/N0.H//g' Drabaniv_expandedfams 
wc -l Drabaniv_expandedfams 
# 1487 Drabaniv_expandedfams

cp ../5spp_analysis/Arabalp_expandedfams .
sed -i 's/N0.H//g' Arabalp_expandedfams 
wc -l Arabalp_expandedfams 
# 461 Arabalp_expandedfams

cp ../5spp_analysis/Rheumnob_expandedfams .
sed -i 's/N0.H//g' Rheumnob_expandedfams 
wc -l Rheumnob_expandedfams 
# 1598

cp ../5spp_analysis/Dryasoct_expandedfams .
sed -i 's/N0.H//g' Dryasoct_expandedfams 
wc -l Dryasoct_expandedfams 
# 1059

cp ../5spp_analysis/Oxydig_expandedfams .
sed -i 's/N0.H//g' Oxydig_expandedfams 
wc -l Oxydig_expandedfams 
# 1230

#---------------------------------------------------
#get a fasta of the total expanded families from the OrthoFinder output - need to use the N0 HOGs from N0.tsv

cd /home/celphin/scratch/Oxyria/CAFE/Total_genomes/5spp_expanded

# see gene names in these orthogroups
grep  ">" /home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/tmp/Draba_nivalis.fa | head
grep  ">" /home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/tmp/Arabis_alpina.fa | head
grep  ">" /home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/tmp/Dryas_octopetala.fa | head
grep  ">" /home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/tmp/Rheum_nobile_H0.fa | head
grep  ">" /home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/tmp/Oxyria_digyna_H1.fa | head
grep  ">" /home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/tmp/Cochlearia_groenlandica.fa | head

#---------------------------------
#Get a list of the specific species genes in these orthogroups (HOGs)
for taxon in Oxydig; do cat "$taxon"_expandedfams | cut -f 1 | while read Orthogroup ; do grep N0.H"$Orthogroup" /home/celphin/scratch/Oxyria/CAFE/Total_genomes/input_data/N0.tsv | awk -F'\t' '{for(i=1; i<=NF; i++) if($i ~ /Oxyria/) print $i}'  >> "$taxon"_total_totalexpanded_geneIDs0.txt ; done ; done
for taxon in Rheumnob; do cat "$taxon"_expandedfams | cut -f 1 | while read Orthogroup ; do grep N0.H"$Orthogroup" /home/celphin/scratch/Oxyria/CAFE/Total_genomes/input_data/N0.tsv | awk -F'\t' '{for(i=1; i<=NF; i++) if($i ~ /Rno/) print $i}'  >> "$taxon"_total_totalexpanded_geneIDs0.txt ; done ; done
for taxon in Dryasoct; do cat "$taxon"_expandedfams | cut -f 1 | while read Orthogroup ; do grep N0.H"$Orthogroup" /home/celphin/scratch/Oxyria/CAFE/Total_genomes/input_data/N0.tsv | awk -F'\t' '{for(i=1; i<=NF; i++) if($i ~ /DoctH0/) print $i}'  >> "$taxon"_total_totalexpanded_geneIDs0.txt ; done ; done
for taxon in Drabaniv; do cat "$taxon"_expandedfams | cut -f 1 | while read Orthogroup ; do grep N0.H"$Orthogroup" /home/celphin/scratch/Oxyria/CAFE/Total_genomes/input_data/N0.tsv | awk -F'\t' '{for(i=1; i<=NF; i++) if($i ~ /\-lg/) print $i}'  >> "$taxon"_total_totalexpanded_geneIDs0.txt ; done ; done
for taxon in Arabalp;do cat "$taxon"_expandedfams | cut -f 1 | while read Orthogroup ; do grep N0.H"$Orthogroup" /home/celphin/scratch/Oxyria/CAFE/Total_genomes/input_data/N0.tsv | awk -F'\t' '{for(i=1; i<=NF; i++) if($i ~ /AALP/) print $i}'  >> "$taxon"_total_totalexpanded_geneIDs0.txt ; done ; done
for taxon in Cochgro;do cat "$taxon"_expandedfams | cut -f 1 | while read Orthogroup ; do grep N0.H"$Orthogroup" /home/celphin/scratch/Oxyria/CAFE/Total_genomes/input_data/N0.tsv | awk -F'\t' '{for(i=1; i<=NF; i++) if($i ~ /g[0-9]+\.t[12]/) print $i}'  >> "$taxon"_total_totalexpanded_geneIDs0.txt ; done ; done

for taxon in Oxydig Rheumnob Dryasoct Drabaniv Arabalp  Cochgro ; do tr ', ' '\n' < "$taxon"_total_totalexpanded_geneIDs0.txt  | awk 'NF' > "$taxon"_total_totalexpanded_geneIDs.txt ; done


#############################################
#Make gene list in orthogroups total contracted
mkdir cd ~/scratch/Oxyria/CAFE/Total_genomes/5spp_contracted
cd ~/scratch/Oxyria/CAFE/Total_genomes/5spp_contracted

cp ../5spp_analysis/Cochgro_contractedfams .
sed -i 's/N0.H//g' Cochgro_contractedfams 
wc -l Cochgro_contractedfams 
#2819 Cochgro_contractedfams

cp ../5spp_analysis/Drabaniv_contractedfams .
sed -i 's/N0.H//g' Drabaniv_contractedfams 
wc -l Drabaniv_contractedfams
# 854 Drabaniv_contractedfams

cp ../5spp_analysis/Arabalp_contractedfams .
sed -i 's/N0.H//g' Arabalp_contractedfams 
wc -l Arabalp_contractedfams 
#4074 Arabalp_contractedfams

cp ../5spp_analysis/Rheumnob_contractedfams .
sed -i 's/N0.H//g' Rheumnob_contractedfams 
wc -l Rheumnob_contractedfams 
# 817

cp ../5spp_analysis/Dryasoct_contractedfams .
sed -i 's/N0.H//g' Dryasoct_contractedfams 
wc -l Dryasoct_contractedfams 
# 1496

cp ../5spp_analysis/Oxydig_contractedfams .
sed -i 's/N0.H//g' Oxydig_contractedfams 
wc -l Oxydig_contractedfams 
# 2556
#---------------------------------------------------
cd /home/celphin/scratch/Oxyria/CAFE/Total_genomes/5spp_contracted

#Get a list of the specific species genes in these orthogroups (HOGs)
for taxon in Oxydig; do cat "$taxon"_contractedfams | cut -f 1 | while read Orthogroup ; do grep N0.H"$Orthogroup" /home/celphin/scratch/Oxyria/CAFE/Total_genomes/input_data/N0.tsv | awk -F'\t' '{for(i=1; i<=NF; i++) if($i ~ /Oxyria/) print $i}'  >> "$taxon"_total_totalcontracted_geneIDs0.txt ; done ; done
for taxon in Rheumnob; do cat "$taxon"_contractedfams | cut -f 1 | while read Orthogroup ; do grep N0.H"$Orthogroup" /home/celphin/scratch/Oxyria/CAFE/Total_genomes/input_data/N0.tsv | awk -F'\t' '{for(i=1; i<=NF; i++) if($i ~ /Rno/) print $i}'  >> "$taxon"_total_totalcontracted_geneIDs0.txt ; done ; done
for taxon in Dryasoct; do cat "$taxon"_contractedfams | cut -f 1 | while read Orthogroup ; do grep N0.H"$Orthogroup" /home/celphin/scratch/Oxyria/CAFE/Total_genomes/input_data/N0.tsv | awk -F'\t' '{for(i=1; i<=NF; i++) if($i ~ /DoctH0/) print $i}'  >> "$taxon"_total_totalcontracted_geneIDs0.txt ; done ; done
for taxon in Drabaniv; do cat "$taxon"_contractedfams | cut -f 1 | while read Orthogroup ; do grep N0.H"$Orthogroup" /home/celphin/scratch/Oxyria/CAFE/Total_genomes/input_data/N0.tsv | awk -F'\t' '{for(i=1; i<=NF; i++) if($i ~ /\-lg/) print $i}'  >> "$taxon"_total_totalcontracted_geneIDs0.txt ; done ; done
for taxon in Arabalp;do cat "$taxon"_contractedfams | cut -f 1 | while read Orthogroup ; do grep N0.H"$Orthogroup" /home/celphin/scratch/Oxyria/CAFE/Total_genomes/input_data/N0.tsv | awk -F'\t' '{for(i=1; i<=NF; i++) if($i ~ /AALP/) print $i}'  >> "$taxon"_total_totalcontracted_geneIDs0.txt ; done ; done
for taxon in Cochgro;do cat "$taxon"_contractedfams | cut -f 1 | while read Orthogroup ; do grep N0.H"$Orthogroup" /home/celphin/scratch/Oxyria/CAFE/Total_genomes/input_data/N0.tsv | awk -F'\t' '{for(i=1; i<=NF; i++) if($i ~ /g[0-9]+\.t[12]/) print $i}'  >> "$taxon"_total_totalcontracted_geneIDs0.txt ; done ; done


for taxon in Oxydig Rheumnob Dryasoct Drabaniv Arabalp Cochgro; do tr ', ' '\n' < "$taxon"_total_totalcontracted_geneIDs0.txt  | awk 'NF' > "$taxon"_total_totalcontracted_geneIDs.txt ; done

######################################################################################################################
#See which  orthogroups expanded and contracted rapidly
cd ~/scratch/Oxyria/CAFE/Total_genomes/5spp_analysis

for taxon in Oxydig Rheumnob Dryasoct Drabaniv Arabalp Cochgro; do grep -f "$taxon"_sig_changes "$taxon"_expandedfams > "$taxon"_expandedfams.sig ; done
for taxon in Oxydig Rheumnob Dryasoct Drabaniv Arabalp Cochgro; do grep -f "$taxon"_sig_changes "$taxon"_contractedfams > "$taxon"_contractedfams.sig ; done

wc -l *_expandedfams*.sig  
  # 36 Arabalp_expandedfams.sig
  # 65 Cochgro_expandedfams.sig
 # 131 Drabaniv_expandedfams.sig
 # 134 Dryasoct_expandedfams.sig
 # 102 Oxydig_expandedfams.sig
 # 130 Rheumnob_expandedfams.sig


wc -l *_contractedfams*.sig
 # 120 Arabalp_contractedfams.sig
  # 72 Cochgro_contractedfams.sig
  # 33 Drabaniv_contractedfams.sig
  # 39 Dryasoct_contractedfams.sig
  # 26 Oxydig_contractedfams.sig
  # 53 Rheumnob_contractedfams.sig


#######################################################################################################################################
#Make gene list in orthogroups expanded rapidly 
cd ~/scratch/Oxyria/CAFE/Total_genomes/5spp_expanded

cp ../5spp_analysis/Oxydig_expandedfams.sig .
sed -i 's/N0.H//g' Oxydig_expandedfams.sig 

cp ../5spp_analysis/Rheumnob_expandedfams.sig .
sed -i 's/N0.H//g' Rheumnob_expandedfams.sig 

cp ../5spp_analysis/Dryasoct_expandedfams.sig .
sed -i 's/N0.H//g' Dryasoct_expandedfams.sig 

cp ../5spp_analysis/Drabaniv_expandedfams.sig .
sed -i 's/N0.H//g' Drabaniv_expandedfams.sig 

cp ../5spp_analysis/Arabalp_expandedfams.sig .
sed -i 's/N0.H//g' Arabalp_expandedfams.sig 

cp ../5spp_analysis/Cochgro_expandedfams.sig .
sed -i 's/N0.H//g' Cochgro_expandedfams.sig 

#---------------------------------------------------
#get a fasta of the rapidly expanded families from the OrthoFinder output

cd /home/celphin/scratch/Oxyria/CAFE/Total_genomes/5spp_expanded

#Get a list of the specific species genes in these orthogroups (HOGs)
for taxon in Oxydig; do cat "$taxon"_expandedfams.sig | cut -f 1 | while read Orthogroup ; do grep N0.H"$Orthogroup" /home/celphin/scratch/Oxyria/CAFE/Total_genomes/input_data/N0.tsv | awk -F'\t' '{for(i=1; i<=NF; i++) if($i ~ /Oxyria/) print $i}'  >> "$taxon"_total_expanded_geneIDs0.txt ; done ; done
for taxon in Rheumnob; do cat "$taxon"_expandedfams.sig | cut -f 1 | while read Orthogroup ; do grep N0.H"$Orthogroup" /home/celphin/scratch/Oxyria/CAFE/Total_genomes/input_data/N0.tsv | awk -F'\t' '{for(i=1; i<=NF; i++) if($i ~ /Rno/) print $i}'  >> "$taxon"_total_expanded_geneIDs0.txt ; done ; done
for taxon in Dryasoct; do cat "$taxon"_expandedfams.sig | cut -f 1 | while read Orthogroup ; do grep N0.H"$Orthogroup" /home/celphin/scratch/Oxyria/CAFE/Total_genomes/input_data/N0.tsv | awk -F'\t' '{for(i=1; i<=NF; i++) if($i ~ /DoctH0/) print $i}'  >> "$taxon"_total_expanded_geneIDs0.txt ; done ; done
for taxon in Drabaniv; do cat "$taxon"_expandedfams.sig | cut -f 1 | while read Orthogroup ; do grep N0.H"$Orthogroup" /home/celphin/scratch/Oxyria/CAFE/Total_genomes/input_data/N0.tsv | awk -F'\t' '{for(i=1; i<=NF; i++) if($i ~ /\-lg/) print $i}'  >> "$taxon"_total_expanded_geneIDs0.txt ; done ; done
for taxon in Arabalp;do cat "$taxon"_expandedfams.sig | cut -f 1 | while read Orthogroup ; do grep N0.H"$Orthogroup" /home/celphin/scratch/Oxyria/CAFE/Total_genomes/input_data/N0.tsv | awk -F'\t' '{for(i=1; i<=NF; i++) if($i ~ /AALP/) print $i}'  >> "$taxon"_total_expanded_geneIDs0.txt ; done ; done
for taxon in Cochgro;do cat "$taxon"_expandedfams.sig | cut -f 1 | while read Orthogroup ; do grep N0.H"$Orthogroup" /home/celphin/scratch/Oxyria/CAFE/Total_genomes/input_data/N0.tsv | awk -F'\t' '{for(i=1; i<=NF; i++) if($i ~ /g[0-9]+\.t[12]/) print $i}'  >> "$taxon"_total_expanded_geneIDs0.txt ; done ; done

for taxon in Oxydig Rheumnob Dryasoct Drabaniv Arabalp  Cochgro ; do tr ', ' '\n' < "$taxon"_total_expanded_geneIDs0.txt  | awk 'NF' > "$taxon"_total_expanded_geneIDs.txt ; done

#############################################
#Make gene list in orthogroups contracted rapidly 
cd ~/scratch/Oxyria/CAFE/Total_genomes/5spp_contracted

cp ../5spp_analysis/Oxydig_contractedfams.sig .
sed -i 's/N0.H//g' Oxydig_contractedfams.sig 

cp ../5spp_analysis/Rheumnob_contractedfams.sig .
sed -i 's/N0.H//g' Rheumnob_contractedfams.sig 

cp ../5spp_analysis/Dryasoct_contractedfams.sig .
sed -i 's/N0.H//g' Dryasoct_contractedfams.sig 

cp ../5spp_analysis/Drabaniv_contractedfams.sig .
sed -i 's/N0.H//g' Drabaniv_contractedfams.sig 

cp ../5spp_analysis/Arabalp_contractedfams.sig .
sed -i 's/N0.H//g' Arabalp_contractedfams.sig 

cp ../5spp_analysis/Cochgro_contractedfams.sig .
sed -i 's/N0.H//g' Cochgro_contractedfams.sig 

#---------------------------------------------------
#get a fasta of the rapidly contracted families from the OrthoFinder output

cd /home/celphin/scratch/Oxyria/CAFE/Total_genomes/5spp_contracted

#Get a list of the specific species genes in these orthogroups (HOGs)
for taxon in Oxydig; do cat "$taxon"_contractedfams.sig | cut -f 1 | while read Orthogroup ; do grep N0.H"$Orthogroup" /home/celphin/scratch/Oxyria/CAFE/Total_genomes/input_data/N0.tsv | awk -F'\t' '{for(i=1; i<=NF; i++) if($i ~ /Oxyria/) print $i}'  >> "$taxon"_total_contracted_geneIDs0.txt ; done ; done
for taxon in Rheumnob; do cat "$taxon"_contractedfams.sig | cut -f 1 | while read Orthogroup ; do grep N0.H"$Orthogroup" /home/celphin/scratch/Oxyria/CAFE/Total_genomes/input_data/N0.tsv | awk -F'\t' '{for(i=1; i<=NF; i++) if($i ~ /Rno/) print $i}'  >> "$taxon"_total_contracted_geneIDs0.txt ; done ; done
for taxon in Dryasoct; do cat "$taxon"_contractedfams.sig | cut -f 1 | while read Orthogroup ; do grep N0.H"$Orthogroup" /home/celphin/scratch/Oxyria/CAFE/Total_genomes/input_data/N0.tsv | awk -F'\t' '{for(i=1; i<=NF; i++) if($i ~ /DoctH0/) print $i}'  >> "$taxon"_total_contracted_geneIDs0.txt ; done ; done
for taxon in Drabaniv; do cat "$taxon"_contractedfams.sig | cut -f 1 | while read Orthogroup ; do grep N0.H"$Orthogroup" /home/celphin/scratch/Oxyria/CAFE/Total_genomes/input_data/N0.tsv | awk -F'\t' '{for(i=1; i<=NF; i++) if($i ~ /\-lg/) print $i}'  >> "$taxon"_total_contracted_geneIDs0.txt ; done ; done
for taxon in Arabalp;do cat "$taxon"_contractedfams.sig | cut -f 1 | while read Orthogroup ; do grep N0.H"$Orthogroup" /home/celphin/scratch/Oxyria/CAFE/Total_genomes/input_data/N0.tsv | awk -F'\t' '{for(i=1; i<=NF; i++) if($i ~ /AALP/) print $i}'  >> "$taxon"_total_contracted_geneIDs0.txt ; done ; done
for taxon in Cochgro;do cat "$taxon"_contractedfams.sig | cut -f 1 | while read Orthogroup ; do grep N0.H"$Orthogroup" /home/celphin/scratch/Oxyria/CAFE/Total_genomes/input_data/N0.tsv | awk -F'\t' '{for(i=1; i<=NF; i++) if($i ~ /g[0-9]+\.t[12]/) print $i}'  >> "$taxon"_total_contracted_geneIDs0.txt ; done ; done

for taxon in Oxydig Rheumnob Dryasoct Drabaniv Arabalp  Cochgro ; do tr ', ' '\n' < "$taxon"_total_contracted_geneIDs0.txt  | awk 'NF' > "$taxon"_total_contracted_geneIDs.txt ; done


#######################################################################################################################################
# look at counts

cd  ~/scratch/Oxyria/CAFE/Total_genomes/5spp_expanded

wc -l *expanded_geneIDs.txt
  # 256 Arabalp_total_expanded_geneIDs.txt
  # 1693 Arabalp_total_totalexpanded_geneIDs.txt
   # 698 Cochgro_total_expanded_geneIDs.txt
  # 8132 Cochgro_total_totalexpanded_geneIDs.txt
  # 1348 Drabaniv_total_expanded_geneIDs.txt
  # 5355 Drabaniv_total_totalexpanded_geneIDs.txt
  # 2251 Dryasoct_total_expanded_geneIDs.txt
  # 5325 Dryasoct_total_totalexpanded_geneIDs.txt
  # 2187 Oxydig_total_expanded_geneIDs.txt
  # 5834 Oxydig_total_totalexpanded_geneIDs.txt
   # 941 Rheumnob_total_expanded_geneIDs.txt
  # 5510 Rheumnob_total_totalexpanded_geneIDs.txt


#----------------
cd  ~/scratch/Oxyria/CAFE/Total_genomes/5spp_contracted

wc -l *contracted_geneIDs.txt
   # 189 Arabalp_total_contracted_geneIDs.txt
  # 1653 Arabalp_total_totalcontracted_geneIDs.txt
   # 112 Cochgro_total_contracted_geneIDs.txt
  # 1793 Cochgro_total_totalcontracted_geneIDs.txt
    # 93 Drabaniv_total_contracted_geneIDs.txt
   # 764 Drabaniv_total_totalcontracted_geneIDs.txt
    # 51 Dryasoct_total_contracted_geneIDs.txt
  # 1014 Dryasoct_total_totalcontracted_geneIDs.txt
    # 19 Oxydig_total_contracted_geneIDs.txt
  # 1221 Oxydig_total_totalcontracted_geneIDs.txt
   # 108 Rheumnob_total_contracted_geneIDs.txt
   # 758 Rheumnob_total_totalcontracted_geneIDs.txt


#####################################################################################################################
#From here: continue at Notes6_Go_enrichment.bash
