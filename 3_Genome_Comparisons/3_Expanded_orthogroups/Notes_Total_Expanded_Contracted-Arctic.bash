#################################################################################################################
# Sept 2024 - Arctic Species expanded gene families
# Tutorial https://github.com/elsemikk/Willisornis_Genome_Assembly/blob/master/4.2_Gene_Family_Expansions.md
#################################################################################################################
tmux new-session -s CAFE1
tmux attach-session -t CAFE1

cd /lustre04/scratch/celphin/Oxyria/CAFE/Total_genomes/analysis

head -n 3 ../output/error_model/Base_clade_results.txt

#get a list of the significantly expanded/contracted/rapid genes in Oxy/Rhuem/Dryas/Draba
#-------------------------------------
# List for HOGs in all spp

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

head -n 3 ../output/error_model/Base_asr.tre

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

#------------------------------------
# High Arctic
#get a list of the significantly expanded/contracted/rapid genes in Oxyria_digyna_H1<6>
printf '%s\n' {1..50} | while read number ; do grep "Oxyria_digyna_H1<6>\*" ../output/error_model_$number/Base_asr.tre | \
sed 's/^.*OG/OG/g' | sed 's/ = .*$//g' > Oxydig_sig_changes_error_model_"$number" ; done 

#get a list of the significantly expanded/contracted/rapid genes in Draba_nivalis<18>
printf '%s\n' {1..50} | while read number ; do grep "Draba_nivalis<18>\*" ../output/error_model_$number/Base_asr.tre | \
sed 's/^.*OG/OG/g' | sed 's/ = .*$//g' > Drabaniv_sig_changes_error_model_"$number" ; done 

#get a list of the significantly expanded/contracted/rapid genes in Dryas_octopetala<10>
printf '%s\n' {1..50} | while read number ; do grep "Dryas_octopetala<10>\*" ../output/error_model_$number/Base_asr.tre | \
sed 's/^.*OG/OG/g' | sed 's/ = .*$//g' > Dryasoct_sig_changes_error_model_"$number" ; done 

#get a list of the significantly expanded/contracted/rapid genes in Cochlearia_groenlandica<14>
printf '%s\n' {1..50} | while read number ; do grep "Cochlearia_groenlandica<14>\*" ../output/error_model_$number/Base_asr.tre | \
sed 's/^.*OG/OG/g' | sed 's/ = .*$//g' > Cochgro_sig_changes_error_model_"$number" ; done 

# Alpine
#get a list of the significantly expanded/contracted/rapid genes in Rheum_nobile_H0<4>
printf '%s\n' {1..50} | while read number ; do grep "Rheum_nobile_H0<4>\*" ../output/error_model_$number/Base_asr.tre | \
sed 's/^.*OG/OG/g' | sed 's/ = .*$//g' > Rheumnob_sig_changes_error_model_"$number" ; done 

#get a list of the significantly expanded/contracted/rapid genes in Arabis_alpina<17>
printf '%s\n' {1..50} | while read number ; do grep "Arabis_alpina<17>\*" ../output/error_model_$number/Base_asr.tre | \
sed 's/^.*OG/OG/g' | sed 's/ = .*$//g' > Arabalp_sig_changes_error_model_"$number" ; done 

# Clade
#get a list of the significantly expanded/contracted/rapid genes in Arabis/Draba clase <38>
printf '%s\n' {1..50} | while read number ; do grep "<38>\*" ../output/error_model_$number/Base_asr.tre | \
sed 's/^.*OG/OG/g' | sed 's/ = .*$//g' > Dra_Ara_sig_changes_error_model_"$number" ; done 

#------------------------
#Data exploration:
wc -l Oxydig_sig_changes_error_model_* | head # 24-29
wc -l Rheumnob_sig_changes_error_model_* | head # 39-46
wc -l Dra_Ara_sig_changes_error_model_* | head # 12-13
wc -l Drabaniv_sig_changes_error_model_* | head # 36-41
wc -l Dryasoct_sig_changes_error_model_* | head # 30-33
wc -l Arabalp_sig_changes_error_model_* | head # 47-54
wc -l Cochgro_sig_changes_error_model_* | head # 50-56

#---------------------
#see what is uniquely found in run 6 (this was the time I had no Hypoch/Rhemel)
printf '%s\n' {1..10} | while read number ; do grep -v -f Oxydig_sig_changes_error_model_"$number" Oxydig_sig_changes_error_model_10 ; done
printf '%s\n' {1..10} | while read number ; do grep -v -f Oxydig_sig_changes_error_model_"$number" Oxydig_sig_changes_error_model_1 ; done
# many differences between runs in Oxyria

##################################################################################################################
#List of expanded and contracted groups (both significant and not)

#Get list of families that have expanded in Oxyria, whether or not significant
# value should be one greater than ID to get right column

printf '%s\n' {1..50} | while read number ; do cut -f 1,7 ../output/error_model_$number/Base_change.tab | \
awk -F'\t' '$2 != 0' | grep -v "-" > Oxydig_expandedfams_error_model_"$number" ; done 

printf '%s\n' {1..50} | while read number ; do cut -f 1,19 ../output/error_model_$number/Base_change.tab | \
awk -F'\t' '$2 != 0' | grep -v "-" > Drabaniv_expandedfams_error_model_"$number" ; done 

printf '%s\n' {1..50} | while read number ; do cut -f 1,11 ../output/error_model_$number/Base_change.tab | \
awk -F'\t' '$2 != 0' | grep -v "-" > Dryasoct_expandedfams_error_model_"$number" ; done 

printf '%s\n' {1..50} | while read number ; do cut -f 1,15 ../output/error_model_$number/Base_change.tab | \
awk -F'\t' '$2 != 0' | grep -v "-" > Cochgro_expandedfams_error_model_"$number" ; done 


printf '%s\n' {1..50} | while read number ; do cut -f 1,18 ../output/error_model_$number/Base_change.tab | \
awk -F'\t' '$2 != 0' | grep -v "-" > Arabalp_expandedfams_error_model_"$number" ; done 

printf '%s\n' {1..50} | while read number ; do cut -f 1,5 ../output/error_model_$number/Base_change.tab | \
awk -F'\t' '$2 != 0' | grep -v "-" > Rheumnob_expandedfams_error_model_"$number" ; done 


printf '%s\n' {1..50} | while read number ; do cut -f 1,39 ../output/error_model_$number/Base_change.tab | \
awk -F'\t' '$2 != 0' | grep -v "-" > Dra_Ara_expandedfams_error_model_"$number" ; done 


wc -l *expandedfams_error_model_1 #remember the header is one line!
  # 232 Arabalp_expandedfams_error_model_1
 # 1125 Cochgro_expandedfams_error_model_1
   # 61 Dra_Ara_expandedfams_error_model_1
  # 567 Drabaniv_expandedfams_error_model_1
  # 281 Dryasoct_expandedfams_error_model_1
  # 398 Oxydig_expandedfams_error_model_1
  # 531 Rheumnob_expandedfams_error_model_1

#--------------------------------------

#Get list of families that have contracted in Oxyria, whether or not significant
# value should be one greater than ID to get right column

printf '%s\n' {1..50} | while read number ; do cut -f 1,7 ../output/error_model_$number/Base_change.tab | \
awk -F'\t' '$2 != 0' | grep "-" > Oxydig_contractedfams_error_model_"$number" ; done 

printf '%s\n' {1..50} | while read number ; do cut -f 1,19 ../output/error_model_$number/Base_change.tab | \
awk -F'\t' '$2 != 0' | grep "-" > Drabaniv_contractedfams_error_model_"$number" ; done 

printf '%s\n' {1..50} | while read number ; do cut -f 1,11 ../output/error_model_$number/Base_change.tab | \
awk -F'\t' '$2 != 0' | grep "-" > Dryasoct_contractedfams_error_model_"$number" ; done 

printf '%s\n' {1..50} | while read number ; do cut -f 1,15 ../output/error_model_$number/Base_change.tab | \
awk -F'\t' '$2 != 0' | grep "-" > Cochgro_contractedfams_error_model_"$number" ; done 


printf '%s\n' {1..50} | while read number ; do cut -f 1,18 ../output/error_model_$number/Base_change.tab | \
awk -F'\t' '$2 != 0' | grep "-" > Arabalp_contractedfams_error_model_"$number" ; done 

printf '%s\n' {1..50} | while read number ; do cut -f 1,5 ../output/error_model_$number/Base_change.tab | \
awk -F'\t' '$2 != 0' | grep "-" > Rheumnob_contractedfams_error_model_"$number" ; done 


printf '%s\n' {1..50} | while read number ; do cut -f 1,39 ../output/error_model_$number/Base_change.tab | \
awk -F'\t' '$2 != 0' | grep "-" > Dra_Ara_contractedfams_error_model_"$number" ; done 


wc -l *contractedfams_error_model_1 #remember the header is one line!
  # 674 Arabalp_contractedfams_error_model_1
  # 589 Cochgro_contractedfams_error_model_1
  # 246 Dra_Ara_contractedfams_error_model_1
  # 227 Drabaniv_contractedfams_error_model_1
  # 344 Dryasoct_contractedfams_error_model_1
  # 446 Oxydig_contractedfams_error_model_1
  # 184 Rheumnob_contractedfams_error_model_1


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
export= open(jn(folderpath,'output_Cochgro'),'w')
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
#Rheum-Oxyria node
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
mkdir ~/scratch/Oxyria/CAFE/Total_genomes/expanded
cd ~/scratch/Oxyria/CAFE/Total_genomes/expanded

cp ../analysis/Cochgro_expandedfams .
sed -i 's/N0.H//g' Cochgro_expandedfams 
wc -l Cochgro_expandedfams
# 1124 Cochgro_expandedfams

cp ../analysis/Drabaniv_expandedfams .
sed -i 's/N0.H//g' Drabaniv_expandedfams 
wc -l Drabaniv_expandedfams 
# 566 Drabaniv_expandedfams

cp ../analysis/Arabalp_expandedfams .
sed -i 's/N0.H//g' Arabalp_expandedfams 
wc -l Arabalp_expandedfams 
# 231 Arabalp_expandedfams

cp ../analysis/Rheumnob_expandedfams .
sed -i 's/N0.H//g' Rheumnob_expandedfams 
wc -l Rheumnob_expandedfams 
# 530

cp ../analysis/Dryasoct_expandedfams .
sed -i 's/N0.H//g' Dryasoct_expandedfams 
wc -l Dryasoct_expandedfams 
# 280

cp ../analysis/Oxydig_expandedfams .
sed -i 's/N0.H//g' Oxydig_expandedfams 
wc -l Oxydig_expandedfams 
# 397

#---------------------------------------------------
#get a fasta of the total expanded families from the OrthoFinder output - need to use the N0 HOGs from N0.tsv

cd /home/celphin/scratch/Oxyria/CAFE/Total_genomes/expanded

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
mkdir ~/scratch/Oxyria/CAFE/Total_genomes/contracted
cd ~/scratch/Oxyria/CAFE/Total_genomes/contracted

cp ../analysis/Cochgro_contractedfams .
sed -i 's/N0.H//g' Cochgro_contractedfams 
wc -l Cochgro_contractedfams 
#588 Cochgro_contractedfams

cp ../analysis/Drabaniv_contractedfams .
sed -i 's/N0.H//g' Drabaniv_contractedfams 
wc -l Drabaniv_contractedfams
# 226 Drabaniv_contractedfams

cp ../analysis/Arabalp_contractedfams .
sed -i 's/N0.H//g' Arabalp_contractedfams 
wc -l Arabalp_contractedfams 
#673 Arabalp_contractedfams

cp ../analysis/Rheumnob_contractedfams .
sed -i 's/N0.H//g' Rheumnob_contractedfams 
wc -l Rheumnob_contractedfams 
# 183

cp ../analysis/Dryasoct_contractedfams .
sed -i 's/N0.H//g' Dryasoct_contractedfams 
wc -l Dryasoct_contractedfams 
# 343

cp ../analysis/Oxydig_contractedfams .
sed -i 's/N0.H//g' Oxydig_contractedfams 
wc -l Oxydig_contractedfams 
# 445
#---------------------------------------------------
cd /home/celphin/scratch/Oxyria/CAFE/Total_genomes/contracted

#Get a list of the specific species genes in these orthogroups (HOGs)
for taxon in Oxydig; do cat "$taxon"_contractedfams | cut -f 1 | while read Orthogroup ; do grep N0.H"$Orthogroup" /home/celphin/scratch/Oxyria/CAFE/Total_genomes/input_data/N0.tsv | awk -F'\t' '{for(i=1; i<=NF; i++) if($i ~ /Oxyria/) print $i}'  >> "$taxon"_total_totalcontracted_geneIDs0.txt ; done ; done
for taxon in Rheumnob; do cat "$taxon"_contractedfams | cut -f 1 | while read Orthogroup ; do grep N0.H"$Orthogroup" /home/celphin/scratch/Oxyria/CAFE/Total_genomes/input_data/N0.tsv | awk -F'\t' '{for(i=1; i<=NF; i++) if($i ~ /Rno/) print $i}'  >> "$taxon"_total_totalcontracted_geneIDs0.txt ; done ; done
for taxon in Dryasoct; do cat "$taxon"_contractedfams | cut -f 1 | while read Orthogroup ; do grep N0.H"$Orthogroup" /home/celphin/scratch/Oxyria/CAFE/Total_genomes/input_data/N0.tsv | awk -F'\t' '{for(i=1; i<=NF; i++) if($i ~ /DoctH0/) print $i}'  >> "$taxon"_total_totalcontracted_geneIDs0.txt ; done ; done
for taxon in Drabaniv; do cat "$taxon"_contractedfams | cut -f 1 | while read Orthogroup ; do grep N0.H"$Orthogroup" /home/celphin/scratch/Oxyria/CAFE/Total_genomes/input_data/N0.tsv | awk -F'\t' '{for(i=1; i<=NF; i++) if($i ~ /\-lg/) print $i}'  >> "$taxon"_total_totalcontracted_geneIDs0.txt ; done ; done
for taxon in Arabalp;do cat "$taxon"_contractedfams | cut -f 1 | while read Orthogroup ; do grep N0.H"$Orthogroup" /home/celphin/scratch/Oxyria/CAFE/Total_genomes/input_data/N0.tsv | awk -F'\t' '{for(i=1; i<=NF; i++) if($i ~ /AALP/) print $i}'  >> "$taxon"_total_totalcontracted_geneIDs0.txt ; done ; done
for taxon in Cochgro;do cat "$taxon"_contractedfams | cut -f 1 | while read Orthogroup ; do grep N0.H"$Orthogroup" /home/celphin/scratch/Oxyria/CAFE/Total_genomes/input_data/N0.tsv | awk -F'\t' '{for(i=1; i<=NF; i++) if($i ~ /g[0-9]+\.t[12]/) print $i}'  >> "$taxon"_total_totalcontracted_geneIDs0.txt ; done ; done

for taxon in Oxydig Rheumnob Dryasoct Drabaniv Arabalp  Cochgro ; do tr ', ' '\n' < "$taxon"_total_totalcontracted_geneIDs0.txt  | awk 'NF' > "$taxon"_total_totalcontracted_geneIDs.txt ; done


######################################################################################################################
#See which  orthogroups expanded and contracted rapidly
cd ~/scratch/Oxyria/CAFE/Total_genomes/analysis
for taxon in Oxydig Rheumnob Dryasoct Drabaniv Arabalp Cochgro; do grep -f "$taxon"_sig_changes "$taxon"_expandedfams > "$taxon"_expandedfams.sig ; done
for taxon in Oxydig Rheumnob Dryasoct Drabaniv Arabalp Cochgro; do grep -f "$taxon"_sig_changes "$taxon"_contractedfams > "$taxon"_contractedfams.sig ; done

wc -l *_expandedfams*.sig  
  # 13 Arabalp_expandedfams.sig
  # 20 Cochgro_expandedfams.sig
  # 24 Drabaniv_expandedfams.sig
  # 13 Dryasoct_expandedfams.sig
  # 15 Oxydig_expandedfams.sig
  # 32 Rheumnob_expandedfams.sig

wc -l *_contractedfams*.sig
  # 33 Arabalp_contractedfams.sig
  # 28 Cochgro_contractedfams.sig
  # 11 Drabaniv_contractedfams.sig
  # 16 Dryasoct_contractedfams.sig
   # 8 Oxydig_contractedfams.sig
   # 5 Rheumnob_contractedfams.sig


#######################################################################################################################################
#Make gene list in orthogroups expanded rapidly 
cd ~/scratch/Oxyria/CAFE/Total_genomes/expanded

cp ../analysis/Oxydig_expandedfams.sig .
sed -i 's/N0.H//g' Oxydig_expandedfams.sig 

cp ../analysis/Rheumnob_expandedfams.sig .
sed -i 's/N0.H//g' Rheumnob_expandedfams.sig 

cp ../analysis/Dryasoct_expandedfams.sig .
sed -i 's/N0.H//g' Dryasoct_expandedfams.sig 

cp ../analysis/Drabaniv_expandedfams.sig .
sed -i 's/N0.H//g' Drabaniv_expandedfams.sig 

cp ../analysis/Arabalp_expandedfams.sig .
sed -i 's/N0.H//g' Arabalp_expandedfams.sig 

cp ../analysis/Cochgro_expandedfams.sig .
sed -i 's/N0.H//g' Cochgro_expandedfams.sig 

#---------------------------------------------------
cd /home/celphin/scratch/Oxyria/CAFE/Total_genomes/expanded

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
cd ~/scratch/Oxyria/CAFE/Total_genomes/contracted

cp ../analysis/Oxydig_contractedfams.sig .
sed -i 's/N0.H//g' Oxydig_contractedfams.sig 

cp ../analysis/Rheumnob_contractedfams.sig .
sed -i 's/N0.H//g' Rheumnob_contractedfams.sig 

cp ../analysis/Dryasoct_contractedfams.sig .
sed -i 's/N0.H//g' Dryasoct_contractedfams.sig 

cp ../analysis/Drabaniv_contractedfams.sig .
sed -i 's/N0.H//g' Drabaniv_contractedfams.sig 

cp ../analysis/Arabalp_contractedfams.sig .
sed -i 's/N0.H//g' Arabalp_contractedfams.sig 

cp ../analysis/Cochgro_contractedfams.sig .
sed -i 's/N0.H//g' Cochgro_contractedfams.sig 

#---------------------------------------------------
#get a fasta of the rapidly contracted families from the OrthoFinder output
cd /home/celphin/scratch/Oxyria/CAFE/Total_genomes/contracted

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

cd  ~/scratch/Oxyria/CAFE/Total_genomes/expanded

wc -l *expanded_geneIDs.txt
    # 94 Arabalp_total_expanded_geneIDs.txt
   # 738 Arabalp_total_totalexpanded_geneIDs.txt
   # 152 Cochgro_total_expanded_geneIDs.txt
  # 3357 Cochgro_total_totalexpanded_geneIDs.txt
   # 176 Drabaniv_total_expanded_geneIDs.txt
  # 1884 Drabaniv_total_totalexpanded_geneIDs.txt
   # 144 Dryasoct_total_expanded_geneIDs.txt
  # 1056 Dryasoct_total_totalexpanded_geneIDs.txt
   # 191 Oxydig_total_expanded_geneIDs.txt
  # 1452 Oxydig_total_totalexpanded_geneIDs.txt
   # 228 Rheumnob_total_expanded_geneIDs.txt
  # 1781 Rheumnob_total_totalexpanded_geneIDs.txt


#----------------
cd  ~/scratch/Oxyria/CAFE/Total_genomes/contracted

wc -l *contracted_geneIDs.txt
   # 64 Arabalp_total_contracted_geneIDs.txt
  # 931 Arabalp_total_totalcontracted_geneIDs.txt
   # 51 Cochgro_total_contracted_geneIDs.txt
  # 844 Cochgro_total_totalcontracted_geneIDs.txt
   # 42 Drabaniv_total_contracted_geneIDs.txt
  # 378 Drabaniv_total_totalcontracted_geneIDs.txt
   # 29 Dryasoct_total_contracted_geneIDs.txt
  # 479 Dryasoct_total_totalcontracted_geneIDs.txt
   # 10 Oxydig_total_contracted_geneIDs.txt
  # 585 Oxydig_total_totalcontracted_geneIDs.txt
    # 5 Rheumnob_total_contracted_geneIDs.txt
  # 268 Rheumnob_total_totalcontracted_geneIDs.txt



#####################################################################################################################
#From here: continue at Notes6_Go_enrichment.bash
cd  ~/scratch/Oxyria/CAFE/Total_genomes/expanded

cat *_expandedfams.sig | sort | uniq -c |sort
      # 3 OG0001637
      # 3 OG0002524

cat Oxydig_expandedfams.sig Drabaniv_expandedfams.sig Dryasoct_expandedfams.sig Cochgro_expandedfams.sig | sort | uniq -c |sort
     # 2 OG0000494
      # 2 OG0001682
      # 2 OG0002229
      # 2 OG0002524
      # 2 OG0002766
      # 2 OG0003181
      # 2 OG0003607
      # 2 OG0004798
      # 2 OG0007299
      # 2 OG0008661



cat *_expandedfams | sort | uniq -c |sort
      # 4 OG0001192
      # 4 OG0003559
      # 4 OG0005805
      # 4 OG0006331
      # 4 OG0007184
      # 4 OG0010794

cat Oxydig_expandedfams Drabaniv_expandedfams Dryasoct_expandedfams Cochgro_expandedfams | sort | uniq -c |sort
      3 OG0000494
      3 OG0000729
      3 OG0001682
      3 OG0001924
      3 OG0002694
      3 OG0002766
      3 OG0003559
      3 OG0003607
      3 OG0004403
      3 OG0004429
      3 OG0005000
      3 OG0005073
      3 OG0005805
      3 OG0006063
      3 OG0006331
      3 OG0007110
      3 OG0007184
      3 OG0008014
      3 OG0009477
      3 OG0010010
      3 OG0010313
      3 OG0010744
      3 OG0010794
      3 OG0011344
      3 OG0011962
      3 OG0012246
      3 OG0012774
      3 OG0013153
      3 OG0013243




# Note they are all HOGs for IDs
grep N0.HOG0001924 /home/celphin/scratch/Oxyria/CAFE/Total_genomes/input_data/N0.tsv
#Oxyria_NCBI_Chr100002225, Oxyria_NCBI_Chr100002224, Oxyria_NCBI_Chr400004801, Oxyria_NCBI_Chr400004804, Oxyria_NCBI_Chr400004808, Oxyria_NCBI_Chr400004805, Oxyria_NCBI_Chr400004802, Oxyria_NCBI_Chr400004799, Oxyria_NCBI_Chr400004806, Oxyria_NCBI_Chr400004809, Oxyria_NCBI_Chr400004810, Oxyria_NCBI_Chr400004803
# DoctH0_Chr400005271, DoctH0_Chr400005059, DoctH0_Chr100009296, DoctH0_Chr100009292, DoctH0_Chr100009297, DoctH0_Chr100009277
# g10183.t1, g12475.t1, g28921.t1   
# maker-scaffold41-snap-gene-1.248-mRNA-1, maker-lg8-snap-gene-50.118-mRNA-1, maker-lg8-exonerate_protein2genome-gene-50.87-mRNA-1, maker-lg8-exonerate_protein2genome-gene-47.124-mRNA-1, snap_masked-lg2-processed-gene-120.69-mRNA-1, maker-lg1-exonerate_protein2genome-gene-100.180-mRNA-1, snap_masked-lg5-processed-gene-0.97-mRNA-1

grep Oxyria_NCBI_Chr100002225 /home/celphin/scratch/Oxyria/CAFE/Total_genomes/input_data/Oxydig_interproscan_edited.tsv
grep Oxyria_NCBI_Chr400004801 /home/celphin/scratch/Oxyria/CAFE/Total_genomes/input_data/Oxydig_interproscan_edited.tsv


#-------------
grep N0.HOG0004403 /home/celphin/scratch/Oxyria/CAFE/Total_genomes/input_data/N0.tsv
# Oxyria_NCBI_Chr500001625, Oxyria_NCBI_Chr400003139, Oxyria_NCBI_Chr400003142, Oxyria_NCBI_Chr800000117 
# DoctH0_Chr800004524, DoctH0_Chr300003839, DoctH0_Chr800004539, DoctH0_Chr800004523, DoctH0_Chr300003840
# RnoG0034402.1, RnoG0004684.1, RnoG0029763.1     RtaG0024705.1   
# g10889.t1, g1459.t1, g20405.t1
# AALP_AA1G058900
# snap_masked-lg7-processed-gene-15.75-mRNA-1, snap_masked-lg7-processed-gene-15.74-mRNA-1, maker-lg5-snap-gene-103.278-mRNA-1, maker-lg7-snap-gene-15.183-mRNA-1, maker-lg5-snap-gene-103.277-mRNA-1

grep Oxyria_NCBI_Chr500001625 /home/celphin/scratch/Oxyria/CAFE/Total_genomes/input_data/Oxydig_interproscan_edited.tsv
grep Oxyria_NCBI_Chr400003139 /home/celphin/scratch/Oxyria/CAFE/Total_genomes/input_data/Oxydig_interproscan_edited.tsv

#####################################

cd  ~/scratch/Oxyria/CAFE/Total_genomes/contracted

cat *_contractedfams.sig | sort | uniq -c |sort



cat Oxydig_contractedfams.sig Drabaniv_contractedfams.sig Dryasoct_contractedfams.sig Cochgro_contractedfams.sig | sort | uniq -c |sort
      2 OG0000499
      2 OG0000872
      2 OG0001631
      2 OG0001801
      2 OG0001804
      2 OG0002440
      2 OG0004125
      2 OG0004638
      3 OG0001477



cat *_contractedfams | sort | uniq -c |sort
      4 OG0000603
      4 OG0001450
      4 OG0001477
      4 OG0001841
      4 OG0002222
      4 OG0002305
      4 OG0002331
      4 OG0002430
      4 OG0002966
      4 OG0003928
      4 OG0003952
      4 OG0004239
      4 OG0004453
      4 OG0004643
      4 OG0004825
      4 OG0004902
      4 OG0005062
      4 OG0005121
      4 OG0005510
      4 OG0005875
      4 OG0006180
      4 OG0006240
      4 OG0006766
      4 OG0006791
      4 OG0007614
      5 OG0001804
      5 OG0002438
      5 OG0002440
      5 OG0010235
      5 OG0010440


cat Oxydig_contractedfams Drabaniv_contractedfams Dryasoct_contractedfams Cochgro_contractedfams | sort | uniq -c |sort
        4 OG0001477
      4 OG0002438
      4 OG0002440
