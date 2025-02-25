#################################################################################################################
# July 2024 - Non-Arctic
# Tutorial https://github.com/elsemikk/Willisornis_Genome_Assembly/blob/master/4.2_Gene_Family_Expansions.md

#################################################################################################################
mkdir /lustre04/scratch/celphin/Oxyria/CAFE/analysis-non-arctic
cd /lustre04/scratch/celphin/Oxyria/CAFE/analysis-non-arctic

head -n 3 ../output/error_model/Base_asr.tre

#get a list of the significantly expanded/contracted/rapid genes in Oxy/Rhuem/Dryas/Draba
# FamilyID        
# Polygunum_aviculare_H0<1>       
# Oxyria_digyna_H1<2>     
# Rheum_tangaticum_H0<3>  
# Rheum_nobile_H0<4>     
# Fagopyrum_tataricum_H1<5>       
# Fagopyrum_escelentum_H2<6>      
# Capsella_rubella<7>   
# Arabidopsis_lyrata<8>    
# Arabidopsis_thaliana<9> 
# Draba_nivalis<10>       
# Arabis_alpina<11>       
# Brassica_oleracea<12>  
# Thlaspi_arvense<13>     
# Rosa_rugosa<14> 
# Argentina_anserina<15>  
# Fragaria_vesca<16>    
# Dryas_octopetala<17>     
# Prunus_persica<18>      
# Malus_sylvestris<19>    
# Pyrus_bretschneideri<20>      
# <21>     <22>    <23>    <24>    <25>    <26>    <27>    <28>    <29>    <30>    <31>    <32>    <33>  <34>     <35>    <36>    <37>    <38>    <39>

#get a list of the significantly expanded/contracted/rapid genes in Polygunum_aviculare_H0<1>
printf '%s\n' {1..50} | while read number ; do grep "Polygunum_aviculare_H0<1>\*" ../output/error_model_$number/Base_asr.tre | \
sed 's/^.*OG/OG/g' | sed 's/ = .*$//g' > Polavi_sig_changes_error_model_"$number" ; done 

#get a list of the significantly expanded/contracted/rapid genes in Fagopyrum_escelentum_H2<6>
printf '%s\n' {1..50} | while read number ; do grep "Fagopyrum_escelentum_H2<6>\*" ../output/error_model_$number/Base_asr.tre | \
sed 's/^.*OG/OG/g' | sed 's/ = .*$//g' > Fagesc_sig_changes_error_model_"$number" ; done 

#get a list of the significantly expanded/contracted/rapid genes in Brassica_oleracea<12>
printf '%s\n' {1..50} | while read number ; do grep "Brassica_oleracea<12>\*" ../output/error_model_$number/Base_asr.tre | \
sed 's/^.*OG/OG/g' | sed 's/ = .*$//g' > Brassica_sig_changes_error_model_"$number" ; done 

#get a list of the significantly expanded/contracted/rapid genes in Rosa_rugosa<14>
printf '%s\n' {1..50} | while read number ; do grep "Rosa_rugosa<14>\*" ../output/error_model_$number/Base_asr.tre | \
sed 's/^.*OG/OG/g' | sed 's/ = .*$//g' > Rosa_sig_changes_error_model_"$number" ; done 

#get a list of the significantly expanded/contracted/rapid genes in Thlaspi_arvense<13>
printf '%s\n' {1..50} | while read number ; do grep "Thlaspi_arvense<13>\*" ../output/error_model_$number/Base_asr.tre | \
sed 's/^.*OG/OG/g' | sed 's/ = .*$//g' > Thlarv_sig_changes_error_model_"$number" ; done 

#------------------------
#Data exploration:
wc -l Polavi_sig_changes_error_model_* # 66-71
wc -l Fagesc_sig_changes_error_model_* # 285-298
wc -l Brassica_sig_changes_error_model_* # 259-273
wc -l Rosa_sig_changes_error_model_* # 174-185
wc -l Thlarv_sig_changes_error_model_* # 161-168


#Data exploration:
wc -l Oxydig_sig_changes_error_model_* # 123-131
wc -l Rheumnob_sig_changes_error_model_* # 169-185
wc -l Rheum_Oxy_sig_changes_error_model_* # 87-93
wc -l Drabaniv_sig_changes_error_model_* # 157-163
wc -l Dryasoct_sig_changes_error_model_* # 170-180
wc -l Arabalp_sig_changes_error_model_* # 174-186

##################################################################################################################
#List of expanded and contracted groups (both significant and not)

#Get list of families that have expanded in Oxyria, whether or not significant

printf '%s\n' {1..50} | while read number ; do cut -f 1,2 ../output/error_model_$number/Base_change.tab | \
awk -F'\t' '$2 != 0' | grep -v "-" > Polavi_expandedfams_error_model_"$number" ; done 

printf '%s\n' {1..50} | while read number ; do cut -f 1,7 ../output/error_model_$number/Base_change.tab | \
awk -F'\t' '$2 != 0' | grep -v "-" > Fagesc_expandedfams_error_model_"$number" ; done 

printf '%s\n' {1..50} | while read number ; do cut -f 1,13 ../output/error_model_$number/Base_change.tab | \
awk -F'\t' '$2 != 0' | grep -v "-" > Brassica_expandedfams_error_model_"$number" ; done 

printf '%s\n' {1..50} | while read number ; do cut -f 1,15 ../output/error_model_$number/Base_change.tab | \
awk -F'\t' '$2 != 0' | grep -v "-" > Rosa_expandedfams_error_model_"$number" ; done 

printf '%s\n' {1..50} | while read number ; do cut -f 1,14 ../output/error_model_$number/Base_change.tab | \
awk -F'\t' '$2 != 0' | grep -v "-" > Thlarv_expandedfams_error_model_"$number" ; done 


wc -l *expandedfams_error_model_1 #remember the header is one line!

  # 6277 Brassica_expandedfams_error_model_1
  # 2727 Fagesc_expandedfams_error_model_1
  # 1598 Polavi_expandedfams_error_model_1
  # 1613 Rosa_expandedfams_error_model_1
   # 988 Thlarv_expandedfams_error_model_1

   # 523 Arabalp_expandedfams_error_model_1
  # 1627 Drabaniv_expandedfams_error_model_1
  # 1258 Dryasoct_expandedfams_error_model_1
  # 1481 Oxydig_expandedfams_error_model_1
  # 1924 Rheumnob_expandedfams_error_model_1
  # 1109 Rheum_Oxy_expandedfams_error_model_1

#--------------------------------------
#Get list of families that have contracted in Oxyria, whether or not significant
printf '%s\n' {1..50} | while read number ; do cut -f 1,3 ../output/error_model_$number/Base_change.tab | \
awk -F'\t' '$2 != 0' | grep "-" > Polavi_contractedfams_error_model_"$number" ; done

printf '%s\n' {1..50} | while read number ; do cut -f 1,5 ../output/error_model_$number/Base_change.tab | \
awk -F'\t' '$2 != 0' | grep "-" > Fagesc_contractedfams_error_model_"$number" ; done 

printf '%s\n' {1..50} | while read number ; do cut -f 1,26 ../output/error_model_$number/Base_change.tab | \
awk -F'\t' '$2 != 0' | grep "-" > Rheum_Oxy_contractedfams_error_model_"$number" ; done 

printf '%s\n' {1..50} | while read number ; do cut -f 1,11 ../output/error_model_$number/Base_change.tab | \
awk -F'\t' '$2 != 0' | grep "-" > Brassica_contractedfams_error_model_"$number" ; done

printf '%s\n' {1..50} | while read number ; do cut -f 1,18 ../output/error_model_$number/Base_change.tab | \
awk -F'\t' '$2 != 0' | grep "-" > Rosa_contractedfams_error_model_"$number" ; done 

printf '%s\n' {1..50} | while read number ; do cut -f 1,12 ../output/error_model_$number/Base_change.tab | \
awk -F'\t' '$2 != 0' | grep "-" > Thlarv_contractedfams_error_model_"$number" ; done 

wc -l *contractedfams_error_model_1 #remember the header is one line!

  # 1017 Brassica_contractedfams_error_model_1
  # 1280 Fagesc_contractedfams_error_model_1
  # 3523 Polavi_contractedfams_error_model_1
   # 972 Rheum_Oxy_contractedfams_error_model_1
  # 1961 Rosa_contractedfams_error_model_1
  # 4437 Thlarv_contractedfams_error_model_1


  # 4437 Arabalp_contractedfams_error_model_1
  # 1017 Drabaniv_contractedfams_error_model_1
  # 1961 Dryasoct_contractedfams_error_model_1
  # 3523 Oxydig_contractedfams_error_model_1
  # 1280 Rheumnob_contractedfams_error_model_1
   # 972 Rheum_Oxy_contractedfams_error_model_1


##################################################################################################################
#Running python script:
module load StdEnv/2020 python/3.9.6 scipy-stack/2021a
python

import os
from os.path import join as jn
import numpy as np

folderpath= os.getcwd()

#############################
#----------------
nameroot='Polavi_sig_changes_error_model_'
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
export= open(jn(folderpath,'output_Polavi'),'w')
export.write(superstring)
export.close()
#---------------------
nameroot='Fagesc_sig_changes_error_model_'
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
export= open(jn(folderpath,'output_Fagesc'),'w')
export.write(superstring)
export.close()

#----------------
nameroot='Brassica_sig_changes_error_model_'
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
export= open(jn(folderpath,'output_Brassica'),'w')
export.write(superstring)
export.close()
#---------------------
nameroot='Rosa_sig_changes_error_model_'
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
export= open(jn(folderpath,'output_Rosa'),'w')
export.write(superstring)
export.close()

#---------------------
nameroot='Thlarv_sig_changes_error_model_'
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
export= open(jn(folderpath,'output_Thlarv'),'w')
export.write(superstring)
export.close()
####################################
#----------------
nameroot='Polavi_expandedfams_error_model_'
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
export= open(jn(folderpath,'output_Polavi_expandedfams'),'w')
export.write(superstring)
export.close()

#---------------------
nameroot='Fagesc_expandedfams_error_model_'
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
export= open(jn(folderpath,'output_Fagesc_expandedfams'),'w')
export.write(superstring)
export.close()

#----------------
nameroot='Brassica_expandedfams_error_model_'
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
export= open(jn(folderpath,'output_Brassica_expandedfams'),'w')
export.write(superstring)
export.close()

#---------------------
nameroot='Rosa_expandedfams_error_model_'
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
export= open(jn(folderpath,'output_Rosa_expandedfams'),'w')
export.write(superstring)
export.close()

#---------------------
nameroot='Thlarv_expandedfams_error_model_'
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
export= open(jn(folderpath,'output_Thlarv_expandedfams'),'w')
export.write(superstring)
export.close()

################################
#----------------
nameroot='Polavi_contractedfams_error_model_'
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
export= open(jn(folderpath,'output_Polavi_contractedfams'),'w')
export.write(superstring)
export.close()

#---------------------
nameroot='Fagesc_contractedfams_error_model_'
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
export= open(jn(folderpath,'output_Fagesc_contractedfams'),'w')
export.write(superstring)
export.close()

#----------------
nameroot='Brassica_contractedfams_error_model_'
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
export= open(jn(folderpath,'output_Brassica_contractedfams'),'w')
export.write(superstring)
export.close()

#---------------------
nameroot='Rosa_contractedfams_error_model_'
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
export= open(jn(folderpath,'output_Rosa_contractedfams'),'w')
export.write(superstring)
export.close()

#---------------------
nameroot='Thlarv_contractedfams_error_model_'
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
export= open(jn(folderpath,'output_Thlarv_contractedfams'),'w')
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
cat output_Polavi_expandedfams | sort | uniq > OGtotaltemp
cat OGtotaltemp | while read line ; do grep -c "$line" output_Polavi_expandedfams >> tempcounts ; done
paste OGtotaltemp tempcounts >> OGcountstemp
cat OGcountstemp | while read orthogroup counts ; do if [ "$counts" -gt 47 ] ; then echo "$orthogroup" >> Polavi_expandedfams ; fi ; done
rm -f OGcountstemp OGtotaltemp tempcounts output


cat output_Polavi_contractedfams | sort | uniq > OGtotaltemp
cat OGtotaltemp | while read line ; do grep -c "$line" output_Polavi_contractedfams >> tempcounts ; done
paste OGtotaltemp tempcounts >> OGcountstemp
cat OGcountstemp | while read orthogroup counts ; do if [ "$counts" -gt 47 ] ; then echo "$orthogroup" >> Polavi_contractedfams ; fi ; done
rm -f OGcountstemp OGtotaltemp tempcounts output


cat output_Polavi | sort | uniq > OGtotaltemp
cat OGtotaltemp | while read line ; do grep -c "$line" output_Polavi >> tempcounts ; done
paste OGtotaltemp tempcounts >> OGcountstemp
cat OGcountstemp | while read orthogroup counts ; do if [ "$counts" -gt 47 ] ; then echo "$orthogroup" >> Polavi_sig_changes ; fi ; done
rm -f OGcountstemp OGtotaltemp tempcounts output

#-------------------------------------------
#Rheum nobile
cat output_Fagesc_expandedfams | sort | uniq > OGtotaltemp
cat OGtotaltemp | while read line ; do grep -c "$line" output_Fagesc_expandedfams >> tempcounts ; done
paste OGtotaltemp tempcounts >> OGcountstemp
cat OGcountstemp | while read orthogroup counts ; do if [ "$counts" -gt 47 ] ; then echo "$orthogroup" >> Fagesc_expandedfams ; fi ; done
rm -f OGcountstemp OGtotaltemp tempcounts output 


cat output_Fagesc_contractedfams | sort | uniq > OGtotaltemp
cat OGtotaltemp | while read line ; do grep -c "$line" output_Fagesc_contractedfams >> tempcounts ; done
paste OGtotaltemp tempcounts >> OGcountstemp
cat OGcountstemp | while read orthogroup counts ; do if [ "$counts" -gt 47 ] ; then echo "$orthogroup" >> Fagesc_contractedfams ; fi ; done
rm -f OGcountstemp OGtotaltemp tempcounts output 


cat output_Fagesc | sort | uniq > OGtotaltemp
cat OGtotaltemp | while read line ; do grep -c "$line" output_Fagesc >> tempcounts ; done
paste OGtotaltemp tempcounts >> OGcountstemp
cat OGcountstemp | while read orthogroup counts ; do if [ "$counts" -gt 47 ] ; then echo "$orthogroup" >>Fagesc_sig_changes ; fi ; done
rm -f OGcountstemp OGtotaltemp tempcounts output

#--------------------------------------------------
# Draba nivalis
cat output_Brassica_expandedfams | sort | uniq > OGtotaltemp
cat OGtotaltemp | while read line ; do grep -c "$line" output_Brassica_expandedfams >> tempcounts ; done
paste OGtotaltemp tempcounts >> OGcountstemp
cat OGcountstemp | while read orthogroup counts ; do if [ "$counts" -gt 47 ] ; then echo "$orthogroup" >> Brassica_expandedfams ; fi ; done
rm -f OGcountstemp OGtotaltemp tempcounts output


cat output_Brassica_contractedfams | sort | uniq > OGtotaltemp
cat OGtotaltemp | while read line ; do grep -c "$line" output_Brassica_contractedfams >> tempcounts ; done
paste OGtotaltemp tempcounts >> OGcountstemp
cat OGcountstemp | while read orthogroup counts ; do if [ "$counts" -gt 47 ] ; then echo "$orthogroup" >> Brassica_contractedfams ; fi ; done
rm -f OGcountstemp OGtotaltemp tempcounts output


cat output_Brassica | sort | uniq > OGtotaltemp
cat OGtotaltemp | while read line ; do grep -c "$line" output_Brassica >> tempcounts ; done
paste OGtotaltemp tempcounts >> OGcountstemp
cat OGcountstemp | while read orthogroup counts ; do if [ "$counts" -gt 47 ] ; then echo "$orthogroup" >> Brassica_sig_changes ; fi ; done
rm -f OGcountstemp OGtotaltemp tempcounts output

#-------------------------------------------
# Dryas octopetala
cat output_Rosa_expandedfams | sort | uniq > OGtotaltemp
cat OGtotaltemp | while read line ; do grep -c "$line" output_Rosa_expandedfams >> tempcounts ; done
paste OGtotaltemp tempcounts >> OGcountstemp
cat OGcountstemp | while read orthogroup counts ; do if [ "$counts" -gt 47 ] ; then echo "$orthogroup" >> Rosa_expandedfams ; fi ; done
rm -f OGcountstemp OGtotaltemp tempcounts output 


cat output_Rosa_contractedfams | sort | uniq > OGtotaltemp
cat OGtotaltemp | while read line ; do grep -c "$line" output_Rosa_contractedfams >> tempcounts ; done
paste OGtotaltemp tempcounts >> OGcountstemp
cat OGcountstemp | while read orthogroup counts ; do if [ "$counts" -gt 47 ] ; then echo "$orthogroup" >> Rosa_contractedfams ; fi ; done
rm -f OGcountstemp OGtotaltemp tempcounts output 


cat output_Rosa | sort | uniq > OGtotaltemp
cat OGtotaltemp | while read line ; do grep -c "$line" output_Rosa >> tempcounts ; done
paste OGtotaltemp tempcounts >> OGcountstemp
cat OGcountstemp | while read orthogroup counts ; do if [ "$counts" -gt 47 ] ; then echo "$orthogroup" >>Rosa_sig_changes ; fi ; done
rm -f OGcountstemp OGtotaltemp tempcounts output

#-------------------------------------------
# Arabis alpina
cat output_Thlarv_expandedfams | sort | uniq > OGtotaltemp
cat OGtotaltemp | while read line ; do grep -c "$line" output_Thlarv_expandedfams >> tempcounts ; done
paste OGtotaltemp tempcounts >> OGcountstemp
cat OGcountstemp | while read orthogroup counts ; do if [ "$counts" -gt 47 ] ; then echo "$orthogroup" >> Thlarv_expandedfams ; fi ; done
rm -f OGcountstemp OGtotaltemp tempcounts output 


cat output_Thlarv_contractedfams | sort | uniq > OGtotaltemp
cat OGtotaltemp | while read line ; do grep -c "$line" output_Thlarv_contractedfams >> tempcounts ; done
paste OGtotaltemp tempcounts >> OGcountstemp
cat OGcountstemp | while read orthogroup counts ; do if [ "$counts" -gt 47 ] ; then echo "$orthogroup" >> Thlarv_contractedfams ; fi ; done
rm -f OGcountstemp OGtotaltemp tempcounts output 


cat output_Thlarv | sort | uniq > OGtotaltemp
cat OGtotaltemp | while read line ; do grep -c "$line" output_Thlarv >> tempcounts ; done
paste OGtotaltemp tempcounts >> OGcountstemp
cat OGcountstemp | while read orthogroup counts ; do if [ "$counts" -gt 47 ] ; then echo "$orthogroup" >>Thlarv_sig_changes ; fi ; done
rm -f OGcountstemp OGtotaltemp tempcounts output


######################################################################################################################
#See which  expanded and contracted orthogroups are also significant
for taxon in Polavi Fagesc Rosa Brassica Thlarv; do grep -f "$taxon"_sig_changes "$taxon"_expandedfams > "$taxon"_expandedfams.sig ; done
for taxon in Polavi Fagesc Rosa Brassica Thlarv; do grep -f "$taxon"_sig_changes "$taxon"_contractedfams > "$taxon"_contractedfams.sig ; done

wc -l *_expandedfams*.sig  
  # 245 Brassica_expandedfams.sig
  # 237 Fagesc_expandedfams.sig
   # 52 Polavi_expandedfams.sig
  # 168 Rosa_expandedfams.sig
   # 58 Thlarv_expandedfams.sig

wc -l *_contractedfams*.sig

  # 56 Brassica_contractedfams.sig
  # 43 Fagesc_contractedfams.sig
  # 17 Polavi_contractedfams.sig
  # 41 Rosa_contractedfams.sig
  # 67 Thlarv_contractedfams.sig


  
 #################################
wc -l *_expandedfams*.sig  
  # 39 Arabalp_expandedfams.sig
 # 128 Drabaniv_expandedfams.sig
 # 142 Dryasoct_expandedfams.sig
 # 106 Oxydig_expandedfams.sig
 # 111 Rheumnob_expandedfams.sig


wc -l *_contractedfams*.sig
# 133 Arabalp_contractedfams.sig
  # 25 Drabaniv_contractedfams.sig
  # 28 Dryasoct_contractedfams.sig
  # 15 Oxydig_contractedfams.sig
  # 53 Rheumnob_contractedfams.sig


#-----------------------------
# Polygonaceae
    # 265 Polavi_expandedfams.sig
  # 610 Fagesc_expandedfams.sig
    # 46 Polavi_contractedfams.sig
  # 67 Fagesc_contractedfams.sig
  # 92 Rheum_Oxy_contractedfams.sig


#######################################################################################################################################
#Make gene list in orthogroups expanded
mkdir ~/scratch/Oxyria/CAFE/expanded-non-arctic
cd ~/scratch/Oxyria/CAFE/expanded-non-arctic

cp ../analysis-non-arctic/Polavi_expandedfams.sig .
sed -i 's/N0.H//g' Polavi_expandedfams.sig 

cp ../analysis-non-arctic/Fagesc_expandedfams.sig .
sed -i 's/N0.H//g' Fagesc_expandedfams.sig 

cp ../analysis-non-arctic/Rosa_expandedfams.sig .
sed -i 's/N0.H//g' Rosa_expandedfams.sig 

cp ../analysis-non-arctic/Brassica_expandedfams.sig .
sed -i 's/N0.H//g' Brassica_expandedfams.sig 

cp ../analysis-non-arctic/Thlarv_expandedfams.sig .
sed -i 's/N0.H//g' Thlarv_expandedfams.sig 

#---------------------------------------------------
#get a fasta of the expanded families from the OrthoFinder output
for taxon in Polavi Fagesc Rosa Brassica Thlarv ; do cat "$taxon"_expandedfams.sig | cut -f 1 | while read Orthogroup ; do cp /home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/orthofinder/Results_Jul03/Orthogroup_Sequences/"$Orthogroup".fa .; done ; done

# see gene names in these orthogroups
grep  ">" /home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/tmp/Polygunum_aviculare_H0.fa | head
grep  ">" /home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/tmp/Fagopyrum_escelentum_H2.fa | head
grep  ">" /home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/tmp/Rosa_rugosa.fa | head
grep  ">" /home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/tmp/Brassica_oleracea.fa | head
grep  ">" /home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/tmp/Thlaspi_arvense.fa | head

cd ~/scratch/Oxyria/CAFE/expanded-non-arctic
#Get a list of the specific species genes in these orthogroups
for taxon in Polavi; do cat "$taxon"_expandedfams.sig | cut -f 1 | while read Orthogroup ; do grep "Polavi" "$Orthogroup".fa >> "$taxon"_expanded_geneIDs.txt ; done ; done
for taxon in Fagesc; do cat "$taxon"_expandedfams.sig | cut -f 1 | while read Orthogroup ; do grep "FEHAP" "$Orthogroup".fa >> "$taxon"_expanded_geneIDs.txt ; done ; done
for taxon in Rosa; do cat "$taxon"_expandedfams.sig | cut -f 1 | while read Orthogroup ; do grep "LOC133" "$Orthogroup".fa >> "$taxon"_expanded_geneIDs.txt ; done ; done
for taxon in Brassica; do cat "$taxon"_expandedfams.sig | cut -f 1 | while read Orthogroup ; do grep -e "106" -e "nad" -e "orf" -e "atp" -e "cox" "$Orthogroup".fa >> "$taxon"_expanded_geneIDs.txt ; done ; done
for taxon in Thlarv; do cat "$taxon"_expandedfams.sig | cut -f 1 | while read Orthogroup ; do grep "TAV2" "$Orthogroup".fa >> "$taxon"_expanded_geneIDs.txt ; done ; done

for taxon in Polavi Fagesc Rosa Brassica Thlarv; do sed -i 's/>//g' "$taxon"_expanded_geneIDs.txt ; done

mkdir orthologues
mv OG*.fa orthologues


#############################################
#Make gene list in orthogroups contracted

mkdir ~/scratch/Oxyria/CAFE/contracted-non-arctic
cd ~/scratch/Oxyria/CAFE/contracted-non-arctic

cp ../analysis-non-arctic/Polavi_contractedfams.sig .
sed -i 's/N0.H//g' Polavi_contractedfams.sig 

cp ../analysis-non-arctic/Fagesc_contractedfams.sig .
sed -i 's/N0.H//g' Fagesc_contractedfams.sig 

cp ../analysis-non-arctic/Rosa_contractedfams.sig .
sed -i 's/N0.H//g' Rosa_contractedfams.sig 

cp ../analysis-non-arctic/Brassica_contractedfams.sig .
sed -i 's/N0.H//g' Brassica_contractedfams.sig 

cp ../analysis-non-arctic/Thlarv_contractedfams.sig .
sed -i 's/N0.H//g' Thlarv_contractedfams.sig 


#---------------------------------------------------
#get a fasta of the contracted families from the OrthoFinder output
for taxon in Polavi Fagesc Rosa Brassica Thlarv ; do cat "$taxon"_contractedfams.sig | cut -f 1 | while read Orthogroup ; do cp /home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/orthofinder/Results_Jul03/Orthogroup_Sequences/"$Orthogroup".fa .; done ; done

# see gene names in these orthogroups
grep  ">" /home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/tmp/Polygunum_aviculare_H0.fa | head
grep  ">" /home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/tmp/Fagopyrum_escelentum_H2.fa | head
grep  ">" /home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/tmp/Rosa_rugosa.fa | head
grep  ">" /home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/tmp/Brassica_oleracea.fa | head
grep  ">" /home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/tmp/Thlaspi_arvense.fa | head

cd ~/scratch/Oxyria/CAFE/contracted-non-arctic
#Get a list of the specific species genes in these orthogroups
for taxon in Polavi; do cat "$taxon"_contractedfams.sig | cut -f 1 | while read Orthogroup ; do grep "Polavi" "$Orthogroup".fa >> "$taxon"_contracted_geneIDs.txt ; done ; done
for taxon in Fagesc; do cat "$taxon"_contractedfams.sig | cut -f 1 | while read Orthogroup ; do grep "FEHAP" "$Orthogroup".fa >> "$taxon"_contracted_geneIDs.txt ; done ; done
for taxon in Rosa; do cat "$taxon"_contractedfams.sig | cut -f 1 | while read Orthogroup ; do grep "LOC133" "$Orthogroup".fa >> "$taxon"_contracted_geneIDs.txt ; done ; done
for taxon in Brassica; do cat "$taxon"_contractedfams.sig | cut -f 1 | while read Orthogroup ; do grep -e "106" -e "nad" -e "orf" -e "atp" -e "cox" "$Orthogroup".fa >> "$taxon"_contracted_geneIDs.txt ; done ; done
for taxon in Thlarv; do cat "$taxon"_contractedfams.sig | cut -f 1 | while read Orthogroup ; do grep "TAV2" "$Orthogroup".fa >> "$taxon"_contracted_geneIDs.txt ; done ; done

for taxon in Polavi Fagesc Rosa Brassica Thlarv; do sed -i 's/>//g' "$taxon"_contracted_geneIDs.txt ; done

mkdir orthologues
mv OG*.fa orthologues
#######################################################################################################################################
#Make intersection files

# Take a look at orthogroups enriched in a non-Arctic species 
# remove these?? that have overlap with Arctic?

###########################################
# try total intersection 
cd  ~/scratch/Oxyria/CAFE/expanded-non-arctic
sort Rosa_expandedfams.sig Polavi_expandedfams.sig Fagesc_expandedfams.sig Brassica_expandedfams.sig Thlarv_expandedfams.sig | uniq -cd > all_spp_expandedfams.sig

for taxon in allspp; do cat "$taxon"_expandedfams.sig | cut -f 1 | while read Orthogroup ; do grep -e "Oxyria" -e "Rno" -e "\-lg" -e "AALP" -e "DoctH0" "$Orthogroup".fa >> "$taxon"_expanded_geneIDs.txt ; done ; done
for taxon in allspp ; do sed -i 's/>//g' "$taxon"_expanded_geneIDs.txt ; done

#----------------------------
cd  ~/scratch/Oxyria/CAFE/contracted-non-arctic
sort Brassica_contractedfams.sig Polavi_contractedfams.sig Fagesc_contractedfams.sig  Thlarv_contractedfams.sig Rosa_contractedfams.sig | uniq -cd > all_spp_contractedfams.sig

      # 2 OG0000166
      # 2 OG0000178
      # 3 OG0000301
      # 2 OG0000433
      # 2 OG0000854
      # 2 OG0000858
      # 2 OG0000859
      # 2 OG0001096
      # 2 OG0001317
      # 2 OG0002083
      # 2 OG0002237
      # 2 OG0002248
      # 2 OG0002344
      # 2 OG0002445
      # 2 OG0009023
      # 2 OG0009910
      # 2 OG0011422

for taxon in allspp ; do cat "$taxon"_contractedfams.sig | cut -f 1 | while read Orthogroup ; do grep -e "Oxyria" -e "Rno" -e "\-lg" -e "AALP" -e "DoctH0" "$Orthogroup".fa >> "$taxon"_contracted_geneIDs.txt ; done ; done
for taxon in allspp ; do sed -i 's/>//g' "$taxon"_contracted_geneIDs.txt ; done

############################################
# Explore
# all - Polavi Fagesc Rosa Brassica Thlarv
cd  ~/scratch/Oxyria/CAFE/expanded-non-arctic
grep -wf Polavi_expandedfams.sig Fagesc_expandedfams.sig Brassica_expandedfams.sig Thlarv_expandedfams.sig Rosa_expandedfams.sig > Polavi_allspp_expandedfams.sig
grep -wf Brassica_expandedfams.sig Polavi_expandedfams.sig Fagesc_expandedfams.sig  Thlarv_expandedfams.sig Rosa_expandedfams.sig > Brassica_allspp_expandedfams.sig
grep -wf Rosa_expandedfams.sig Polavi_expandedfams.sig Fagesc_expandedfams.sig Brassica_expandedfams.sig Thlarv_expandedfams.sig  > Rosa_allspp_expandedfams.sig


#--------------------------------
cd  ~/scratch/Oxyria/CAFE/contracted-non-arctic
grep -f Polavi_contractedfams.sig Fagesc_contractedfams.sig Brassica_contractedfams.sig Thlarv_contractedfams.sig Rosa_contractedfams.sig > Polavi_allspp_contractedfams.sig
grep -f Rosa_contractedfams.sig Polavi_contractedfams.sig Fagesc_contractedfams.sig Brassica_contractedfams.sig Thlarv_contractedfams.sig  > Rosa_allspp_contractedfams.sig
grep -f Brassica_contractedfams.sig Polavi_contractedfams.sig Fagesc_contractedfams.sig  Thlarv_contractedfams.sig Rosa_contractedfams.sig > Brassica_allspp_contractedfams.sig

#------------------------------
# more than 2 duplicates
sort Rosa_expandedfams.sig Polavi_expandedfams.sig Fagesc_expandedfams.sig Brassica_expandedfams.sig Thlarv_expandedfams.sig |\
awk '{ a[$0]++ } END{ for(x in a) if(a[x]>2) print a[x], x }'
# 3 OG0000703  * Draba, Rheum, Dryas share
# 3 OG0002099  * Draba, Rheum, Dryas share

# non -arctic share more 
# 3 OG0001662
# 3 OG0000433
# 3 OG0000163
# 3 OG0001091
# 3 OG0001713
# 3 OG0001449
# 3 OG0000301  * Rosa, Brassica, Polavi
# 3 OG0002344
# 3 OG0000858

#-----------------------
# Summary 
# 444 orthogroups across all
# 2 overlap with 3 species
# 80 overlap with one other species

cd  ~/scratch/Oxyria/CAFE/expanded-non-arctic
wc -l *.sig
  # 129 all_spp_expandedfams.sig
   # 75 Brassica_allspp_expandedfams.sig
  # 245 Brassica_expandedfams.sig
  # 237 Fagesc_expandedfams.sig
   # 34 Polavi_allspp_expandedfams.sig
   # 52 Polavi_expandedfams.sig
   # 74 Rosa_allspp_expandedfams.sig
  # 168 Rosa_expandedfams.sig
   # 58 Thlarv_expandedfams.sig
   
   
# 30%, 65%, 44% of genes overlap

#--------------------------------
cd  ~/scratch/Oxyria/CAFE/contracted-non-arctic
wc -l *.sig

  # 17 all_spp_contractedfams.sig
   # 7 Brassica_allspp_contractedfams.sig
  # 56 Brassica_contractedfams.sig
  # 43 Fagesc_contractedfams.sig
   # 9 Polavi_allspp_contractedfams.sig
  # 17 Polavi_contractedfams.sig
   # 8 Rosa_allspp_contractedfams.sig
  # 41 Rosa_contractedfams.sig
  # 67 Thlarv_contractedfams.sig

# Combine Brassica and Thlarv


########################
# Check between Arctic and non-Arctic

sort ~/scratch/Oxyria/CAFE/expanded-non-arctic/all_spp_expandedfams.sig ~/scratch/Oxyria/CAFE/expanded/all_spp_expandedfams.sig | awk   '{print $2}' | sort | uniq -cd > Arctic_non-Arctic_expandedfams.sig

      # 2 OG0000149
      # 2 OG0000308
      # 2 OG0000620
      # 2 OG0000627
      # 2 OG0001065
      # 2 OG0001091
      # 2 OG0001517
      # 2 OG0001912
      # 2 OG0006410
      # 2 OG0009771
      # 2 OG0015176

sort ~/scratch/Oxyria/CAFE/contracted-non-arctic/all_spp_contractedfams.sig ~/scratch/Oxyria/CAFE/contracted/all_spp_contractedfams.sig | awk   '{print $2}' | sort | uniq -cd > Arctic_non-Arctic_contractedfams.sig

      # 2 OG0002344
      # 2 OG0009023
      # 2 OG0009910

for taxon in Oxydig ; do cat "$taxon"_contractedfams.sig | cut -f 1 | while read Orthogroup ; do grep -e "Oxyria" -e "Rno" -e "\-lg" -e "AALP" -e "DoctH0" "$Orthogroup".fa >> "$taxon"_contracted_geneIDs.txt ; done ; done
for taxon in allspp ; do sed -i 's/>//g' "$taxon"_contracted_geneIDs.txt ; done

#####################################################################################################################
#From here: continue at Notes5B_May2024_Go_enrichment.bash



