#################################################################################################################
# Brassicaceae expanded and contract gene families
# Sept 2024
#################################################################################################################

tmux new-session -s CAFE1
tmux attach-session -t CAFE1

cd /home/celphin/scratch/Oxyria/CAFE/Brassicaceae_data/analysis

#Check what clade number the taxa of interest are
more ../output/error_model/Base_clade_results.txt

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

more /lustre04/scratch/celphin/Oxyria/CAFE/Brassicaceae_data/output/error_model/Base_asr.tre

 # TREE N0.HOG0000037 = ((Brassica_oleracea<1>_1:34.4543,Thlaspi_arvense<2>_2:34.4543)<8>_1:3.54574,
# (Cochlearia_groenlandica<3>_1:35.7726,(Arabidopsis_lyrata<4>_1:31.3765,(Draba_nivalis<5>_1:23.1667,
# Arabis_alpina<6>_1:23.1667)<11>_1:8.20977)<10>_1:4.39608)<9>_1:2.22743)<7>_1;

#--------------------------------
#get a list of the significantly expanded/contracted/rapid genes in Draba/Arabis clade <11>
printf '%s\n' {1..50} | while read number ; do grep "<11>\*" ../output/error_model_$number/Base_asr.tre | \
sed 's/^.*OG/OG/g' | sed 's/ = .*$//g' > Dra_Ara_sig_changes_error_model_"$number" ; done 

#get a list of the significantly expanded/contracted/rapid genes in Draba_nivalis<5>
printf '%s\n' {1..50} | while read number ; do grep "Draba_nivalis<5>\*" ../output/error_model_$number/Base_asr.tre | \
sed 's/^.*OG/OG/g' | sed 's/ = .*$//g' > Drabaniv_sig_changes_error_model_"$number" ; done 

#get a list of the significantly expanded/contracted/rapid genes in Cochlearia_groenlandica<3>
printf '%s\n' {1..50} | while read number ; do grep "Cochlearia_groenlandica<3>\*" ../output/error_model_$number/Base_asr.tre | \
sed 's/^.*OG/OG/g' | sed 's/ = .*$//g' > Cochgro_sig_changes_error_model_"$number" ; done 

#get a list of the significantly expanded/contracted/rapid genes in Arabis_alpina<6>
printf '%s\n' {1..50} | while read number ; do grep "Arabis_alpina<6>\*" ../output/error_model_$number/Base_asr.tre | \
sed 's/^.*OG/OG/g' | sed 's/ = .*$//g' > Arabalp_sig_changes_error_model_"$number" ; done 

#------------------------
#Data exploration:
wc -l Dra_Ara_sig_changes_error_model_* | head # 7-8
wc -l Drabaniv_sig_changes_error_model_* | head # 44-55
wc -l Cochgro_sig_changes_error_model_* | head # 23-25
wc -l Arabalp_sig_changes_error_model_* | head # 20-25

#---------------------
#see what is uniquely found in run 6 (this was the time I had no Hypoch/Rhemel)
printf '%s\n' {1..10} | while read number ; do grep -v -f Drabaniv_sig_changes_error_model_"$number" Drabaniv_sig_changes_error_model_10 ; done |wc -l
2
printf '%s\n' {1..10} | while read number ; do grep -v -f Drabaniv_sig_changes_error_model_"$number" Drabaniv_sig_changes_error_model_5 ; done |wc -l
10
printf '%s\n' {1..10} | while read number ; do grep -v -f Drabaniv_sig_changes_error_model_"$number" Drabaniv_sig_changes_error_model_1 ; done |wc -l
2

##################################################################################################################
#List of expanded and contracted groups (both significant and not)

#Get list of families that have expanded in each group, whether or not significant
# value should be one greater than ID to get right column

printf '%s\n' {1..50} | while read number ; do cut -f 1,12 ../output/error_model_$number/Base_change.tab | \
awk -F'\t' '$2 != 0' | grep -v "-" > Dra_Ara_expandedfams_error_model_"$number" ; done 

printf '%s\n' {1..50} | while read number ; do cut -f 1,6 ../output/error_model_$number/Base_change.tab | \
awk -F'\t' '$2 != 0' | grep -v "-" > Drabaniv_expandedfams_error_model_"$number" ; done 

printf '%s\n' {1..50} | while read number ; do cut -f 1,4 ../output/error_model_$number/Base_change.tab | \
awk -F'\t' '$2 != 0' | grep -v "-" > Cochgro_expandedfams_error_model_"$number" ; done 

printf '%s\n' {1..50} | while read number ; do cut -f 1,7 ../output/error_model_$number/Base_change.tab | \
awk -F'\t' '$2 != 0' | grep -v "-" > Arabalp_expandedfams_error_model_"$number" ; done 


wc -l *expandedfams_error_model_1 #remember the header is one line!

# This run Brassicaceae only
  # 363 Arabalp_expandedfams_error_model_1
 # 1183 Cochgro_expandedfams_error_model_1
   # 71 Dra_Ara_expandedfams_error_model_1
  # 810 Drabaniv_expandedfams_error_model_1

#--------------------------------------
#Get list of families that have contracted in Oxyria, whether or not significant

printf '%s\n' {1..50} | while read number ; do cut -f 1,12 ../output/error_model_$number/Base_change.tab | \
awk -F'\t' '$2 != 0' | grep "-" > Dra_Ara_contractedfams_error_model_"$number" ; done 

printf '%s\n' {1..50} | while read number ; do cut -f 1,6 ../output/error_model_$number/Base_change.tab | \
awk -F'\t' '$2 != 0' | grep "-" > Drabaniv_contractedfams_error_model_"$number" ; done

printf '%s\n' {1..50} | while read number ; do cut -f 1,4 ../output/error_model_$number/Base_change.tab | \
awk -F'\t' '$2 != 0' | grep "-" > Cochgro_contractedfams_error_model_"$number" ; done 

printf '%s\n' {1..50} | while read number ; do cut -f 1,7 ../output/error_model_$number/Base_change.tab | \
awk -F'\t' '$2 != 0' | grep "-" > Arabalp_contractedfams_error_model_"$number" ; done 

wc -l *contractedfams_error_model_1 #remember the header is one line!

# This run Brassicaceae only
  # 189 Arabalp_contractedfams_error_model_1
  # 260 Cochgro_contractedfams_error_model_1
   # 75 Dra_Ara_contractedfams_error_model_1
   # 92 Drabaniv_contractedfams_error_model_1


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

exit()

##################################################################################################################
#Clean up:
mkdir sig_genes; mv *sig_changes* sig_genes
mkdir expandedfams_error_models; mv *expandedfams_error_model* expandedfams_error_models
mkdir contractedfams_error_models; mv *contractedfams_error_model* contractedfams_error_models


##############################################################
#Draba-Arabis node
cat output_Dra_Ara_expandedfams | sort | uniq > OGtotaltemp
cat OGtotaltemp | while read line ; do grep -c "$line" output_Dra_Ara_expandedfams >> tempcounts ; done
paste OGtotaltemp tempcounts >> OGcountstemp
cat OGcountstemp | while read orthogroup counts ; do if [ "$counts" -gt 47 ] ; then echo "$orthogroup" >> Dra_Ara_expandedfams ; fi ; done
rm -f OGcountstemp OGtotaltemp tempcounts output 

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
# Cochlearia groenlandica
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


#######################################################################################################################################
#Make gene list in orthogroups total expanded
mkdir ~/scratch/Oxyria/CAFE/Brassicaceae_data/expanded
cd ~/scratch/Oxyria/CAFE/Brassicaceae_data/expanded

cp ../analysis/Cochgro_expandedfams .
sed -i 's/N0.H//g' Cochgro_expandedfams 
wc -l Cochgro_expandedfams
# 1182 Cochgro_expandedfams

cp ../analysis/Drabaniv_expandedfams .
sed -i 's/N0.H//g' Drabaniv_expandedfams 
wc -l Drabaniv_expandedfams 
# 809 Drabaniv_expandedfams

cp ../analysis/Arabalp_expandedfams .
sed -i 's/N0.H//g' Arabalp_expandedfams 
wc -l Arabalp_expandedfams 
# 362 Arabalp_expandedfams

#---------------------------------------------------
#get a fasta of the total expanded families from the OrthoFinder output
for taxon in Cochgro Drabaniv Arabalp ; do cat "$taxon"_expandedfams | cut -f 1 | while read Orthogroup ; do cp /home/celphin/scratch/Oxyria/GeneSpace/Brassicaceae_genomes/orthofinder/Results_Oct09/Orthogroup_Sequences/"$Orthogroup".fa .; done ; done

# see gene names in these orthogroups
grep  ">" /home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/tmp/Draba_nivalis.fa | head
grep  ">" /home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/tmp/Arabis_alpina.fa | head
grep  ">" /home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/tmp/Cochlearia_groenlandica.fa | head

cd ~/scratch/Oxyria/CAFE/Brassicaceae_data/expanded
#Get a list of the specific species genes in these orthogroups
for taxon in Cochgro; do cat "$taxon"_expandedfams | cut -f 1 | while read Orthogroup ; do grep ">g" "$Orthogroup".fa >> "$taxon"_totalexpanded_geneIDs.txt ; done ; done
for taxon in Drabaniv; do cat "$taxon"_expandedfams | cut -f 1 | while read Orthogroup ; do grep "\-lg" "$Orthogroup".fa >> "$taxon"_totalexpanded_geneIDs.txt ; done ; done
for taxon in Arabalp; do cat "$taxon"_expandedfams | cut -f 1 | while read Orthogroup ; do grep "AALP" "$Orthogroup".fa >> "$taxon"_totalexpanded_geneIDs.txt ; done ; done

for taxon in Cochgro Drabaniv Arabalp; do sed -i 's/>//g' "$taxon"_totalexpanded_geneIDs.txt ; done

#############################################
#Make gene list in orthogroups total contracted
mkdir cd ~/scratch/Oxyria/CAFE/Brassicaceae_data/contracted
cd ~/scratch/Oxyria/CAFE/Brassicaceae_data/contracted

cp ../analysis/Cochgro_contractedfams .
sed -i 's/N0.H//g' Cochgro_contractedfams 
wc -l Cochgro_contractedfams 
#259 Cochgro_contractedfams

cp ../analysis/Drabaniv_contractedfams .
sed -i 's/N0.H//g' Drabaniv_contractedfams 
wc -l Drabaniv_contractedfams
# 91 Drabaniv_contractedfams

cp ../analysis/Arabalp_contractedfams .
sed -i 's/N0.H//g' Arabalp_contractedfams 
wc -l Arabalp_contractedfams 
#188 Arabalp_contractedfams

#---------------------------------------------------
#get a fasta of the contracted families from the OrthoFinder output
for taxon in Cochgro Drabaniv Arabalp ; do cat "$taxon"_contractedfams | cut -f 1 | while read Orthogroup ; do cp /home/celphin/scratch/Oxyria/GeneSpace/Brassicaceae_genomes/orthofinder/Results_Oct09/Orthogroup_Sequences/"$Orthogroup".fa .; done ; done

# see gene names in these orthogroups
grep  ">" /home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/tmp/Draba_nivalis.fa | head
grep  ">" /home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/tmp/Arabis_alpina.fa | head
grep  ">" /home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/tmp/Cochlearia_groenlandica.fa | head

cd ~/scratch/Oxyria/CAFE/Brassicaceae_data/contracted
#Get a list of the specific species genes in these orthogroups
for taxon in Cochgro; do cat "$taxon"_contractedfams | cut -f 1 | while read Orthogroup ; do grep ">g" "$Orthogroup".fa >> "$taxon"_totalcontracted_geneIDs.txt ; done ; done
for taxon in Drabaniv; do cat "$taxon"_contractedfams | cut -f 1 | while read Orthogroup ; do grep "\-lg" "$Orthogroup".fa >> "$taxon"_totalcontracted_geneIDs.txt ; done ; done
for taxon in Arabalp; do cat "$taxon"_contractedfams | cut -f 1 | while read Orthogroup ; do grep "AALP" "$Orthogroup".fa >> "$taxon"_totalcontracted_geneIDs.txt ; done ; done

for taxon in Cochgro Drabaniv Arabalp; do sed -i 's/>//g' "$taxon"_totalcontracted_geneIDs.txt ; done

######################################################################################################################
#See which  expanded and contracted orthogroups rapidly expanded
cd /home/celphin/scratch/Oxyria/CAFE/Brassicaceae_data/analysis
for taxon in Cochgro Drabaniv Arabalp Dra_Ara; do grep -f "$taxon"_sig_changes "$taxon"_expandedfams > "$taxon"_expandedfams.sig ; done
for taxon in Cochgro Drabaniv Arabalp Dra_Ara; do grep -f "$taxon"_sig_changes "$taxon"_contractedfams > "$taxon"_contractedfams.sig ; done

wc -l *_expandedfams*.sig  
  # 12 Arabalp_expandedfams.sig
  # 15 Cochgro_expandedfams.sig
   # 6 Dra_Ara_expandedfams.sig
  # 41 Drabaniv_expandedfams.sig

wc -l *_contractedfams*.sig
  # 7 Arabalp_contractedfams.sig
  # 7 Cochgro_contractedfams.sig
  # 0 Dra_Ara_contractedfams.sig
  # 1 Drabaniv_contractedfams.sig

#######################################################################################################################################
#Make gene list in orthogroups rapidly expanded

cd ~/scratch/Oxyria/CAFE/Brassicaceae_data/expanded

cp ../analysis/Cochgro_expandedfams.sig .
sed -i 's/N0.H//g' Cochgro_expandedfams.sig 

cp ../analysis/Drabaniv_expandedfams.sig .
sed -i 's/N0.H//g' Drabaniv_expandedfams.sig 

cp ../analysis/Arabalp_expandedfams.sig .
sed -i 's/N0.H//g' Arabalp_expandedfams.sig 

#---------------------------------------------------
#get a fasta of the expanded families from the OrthoFinder output
for taxon in Cochgro Drabaniv Arabalp ; do cat "$taxon"_expandedfams.sig | cut -f 1 | while read Orthogroup ; do cp /home/celphin/scratch/Oxyria/GeneSpace/Brassicaceae_genomes/orthofinder/Results_Oct09/Orthogroup_Sequences/"$Orthogroup".fa .; done ; done

# see gene names in these orthogroups
grep  ">" /home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/tmp/Draba_nivalis.fa | head
grep  ">" /home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/tmp/Arabis_alpina.fa | head
grep  ">" /home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/tmp/Cochlearia_groenlandica.fa | head

cd ~/scratch/Oxyria/CAFE/Brassicaceae_data/expanded
#Get a list of the specific species genes in these orthogroups
for taxon in Cochgro; do cat "$taxon"_expandedfams.sig | cut -f 1 | while read Orthogroup ; do grep ">g" "$Orthogroup".fa >> "$taxon"_expanded_geneIDs.txt ; done ; done
for taxon in Drabaniv; do cat "$taxon"_expandedfams.sig | cut -f 1 | while read Orthogroup ; do grep "\-lg" "$Orthogroup".fa >> "$taxon"_expanded_geneIDs.txt ; done ; done
for taxon in Arabalp; do cat "$taxon"_expandedfams.sig | cut -f 1 | while read Orthogroup ; do grep "AALP" "$Orthogroup".fa >> "$taxon"_expanded_geneIDs.txt ; done ; done

for taxon in Cochgro Drabaniv Arabalp; do sed -i 's/>//g' "$taxon"_expanded_geneIDs.txt ; done

#############################################
#Make gene list in orthogroups rapidly contracted

cd ~/scratch/Oxyria/CAFE/Brassicaceae_data/contracted

cp ../analysis/Cochgro_contractedfams.sig .
sed -i 's/N0.H//g' Cochgro_contractedfams.sig 
wc -l Cochgro_contractedfams.sig 
#7 Cochgro_contractedfams.sig

cp ../analysis/Drabaniv_contractedfams.sig .
sed -i 's/N0.H//g' Drabaniv_contractedfams.sig 
wc -l Drabaniv_contractedfams.sig 
#1 Drabaniv_contractedfams.sig

cp ../analysis/Arabalp_contractedfams.sig .
sed -i 's/N0.H//g' Arabalp_contractedfams.sig 
wc -l Arabalp_contractedfams.sig 
#7 Arabalp_contractedfams.sig

#---------------------------------------------------
#get a fasta of the contracted families from the OrthoFinder output
for taxon in Cochgro Drabaniv Arabalp ; do cat "$taxon"_contractedfams.sig | cut -f 1 | while read Orthogroup ; do cp /home/celphin/scratch/Oxyria/GeneSpace/Brassicaceae_genomes/orthofinder/Results_Oct09/Orthogroup_Sequences/"$Orthogroup".fa .; done ; done

# see gene names in these orthogroups
grep  ">" /home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/tmp/Draba_nivalis.fa | head
grep  ">" /home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/tmp/Arabis_alpina.fa | head
grep  ">" /home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/tmp/Cochlearia_groenlandica.fa | head

cd ~/scratch/Oxyria/CAFE/Brassicaceae_data/contracted
#Get a list of the specific species genes in these orthogroups
for taxon in Cochgro; do cat "$taxon"_contractedfams.sig | cut -f 1 | while read Orthogroup ; do grep ">g" "$Orthogroup".fa >> "$taxon"_contracted_geneIDs.txt ; done ; done
for taxon in Drabaniv; do cat "$taxon"_contractedfams.sig | cut -f 1 | while read Orthogroup ; do grep "\-lg" "$Orthogroup".fa >> "$taxon"_contracted_geneIDs.txt ; done ; done
for taxon in Arabalp; do cat "$taxon"_contractedfams.sig | cut -f 1 | while read Orthogroup ; do grep "AALP" "$Orthogroup".fa >> "$taxon"_contracted_geneIDs.txt ; done ; done

for taxon in Cochgro Drabaniv Arabalp; do sed -i 's/>//g' "$taxon"_contracted_geneIDs.txt ; done


###############################
cd  /home/celphin/scratch/Oxyria/CAFE/Brassicaceae_data/expanded

wc -l *geneIDs.txt
     # 7 Arabalp_expanded_geneIDs.txt
   # 636 Arabalp_totalexpanded_geneIDs.txt
    # 37 Cochgro_expanded_geneIDs.txt
  # 2574 Cochgro_totalexpanded_geneIDs.txt
    # 64 Drabaniv_expanded_geneIDs.txt


#----------------
cd  /home/celphin/scratch/Oxyria/CAFE/Brassicaceae_data/contracted

wc -l *geneIDs.txt
    # 7 Arabalp_contracted_geneIDs.txt
  # 187 Arabalp_totalcontracted_geneIDs.txt
    # 8 Cochgro_contracted_geneIDs.txt
  # 378 Cochgro_totalcontracted_geneIDs.txt
    # 1 Drabaniv_contracted_geneIDs.txt
  # 111 Drabaniv_totalcontracted_geneIDs.txt



#####################################################################################################################
#From here: continue at Notes6_May2024_Go_enrichment.bash
