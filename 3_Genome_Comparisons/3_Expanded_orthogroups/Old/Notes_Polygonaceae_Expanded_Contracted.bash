#################################################################################################################
# Polygonaceae expanded and contracted gene families
# Sept 2024
#################################################################################################################

tmux new-session -s CAFE
tmux attach-session -t CAFE

cd /home/celphin/scratch/Oxyria/CAFE/Polygonaceae_data/analysis

#Check what clade number the taxa of interest are
more ../output/error_model/Base_clade_results.txt

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


more /lustre04/scratch/celphin/Oxyria/CAFE/Polygonaceae_data/output/error_model/Base_asr.tre

# TREE N0.HOG0000018 = ((Fagopyrum_tataricum_H1<1>*_26:8.47061,Fagopyrum_escelentum_H2<2>*_1:8.4706
# 1)<8>*_5:17.765,(Polygunum_aviculare_H0<3>_1:20.0804,(Oxyria_digyna_H1<4>_1:16.0505,(Rheum_nobile_H
# 0<5>_0:6.703,Rheum_tangaticum_H0<6>_1:6.703)<11>_1:9.34751)<10>_1:4.02989)<9>_1:6.15523)<7>_2;

#-----------------------------------

#get a list of the significantly expanded/contracted/rapid genes in Oxyria_digyna_H1<4>
printf '%s\n' {1..50} | while read number ; do grep "Oxyria_digyna_H1<4>\*" ../output/error_model_$number/Base_asr.tre | \
sed 's/^.*OG/OG/g' | sed 's/ = .*$//g' > Oxydig_sig_changes_error_model_"$number" ; done 

#get a list of the significantly expanded/contracted/rapid genes in Rheum_nobile_H0<5>
printf '%s\n' {1..50} | while read number ; do grep "Rheum_nobile_H0<5>\*" ../output/error_model_$number/Base_asr.tre | \
sed 's/^.*OG/OG/g' | sed 's/ = .*$//g' > Rheumnob_sig_changes_error_model_"$number" ; done 

#------------------------
#Data exploration:
wc -l Oxydig_sig_changes_error_model_* # 102-110
wc -l Rheumnob_sig_changes_error_model_* # 247-327

#---------------------
#see what is uniquely found in run 6 (this was the time I had no Hypoch/Rhemel)
printf '%s\n' {1..10} | while read number ; do grep -v -f Oxydig_sig_changes_error_model_"$number" Oxydig_sig_changes_error_model_10 ; done |wc -l
# 33
printf '%s\n' {1..10} | while read number ; do grep -v -f Oxydig_sig_changes_error_model_"$number" Oxydig_sig_changes_error_model_1 ; done |wc -l
# 0

##################################################################################################################
#List of expanded and contracted groups (both significant and not)

#Get list of families that have expanded in Oxyria, whether or not significant
cd /home/celphin/scratch/Oxyria/CAFE/Polygonaceae_data/analysis

printf '%s\n' {1..50} | while read number ; do cut -f 1,5 ../output/error_model_$number/Base_change.tab | \
awk -F'\t' '$2 != 0' | grep -v "-" > Oxydig_expandedfams_error_model_"$number" ; done 

printf '%s\n' {1..50} | while read number ; do cut -f 1,6 ../output/error_model_$number/Base_change.tab | \
awk -F'\t' '$2 != 0' | grep -v "-" > Rheumnob_expandedfams_error_model_"$number" ; done 

wc -l *expandedfams_error_model_1 #remember the header is one line!

  # 597 Oxydig_expandedfams_error_model_1
 # 1037 Rheumnob_expandedfams_error_model_1


#--------------------------------------
#Get list of families that have contracted in Oxyria, whether or not significant
# value should be one greater than ID to get right column

printf '%s\n' {1..50} | while read number ; do cut -f 1,5 ../output/error_model_$number/Base_change.tab | \
awk -F'\t' '$2 != 0' | grep "-" > Oxydig_contractedfams_error_model_"$number" ; done

printf '%s\n' {1..50} | while read number ; do cut -f 1,6 ../output/error_model_$number/Base_change.tab | \
awk -F'\t' '$2 != 0' | grep "-" > Rheumnob_contractedfams_error_model_"$number" ; done 

wc -l *contractedfams_error_model_1 #remember the header is one line!

 # 367 Oxydig_contractedfams_error_model_1
 # 192 Rheumnob_contractedfams_error_model_1

##################################################################################################################
#Running python script:
cd /home/celphin/scratch/Oxyria/CAFE/Polygonaceae_data/analysis

module load StdEnv/2020 python/3.9.6 scipy-stack/2021a
python

import os
from os.path import join as jn
import numpy as np

folderpath= os.getcwd()

#############################

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


####################################

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

################################

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

#######################################################################################################################################
#Make gene list of total orthogroups expanded
cd ~/scratch/Oxyria/CAFE/Polygonaceae_data/expanded

cp ../analysis/Oxydig_expandedfams .
sed -i 's/N0.H//g' Oxydig_expandedfams 
wc -l Oxydig_expandedfams 
#596 Oxydig_expandedfams

cp ../analysis/Rheumnob_expandedfams .
sed -i 's/N0.H//g' Rheumnob_expandedfams 
wc -l Rheumnob_expandedfams 
#1036 Rheumnob_expandedfams

#---------------------------------------------------
cd ~/scratch/Oxyria/CAFE/Polygonaceae_data/expanded

#get a fasta of the expanded families from the OrthoFinder output
for taxon in Oxydig Rheumnob ; do cat "$taxon"_expandedfams | cut -f 1 | while read Orthogroup ; do cp /home/celphin/scratch/Oxyria/GeneSpace/Polygonaceae_genomes/orthofinder/Results_Oct07/Orthogroup_Sequences/"$Orthogroup".fa .; done ; done

# see gene names in these orthogroups
grep  ">" /home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/tmp/Rheum_nobile_H0.fa | head
grep  ">" /home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/tmp/Oxyria_digyna_H1.fa | head

cd ~/scratch/Oxyria/CAFE/Polygonaceae_data/expanded
#Get a list of the specific species genes in these orthogroups
for taxon in Oxydig; do cat "$taxon"_expandedfams | cut -f 1 | while read Orthogroup ; do grep "Oxyria" "$Orthogroup".fa >> "$taxon"_totalexpanded_geneIDs.txt ; done ; done
for taxon in Rheumnob; do cat "$taxon"_expandedfams | cut -f 1 | while read Orthogroup ; do grep "Rno" "$Orthogroup".fa >> "$taxon"_totalexpanded_geneIDs.txt ; done ; done

for taxon in Oxydig Rheumnob ; do sed -i 's/>//g' "$taxon"_totalexpanded_geneIDs.txt ; done

#############################################
#Make gene list in total orthogroups contracted

cd ~/scratch/Oxyria/CAFE/Polygonaceae_data/contracted

cp ../analysis/Oxydig_contractedfams .
sed -i 's/N0.H//g' Oxydig_contractedfams 
wc -l Oxydig_contractedfams 

cp ../analysis/Rheumnob_contractedfams .
sed -i 's/N0.H//g' Rheumnob_contractedfams 
wc -l Rheumnob_contractedfams

#---------------------------------------------------
#get a fasta of the contracted families from the OrthoFinder output
for taxon in Oxydig Rheumnob ; do cat "$taxon"_contractedfams | cut -f 1 | while read Orthogroup ; do cp /home/celphin/scratch/Oxyria/GeneSpace/Polygonaceae_genomes/orthofinder/Results_Oct07/Orthogroup_Sequences/"$Orthogroup".fa .; done ; done

# see gene names in these orthogroups
grep  ">" /home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/tmp/Rheum_nobile_H0.fa | head
grep  ">" /home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/tmp/Oxyria_digyna_H1.fa | head

cd ~/scratch/Oxyria/CAFE/Polygonaceae_data/contracted
#Get a list of the specific species genes in these orthogroups
for taxon in Oxydig; do cat "$taxon"_contractedfams | cut -f 1 | while read Orthogroup ; do grep "Oxyria" "$Orthogroup".fa >> "$taxon"_totalcontracted_geneIDs.txt ; done ; done
for taxon in Rheumnob; do cat "$taxon"_contractedfams | cut -f 1 | while read Orthogroup ; do grep "Rno" "$Orthogroup".fa >> "$taxon"_totalcontracted_geneIDs.txt ; done ; done

for taxon in Oxydig Rheumnob ; do sed -i 's/>//g' "$taxon"_totalcontracted_geneIDs.txt ; done


######################################################################################################################
#See which  expanded and contracted orthogroups are also significant for rapid expansion/contraction
cd /home/celphin/scratch/Oxyria/CAFE/Polygonaceae_data/analysis
for taxon in Oxydig Rheumnob ; do grep -f "$taxon"_sig_changes "$taxon"_expandedfams > "$taxon"_expandedfams.sig ; done
for taxon in Oxydig Rheumnob ; do grep -f "$taxon"_sig_changes "$taxon"_contractedfams > "$taxon"_contractedfams.sig ; done

wc -l *_expandedfams*.sig  
  # 92 Oxydig_expandedfams.sig
 # 241 Rheumnob_expandedfams.sig


wc -l *_contractedfams*.sig
  # 9 Oxydig_contractedfams.sig
  # 6 Rheumnob_contractedfams.sig


#######################################################################################################################################
#Make gene list of rapid orthogroups expanded
cd ~/scratch/Oxyria/CAFE/Polygonaceae_data/expanded

cp ../analysis/Oxydig_expandedfams.sig .
sed -i 's/N0.H//g' Oxydig_expandedfams.sig 

cp ../analysis/Rheumnob_expandedfams.sig .
sed -i 's/N0.H//g' Rheumnob_expandedfams.sig 

#---------------------------------------------------
cd ~/scratch/Oxyria/CAFE/Polygonaceae_data/expanded

#get a fasta of the expanded families from the OrthoFinder output
for taxon in Oxydig Rheumnob ; do cat "$taxon"_expandedfams.sig | cut -f 1 | while read Orthogroup ; do cp /home/celphin/scratch/Oxyria/GeneSpace/Polygonaceae_genomes/orthofinder/Results_Oct07/Orthogroup_Sequences/"$Orthogroup".fa .; done ; done

# see gene names in these orthogroups
grep  ">" /home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/tmp/Rheum_nobile_H0.fa | head
grep  ">" /home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/tmp/Oxyria_digyna_H1.fa | head

cd ~/scratch/Oxyria/CAFE/Polygonaceae_data/expanded
#Get a list of the specific species genes in these orthogroups
for taxon in Oxydig; do cat "$taxon"_expandedfams.sig | cut -f 1 | while read Orthogroup ; do grep "Oxyria" "$Orthogroup".fa >> "$taxon"_expanded_geneIDs.txt ; done ; done
for taxon in Rheumnob; do cat "$taxon"_expandedfams.sig | cut -f 1 | while read Orthogroup ; do grep "Rno" "$Orthogroup".fa >> "$taxon"_expanded_geneIDs.txt ; done ; done

for taxon in Oxydig Rheumnob ; do sed -i 's/>//g' "$taxon"_expanded_geneIDs.txt ; done

#############################################
#Make gene list of rapid orthogroups contracted

cd ~/scratch/Oxyria/CAFE/Polygonaceae_data/contracted

cp ../analysis/Oxydig_contractedfams.sig .
sed -i 's/N0.H//g' Oxydig_contractedfams.sig 

cp ../analysis/Rheumnob_contractedfams.sig .
sed -i 's/N0.H//g' Rheumnob_contractedfams.sig 


#---------------------------------------------------
#get a fasta of the contracted families from the OrthoFinder output
for taxon in Oxydig Rheumnob ; do cat "$taxon"_contractedfams.sig | cut -f 1 | while read Orthogroup ; do cp /home/celphin/scratch/Oxyria/GeneSpace/Polygonaceae_genomes/orthofinder/Results_Oct07/Orthogroup_Sequences/"$Orthogroup".fa .; done ; done

# see gene names in these orthogroups
grep  ">" /home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/tmp/Rheum_nobile_H0.fa | head
grep  ">" /home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/tmp/Oxyria_digyna_H1.fa | head

cd ~/scratch/Oxyria/CAFE/Polygonaceae_data/contracted
#Get a list of the specific species genes in these orthogroups
for taxon in Oxydig; do cat "$taxon"_contractedfams.sig | cut -f 1 | while read Orthogroup ; do grep "Oxyria" "$Orthogroup".fa >> "$taxon"_contracted_geneIDs.txt ; done ; done
for taxon in Rheumnob; do cat "$taxon"_contractedfams.sig | cut -f 1 | while read Orthogroup ; do grep "Rno" "$Orthogroup".fa >> "$taxon"_contracted_geneIDs.txt ; done ; done

for taxon in Oxydig Rheumnob ; do sed -i 's/>//g' "$taxon"_contracted_geneIDs.txt ; done


###############################
cd  ~/scratch/Oxyria/CAFE/Polygonaceae_data/expanded

wc -l *geneIDs.txt
  151 Oxydig_expanded_geneIDs.txt
  802 Oxydig_totalexpanded_geneIDs.txt
  408 Rheumnob_expanded_geneIDs.txt
 1391 Rheumnob_totalexpanded_geneIDs.txt


#----------------
cd  ~/scratch/Oxyria/CAFE/Polygonaceae_data/contracted

wc -l *geneIDs.txt
   59 Oxydig_contracted_geneIDs.txt
  675 Oxydig_totalcontracted_geneIDs.txt
    9 Rheumnob_contracted_geneIDs.txt
  317 Rheumnob_totalcontracted_geneIDs.txt




#####################################################################################################################
#From here: continue at Notes6_May2024_Go_enrichment.bash













