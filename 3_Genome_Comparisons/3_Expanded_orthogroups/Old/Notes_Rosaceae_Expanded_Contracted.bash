#################################################################################################################
# Rosaceae expanded gene families
# Sept 2024
# Tutorial https://github.com/elsemikk/Willisornis_Genome_Assembly/blob/master/4.2_Gene_Family_Expansions.md
#################################################################################################################

tmux new-session -s CAFE2
tmux attach-session -t CAFE2

cd /lustre04/scratch/celphin/Oxyria/CAFE/Rosaceae_data/analysis

more ../output/error_model/Base_clade_results.txt

#get a list of the significantly expanded/contracted/rapid genes in Oxy/Rhuem/Dryas/Draba
#Taxon_ID                    Increase  Decrease
                       # <12>    284     44
       # Dryas_octopetala<1>     378     795
   # Pyrus_bretschneideri<5>     1315    106
       # Malus_sylvestris<6>     1226    153
                       # <9>     15      44
                       # <10>    225     294
            # Rosa_rugosa<2>     1063    161
                       # <11>    26      176
     # Argentina_anserina<3>     333     460
         # Fragaria_vesca<4>     561     218
                       # <13>    6747    102
         # Prunus_persica<7>     345     679


more ../output/error_model/Base_asr.tre 
##################################################################################################################
# look at only top 3 runs with lowest error and move to new folder

cd /home/celphin/scratch/Oxyria/CAFE/Rosaceae_data/analysis

# get epsilon values for each of the 50 runs
grep "Epsilon:" ../output/error_model_*/Base_results.txt | sed 's/\.\.\/output\/error_model\_//g' |  \
sed 's/\/Base_results.txt\:Epsilon\:/ /g' | \
sort -k 2 > Dryasoct_epsilonval_error_models 

Dryasoct_low_error_models=$(head -n 3 Dryasoct_epsilonval_error_models | awk  -F $' ' '{print $1}')

#get a list of the significantly expanded/contracted/rapid genes in Dryas_octopetala<1>
printf '%s\n' $Dryasoct_low_error_models | while read number ; do grep "y" ../output/error_model_$number/Base_family_results.txt | \
sed 's/^.*OG/OG/g' | sed 's/ = .*$//g' > Dryasoct_sig_changes_error_model_"$number" ; done 

#------------------------------
# Old
#get a list of the significantly expanded/contracted/rapid genes in Dryas_octopetala<1>
printf '%s\n' {1..50} | while read number ; do grep "Dryas_octopetala<1>\*" ../output/error_model_$number/Base_asr.tre | \
sed 's/^.*OG/OG/g' | sed 's/ = .*$//g' > Dryasoct_sig_changes_error_model_"$number" ; done 

#------------------------
#Data exploration:
wc -l Dryasoct_sig_changes_error_model_* # 106-112
# 106

#---------------------
#see what is uniquely found in run 6 
printf '%s\n' $Dryasoct_low_error_models | while read number ; do grep -v -f Dryasoct_sig_changes_error_model_"$number" Dryasoct_sig_changes_error_model_15; done |wc -l
printf '%s\n' $Dryasoct_low_error_models | while read number ; do grep -v -f Dryasoct_sig_changes_error_model_"$number" Dryasoct_sig_changes_error_model_1 ; done

##################################################################################################################
#List of expanded and contracted groups (both significant and not)

#Get list of families that have expanded in Oxyria, whether or not significant
printf '%s\n' {1..50} | while read number ; do cut -f 1,2 ../output/error_model_$number/Base_change.tab | \
awk -F'\t' '$2 != 0' | grep -v "-" > Dryasoct_expandedfams_error_model_"$number" ; done 

wc -l *expandedfams_error_model_1 #remember the header is one line!
# 379 Dryasoct_expandedfams_error_model_1

#--------------------------------------
#Get list of families that have contracted in Oxyria, whether or not significant

printf '%s\n' {1..50} | while read number ; do cut -f 1,2 ../output/error_model_$number/Base_change.tab | \
awk -F'\t' '$2 != 0' | grep "-" > Dryasoct_contractedfams_error_model_"$number" ; done 

wc -l *contractedfams_error_model_1 #remember the header is one line!
# 795 Dryasoct_contractedfams_error_model_1

##################################################################################################################
#Running python script:
module load StdEnv/2020 python/3.9.6 scipy-stack/2021a
python

import os
from os.path import join as jn
import numpy as np

folderpath= os.getcwd()

#############################

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

######################
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

##############################
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


exit()

##################################################################################################################
#Clean up:
# rm -r contractedfams_error_models
# rm -r expandedfams_error_models
# rm -r sig_genes

mkdir sig_genes; mv *sig_changes* sig_genes
mkdir expandedfams_error_models; mv *expandedfams_error_model* expandedfams_error_models
mkdir contractedfams_error_models; mv *contractedfams_error_model* contractedfams_error_models

##################################################################################################################

# Dryas octopetala
cat output_Dryasoct_expandedfams | sort | uniq > OGtotaltemp
cat OGtotaltemp | while read line ; do grep -c "$line" output_Dryasoct_expandedfams >> tempcounts ; done
paste OGtotaltemp tempcounts >> OGcountstemp
cat OGcountstemp | while read orthogroup counts ; do if [ "$counts" -gt 3 ] ; then echo "$orthogroup" >> Dryasoct_expandedfams ; fi ; done
rm -f OGcountstemp OGtotaltemp tempcounts output 


cat output_Dryasoct_contractedfams | sort | uniq > OGtotaltemp
cat OGtotaltemp | while read line ; do grep -c "$line" output_Dryasoct_contractedfams >> tempcounts ; done
paste OGtotaltemp tempcounts >> OGcountstemp
cat OGcountstemp | while read orthogroup counts ; do if [ "$counts" -gt 3 ] ; then echo "$orthogroup" >> Dryasoct_contractedfams ; fi ; done
rm -f OGcountstemp OGtotaltemp tempcounts output 


cat output_Dryasoct | sort | uniq > OGtotaltemp
cat OGtotaltemp | while read line ; do grep -c "$line" output_Dryasoct >> tempcounts ; done
paste OGtotaltemp tempcounts >> OGcountstemp
cat OGcountstemp | while read orthogroup counts ; do if [ "$counts" -gt 3 ] ; then echo "$orthogroup" >>Dryasoct_sig_changes ; fi ; done
rm -f OGcountstemp OGtotaltemp tempcounts output


#######################################################################################################################################
#Make gene list in orthogroups total expanded
mkdir ~/scratch/Oxyria/CAFE/Rosaceae_data/expanded
cd ~/scratch/Oxyria/CAFE/Rosaceae_data/expanded

cp ../analysis/Dryasoct_expandedfams  .
sed -i 's/N0.H//g' Dryasoct_expandedfams 
wc -l Dryasoct_expandedfams
#378 Dryasoct_expandedfams

#---------------------------------------------------
#get a fasta of the expanded families from the OrthoFinder output
for taxon in Dryasoct ; do cat "$taxon"_expandedfams | cut -f 1 | while read Orthogroup ; do cp -v /home/celphin/scratch/Oxyria/GeneSpace/Rosaceae_genomes/orthofinder/Results_Aug16/Orthogroup_Sequences/"$Orthogroup".fa .; done ; done

# see gene names in these orthogroups
grep  ">" /home/celphin/scratch/Oxyria/GeneSpace/Rosaceae_genomes/tmp/Dryas_octopetala.fa | head

cd ~/scratch/Oxyria/CAFE/Rosaceae_data/expanded
#Get a list of the specific species genes in these orthogroups
for taxon in Dryasoct; do cat "$taxon"_expandedfams | cut -f 1 | while read Orthogroup ; do grep "DoctH0" "$Orthogroup".fa >> "$taxon"_totalexpanded_geneIDs.txt ; done ; done
for taxon in Dryasoct ; do sed -i 's/>//g' "$taxon"_totalexpanded_geneIDs.txt ; done

#############################################
#Make gene list in orthogroups total contracted

mkdir ~/scratch/Oxyria/CAFE/Rosaceae_data/contracted
cd ~/scratch/Oxyria/CAFE/Rosaceae_data/contracted

cp ../analysis/Dryasoct_contractedfams .
sed -i 's/N0.H//g' Dryasoct_contractedfams
wc -l Dryasoct_contractedfams
#---------------------------------------------------
#get a fasta of the contracted families from the OrthoFinder output
for taxon in Dryasoct ; do cat "$taxon"_contractedfams | cut -f 1 | while read Orthogroup ; do cp -v /home/celphin/scratch/Oxyria/GeneSpace/Rosaceae_genomes/orthofinder/Results_Aug16/Orthogroup_Sequences/"$Orthogroup".fa .; done ; done

# see gene names in these orthogroups
grep  ">" /home/celphin/scratch/Oxyria/GeneSpace/Rosaceae_genomes/tmp/Dryas_octopetala.fa | head

#---------------------------
cd ~/scratch/Oxyria/CAFE/Rosaceae_data/contracted
#Get a list of the specific species genes in these orthogroups
for taxon in Dryasoct; do cat "$taxon"_contractedfams | cut -f 1 | while read Orthogroup ; do grep "DoctH0" "$Orthogroup".fa >> "$taxon"_totalcontracted_geneIDs.txt ; done ; done

for taxon in Dryasoct ; do sed -i 's/>//g' "$taxon"_totalcontracted_geneIDs.txt ; done


######################################################################################################################
#See which  expanded and contracted orthogroups are also significant
for taxon in Dryasoct ; do grep -f "$taxon"_sig_changes "$taxon"_expandedfams > "$taxon"_expandedfams.sig ; done
for taxon in Dryasoct ; do grep -f "$taxon"_sig_changes "$taxon"_contractedfams > "$taxon"_contractedfams.sig ; done

wc -l *_expandedfams*.sig  
 # 29 Dryasoct_expandedfams.sig

wc -l *_contractedfams*.sig
  # 74 Dryasoct_contractedfams.sig


#######################################################################################################################################
#Make gene list in orthogroups sig rapidly expanded
mkdir ~/scratch/Oxyria/CAFE/Rosaceae_data/expanded
cd ~/scratch/Oxyria/CAFE/Rosaceae_data/expanded

cp ../analysis/Dryasoct_expandedfams.sig .
sed -i 's/N0.H//g' Dryasoct_expandedfams.sig 

#---------------------------------------------------
#get a fasta of the expanded families from the OrthoFinder output
#for taxon in Dryasoct ; do cat "$taxon"_expandedfams.sig | cut -f 1 | while read Orthogroup ; do cp /home/celphin/scratch/Oxyria/GeneSpace/Rosaceae_genomes/orthofinder/Results_Aug16/Orthogroup_Sequences/"$Orthogroup".fa .; done ; done

# see gene names in these orthogroups
grep  ">" /home/celphin/scratch/Oxyria/GeneSpace/Rosaceae_genomes/tmp/Dryas_octopetala.fa | head


cd ~/scratch/Oxyria/CAFE/Rosaceae_data/expanded
#Get a list of the specific species genes in these orthogroups
for taxon in Dryasoct; do cat "$taxon"_expandedfams.sig | cut -f 1 | while read Orthogroup ; do grep "DoctH0" "$Orthogroup".fa >> "$taxon"_expanded_geneIDs.txt ; done ; done
for taxon in Dryasoct ; do sed -i 's/>//g' "$taxon"_expanded_geneIDs.txt ; done

#############################################
#Make gene list in orthogroups sig rapidly contracted

mkdir ~/scratch/Oxyria/CAFE/Rosaceae_data/contracted
cd ~/scratch/Oxyria/CAFE/Rosaceae_data/contracted

cp ../analysis/Dryasoct_contractedfams.sig .
sed -i 's/N0.H//g' Dryasoct_contractedfams.sig 

#---------------------------------------------------
#get a fasta of the contracted families from the OrthoFinder output
#for taxon in Dryasoct ; do cat "$taxon"_contractedfams.sig | cut -f 1 | while read Orthogroup ; do cp /home/celphin/scratch/Oxyria/GeneSpace/Rosaceae_genomes/orthofinder/Results_Aug16/Orthogroup_Sequences/"$Orthogroup".fa .; done ; done

# see gene names in these orthogroups
grep  ">" /home/celphin/scratch/Oxyria/GeneSpace/Rosaceae_genomes/tmp/Dryas_octopetala.fa | head

#---------------------------
cd ~/scratch/Oxyria/CAFE/Rosaceae_data/contracted
#Get a list of the specific species genes in these orthogroups
for taxon in Dryasoct; do cat "$taxon"_contractedfams.sig | cut -f 1 | while read Orthogroup ; do grep "DoctH0" "$Orthogroup".fa >> "$taxon"_contracted_geneIDs.txt ; done ; done

for taxon in Dryasoct ; do sed -i 's/>//g' "$taxon"_contracted_geneIDs.txt ; done

###############################
cd  ~/scratch/Oxyria/CAFE/Rosaceae_data/expanded

wc -l *geneIDs.txt
  # 122 Dryasoct_expanded_geneIDs.txt
  # 745 Dryasoct_totalexpanded_geneIDs.txt


#----------------
cd  ~/scratch/Oxyria/CAFE/Rosaceae_data/contracted

wc -l *geneIDs.txt
  # 256 Dryasoct_contracted_geneIDs.txt
 # 1636 Dryasoct_totalcontracted_geneIDs.txt




#####################################################################################################################
#From here: continue at Notes6_May2024_Go_enrichment.bash


