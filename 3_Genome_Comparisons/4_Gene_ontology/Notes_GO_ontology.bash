###########################################################
# Sept 2024
# following this: https://github.com/elsemikk/Willisornis_Genome_Assembly/blob/master/4.2_Gene_Family_Expansions.md
# GO analysis section
###########################################################

#######################################################################################################################################
#In home, install ermine:
cd ~
wget https://home.pavlab.msl.ubc.ca/ermineJ/distributions/ermineJ-3.2-generic-bundle.zip
unzip ermineJ-3.2-generic-bundle.zip
module load StdEnv/2020  python/3.10.2 scipy-stack/2021a diamond/2.1.7 fastme/2.1.6.2 gcc/9.3.0 mcl/14.137 java/13.0.2 r/4.2.2 glpk/5.0
cd ermineJ-3.2/bin/
chmod +x ermineJ.sh

#######################################################################################################################################
#Get gene ontology file (gene set file)
mkdir ~/ermineJ.data
cd ~/ermineJ.data
wget https://release.geneontology.org/2024-04-24/ontology/go.obo

mkdir ~/scratch/Oxyria/CAFE/enrichment_analysis
cd ~/scratch/Oxyria/CAFE/enrichment_analysis
mkdir Interproscan_files

# copy over Interproscan data
cd /home/celphin/scratch/Oxyria/CAFE/enrichment_analysis/Interproscan_files
cp -v ~/scratch/Oxyria/Interproscan/*_interproscan_output.tsv . 
cp -v /home/celphin/scratch/Oxyria/Genomes_Annotations/*/*.tsv .

#------------------
# Join Interpro files but keep file name - can just copy same file as done below over from genome synteny
# Create or clear the output file
> Total_interproscan_output0.tsv

# Loop through each file and append its contents to the output file
for file in *_interproscan_output.tsv; do
    awk -v filename="$(basename "$file")" '{print filename "\t" $0}' "$file" >> Total_interproscan_output0.tsv
done

# format 
awk -v FS="\t" '{print $1 "\t" $2 "\t" $5 "\t" $13 "\t" $14 "\t" $15}' Total_interproscan_output0.tsv | wc -l 
# 1669585

awk -v FS="\t" '{print $1 "\t" $2 "\t" $5 "\t" $13 "\t" $14 "\t" $15}' Total_interproscan_output0.tsv | sort | uniq
# 1112967

awk -v FS="\t" '{print $1 "\t" $2 "\t" $5 "\t" $13 "\t" $14 "\t" $15}' Total_interproscan_output0.tsv | sort | uniq |wc -l
# 669657

awk -v FS="\t" '{print $1 "\t" $2 "\t" $5 "\t" $13 "\t" $14 "\t" $15}' Total_interproscan_output0.tsv | sort | uniq > Total_interproscan_output_edited.tsv

grep -v $'\t''-'$'\t''-'$'\t' Total_interproscan_output_edited.tsv > Total_interproscan_output_edited1.tsv

sed 's/|/,/g' Total_interproscan_output_edited1.tsv | sort -u > Total_interproscan_output_edited2.tsv

#############################################
# Prepare gene score file for gene families found in all four Arctic spp but not other species
# Oxydig Rheumnob Arabalp Dryasoct Drabaniv Cochgro 

# copy over important files
cd /home/celphin/scratch/Oxyria/CAFE/enrichment_analysis/
cp /home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/orthofinder/Results_Aug18/Phylogenetic_Hierarchical_Orthogroups/N0.tsv .
cp /home/celphin/scratch/Oxyria/synteny_quantity/Total_interproscan_output_edited3.tsv .

# in R
module load StdEnv/2020 r/4.2.2 
R

#install.packages("ape")
# https://rdrr.io/cran/ape/man/chronos.html
# https://phylobotanist.blogspot.com/2016/12/how-to-set-multiple-calibration-points.html

library(ape)
library(data.table)

# Read in the HOG data
hog <- fread('./N0.tsv')
hog[, OG := NULL]
hog[, `Gene Tree Parent Clade` := NULL]
hog <- melt(hog, id.vars='HOG', variable.name='species', value.name='pid')
hog <- hog[pid != '']
hog$n <- sapply(hog$pid, function(x) length(strsplit(x, ', ')[[1]]))

# Define the list of species you're interested in
species_of_interest <- c("Oxyria_digyna_H1",  "Dryas_octopetala", "Draba_nivalis", "Cochlearia_groenlandica")

# Filter hog data for the specified species
filtered_hog <- hog[species %in% species_of_interest]

# Identify HOGs that are present in any species outside the list
all_hogs <- unique(hog$HOG)
hogs_in_other_species <- unique(hog[!species %in% species_of_interest, HOG])

# Filter out HOGs found in species not of interest
hogs_exclusive <- filtered_hog[!HOG %in% hogs_in_other_species]

# Identify HOGs that are present in more than one species of interest
hogs_multiple_species <- hogs_exclusive[, .(species_count = uniqueN(species)), by = HOG][species_count > 1]$HOG

# Final filtering to keep only those HOGs
final_hog <- hogs_exclusive[HOG %in% hogs_multiple_species]

# Create a count table
counts <- dcast(final_hog, HOG ~ species, value.var='n', fill=0)
counts[, Desc := 'n/a']
setcolorder(counts, 'Desc')

# Write out the filtered file
fwrite(counts, 'Arctic_specific_hogs_filtered.tsv', sep='\t')

# Display the filtered counts
print(counts)
nrow(counts)
# 144

       # Desc           HOG Cochlearia_groenlandica Draba_nivalis
     # <char>        <char>                   <int>         <int>
  # 1:    n/a N0.HOG0000050                       0             1
  # 2:    n/a N0.HOG0000057                       0             1
  # 3:    n/a N0.HOG0000139                       0             0
  # 4:    n/a N0.HOG0000209                       0             0
  # 5:    n/a N0.HOG0000253                       0             0
 # ---
# 140:    n/a N0.HOG0033962                       0             0
# 141:    n/a N0.HOG0034042                       0             0
# 142:    n/a N0.HOG0034049                       0             0
# 143:    n/a N0.HOG0034112                       0             0
# 144:    n/a N0.HOG0034227                       0             0
     # Dryas_octopetala Oxyria_digyna_H1
                # <int>            <int>
  # 1:                1                0
  # 2:                0                1
  # 3:               15                1
  # 4:                2                1
  # 5:                2                1
 # ---
# 140:                1                1
# 141:                1                1
# 142:                1                1
# 143:                1                1
# 144:                1                1

##########################
# NEED TO GET GENE IDs for Arctic specific
cd /home/celphin/scratch/Oxyria/CAFE/enrichment_analysis/

# remove the header line
nano Arctic_specific_hogs_filtered.tsv

#Get a list of the specific species genes in these orthogroups (HOGs)
for taxon in Oxydig; do
    cat Arctic_specific_hogs_filtered.tsv | cut -f 2 | while read Orthogroup; do
        grep "$Orthogroup" /home/celphin/scratch/Oxyria/CAFE/Total_genomes/input_data/N0.tsv | \
        awk -F'\t' '{for(i=1; i<=NF; i++) if($i ~ /Oxyria/) print $i}' >> "${taxon}_Arctic_specific_geneIDs0.txt"
    done
done

for taxon in Dryasoct; do
    cat Arctic_specific_hogs_filtered.tsv | cut -f 2 | while read Orthogroup; do
        grep "$Orthogroup" /home/celphin/scratch/Oxyria/CAFE/Total_genomes/input_data/N0.tsv | \
        awk -F'\t' '{for(i=1; i<=NF; i++) if($i ~ /DoctH0/) print $i}' >> "${taxon}_Arctic_specific_geneIDs0.txt"
    done
done

for taxon in Drabaniv; do
    cat Arctic_specific_hogs_filtered.tsv | cut -f 2 | while read Orthogroup; do
        grep "$Orthogroup" /home/celphin/scratch/Oxyria/CAFE/Total_genomes/input_data/N0.tsv | \
        awk -F'\t' '{for(i=1; i<=NF; i++) if($i ~ /\-lg/) print $i}' >> "${taxon}_Arctic_specific_geneIDs0.txt"
    done
done

for taxon in Cochgro; do
    cat Arctic_specific_hogs_filtered.tsv | cut -f 2 | while read Orthogroup; do
        grep "$Orthogroup" /home/celphin/scratch/Oxyria/CAFE/Total_genomes/input_data/N0.tsv | \
        awk -F'\t' '{for(i=1; i<=NF; i++) if($i ~ /g[0-9]+\.t[12]/) print $i}' >> "${taxon}_Arctic_specific_geneIDs0.txt"
    done
done

for taxon in Oxydig Dryasoct Drabaniv Cochgro ; do tr ', ' '\n' < "$taxon"_Arctic_specific_geneIDs0.txt  | awk 'NF' | sort -u > "$taxon"_Arctic_specific_geneIDs.txt ; done

wc -l *_Arctic_specific_geneIDs0.txt
   # 82 Cochgro_Arctic_specific_geneIDs0.txt
   # 75 Drabaniv_Arctic_specific_geneIDs0.txt
   # 64 Dryasoct_Arctic_specific_geneIDs0.txt
   # 67 Oxydig_Arctic_specific_geneIDs0.txt

wc -l *_Arctic_specific_geneIDs.txt
 # 180 Cochgro_Arctic_specific_geneIDs.txt
  # 245 Drabaniv_Arctic_specific_geneIDs.txt
  # 168 Dryasoct_Arctic_specific_geneIDs.txt
  # 211 Oxydig_Arctic_specific_geneIDs.txt


#######################################################################################################################################
# edit combined Interproscan file and then split by species 
wc -l Total_interproscan_output_edited3.tsv
# 109117 Total_interproscan_output_edited3.tsv

cd /home/celphin/scratch/Oxyria/CAFE/enrichment_analysis/

module load  StdEnv/2020 r/4.2.2
export R_LIBS_USER=/home/celphin/R/x86_64-pc-linux-gnu-library/4.1.0/

R
library(dplyr)
library(tidyr)

#-------------------------
# load GO ont data
# formatted Interproscan to have no duplicates of genes - one row per gene

path="/home/celphin/scratch/Oxyria/CAFE/enrichment_analysis/"
Gene_ont_file <- "Total_interproscan_output_edited3.tsv"
gene_ont <- read.delim(paste0(path,"/", Gene_ont_file), header = TRUE, sep = "\t", na.strings = "-", colClasses = c("character", "character", "character", "character"))

nrow(gene_ont)
# [1] 109 116

colnames(gene_ont) <- c("spp", "gene", "INTPRO", "descrip", "GOterm")
length(unique(gene_ont$INTPRO))
# [1] 15 434

#-----------------------------
# subset by spp
unique(gene_ont$spp)
# [1] "Arabis_alpina_interproscan_output.tsv"
# [2] "Cochlearia_groenlandica_interproscan_output.tsv"
# [3] "Draba_nivalis_interproscan_output.tsv"
# [4] "Dryas_octopetala_interproscan_output.tsv"
# [5] "Oxyria_digyna_H1_interproscan_output.tsv"
# [6] "Rheum_nobile_H0_interproscan_output.tsv"

# Oxyria
gene_ont_Oxydig <- gene_ont[which(gene_ont$spp=="Oxyria_digyna_H1_interproscan_output.tsv"),]
# format for ermineJ
# https://erminej.msl.ubc.ca/help/input-files/
# https://erminej.msl.ubc.ca/help/input-files/gene-annotations/
gene_ont_Oxydig1 <- gene_ont_Oxydig %>%
  mutate(gene2 = gene) %>%
  select(gene, gene2, descrip, GOterm)
# write output for each spp
write.table(gene_ont_Oxydig1, "Oxydig_GO_mappings.ermineJ.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gene_ont_Oxydig, "Oxydig_interproscan_edited.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

#------------------------------------------
# Rheumnob
gene_ont_Rheumnob <- gene_ont[which(gene_ont$spp=="Rheum_nobile_H0_interproscan_output.tsv"),]
# format for ermineJ
# https://erminej.msl.ubc.ca/help/input-files/
# https://erminej.msl.ubc.ca/help/input-files/gene-annotations/
gene_ont_Rheumnob1 <- gene_ont_Rheumnob %>%
  mutate(gene2 = gene) %>%
  select(gene, gene2, descrip, GOterm)
# write output for each spp
write.table(gene_ont_Rheumnob1, "Rheumnob_GO_mappings.ermineJ.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gene_ont_Rheumnob, "Rheumnob_interproscan_edited.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

#------------------------------------------
# Drabaniv
gene_ont_Drabaniv <- gene_ont[which(gene_ont$spp=="Draba_nivalis_interproscan_output.tsv"),]
# format for ermineJ
# https://erminej.msl.ubc.ca/help/input-files/
# https://erminej.msl.ubc.ca/help/input-files/gene-annotations/
gene_ont_Drabaniv1 <- gene_ont_Drabaniv %>%
  mutate(gene2 = gene) %>%
  select(gene, gene2, descrip, GOterm)
# write output for each spp
write.table(gene_ont_Drabaniv1, "Drabaniv_GO_mappings.ermineJ.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gene_ont_Drabaniv, "Drabaniv_interproscan_edited.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

#------------------------------------------
# Dryasoct
gene_ont_Dryasoct <- gene_ont[which(gene_ont$spp=="Dryas_octopetala_interproscan_output.tsv"),]
# format for ermineJ
# https://erminej.msl.ubc.ca/help/input-files/
# https://erminej.msl.ubc.ca/help/input-files/gene-annotations/
gene_ont_Dryasoct1 <- gene_ont_Dryasoct %>%
  mutate(gene2 = gene) %>%
  select(gene, gene2, descrip, GOterm)
# write output for each spp
write.table(gene_ont_Dryasoct1, "Dryasoct_GO_mappings.ermineJ.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gene_ont_Dryasoct, "Dryasoct_interproscan_edited.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

#------------------------------------------
# Arabalp
gene_ont_Arabalp <- gene_ont[which(gene_ont$spp=="Arabis_alpina_interproscan_output.tsv"),]
# format for ermineJ
# https://erminej.msl.ubc.ca/help/input-files/
# https://erminej.msl.ubc.ca/help/input-files/gene-annotations/
gene_ont_Arabalp1 <- gene_ont_Arabalp %>%
  mutate(gene2 = gene) %>%
  select(gene, gene2, descrip, GOterm)
# write output for each spp
write.table(gene_ont_Arabalp1, "Arabalp_GO_mappings.ermineJ.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gene_ont_Arabalp, "Arabalp_interproscan_edited.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

#------------------------------------------
# Cochgro
gene_ont_Cochgro <- gene_ont[which(gene_ont$spp=="Cochlearia_groenlandica_interproscan_output.tsv"),]
# format for ermineJ
# https://erminej.msl.ubc.ca/help/input-files/
# https://erminej.msl.ubc.ca/help/input-files/gene-annotations/
gene_ont_Cochgro1 <- gene_ont_Cochgro %>%
  mutate(gene2 = gene) %>%
  select(gene, gene2, descrip, GOterm)
# write output for each spp
write.table(gene_ont_Cochgro1, "Cochgro_GO_mappings.ermineJ.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gene_ont_Cochgro, "Cochgro_interproscan_edited.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

q()
n

#######################################################################################################################################
# Prepare gene score files
# https://erminej.msl.ubc.ca/help/input-files/gene-scores/

cd /home/celphin/scratch/Oxyria/CAFE/enrichment_analysis/

# get baseline file
# all orthogroups
cp /home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/orthofinder/Results_Aug18/Orthogroups/Orthogroups.tsv .

# subset for all hogs found in all spp - total_genomes
cp /home/celphin/scratch/Oxyria/CAFE/Total_genomes/input_data/hog_gene_counts.tsv  ./Total_hogs_allspp_filtered.tsv
# NOTE: will only run this total set for the total families counts

#Get a list of the specific species genes in these orthogroups

for taxon in Dryasoct; do
    cat Total_hogs_allspp_filtered.tsv | cut -f 2 | while read Orthogroup; do
        grep "$Orthogroup" /home/celphin/scratch/Oxyria/CAFE/Total_genomes/input_data/N0.tsv | \
        awk -F'\t' '{for(i=1; i<=NF; i++) if($i ~ /DoctH0/) print $i}' >> "${taxon}_Total_geneIDs0.txt"
    done
done

for taxon in Oxydig; do
    cat Total_hogs_allspp_filtered.tsv | cut -f 2 | while read Orthogroup; do
        grep "$Orthogroup" /home/celphin/scratch/Oxyria/CAFE/Total_genomes/input_data/N0.tsv | \
        awk -F'\t' '{for(i=1; i<=NF; i++) if($i ~ /Oxyria/) print $i}' >> "${taxon}_Total_geneIDs0.txt"
    done
done

for taxon in Drabaniv; do
    cat Total_hogs_allspp_filtered.tsv | cut -f 2 | while read Orthogroup; do
        grep "$Orthogroup" /home/celphin/scratch/Oxyria/CAFE/Total_genomes/input_data/N0.tsv | \
        awk -F'\t' '{for(i=1; i<=NF; i++) if($i ~ /\-lg/) print $i}' >> "${taxon}_Total_geneIDs0.txt"
    done
done

for taxon in Cochgro; do
    cat Total_hogs_allspp_filtered.tsv | cut -f 2 | while read Orthogroup; do
        grep "$Orthogroup" /home/celphin/scratch/Oxyria/CAFE/Total_genomes/input_data/N0.tsv | \
        awk -F'\t' '{for(i=1; i<=NF; i++) if($i ~ /g[0-9]+\.t[12]/) print $i}' >> "${taxon}_Total_geneIDs0.txt"
    done
done

for taxon in Rheumnob; do
    cat Total_hogs_allspp_filtered.tsv | cut -f 2 | while read Orthogroup; do
        grep "$Orthogroup" /home/celphin/scratch/Oxyria/CAFE/Total_genomes/input_data/N0.tsv | \
        awk -F'\t' '{for(i=1; i<=NF; i++) if($i ~ /Rno/) print $i}' >> "${taxon}_Total_geneIDs0.txt"
    done
done

for taxon in Arabalp; do
    cat Total_hogs_allspp_filtered.tsv | cut -f 2 | while read Orthogroup; do
        grep "$Orthogroup" /home/celphin/scratch/Oxyria/CAFE/Total_genomes/input_data/N0.tsv | \
        awk -F'\t' '{for(i=1; i<=NF; i++) if($i ~ /AALP/) print $i}' >> "${taxon}_Total_geneIDs0.txt"
    done
done

for taxon in Oxydig Dryasoct Drabaniv Cochgro Arabalp Rheumnob ; do tr ', ' '\n' < "${taxon}_Total_geneIDs0.txt"  | awk 'NF' > "${taxon}_Total_geneIDs.txt" ; done

# These files are actually more genes then annotated list
###################################

# Collected rapidly expanded and contracted genes sets from various folders
mkdir ~/scratch/Oxyria/CAFE/contracted/
mkdir ~/scratch/Oxyria/CAFE/expanded/

# cp /lustre04/scratch/celphin/Oxyria/CAFE/Polygonaceae_data/expanded/*_expanded_geneIDs.txt ~/scratch/Oxyria/CAFE/expanded
# cp /lustre04/scratch/celphin/Oxyria/CAFE/Polygonaceae_data/contracted/*_contracted_geneIDs.txt ~/scratch/Oxyria/CAFE/contracted

# cp /lustre04/scratch/celphin/Oxyria/CAFE/Brassicaceae_data/expanded/*_expanded_geneIDs.txt ~/scratch/Oxyria/CAFE/expanded
# cp /lustre04/scratch/celphin/Oxyria/CAFE/Brassicaceae_data/contracted/*_contracted_geneIDs.txt ~/scratch/Oxyria/CAFE/contracted

# cp /lustre04/scratch/celphin/Oxyria/CAFE/Rosaceae_data/expanded/*_expanded_geneIDs.txt ~/scratch/Oxyria/CAFE/expanded
# cp /lustre04/scratch/celphin/Oxyria/CAFE/Rosaceae_data/contracted/*_contracted_geneIDs.txt ~/scratch/Oxyria/CAFE/contracted

cp /lustre04/scratch/celphin/Oxyria/CAFE/Total_genomes/expanded/*_total_expanded_geneIDs.txt ~/scratch/Oxyria/CAFE/expanded
cp /lustre04/scratch/celphin/Oxyria/CAFE/Total_genomes/contracted/*_total_contracted_geneIDs.txt ~/scratch/Oxyria/CAFE/contracted

#----------------------------
# #Prepare gene score file, expanded within families

# for taxon in Oxydig Rheumnob Arabalp Dryasoct Drabaniv Cochgro  ; \
# do awk 'BEGIN{FS="\t"}{print $2,"0"}' ~/scratch/Oxyria/CAFE/enrichment_analysis/"$taxon"_interproscan_edited.tsv | sort -u > "$taxon"_expanded_genesets ; done
# wc -l Oxydig_expanded_genesets
# # 20406 Oxydig_expanded_genesets

# for taxon in Oxydig Rheumnob Arabalp Dryasoct Drabaniv Cochgro ; \
# do cat ~/scratch/Oxyria/CAFE/expanded/"$taxon"_expanded_geneIDs.txt | \
# sed 's/ .*$//g' | while read gene ; do sed -i "s/$gene 0/$gene 1/g" "$taxon"_expanded_genesets ; done ; done

# for taxon in Oxydig Rheumnob Arabalp Dryasoct Drabaniv Cochgro ; \
# do sed -i 's/ /\t/g' "$taxon"_expanded_genesets; done

# #--------------------------------
# #Prepare gene score file, contracted within families
# for taxon in Oxydig Rheumnob Arabalp Dryasoct Drabaniv Cochgro  ; \
# do awk 'BEGIN{FS="\t"}{print $2,"0"}' ~/scratch/Oxyria/CAFE/enrichment_analysis/"$taxon"_interproscan_edited.tsv | sort -u > "$taxon"_contracted_genesets ; done

# for taxon in Oxydig Rheumnob Arabalp Dryasoct Drabaniv Cochgro ; \
# do cat ~/scratch/Oxyria/CAFE/contracted/"$taxon"_contracted_geneIDs.txt | \
# sed 's/ .*$//g' | while read gene ; do sed -i "s/$gene 0/$gene 1/g" "$taxon"_contracted_genesets ; done ; done

# for taxon in Oxydig Rheumnob Arabalp Dryasoct Drabaniv Cochgro ; \
# do sed -i 's/ /\t/g' "$taxon"_contracted_genesets; done

#------------------------------
# Prepare gene score file across families (total) - expanded

for taxon in Oxydig Rheumnob Arabalp Dryasoct Drabaniv Cochgro ; \
do awk 'BEGIN{FS="\t"}{print $2,"0"}' ~/scratch/Oxyria/CAFE/enrichment_analysis/"$taxon"_interproscan_edited.tsv | sort -u > "$taxon"_total_expanded_genesets ; done
wc -l *_total_expanded_genesets
  # 12568 Arabalp_total_expanded_genesets
  # 17887 Cochgro_total_expanded_genesets
  # 19417 Drabaniv_total_expanded_genesets
  # 19364 Dryasoct_total_expanded_genesets
  # 20406 Oxydig_total_expanded_genesets
  # 19480 Rheumnob_total_expanded_genesets

wc -l *_total_expanded_genesets # for genes in all orthogroups found in all spp
  # 20389 Arabalp_total_expanded_genesets
  # 29541 Cochgro_total_expanded_genesets
  # 30455 Drabaniv_total_expanded_genesets
  # 37338 Dryasoct_total_expanded_genesets
  # 31594 Oxydig_total_expanded_genesets
  # 32596 Rheumnob_total_expanded_genesets

for taxon in Oxydig Rheumnob Arabalp Dryasoct Drabaniv Cochgro ; \
do cat ~/scratch/Oxyria/CAFE/expanded/"$taxon"_total_expanded_geneIDs.txt | \
sed 's/ .*$//g' | while read gene ; do sed -i "s/$gene 0/$gene 1/g" "$taxon"_total_expanded_genesets ; done ; done

for taxon in Oxydig Rheumnob Arabalp Dryasoct Drabaniv Cochgro ; \
do sed -i 's/ /\t/g' "$taxon"_total_expanded_genesets; done

#--------------------------------
# Prepare gene score file across families (total) - contracted
for taxon in Oxydig Rheumnob Arabalp Dryasoct Drabaniv Cochgro ; \
do awk 'BEGIN{FS="\t"}{print $2,"0"}' ~/scratch/Oxyria/CAFE/enrichment_analysis/"$taxon"_interproscan_edited.tsv | sort -u > "$taxon"_total_contracted_genesets ; done

for taxon in Oxydig Rheumnob Arabalp Dryasoct Drabaniv Cochgro ; \
do cat ~/scratch/Oxyria/CAFE/contracted/"$taxon"_total_contracted_geneIDs.txt | \
sed 's/ .*$//g' | while read gene ; do sed -i "s/$gene 0/$gene 1/g" "$taxon"_total_contracted_genesets ; done ; done

for taxon in Oxydig Rheumnob Arabalp Dryasoct Drabaniv Cochgro ; \
do sed -i 's/ /\t/g' "$taxon"_total_contracted_genesets; done

#---------------------------------
# Prepare gene score file for Arctic specific geneIDs
for taxon in Oxydig Dryasoct Drabaniv Cochgro ; \
do awk 'BEGIN{FS="\t"}{print $2,"0"}' ~/scratch/Oxyria/CAFE/enrichment_analysis/"$taxon"_interproscan_edited.tsv | sort -u > "$taxon"_Arctic_specific_genesets ; done
wc -l *_Arctic_specific_genesets

for taxon in Oxydig Dryasoct Drabaniv Cochgro ; \
do cat ~/scratch/Oxyria/CAFE/enrichment_analysis/"$taxon"_Arctic_specific_geneIDs.txt | \
sed 's/ .*$//g' | while read gene ; do sed -i "s/$gene 0/$gene 1/g" "$taxon"_Arctic_specific_genesets ; done ; done

for taxon in Oxydig Dryasoct Drabaniv Cochgro ; \
do sed -i 's/ /\t/g' "$taxon"_Arctic_specific_genesets; done


########################
# Run again for not rapidly expanded families but all expanded and all contracted gene families

# Collected expanded and contracted genes sets from various folders

cd ~/scratch/Oxyria/CAFE/enrichment_analysis

# cp /lustre04/scratch/celphin/Oxyria/CAFE/Polygonaceae_data/expanded/*_totalexpanded_geneIDs.txt ~/scratch/Oxyria/CAFE/expanded
# cp /lustre04/scratch/celphin/Oxyria/CAFE/Polygonaceae_data/contracted/*_totalcontracted_geneIDs.txt ~/scratch/Oxyria/CAFE/contracted

# cp /lustre04/scratch/celphin/Oxyria/CAFE/Brassicaceae_data/expanded/*_totalexpanded_geneIDs.txt ~/scratch/Oxyria/CAFE/expanded
# cp /lustre04/scratch/celphin/Oxyria/CAFE/Brassicaceae_data/contracted/*_totalcontracted_geneIDs.txt ~/scratch/Oxyria/CAFE/contracted

# cp /lustre04/scratch/celphin/Oxyria/CAFE/Rosaceae_data/expanded/*_totalexpanded_geneIDs.txt ~/scratch/Oxyria/CAFE/expanded
# cp /lustre04/scratch/celphin/Oxyria/CAFE/Rosaceae_data/contracted/*_totalcontracted_geneIDs.txt ~/scratch/Oxyria/CAFE/contracted

cp /lustre04/scratch/celphin/Oxyria/CAFE/Total_genomes/expanded/*_total_totalexpanded_geneIDs.txt ~/scratch/Oxyria/CAFE/expanded
cp /lustre04/scratch/celphin/Oxyria/CAFE/Total_genomes/contracted/*_total_totalcontracted_geneIDs.txt ~/scratch/Oxyria/CAFE/contracted

# #----------------------------
# #Prepare gene score file, expanded within families
# for taxon in Oxydig Rheumnob Arabalp Dryasoct Drabaniv Cochgro ; \
# do awk 'BEGIN{FS="\t"}{print $2,"0"}' ~/scratch/Oxyria/CAFE/enrichment_analysis/"$taxon"_interproscan_edited.tsv | sort -u > "$taxon"_totalexpanded_genesets ; done

# for taxon in Oxydig Rheumnob Arabalp Dryasoct Drabaniv Cochgro ; \
# do cat ~/scratch/Oxyria/CAFE/expanded/"$taxon"_totalexpanded_geneIDs.txt | \
# sed 's/ .*$//g' | while read gene ; do sed -i "s/$gene 0/$gene 1/g" "$taxon"_totalexpanded_genesets ; done ; done

# for taxon in Oxydig Rheumnob Arabalp Dryasoct Drabaniv Cochgro ; \
# do sed -i 's/ /\t/g' "$taxon"_totalexpanded_genesets; done

# #--------------------------------
# #Prepare gene score file, contracted within families
# for taxon in Oxydig Rheumnob Arabalp Dryasoct Drabaniv Cochgro  ; \
# do awk 'BEGIN{FS="\t"}{print $2,"0"}' ~/scratch/Oxyria/CAFE/enrichment_analysis/"$taxon"_interproscan_edited.tsv | sort -u > "$taxon"_totalcontracted_genesets ; done

# for taxon in Oxydig Rheumnob Arabalp Dryasoct Drabaniv Cochgro ; \
# do cat ~/scratch/Oxyria/CAFE/contracted/"$taxon"_totalcontracted_geneIDs.txt | \
# sed 's/ .*$//g' | while read gene ; do sed -i "s/$gene 0/$gene 1/g" "$taxon"_totalcontracted_genesets ; done ; done

# for taxon in Oxydig Rheumnob Arabalp Dryasoct Drabaniv Cochgro ; \
# do sed -i 's/ /\t/g' "$taxon"_totalcontracted_genesets; done

#----------------------------------
# Prepare gene score file across families (total) - expanded

for taxon in Oxydig Rheumnob Arabalp Dryasoct Drabaniv Cochgro ; \
do awk 'BEGIN{FS="\t"}{print $2,"0"}' ~/scratch/Oxyria/CAFE/enrichment_analysis/"$taxon"_interproscan_edited.tsv | sort -u > "$taxon"_total_totalexpanded_genesets ; done

for taxon in Oxydig Rheumnob Arabalp Dryasoct Drabaniv Cochgro ; \
do cat ~/scratch/Oxyria/CAFE/expanded/"$taxon"_total_totalexpanded_geneIDs.txt | \
sed 's/ .*$//g' | while read gene ; do sed -i "s/$gene 0/$gene 1/g" "$taxon"_total_totalexpanded_genesets ; done ; done

for taxon in Oxydig Rheumnob Arabalp Dryasoct Drabaniv Cochgro ; \
do sed -i 's/ /\t/g' "$taxon"_total_totalexpanded_genesets; done

#--------------------------------
# Prepare gene score file across families (total) - contracted

for taxon in Oxydig Rheumnob Arabalp Dryasoct Drabaniv Cochgro ; \
do awk 'BEGIN{FS="\t"}{print $2,"0"}' ~/scratch/Oxyria/CAFE/enrichment_analysis/"$taxon"_interproscan_edited.tsv | sort -u > "$taxon"_total_totalcontracted_genesets ; done

for taxon in Oxydig Rheumnob Arabalp Dryasoct Drabaniv Cochgro ; \
do cat ~/scratch/Oxyria/CAFE/contracted/"$taxon"_total_totalcontracted_geneIDs.txt | \
sed 's/ .*$//g' | while read gene ; do sed -i "s/$gene 0/$gene 1/g" "$taxon"_total_totalcontracted_genesets ; done ; done

for taxon in Oxydig Rheumnob Arabalp Dryasoct Drabaniv Cochgro ; \
do sed -i 's/ /\t/g' "$taxon"_total_totalcontracted_genesets; done

####################################################
# get counts
cd /home/celphin/scratch/Oxyria/CAFE/enrichment_analysis/
wc -l *_genesets
  12568 Arabalp_Arctic_specific_genesets
  17887 Cochgro_Arctic_specific_genesets
  19417 Drabaniv_Arctic_specific_genesets
  19364 Dryasoct_Arctic_specific_genesets
  20406 Oxydig_Arctic_specific_genesets
  19480 Rheumnob_Arctic_specific_genesets


wc -l *_geneIDs.txt
  180 Cochgro_Arctic_specific_geneIDs.txt
  245 Drabaniv_Arctic_specific_geneIDs.txt
  168 Dryasoct_Arctic_specific_geneIDs.txt
  211 Oxydig_Arctic_specific_geneIDs.txt


cd ~/scratch/Oxyria/CAFE/expanded
wc -l *_geneIDs.txt
    94 Arabalp_total_expanded_geneIDs.txt
   738 Arabalp_total_totalexpanded_geneIDs.txt
   152 Cochgro_total_expanded_geneIDs.txt
  3357 Cochgro_total_totalexpanded_geneIDs.txt
   176 Drabaniv_total_expanded_geneIDs.txt
  1884 Drabaniv_total_totalexpanded_geneIDs.txt
   144 Dryasoct_total_expanded_geneIDs.txt
  1056 Dryasoct_total_totalexpanded_geneIDs.txt
   191 Oxydig_total_expanded_geneIDs.txt
  1452 Oxydig_total_totalexpanded_geneIDs.txt
   228 Rheumnob_total_expanded_geneIDs.txt
  1781 Rheumnob_total_totalexpanded_geneIDs.txt


cd ~/scratch/Oxyria/CAFE/contracted
wc -l *_geneIDs.txt
   64 Arabalp_total_contracted_geneIDs.txt
  931 Arabalp_total_totalcontracted_geneIDs.txt
   51 Cochgro_total_contracted_geneIDs.txt
  844 Cochgro_total_totalcontracted_geneIDs.txt
   42 Drabaniv_total_contracted_geneIDs.txt
  378 Drabaniv_total_totalcontracted_geneIDs.txt
   29 Dryasoct_total_contracted_geneIDs.txt
  479 Dryasoct_total_totalcontracted_geneIDs.txt
   10 Oxydig_total_contracted_geneIDs.txt
  585 Oxydig_total_totalcontracted_geneIDs.txt
    5 Rheumnob_total_contracted_geneIDs.txt
  268 Rheumnob_total_totalcontracted_geneIDs.txt


##############################################
# ermineJ

# https://erminej.msl.ubc.ca/help/tutorials/running-an-analysis-ora/

# As of ErmineJ 3, when using the ‘ORA’ method you have the option to use a simple “hit list” of genes,
# rather than preparing a score file yourself (a “quick list”). Caution: If you use this feature, 
# the “non-hits” will be all the rest of the genes listed in your annotation file. That might not 
# be appropriate if the annotation file includes genes that were not assayed in your experiment. 
# This is most likely to be a problem if your annotation file is a list of all the genes in the genome

# Note I should switch total to just be the orthogroups shared by all 


#######################################################################################################################################
tmux new-session -s Enrichment
tmux attach-session -t Enrichment

salloc -c1 --time 3:00:00 --mem 120000m --account def-rieseber

cd ~/scratch/Oxyria/CAFE/enrichment_analysis

#Expanded : Oxyria
ERMINEJ_HOME=/home/celphin/ermineJ-3.2
export JAVA_HOME=/cvmfs/soft.computecanada.ca/easybuild/software/2020/Core/java/13.0.2/

module load java/13.0.2
#-------------------
# # families separate rapidly expanded/contracted
# for taxon in Oxydig Cochgro Dryasoct Drabaniv ; do $ERMINEJ_HOME/bin/ermineJ.sh \
# -a "$taxon"_GO_mappings.ermineJ.txt \
# -s "$taxon"_expanded_genesets \
# -c /home/celphin/ermineJ.data/go.obo \
# --genesOut -aspects BCM  \
# -o "$taxon"_rapidly_expanded_genesets.ermine.results -y 5 -b ; done

# for taxon in Oxydig Cochgro Dryasoct Drabaniv ; do $ERMINEJ_HOME/bin/ermineJ.sh \
# -a "$taxon"_GO_mappings.ermineJ.txt \
# -s "$taxon"_contracted_genesets \
# -c /home/celphin/ermineJ.data/go.obo \
# --genesOut -aspects BCM  \
# -o "$taxon"_rapidly_contracted_genesets.ermine.results -y 5 -b ; done

#-------------------
# total families combined - rapidly expanded/contracted

for taxon in Oxydig Cochgro Dryasoct Drabaniv ; do $ERMINEJ_HOME/bin/ermineJ.sh \
-a "$taxon"_GO_mappings.ermineJ.txt \
-s "$taxon"_total_expanded_genesets \
-c /home/celphin/ermineJ.data/go.obo \
--genesOut -aspects BCM \
-o "$taxon"_totalfam_rapidly_expanded_genesets.ermine.results -y 5 -b ; done

for taxon in Oxydig Cochgro Dryasoct Drabaniv ; do $ERMINEJ_HOME/bin/ermineJ.sh \
-a "$taxon"_GO_mappings.ermineJ.txt \
-s "$taxon"_total_contracted_genesets \
-c /home/celphin/ermineJ.data/go.obo \
--genesOut -aspects BCM  \
-o "$taxon"_totalfam_rapidly_contracted_genesets.ermine.results -y 5 -b ; done

#-------------------
# for families separately for all expanded and contracted

# for taxon in Oxydig Cochgro Dryasoct Drabaniv ; do $ERMINEJ_HOME/bin/ermineJ.sh \
# -a "$taxon"_GO_mappings.ermineJ.txt \
# -s "$taxon"_totalexpanded_genesets \
# -c /home/celphin/ermineJ.data/go.obo \
# --genesOut -aspects BCM  \
# -o "$taxon"_all_expanded_genesets.ermine.results -y 5 -b ; done

# for taxon in Oxydig Cochgro Dryasoct Drabaniv ; do $ERMINEJ_HOME/bin/ermineJ.sh \
# -a "$taxon"_GO_mappings.ermineJ.txt \
# -s "$taxon"_totalcontracted_genesets \
# -c /home/celphin/ermineJ.data/go.obo \
# --genesOut -aspects BCM  \
# -o "$taxon"_all_contracted_genesets.ermine.results -y 5 -b ; done

#-------------------
# for families combined for all expanded and contracted

for taxon in Oxydig Cochgro Dryasoct Drabaniv ; do $ERMINEJ_HOME/bin/ermineJ.sh \
-a "$taxon"_GO_mappings.ermineJ.txt \
-s "$taxon"_total_totalexpanded_genesets \
-c /home/celphin/ermineJ.data/go.obo \
--genesOut -aspects BCM  \
-o "$taxon"_totalfam_all_expanded_genesets.ermine.results -y 5 -b ; done

for taxon in Oxydig Cochgro Dryasoct Drabaniv ; do $ERMINEJ_HOME/bin/ermineJ.sh \
-a "$taxon"_GO_mappings.ermineJ.txt \
-s "$taxon"_total_totalcontracted_genesets \
-c /home/celphin/ermineJ.data/go.obo \
--genesOut -aspects BCM  \
-o "$taxon"_totalfam_all_contracted_genesets.ermine.results -y 5 -b ; done

#-------------------------
# Arctic specific

for taxon in Oxydig Cochgro Dryasoct Drabaniv ; do $ERMINEJ_HOME/bin/ermineJ.sh \
-a "$taxon"_GO_mappings.ermineJ.txt \
-s "$taxon"_Arctic_specific_genesets \
-c /home/celphin/ermineJ.data/go.obo \
--genesOut -aspects BCM  \
-o "$taxon"_Arctic_specific_genesets.ermine.results -y 5 -b ; done


#####################################
# join files for all spp

cat *_GO_mappings.ermineJ.txt > all_spp_GO_mappings.ermineJ.txt
cat *_total_expanded_genesets  > all_spp_total_expanded_genesets
cat *_total_contracted_genesets  > all_spp_total_contracted_genesets
cat *_Arctic_specific_genesets  > all_spp_Arctic_specific_genesets

cat *_total_totalexpanded_genesets > all_spp_total_totalexpanded_genesets
cat *_total_totalcontracted_genesets   > all_spp_total_totalcontracted_genesets


#-------------------
# total families combined - rapidly expanded/contracted

for taxon in all_spp ; do $ERMINEJ_HOME/bin/ermineJ.sh \
-a "$taxon"_GO_mappings.ermineJ.txt \
-s "$taxon"_total_expanded_genesets \
-c /home/celphin/ermineJ.data/go.obo \
--genesOut -aspects BCM \
-o "$taxon"_totalfam_rapidly_expanded_genesets.ermine.results -y 5 -b ; done

for taxon in all_spp ; do $ERMINEJ_HOME/bin/ermineJ.sh \
-a "$taxon"_GO_mappings.ermineJ.txt \
-s "$taxon"_total_contracted_genesets \
-c /home/celphin/ermineJ.data/go.obo \
--genesOut -aspects BCM  \
-o "$taxon"_totalfam_rapidly_contracted_genesets.ermine.results -y 5 -b ; done

#-----------------
# for families combined for all expanded and contracted

for taxon in all_spp ; do $ERMINEJ_HOME/bin/ermineJ.sh \
-a "$taxon"_GO_mappings.ermineJ.txt \
-s "$taxon"_total_totalexpanded_genesets \
-c /home/celphin/ermineJ.data/go.obo \
--genesOut -aspects BCM  \
-o "$taxon"_totalfam_all_expanded_genesets.ermine.results -y 5 -b ; done

for taxon in all_spp ; do $ERMINEJ_HOME/bin/ermineJ.sh \
-a "$taxon"_GO_mappings.ermineJ.txt \
-s "$taxon"_total_totalcontracted_genesets \
-c /home/celphin/ermineJ.data/go.obo \
--genesOut -aspects BCM  \
-o "$taxon"_totalfam_all_contracted_genesets.ermine.results -y 5 -b ; done

#-------------------------
# Arctic specific

for taxon in all_spp ; do $ERMINEJ_HOME/bin/ermineJ.sh \
-a "$taxon"_GO_mappings.ermineJ.txt \
-s "$taxon"_Arctic_specific_genesets \
-c /home/celphin/ermineJ.data/go.obo \
--genesOut -aspects BCM  \
-o "$taxon"_Arctic_specific_genesets.ermine.results -y 5 -b ; done
#446

###################################
# Spp combined - gene IDs for shared gene families

cd  ~/scratch/Oxyria/CAFE/Total_genomes/expanded

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

cat Oxydig_expandedfams Drabaniv_expandedfams Dryasoct_expandedfams Cochgro_expandedfams | sort | uniq -c |sort
     # 3 OG0000494

