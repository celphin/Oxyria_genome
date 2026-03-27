#########################
# Running Gene Ontology enrichment on guidance results
# Mar 2026
############################

# Narval1
tmux new-session -s total
tmux attach-session -t total

cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees/Arctic_trees

# Extract full rows with all 4 species
awk '{split($3,species,","); if(length(species)==4) print $0}' \
OG_paralogs_species_counts_sorted_guidance.tsv \
> OG_4species_rows.tsv

# Extract full rows with all 4 species
awk '{split($3,species,","); if(length(species)==4) print $0}' \
OG_single_copy_species_counts_sorted_guidance.tsv \
> single_copy_4species_rows.tsv


#---------------------------
# Copy files to GO enrichment folder
cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/Gene_ontology

cp /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees/Arctic_trees/single_copy_4species_rows.tsv . 
cp /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees/Arctic_trees/OG_4species_rows.tsv .
# single_copy_4species_rows.tsv
# OG_4species_rows.tsv



#------------------------
# Look for overlap with GO enrichment results

# Step 1: Extract OG IDs
cut -f1 OG_4species_rows.tsv > OG_4species_IDs.txt

# Step 2: Get genes with p-value < 0.05 for those OGs
awk -F'\t' 'NR==FNR {og[$1]; next} $1 in og && $3 < 0.05 {print $1, $2, $3}' OFS='\t' OG_4species_IDs.txt baseline_guidance_all_genes.tsv > genes_filtered.txt

# Step 3: Join with InterProScan data and assign species
awk -F'\t' 'NR==FNR {g[$2]=$1"\t"$3; next} 
     $2 in g {
         split(g[$2], ogp, "\t")
         og = ogp[1]
         pval = ogp[2]
         gene = $2
         
         # assign species
         species=""
         if (gene ~ /^DOCTH0_CHR/) species="Dryas"
         else if (gene ~ /^OXYRIA_NCBI_CHR/) species="Oxyria"
         else if (gene ~ /^MAKER_/ || gene ~ /^SNAP_MASKED_/) species="Draba"
         else if (gene ~ /^G[0-9]+_T[0-9]+/) species="Coch"
         else next

         # print OG, p-value, gene, species, description, GO
         print og, pval, gene, species, $10, $11
     }' OFS='\t' genes_filtered.txt InterProscan_guidance_all_genes.tsv \
     > genes_GO_description_4species_filtered.tsv

more genes_GO_description_4species_filtered.tsv

# extract GO terms
cut -f6 genes_GO_description_4species_filtered.tsv > GO_terms_4spp.txt
# split
cat GO_terms_4spp.txt | tr ',' '\n' > GO_terms_4spp_split.txt
# count occurances
sort GO_terms_4spp_split.txt | uniq -c | sort -nr > GO_terms_4spp_counts.txt

 # 176 GO:0005515
    # 165 NA
    # 138 GO:0005524
     # 69 GO:0016020
     # 65 GO:0006468
     # 65 GO:0004672
     # 36 GO:0055085
     # 35 GO:0008270
     # 31 GO:0003824
     # 29 GO:0003677
     # 28 GO:0003676
     # 27 GO:0003723
     # 25 GO:0016887
     # 23 GO:0016757
     # 20 GO:0006888
     # 20 GO:0005975
     # 20 GO:0005634
     # 18 GO:0006508
     # 16 GO:0140658
     # 16 GO:0016491
     # 13 GO:0006886
     # 13 GO:0006486
     # 12 GO:0050660
     # 12 GO:0030127
     # 12 GO:0022857
     # 12 GO:0009058
     # 12 GO:0008017
     # 12 GO:0006355
     # 12 GO:0006281
     # 12 GO:0006260
     # 12 GO:0005737
     # 12 GO:0004553
      # 9 GO:0046872
      # 8 GO:0045454
      # 8 GO:0042555
      # 8 GO:0032508
      # 8 GO:0030170
      # 8 GO:0016787
      # 8 GO:0016740
      # 8 GO:0016567
      # 8 GO:0015031
      # 8 GO:0008236

# High-frequency GO terms:
# GO:0005515 (protein binding) – 176 times
# GO:0005524 (ATP binding) – 138 times
# GO:0016020 (membrane) – 69 times
# GO:0006468 (protein phosphorylation) – 65 times
# GO:0004672 (protein kinase activity) – 65 times
# NA entries: 165 — these genes don’t have a GO term assigned in your InterProScan output.
# Mid-frequency terms reflect diverse functions: nucleic acid binding, transport, enzymatic activity, etc.

# Potentially interesting/adaptation-related:
# protein kinase activity / protein phosphorylation (signaling adaptation)
# transmembrane transport (ion/homeostasis adaptation)
# DNA binding / RNA binding (regulatory adaptation)
# zinc ion binding (if transcription factors are involved)



#########################
# Run GO enrichment for combined all species

# Join files in R 


module load StdEnv/2020 r/4.2.2

export R_LIBS_USER=$HOME/R/x86_64-pc-linux-gnu-library/4.2
mkdir -p $R_LIBS_USER

R

library(dplyr)
library(tidyr)
library(stringr)

#-------------------------
# load GO ont data
# formatted Interproscan to have no duplicates of genes - one row per gene

path="/home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/Gene_ontology/"
Gene_ont_file <- "Total_interproscan_output_edited3.tsv"
gene_ont <- read.delim(paste0(path,"/", Gene_ont_file), header = TRUE, sep = "\t", na.strings = "-", colClasses = c("character", "character", "character", "character"))

nrow(gene_ont)
# [1] 109 116

colnames(gene_ont) <- c("spp", "gene", "INTPRO", "descrip", "GOterm")
length(unique(gene_ont$INTPRO))
# [1] 15 433

unique(gene_ont$spp)

#-------------------------------
# Load all genes list

path="/home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/Gene_ontology/"
all_genes <- read.delim(paste0(path,"/baseline_guidance_all_genes.tsv"), header = TRUE, sep = "\t")

nrow(all_genes)
# 626 371

colnames(all_genes) 

head(all_genes)

all_genes$orthogroup <- sub("_tree\\.txt_unique\\.nxh\\.ABSREL\\.json$", "", all_genes$file)

#-------------------------
# Fix capitalization change 
# Also fix . to _ and - to _

all_genes <- all_genes %>%
  mutate(gene = gene %>%
           toupper() %>%
           str_replace_all("[.-]", "_"))

gene_ont <- gene_ont %>%
  mutate(gene = gene %>%
           toupper() %>%
           str_replace_all("[.-]", "_"))


#----------------------------
# Join sig_genes with GO info

merged_data <- all_genes %>% left_join(gene_ont, by = "gene")

head(merged_data, 10)
colnames(merged_data)

write.table(merged_data, "InterProscan_guidance_all_genes.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# remove descrip column
merged_data_nodescrip <- merged_data[,-6]
head(merged_data_nodescrip, 100)

# remove rows with no GO terms
merged_data_filtered <- merged_data %>%
  filter(!is.na(GOterm))

head(merged_data_filtered, 10)

gene_ont <- merged_data_filtered

#-------------------------------------------
# Need to compare to frequencies in whole genome - ErmineJ
# Need to get list of GO terms that were in original ABSREL test

#subset Interproscan info by spp - to get null model with all genes

unique(gene_ont$spp)
# [2] "Cochlearia_groenlandica_interproscan_output.tsv"
# [3] "Draba_nivalis_interproscan_output.tsv"
# [4] "Dryas_octopetala_interproscan_output.tsv"
# [5] "Oxyria_digyna_H1_interproscan_output.tsv"

###########################
# format for ermineJ
# https://erminej.msl.ubc.ca/help/input-files/
# https://erminej.msl.ubc.ca/help/input-files/gene-annotations/

# Make mapping file and genesets

# Combined?
gene_ont_Oxydig <- gene_ont[which(gene_ont$spp=="Oxyria_digyna_H1_interproscan_output.tsv"),]
write.table(gene_ont_Oxydig, "Oxydig_interproscan_guidance.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

gene_ont_Oxydig1 <- gene_ont_Oxydig %>%
  mutate(gene2 = gene) %>%
  select(gene, gene2, descrip, GOterm)
write.table(gene_ont_Oxydig1, "Oxydig_GO_mappings_guidance.ermineJ.txt", sep = "\t", quote = FALSE, row.names = FALSE)

gene_ont_Oxydig2 <- gene_ont_Oxydig %>%
  select(gene, corrected_p)
write.table(gene_ont_Oxydig2, "Oxydig_genescores_guidance", sep = "\t", quote = FALSE, col.names=TRUE, row.names = FALSE)


##################################
# Remove gene families that have paralogs in Arctic species

# Combined?
single_gene_ont_Oxydig <- gene_ont_Oxydig %>%
  group_by(orthogroup) %>%
  filter(n_distinct(gene) == 1) %>%
  ungroup()

single_gene_ont_Oxydig1 <- single_gene_ont_Oxydig %>%
  mutate(gene2 = gene) %>%
  select(gene, gene2, descrip, GOterm)
write.table(single_gene_ont_Oxydig1, "single_Oxydig_GO_mappings_guidance.ermineJ.txt", sep = "\t", quote = FALSE, row.names = FALSE)

single_gene_ont_Oxydig2 <- single_gene_ont_Oxydig %>%
  select(gene, corrected_p)
write.table(single_gene_ont_Oxydig2, "single_Oxydig_genescores_guidance", sep = "\t", quote = FALSE, col.names=TRUE, row.names = FALSE)


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
# narval1
tmux new-session -s Enrichment
tmux attach-session -t Enrichment

# https://erminej.msl.ubc.ca/help/tutorials/running-an-analysis-ora/

salloc -c1 --time 3:00:00 --mem 120000m --account def-rieseber

cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/Gene_ontology

ERMINEJ_HOME=/home/celphin/ermineJ-3.2
export JAVA_HOME=/cvmfs/soft.computecanada.ca/easybuild/software/2020/Core/java/13.0.2/

module load StdEnv/2020 java/13.0.2

#-------------------
# Scores are p-values

for taxon in Combined ; do $ERMINEJ_HOME/bin/ermineJ.sh \
-a "$taxon"_GO_mappings_guidance.ermineJ.txt \
-s "$taxon"_genescores_guidance \
-c /home/celphin/ermineJ.data/go.obo \
--genesOut --logTrans -n ORA \
-o "$taxon"_ORA_guidance.ermine.results -y 5 --mtc FDR ; done


#------------------------
# Run again removing any genes that are paralogs in Arctic species
for taxon in Combined ; do $ERMINEJ_HOME/bin/ermineJ.sh \
-a single_"$taxon"_GO_mappings_guidance.ermineJ.txt \
-s single_"$taxon"_genescores_guidance \
-c /home/celphin/ermineJ.data/go.obo \
--genesOut --logTrans -n ORA \
-o single_"$taxon"_ORA_guidance.ermine.results -y 5 --mtc FDR ; done


###########################################
# Explore data

# Extract GO ID 
more Combined_ORA_guidance.ermine.results

# Broad themes
# Oxyria glycosylation
# Dryas DNA repair
# Draba defense Response
# Coch Oxidioreductase 

#-----------------------------

awk -F'\t' '$1=="!" {print $2,$3,$6,$7}' OFS='\t' Combined_ORA_guidance.ermine.results \
| sort -t $'\t' -k4,4g | head -n 30

