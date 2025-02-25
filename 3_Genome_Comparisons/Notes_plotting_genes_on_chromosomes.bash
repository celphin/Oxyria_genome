############################
# Plot gene or repeat densities
# https://genviz.org/module-02-r/0002/03/03/ggplot2_exercises/
# Jan 2025
##############################

# Narval2
tmux new-session -s gene_plots
tmux attach-session -t gene_plots

mkdir /lustre04/scratch/celphin/Oxyria/gene_plots
cd /lustre04/scratch/celphin/Oxyria/gene_plots

# get chromosome lengths

# remove sequneces that are too short
module load StdEnv/2023 seqkit/2.5.1
seqkit seq -m 10000000 ./Dryas/Dry-octo-H2_DoctH0_Main.fasta > DryOcto_chr.fasta
seqkit seq -m 10000000 Oxyria_Main.fasta > Oxyria_digyna_chr.fasta

module load  StdEnv/2020 bioawk/1.0
bioawk -c fastx '{print $name "\t" length($seq)}' DryOcto_chr.fasta  > DryOcto_chr_sizes.txt
bioawk -c fastx '{print $name "\t" length($seq)}' Oxyria_digyna_chr.fasta > Oxyria_digyna_chr_sizes.txt

#-----------------------

module load StdEnv/2023
module load r/4.4.0
 
R

library(dplyr)
library(ggplot2) 
library(tidyverse)
#library(statebins)

# import a text file with gene positions
# Dryas
Dry_genes0 <- read.table("/lustre04/scratch/celphin/Oxyria/GeneSpace/Total_genomes/genomes/Dryas_octopetala/Dryas_octopetala.gff3",sep="\t",header=F)
Dry_chr_sizes0 <- read.table("/lustre04/scratch/celphin/Oxyria/EDTA/DryOcto_chr_sizes.txt",sep="\t",header=F)
Dry_TE_repeats0 <- read.table("/lustre04/scratch/celphin/Oxyria/EDTA/DoctH0_Main.fasta.mod.EDTA.TEanno.gff3",sep="\t",header=F)
Dry_wgd0 <- read.table("/lustre04/scratch/celphin/Oxyria/DupGen_finder/output/7.wgd.pairs",sep="\t",header=T)
Dry_tandem0 <- read.table("/lustre04/scratch/celphin/Oxyria/DupGen_finder/output/7.tandem.pairs",sep="\t",header=T)
Dry_proximal0 <- read.table("/lustre04/scratch/celphin/Oxyria/DupGen_finder/output/7.proximal.pairs",sep="\t",header=T)
Dry_transposed0 <- read.table("/lustre04/scratch/celphin/Oxyria/DupGen_finder/output/7.transposed.pairs",sep="\t",header=T)
Dry_dispersed0 <- read.table("/lustre04/scratch/celphin/Oxyria/DupGen_finder/output/7.dispersed.pairs",sep="\t",header=T)

# Oxyria
Oxy_genes0 <- read.table("/lustre04/scratch/celphin/Oxyria/GeneSpace/Total_genomes/genomes/Oxyria_digyna_H1/Oxyria_digyna_H1.gff3",sep="\t",header=F)
Oxy_chr_sizes0 <- read.table("/lustre04/scratch/celphin/Oxyria/EDTA/Oxyria_digyna_chr_sizes.txt",sep="\t",header=F)
Oxy_TE_repeats1 <- read.table("/lustre04/scratch/celphin/Oxyria/EDTA/Oxyria_digyna.fasta.mod.EDTA.TEanno.gff3",sep="\t",header=F)
Oxy_TE_repeats0 <- read.table("/lustre04/scratch/celphin/Oxyria/EDTA/Oxyria_Main.fasta.mod.EDTA.intact.gff3",sep="\t",header=F)

Oxy_wgd0 <- read.table("/lustre04/scratch/celphin/Oxyria/DupGen_finder/output/12.wgd.pairs",sep="\t",header=T)
Oxy_tandem0 <- read.table("/lustre04/scratch/celphin/Oxyria/DupGen_finder/output/12.tandem.pairs",sep="\t",header=T)
Oxy_proximal0 <- read.table("/lustre04/scratch/celphin/Oxyria/DupGen_finder/output/12.proximal.pairs",sep="\t",header=T)
Oxy_transposed0 <- read.table("/lustre04/scratch/celphin/Oxyria/DupGen_finder/output/12.transposed.pairs",sep="\t",header=T)
Oxy_dispersed0 <- read.table("/lustre04/scratch/celphin/Oxyria/DupGen_finder/output/12.dispersed.pairs",sep="\t",header=T)

# DupGen seq IDs
SequenceIDs <- read.table("/lustre04/scratch/celphin/Oxyria/DupGen_finder/data/SequenceIDs.txt",sep=":",header=F)

# Interproscan data for all species
Gene_ont_file <- "/lustre04/scratch/celphin/Oxyria/synteny_quantity/Total_interproscan_output_edited3.tsv"
gene_ont <- read.delim(Gene_ont_file, header = TRUE, sep = "\t", na.strings = "-", colClasses = c("character", "character", "character", "character"))

# DMR gene lists
Dryas_gene_ont <- read.delim("/lustre04/scratch/celphin/Dryas/GO_enrichment/interproscan_dryas_full3.tsv", header = TRUE, sep = "\t", na.strings = "-")
DMR_DEG <- read.delim("/lustre04/scratch/celphin/Dryas/GO_enrichment/genes_RNA_MethylkitDMR_merged_data.tsv", header = TRUE, sep = "\t")

#----------------------
# DMRs per site

# "ALAS_W_C"
DMR_Alaska_W_C <- DMR_DEG[which(DMR_DEG$site=="ALAS_W_C" & DMR_DEG$perdiff== 10 & DMR_DEG$random=="non-rand" & DMR_DEG$context=="CpG"),]
Alaska_W_C <- as.data.frame(cbind(DMR_Alaska_W_C$Gene, DMR_Alaska_W_C$qvalue))
Alaska_W_C <- distinct(Alaska_W_C)
write.table(Alaska_W_C, "Alaska_W_C_genes.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
#142

# "LAT_W_C"
DMR_Sweden_W_C <- DMR_DEG[which(DMR_DEG$site=="LAT_W_C" & DMR_DEG$perdiff== 10 & DMR_DEG$random=="non-rand" & DMR_DEG$context=="CpG"),]
Sweden_W_C <- as.data.frame(cbind(DMR_Sweden_W_C$Gene, DMR_Sweden_W_C$qvalue))
Sweden_W_C <- distinct(Sweden_W_C)
write.table(Sweden_W_C, "Sweden_W_C_genes.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
# 808

# "CASS_W_C"
DMR_Nunavut_W_C <- DMR_DEG[which(DMR_DEG$site=="CASS_W_C" & DMR_DEG$perdiff== 10 & DMR_DEG$random=="non-rand" & DMR_DEG$context=="CpG"),]
Nunavut_W_C <- as.data.frame(cbind(DMR_Nunavut_W_C$Gene, DMR_Nunavut_W_C$qvalue))
Nunavut_W_C <- distinct(Nunavut_W_C)
write.table(Nunavut_W_C, "Nunavut_W_C_genes.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
#78

# "SVAL_W_C"
DMR_Svalbard_W_C <- DMR_DEG[which(DMR_DEG$site=="SVAL_W_C" & DMR_DEG$perdiff== 10 & DMR_DEG$random=="non-rand" & DMR_DEG$context=="CpG"),]
Svalbard_W_C <- as.data.frame(cbind(DMR_Svalbard_W_C$Gene, DMR_Svalbard_W_C$qvalue))
Svalbard_W_C <- distinct(Svalbard_W_C)
write.table(Svalbard_W_C, "Svalbard_W_C_genes.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
#103

# "Pheno"
DMR_Mat_Sen <- DMR_DEG[which(DMR_DEG$site=="Pheno" & DMR_DEG$perdiff== 10 & DMR_DEG$random=="non-rand" & DMR_DEG$context=="CpG"),]
Mat_Sen <- as.data.frame(cbind(DMR_Mat_Sen$Gene, DMR_Mat_Sen$qvalue))
Mat_Sen <- distinct(Mat_Sen)
write.table(Mat_Sen, "Mat_Sen_genes.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
# 8584

# "HL"
DMR_Wild_Lat_L_H <- DMR_DEG[which(DMR_DEG$site=="HL" & DMR_DEG$perdiff== 10 & DMR_DEG$random=="non-rand" & DMR_DEG$context=="CpG"),]
Wild_Lat_L_H <- as.data.frame(cbind(DMR_Wild_Lat_L_H$Gene, DMR_Wild_Lat_L_H$qvalue))
Wild_Lat_L_H <- distinct(Wild_Lat_L_H)
write.table(Wild_Lat_L_H, "Wild_Lat_L_H_genes.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
#7621

# "SE_W_C"
DMR_SE_W_C <- DMR_DEG[which(DMR_DEG$site=="SE_W_C" & DMR_DEG$perdiff== 10 & DMR_DEG$random=="non-rand" & DMR_DEG$context=="CpG"),]
SE_W_C <- as.data.frame(cbind(DMR_SE_W_C$Gene, DMR_SE_W_C$qvalue))
SE_W_C <- distinct(SE_W_C)
write.table(SE_W_C, "SE_W_C_genes.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
#1434

# "SE_HL"
DMR_SE_L_H <- DMR_DEG[which(DMR_DEG$site=="SE_HL" & DMR_DEG$perdiff== 10 & DMR_DEG$random=="non-rand" & DMR_DEG$context=="CpG"),]
SE_L_H <- as.data.frame(cbind(DMR_SE_L_H$Gene, DMR_SE_L_H$qvalue))
SE_L_H <- distinct(SE_L_H)
write.table(SE_L_H, "SE_L_H_genes.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
# 388

#-------------------
# Extract RNA info separately

# "LAT_W_C"
LAT_RNA <- DMR_DEG[which(DMR_DEG$RNAsite=="LAT_W_C"),]
LAT_RNA_gene <- as.data.frame(cbind(LAT_RNA$Gene, LAT_RNA$PValue))
LAT_RNA_gene <- distinct(LAT_RNA_gene)
write.table(LAT_RNA_gene, "Sweden_RNA_genes.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
#59

# "ALAS_W_C"
ALAS_RNA <- DMR_DEG[which(DMR_DEG$RNAsite=="ALAS_W_C"),]
ALAS_RNA_gene <- as.data.frame(cbind(ALAS_RNA$Gene, ALAS_RNA$PValue))
ALAS_RNA_gene <- distinct(ALAS_RNA_gene)
write.table(ALAS_RNA_gene , "Alaska_RNA_genes.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
#19

# "ALEX_W_C"
ALEX_RNA <- DMR_DEG[which(DMR_DEG$RNAsite=="ALEX_W_C"),]
ALEX_RNA_gene <- as.data.frame(cbind(ALEX_RNA$Gene, ALEX_RNA$PValue))
ALEX_RNA_gene <- distinct(ALEX_RNA_gene)
write.table(ALEX_RNA_gene , "Nunavut_RNA_genes.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
#21

# "NORW_W_C"
NORW_RNA <- DMR_DEG[which(DMR_DEG$RNAsite=="NORW_W_C"),]
NORW_RNA_gene <- as.data.frame(cbind(NORW_RNA$Gene, NORW_RNA$PValue))
NORW_RNA_gene <- distinct(NORW_RNA_gene)
write.table(NORW_RNA_gene , "Norway_RNA_genes.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
#6

# "SE_W_C"
SE_RNA <- DMR_DEG[which(DMR_DEG$RNAsite=="SE_W_C"),]
SE_RNA_gene <- as.data.frame(cbind(SE_RNA$Gene, SE_RNA$PValue))
SE_RNA_gene <- distinct(SE_RNA_gene)
write.table(SE_RNA_gene , "Seedling_RNA_genes.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
#45



#############################

Spp_genes0 <- Oxy_genes0
Spp_TE_repeats0 <- Oxy_TE_repeats0
Spp="Oxyria_digyna"


Spp_genes0 <- Dry_genes0
Spp_TE_repeats0 <- Dry_TE_repeats0
Spp="Dryas_octopetala_H0"

#-------------------------
# Plot genes
#---------------------------
scaffold_lengths <- Spp_genes0 %>%
  as.data.frame() %>%
  group_by(V1) %>%
  summarise(length = max(V5)) %>%
  ungroup()
  
threshold <- 1e7
long_scaffolds <- scaffold_lengths %>%
  filter(length > threshold) %>%
  pull(V1)

filtered_genes <- Spp_genes0[which(Spp_genes0$V1 %in% long_scaffolds & Spp_genes0$V3 == "gene"),]

# Edit so columns are: chr, position (no end or gene name required)
Spp_genes <- as.data.frame(cbind(filtered_genes$V1, filtered_genes$V4))
colnames(Spp_genes) <- c("chr", "pos")
Spp_genes$pos <- as.numeric(Spp_genes$pos)

# make a histogram plot of genes over the provided chromosomes 
plottedSppGenes <- ggplot(Spp_genes) + 
	geom_histogram(aes(x=pos),binwidth=1000000) + 
	facet_wrap(~chr,ncol=1) + 
	xlab("Genomic position (bins 1 Mb)") + 
	ggplot2::theme_classic() +
	ylab("Number of genes")

# save it to an image
png(paste0("./plots/", Spp, "_gene_density.png"),width=700,height=1500)
print(plottedSppGenes)
dev.off()

pdf(paste0("./plots/", Spp, "_gene_density.pdf"))
print(plottedSppGenes)
dev.off()

#---------------------------------
# run through all the repeat types
#--------------------------------
# Step 1: Filter out repeats on short scaffolds
scaffold_lengths <- Spp_TE_repeats0 %>%
  as.data.frame() %>%
  group_by(V1) %>%
  summarise(length = max(V5)) %>%
  ungroup()

threshold <- 1e7
long_scaffolds <- scaffold_lengths %>%
  filter(length > threshold) %>%
  pull(V1)

# Step 2: Get unique repeat types
unique_repeat_types <- unique(Spp_TE_repeats0$V3)

# Step 3: Loop through each repeat type and create plots
for (repeat_type in unique_repeat_types) {
  
  # Filter for the current repeat type
  filtered_repeats <- Spp_TE_repeats0 %>%
    filter(V1 %in% long_scaffolds & V3 == repeat_type)

  # Check if the filtered data is not empty
  if (nrow(filtered_repeats) > 0) {
    # Prepare the data for plotting
    Spp_repeats <- as.data.frame(cbind(filtered_repeats$V1, filtered_repeats$V4))
    colnames(Spp_repeats) <- c("chr", "pos")
    Spp_repeats$pos <- as.numeric(Spp_repeats$pos)

    # Create the histogram plot
    plottedSpp_repeats <- ggplot(Spp_repeats) + 
      geom_histogram(aes(x = pos), binwidth = 1000000) + 
      facet_wrap(~ chr, ncol = 1) + 
      xlab("Genomic position (bins 1 Mb)") + 
      theme_classic() +
      ylab("Number of repeats") +
      ggtitle(paste("Histogram of", repeat_type))

    # Step 4: Save the plot to an image file
    filename <- paste0("./plots/", Spp, "_repeats_", gsub(" ", "_", repeat_type), "_density.png")
	filename1 <- paste0("./plots/", Spp, "_repeats_", gsub(" ", "_", repeat_type), "_density.pdf")
    png(filename, width = 700, height = 1500)
    print(plottedSpp_repeats)
    dev.off()
	pdf(filename1)
	print(plottedSpp_repeats)
	dev.off()
  } else {
    message(paste("No data on more than 1 chromosome for repeat type:", repeat_type))
  }
}

# Oxyria - No data for repeat type: PIF_Harbinger_TIR_transposon

#----------------------------
# Gene duplicates
#--------------------------
Spp="Dryas_octopetala_H0"
wgddata <- Dry_wgd0
tanddata <- Dry_tandem0
proxdata <- Dry_proximal0
transdata <- Dry_transposed0 
dispdata <- Dry_dispersed0

Spp="Oxyria_digyna"
wgddata <- Oxy_wgd0
tanddata <- Oxy_tandem0
proxdata <- Oxy_proximal0
transdata <- Oxy_transposed0 
dispdata <- Oxy_dispersed0

#---------------------------
#WGD
 # Step 1: Split Location columns to extract chromosome and position
location_data <- wgddata %>%
  select(Location, Location.1) %>%
  pivot_longer(cols = everything(), names_to = "Type", values_to = "Location") %>%
  separate(Location, into = c("chr", "pos"), sep = ":", convert = TRUE)

scaffold_lengths <- location_data %>%
  as.data.frame() %>%
  group_by(chr) %>%
  summarise(length = max(pos)) %>%
  ungroup()

threshold <- 1e7
long_scaffolds <- scaffold_lengths %>%
  filter(length > threshold) %>%
  pull(chr)

  # Filter for the current repeat type
  filtered_data <- location_data %>%
    filter(chr %in% long_scaffolds)

location_data <- filtered_data

# Convert pos to numeric
location_data$pos <- as.numeric(location_data$pos)

# Step 2: Create the histogram plot
plotted_locations <- ggplot(location_data) + 
  geom_histogram(aes(x = pos), binwidth = 1000000) + 
  facet_wrap(~chr, ncol = 1) + 
  xlab("Genomic position (bins 1 Mb)") + 
  theme_classic() +
  ylab("Number of locations") +
  ggtitle("Histogram of Locations")

# Step 3: Save the plot to an image file
png(paste0(Spp, "_wgd_location_histogram.png"), width = 700, height = 1500)
print(plotted_locations)
dev.off()
pdf(paste0(Spp, "_wgd_location_histogram.pdf"))
print(plotted_locations)
dev.off()
#---------------------------
#Tandem
 # Step 1: Split Location columns to extract chromosome and position
location_data <- tanddata %>%
  select(Location, Location.1) %>%
  pivot_longer(cols = everything(), names_to = "Type", values_to = "Location") %>%
  separate(Location, into = c("chr", "pos"), sep = ":", convert = TRUE)

scaffold_lengths <- location_data %>%
  as.data.frame() %>%
  group_by(chr) %>%
  summarise(length = max(pos)) %>%
  ungroup()

threshold <- 1e7
long_scaffolds <- scaffold_lengths %>%
  filter(length > threshold) %>%
  pull(chr)

  # Filter for the current repeat type
  filtered_data <- location_data %>%
    filter(chr %in% long_scaffolds)

location_data <- filtered_data

# Convert pos to numeric
location_data$pos <- as.numeric(location_data$pos)

# Step 2: Create the histogram plot
plotted_locations <- ggplot(location_data) + 
  geom_histogram(aes(x = pos), binwidth = 1000000) + 
  facet_wrap(~chr, ncol = 1) + 
  xlab("Genomic position (bins 1 Mb)") + 
  theme_classic() +
  ylab("Number of locations") +
  ggtitle("Histogram of Locations")

# Step 3: Save the plot to an image file
png(paste0("./plots/", Spp, "_tandem_location_histogram.png"), width = 700, height = 1500)
print(plotted_locations)
dev.off()
pdf(paste0("./plots/", Spp, "_tandem_location_histogram.pdf"))
print(plotted_locations)
dev.off()
#---------------------------
# Proximal
 # Step 1: Split Location columns to extract chromosome and position
location_data <- proxdata %>%
  select(Location, Location.1) %>%
  pivot_longer(cols = everything(), names_to = "Type", values_to = "Location") %>%
  separate(Location, into = c("chr", "pos"), sep = ":", convert = TRUE)

scaffold_lengths <- location_data %>%
  as.data.frame() %>%
  group_by(chr) %>%
  summarise(length = max(pos)) %>%
  ungroup()

threshold <- 1e7
long_scaffolds <- scaffold_lengths %>%
  filter(length > threshold) %>%
  pull(chr)

  # Filter for the current repeat type
  filtered_data <- location_data %>%
    filter(chr %in% long_scaffolds)

location_data <- filtered_data

# Convert pos to numeric
location_data$pos <- as.numeric(location_data$pos)

# Step 2: Create the histogram plot
plotted_locations <- ggplot(location_data) + 
  geom_histogram(aes(x = pos), binwidth = 1000000) + 
  facet_wrap(~chr, ncol = 1) + 
  xlab("Genomic position (bins 1 Mb)") + 
  theme_classic() +
  ylab("Number of locations") +
  ggtitle("Histogram of Locations")

# Step 3: Save the plot to an image file
png(paste0("./plots/", Spp, "_proximal_location_histogram.png"), width = 700, height = 1500)
print(plotted_locations)
dev.off()
pdf(paste0("./plots/", Spp, "_proximal_location_histogram.pdf"))
print(plotted_locations)
dev.off()
#---------------------------
# Transposed
 # Step 1: Split Location columns to extract chromosome and position
location_data <- transdata %>%
  select(Location, Location.1) %>%
  pivot_longer(cols = everything(), names_to = "Type", values_to = "Location") %>%
  separate(Location, into = c("chr", "pos"), sep = ":", convert = TRUE)

scaffold_lengths <- location_data %>%
  as.data.frame() %>%
  group_by(chr) %>%
  summarise(length = max(pos)) %>%
  ungroup()

threshold <- 1e7
long_scaffolds <- scaffold_lengths %>%
  filter(length > threshold) %>%
  pull(chr)

  # Filter for the current repeat type
  filtered_data <- location_data %>%
    filter(chr %in% long_scaffolds)

location_data <- filtered_data

# Convert pos to numeric
location_data$pos <- as.numeric(location_data$pos)

# Step 2: Create the histogram plot
plotted_locations <- ggplot(location_data) + 
  geom_histogram(aes(x = pos), binwidth = 1000000) + 
  facet_wrap(~chr, ncol = 1) + 
  xlab("Genomic position (bins 1 Mb)") + 
  theme_classic() +
  ylab("Number of locations") +
  ggtitle("Histogram of Locations")

# Step 3: Save the plot to an image file
png(paste0("./plots/", Spp, "_transposed_location_histogram.png"), width = 700, height = 1500)
print(plotted_locations)
dev.off()
pdf(paste0("./plots/", Spp, "_transposed_location_histogram.pdf"))
print(plotted_locations)
dev.off()
#---------------------------
# Dispersed
 # Step 1: Split Location columns to extract chromosome and position
location_data <- dispdata %>%
  select(Location, Location.1) %>%
  pivot_longer(cols = everything(), names_to = "Type", values_to = "Location") %>%
  separate(Location, into = c("chr", "pos"), sep = ":", convert = TRUE)

scaffold_lengths <- location_data %>%
  as.data.frame() %>%
  group_by(chr) %>%
  summarise(length = max(pos)) %>%
  ungroup()

threshold <- 1e7
long_scaffolds <- scaffold_lengths %>%
  filter(length > threshold) %>%
  pull(chr)

  # Filter for the current repeat type
  filtered_data <- location_data %>%
    filter(chr %in% long_scaffolds)

location_data <- filtered_data

# Convert pos to numeric
location_data$pos <- as.numeric(location_data$pos)

# Step 2: Create the histogram plot
plotted_locations <- ggplot(location_data) + 
  geom_histogram(aes(x = pos), binwidth = 1000000) + 
  facet_wrap(~chr, ncol = 1) + 
  xlab("Genomic position (bins 1 Mb)") + 
  theme_classic() +
  ylab("Number of locations") +
  ggtitle("Histogram of Locations")

# Step 3: Save the plot to an image file
png(paste0("./plots/", Spp, "_dispersed_location_histogram.png"), width = 700, height = 1500)
print(plotted_locations)
dev.off()
pdf(paste0("./plots/", Spp, "_dispersed_location_histogram.pdf"))
print(plotted_locations)
dev.off()

#############################
# Join gene duplicates with Sequence IDs
colnames(SequenceIDs) <- c("Duplicate.1", "gene")
# join with Interproscan data
colnames(gene_ont) <- c("spp", "gene", "INTPRO", "descrip", "GOterm")
length(unique(gene_ont$INTPRO))
unique(gene_ont$spp)
# [1] "Arabis_alpina_interproscan_output.tsv"
# [2] "Cochlearia_groenlandica_interproscan_output.tsv"
# [3] "Draba_nivalis_interproscan_output.tsv"
# [4] "Dryas_octopetala_interproscan_output.tsv"
# [5] "Oxyria_digyna_H1_interproscan_output.tsv"
# [6] "Rheum_nobile_H0_interproscan_output.tsv"

#-----------------------

Spp_tand_genes<- dplyr::left_join(tanddata, SequenceIDs, by="Duplicate.1")
gene_ont_Spp <- gene_ont[which(gene_ont$spp=="Dryas_octopetala_interproscan_output.tsv"),]
gene_ont_Spp$gene <- as.factor(gene_ont_Spp$gene)
Spp_tand_genes$gene <- gsub(" ", "", Spp_tand_genes$gene)
Spp_tand_genes_ont <- dplyr::left_join(Spp_tand_genes, gene_ont_Spp, by="gene")

unique(Spp_tand_genes_ont$GOterm)
unique(Spp_tand_genes_ont$descrip)

# Count occurrences of each unique GOterm
go_counts <- table(Spp_tand_genes_ont$descrip)

# Convert to a data frame for easier viewing
go_counts_df <- as.data.frame(go_counts)

# Rename columns
colnames(go_counts_df) <- c("GOterm", "Count")

# Order by Count in descending order and get the top 10
top_go_counts <- go_counts_df %>%
  arrange(desc(Count)) %>%
  head(10)

# Print the top 10
print(head(top_go_counts))
#  conserved site, E-class, group I,Cytochrome P450,Cytochrome P450 superfamily
# conserved site,UDP-glucuronosyl/UDP-glucosyltransferase,UDP-glycosyltransferase family
# family 28,Glycoside hydrolase,Parallel beta-helix repeat,Pectin lyase fold,Pectin lyase fold/virulence factor
# UDP-glucuronosyl/UDP-glucosyltransferase
# Small auxin-up RNA


#----------------------
Spp_wgd_genes<- dplyr::left_join(wgddata, SequenceIDs, by="Duplicate.1")
gene_ont_Spp <- gene_ont[which(gene_ont$spp=="Dryas_octopetala_interproscan_output.tsv"),]
gene_ont_Spp$gene <- as.factor(gene_ont_Spp$gene)
Spp_wgd_genes$gene <- gsub(" ", "", Spp_wgd_genes$gene)
Spp_wgd_genes_ont <- dplyr::left_join(Spp_wgd_genes, gene_ont_Spp, by="gene")

unique(Spp_wgd_genes_ont$GOterm)
unique(Spp_wgd_genes_ont$descrip)

# Count occurrences of each unique GOterm
go_counts <- table(Spp_wgd_genes_ont$descrip)

# Convert to a data frame for easier viewing
go_counts_df <- as.data.frame(go_counts)

# Rename columns
colnames(go_counts_df) <- c("GOterm", "Count")

# Order by Count in descending order and get the top 10
top_go_counts <- go_counts_df %>%
  arrange(desc(Count)) %>%
  head(10)

# Print the top 10
print(head(top_go_counts))
# DNA integration 
# https://www.ebi.ac.uk/QuickGO/term/GO:0015074
# Integrase,Integrase zinc-binding domain,Retrotransposon gag domain,Reverse transcriptase,Reverse transcriptase domain,Reverse transcriptase/Diguanylate cyclase domain,Ribonuclease H superfamily,Ribonuclease H-like superfamily


####################################
# make a density plot of genes over the provided chromosomes 

library(dplyr)
library(tidyr)
library(ggplot2)
library(grid)

# Prepare gene data
gene_data <- Spp_genes %>%
  mutate(Type = "Gene")

# Prepare repeat data
repeat_data_list <- lapply(unique_repeat_types, function(repeat_type) {
  filtered_repeats <- Spp_TE_repeats0 %>%
    filter(V3 == repeat_type) %>%
    select(chr = V1, pos = V4) %>%
    mutate(Type = repeat_type)
  return(filtered_repeats)
})

# Combine all repeat datasets
repeat_data <- bind_rows(repeat_data_list)

# Define the minimum chromosome length
min_length <- 1e7  # Set your desired minimum length here

# Combine gene and repeat data
combined_data <- bind_rows(
  gene_data,
  repeat_data
)

# Set up a plotting area for each chromosome and type
plot_data <- combined_data %>%
  group_by(chr) %>%
  summarise(Start = min(pos), End = max(pos), .groups = 'drop') %>%
  ungroup() %>%
  mutate(Length = End - Start) %>%
  filter(Length >= min_length)

# Filter combined data based on calculated chromosome lengths
filtered_combined_data <- combined_data %>%
  inner_join(plot_data, by = "chr")

# Create a data frame with y positions dynamically
unique_types=3
n_types <- length(unique_types)  # Count unique types

# Generate y positions based on the number of unique types
y_positions <- data.frame(
  Type = unique_types,
  ymin = seq(0, (n_types - 1) * 0.1, by = 0.1),  # Each type spaced by 0.1
  ymax = seq(0.1, n_types * 0.1, by = 0.1)        # Each type spaced by 0.1
)

# Convert chr to a factor or numeric if it's not already
filtered_combined_data$chr <- as.factor(filtered_combined_data$chr)
plot_data$chr <- as.factor(plot_data$chr)

# Create plot
png("./plots/Total_gene_repeat_density_plot.png", width = 1500, height = 900)
print(ggplot() +
  geom_rect(data = plot_data,
            aes(xmin = Start, xmax = End, ymin = as.numeric(chr)-0.5, ymax = as.numeric(chr)+0.5),
            colour="black", fill = NA, alpha = 0.2) +

  # Adjust the y aesthetics to use factor levels for chr
  geom_segment(data = filtered_combined_data %>% filter(Type == "Gene"),
               aes(x = pos, xend = pos, 
                   y = as.numeric(chr)-0.5, yend = as.numeric(chr)-0.5+0.33),
               color = "tomato3", alpha = 0.05, size = 1) +

  geom_segment(data = filtered_combined_data %>% filter(Type == "LTR_retrotransposon"),
               aes(x = pos, xend = pos, 
                   y = as.numeric(chr)-0.5+0.33, yend = as.numeric(chr)-0.5 +0.67),
               color = "blue", alpha = 0.05, size = 1) +

  geom_segment(data = filtered_combined_data %>% filter(Type == "repeat_region"),
               aes(x = pos, xend = pos, 
                   y = as.numeric(chr)-0.5+0.67, yend = as.numeric(chr)+0.5),
               color = "green", alpha = 0.05, size = 1) +

  scale_y_reverse() +
  labs(x = "Genomic Position (Mbp)", y = "Chromosome") +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")
)
dev.off()

# Create plot
pdf("./plots/Total_gene_repeat_density_plot.pdf")
print(ggplot() +
  geom_rect(data = plot_data,
            aes(xmin = Start, xmax = End, ymin = as.numeric(chr)-0.5, ymax = as.numeric(chr)+0.5),
            colour="black", fill = NA, alpha = 0.2) +

  # Adjust the y aesthetics to use factor levels for chr
  geom_segment(data = filtered_combined_data %>% filter(Type == "Gene"),
               aes(x = pos, xend = pos, 
                   y = as.numeric(chr)-0.5, yend = as.numeric(chr)-0.5+0.33),
               color = "tomato3", alpha = 0.05, size = 0.1) +

  geom_segment(data = filtered_combined_data %>% filter(Type == "LTR_retrotransposon"),
               aes(x = pos, xend = pos, 
                   y = as.numeric(chr)-0.5+0.33, yend = as.numeric(chr)-0.5 +0.67),
               color = "blue", alpha = 0.05, size = 0.1) +

  geom_segment(data = filtered_combined_data %>% filter(Type == "repeat_region"),
               aes(x = pos, xend = pos, 
                   y = as.numeric(chr)-0.5+0.67, yend = as.numeric(chr)+0.5),
               color = "green", alpha = 0.05, size = 0.1) +

  scale_y_reverse() +
  labs(x = "Genomic Position (Mbp)", y = "Chromosome") +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")
)
dev.off()

###############################
# Make some similar plots of Interproscan data 

# join with chr and start Position
library(dplyr)
library(stringr)

Spp_genes0 <- Spp_genes0 %>%
  mutate(gene = str_extract(V9, "(?<=ID=)[^;]+"))

gene_ont_Spp0 <- left_join(gene_ont_Spp, Spp_genes0, by="gene")

# test
gene_ont_Spp0$descrip[grep("Homeobox domain", gene_ont_Spp0$descrip)]
gene_ont_Spp0$descrip[grep("MADS", gene_ont_Spp0$descrip)]
gene_ont_Spp0$descrip[grep("transcription factor", gene_ont_Spp0$descrip)]


# extract data
geneont_data1 <- gene_ont_Spp0 %>%
  filter(grepl("Homeobox domain", descrip))%>%
  select(chr = V1, pos = V4) %>%
  mutate(Type = "Homeobox")

geneont_data2 <- gene_ont_Spp0 %>%
  filter(grepl("MADS", descrip))%>%
  select(chr = V1, pos = V4) %>%
  mutate(Type = "MADS")

geneont_data3 <- gene_ont_Spp0 %>%
  filter(grepl("transcription factor", descrip))%>%
  select(chr = V1, pos = V4) %>%
  mutate(Type = "transcription factor")

# Define the minimum chromosome length
min_length <- 1e7  # Set your desired minimum length here

# Combine gene and repeat data
combined_data <- bind_rows(
  geneont_data1,
  geneont_data2,
  geneont_data3
)

unique_combined_data <- combined_data %>%
  distinct()

# Set up a plotting area for each chromosome and type
plot_data <- unique_combined_data %>%
  group_by(chr) %>%
  summarise(Start = min(pos), End = max(pos), .groups = 'drop') %>%
  ungroup() %>%
  mutate(Length = End - Start) %>%
  filter(Length >= min_length)

# Filter combined data based on calculated chromosome lengths
filtered_combined_data <- unique_combined_data %>%
  inner_join(plot_data, by = "chr")

# Create a data frame with y positions dynamically
n_types <- length(unique_types)  # Count unique types

# Generate y positions based on the number of unique types
y_positions <- data.frame(
  Type = unique_types,
  ymin = seq(0, (n_types - 1) * 0.1, by = 0.1),  # Each type spaced by 0.1
  ymax = seq(0.1, n_types * 0.1, by = 0.1)        # Each type spaced by 0.1
)

# Convert chr to a factor or numeric if it's not already
filtered_combined_data$chr <- as.factor(filtered_combined_data$chr)
plot_data$chr <- as.factor(plot_data$chr)

# Create plot
png("./plots/Homeobox_pos_plot.png", width = 1200, height = 700)
print(ggplot() +
  geom_rect(data = plot_data,
            aes(xmin = 0, xmax = End+1000, ymin = as.numeric(chr)-0.5, ymax = as.numeric(chr)+0.5),
            colour="black", fill = NA, alpha = 0.2) +

  # Adjust the y aesthetics to use factor levels for chr
  geom_segment(data = filtered_combined_data %>% filter(Type == "Homeobox"),
               aes(x = pos, xend = pos, 
                   y = as.numeric(chr)-0.5, yend = as.numeric(chr)-0.5+0.33),
               color = "tomato3", alpha = 0.5, size = 2) +

  geom_segment(data = filtered_combined_data %>% filter(Type == "MADS"),
               aes(x = pos, xend = pos, 
                   y = as.numeric(chr)-0.5+0.33, yend = as.numeric(chr)-0.5 +0.67),
               color = "blue", alpha = 0.5, size = 2) +

  geom_segment(data = filtered_combined_data %>% filter(Type == "transcription factor"),
               aes(x = pos, xend = pos, 
                   y = as.numeric(chr)-0.5+0.67, yend = as.numeric(chr)+0.5),
               color = "green", alpha = 0.5, size = 2) +

  scale_y_reverse() +
  labs(x = "Genomic Position (Mbp)", y = "Chromosome") +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")
)
dev.off()

# Create plot
pdf("./plots/Homeobox_pos_plot.pdf", width = 1200, height = 700)
print(ggplot() +
  geom_rect(data = plot_data,
            aes(xmin = 0, xmax = End+1000, ymin = as.numeric(chr)-0.5, ymax = as.numeric(chr)+0.5),
            colour="black", fill = NA, alpha = 0.2) +

  # Adjust the y aesthetics to use factor levels for chr
  geom_segment(data = filtered_combined_data %>% filter(Type == "Homeobox"),
               aes(x = pos, xend = pos, 
                   y = as.numeric(chr)-0.5, yend = as.numeric(chr)-0.5+0.33),
               color = "tomato3", alpha = 0.5, size = 2) +

  geom_segment(data = filtered_combined_data %>% filter(Type == "MADS"),
               aes(x = pos, xend = pos, 
                   y = as.numeric(chr)-0.5+0.33, yend = as.numeric(chr)-0.5 +0.67),
               color = "blue", alpha = 0.5, size = 2) +

  geom_segment(data = filtered_combined_data %>% filter(Type == "transcription factor"),
               aes(x = pos, xend = pos, 
                   y = as.numeric(chr)-0.5+0.67, yend = as.numeric(chr)+0.5),
               color = "green", alpha = 0.5, size = 2) +

  scale_y_reverse() +
  labs(x = "Genomic Position (Mbp)", y = "Chromosome") +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")
)
dev.off()

###############################
# Look at types of TFs 

# join with chr and start Position
library(dplyr)
library(stringr)

Spp_genes0 <- Spp_genes0 %>%
  mutate(gene = str_extract(V9, "(?<=ID=)[^;]+"))

gene_ont_Spp0 <- left_join(gene_ont_Spp, Spp_genes0, by="gene")

# test
gene_ont_Spp0$descrip[grep("transcription factor", gene_ont_Spp0$descrip)]

TF_1 = "RF2" #"red" "top"
TF_2 = "bZIP" #"pink" "mid"
TF_3 = "GATA" #"orange" "bottom"
TF_4 = "Myc-type" #"yellow" "top"
TF_5 = "AP2/ERF" #"green" "mid"
TF_6 = "Heat shock" #"darkgreen" "bottom"
TF_7 = "SANT/Myb" #"blue" "top"
TF_8 = "PIF1" #"black" "mid"
TF_9 = "WRKY" #"cyan" "bottom"

#------------------
gene_ont_Spp0$descrip[grep(TF_1, gene_ont_Spp0$descrip)]
gene_ont_Spp0$descrip[grep(TF_2, gene_ont_Spp0$descrip)]
gene_ont_Spp0$descrip[grep(TF_3, gene_ont_Spp0$descrip)]

# extract data
geneont_data1 <- gene_ont_Spp0 %>%
  filter(grepl(TF_1, descrip))%>%
  select(chr = V1, pos = V4) %>%
  mutate(Type = TF_1)

geneont_data2 <- gene_ont_Spp0 %>%
  filter(grepl(TF_2, descrip))%>%
  select(chr = V1, pos = V4) %>%
  mutate(Type = TF_2)

geneont_data3 <- gene_ont_Spp0 %>%
  filter(grepl(TF_3, descrip))%>%
  select(chr = V1, pos = V4) %>%
  mutate(Type = TF_3)

geneont_data4 <- gene_ont_Spp0 %>%
  filter(grepl(TF_4, descrip))%>%
  select(chr = V1, pos = V4) %>%
  mutate(Type = TF_4)
  
geneont_data5 <- gene_ont_Spp0 %>%
  filter(grepl(TF_5, descrip))%>%
  select(chr = V1, pos = V4) %>%
  mutate(Type = TF_5)
  
geneont_data6 <- gene_ont_Spp0 %>%
  filter(grepl(TF_6, descrip))%>%
  select(chr = V1, pos = V4) %>%
  mutate(Type = TF_6)
  
geneont_data7 <- gene_ont_Spp0 %>%
  filter(grepl(TF_7, descrip))%>%
  select(chr = V1, pos = V4) %>%
  mutate(Type = TF_7)
  
  geneont_data8 <- gene_ont_Spp0 %>%
  filter(grepl(TF_8, descrip))%>%
  select(chr = V1, pos = V4) %>%
  mutate(Type = TF_8)

geneont_data9 <- gene_ont_Spp0 %>%
  filter(grepl(TF_9, descrip))%>%
  select(chr = V1, pos = V4) %>%
  mutate(Type = TF_9)

# Define the minimum chromosome length
min_length <- 1e7  # Set your desired minimum length here

# Combine gene and repeat data
combined_data <- bind_rows(
  geneont_data1,
  geneont_data2,
  geneont_data3,
  geneont_data4,
  geneont_data5,
  geneont_data6,
  geneont_data7,
  geneont_data8,
  geneont_data9
)

unique_combined_data <- combined_data %>%
  distinct()

# Set up a plotting area for each chromosome and type
plot_data <- unique_combined_data %>%
  group_by(chr) %>%
  summarise(Start = min(pos), End = max(pos), .groups = 'drop') %>%
  ungroup() %>%
  mutate(Length = End - Start) %>%
  filter(Length >= min_length)

# Filter combined data based on calculated chromosome lengths
filtered_combined_data <- unique_combined_data %>%
  inner_join(plot_data, by = "chr")
unique_types=9
# Create a data frame with y positions dynamically
n_types <- length(unique_types)  # Count unique types

# Generate y positions based on the number of unique types
y_positions <- data.frame(
  Type = unique_types,
  ymin = seq(0, (n_types - 1) * 0.1, by = 0.1),  # Each type spaced by 0.1
  ymax = seq(0.1, n_types * 0.1, by = 0.1)        # Each type spaced by 0.1
)

# Convert chr to a factor or numeric if it's not already
filtered_combined_data$chr <- as.factor(filtered_combined_data$chr)
plot_data$chr <- as.factor(plot_data$chr)

# Create plot
png("./plots/TF_types_pos_plot.png", width = 1200, height = 700)
print(ggplot() +
  geom_rect(data = plot_data,
            aes(xmin = 0, xmax = End+1000, ymin = as.numeric(chr)-0.5, ymax = as.numeric(chr)+0.5),
            colour="black", fill = NA, alpha = 0.2) +

  # Adjust the y aesthetics to use factor levels for chr
  geom_segment(data = filtered_combined_data %>% filter(Type == TF_1),
               aes(x = pos, xend = pos, 
                   y = as.numeric(chr)-0.5, yend = as.numeric(chr)-0.5+0.33),
               color = "red", alpha = 0.5, size = 2) +

  geom_segment(data = filtered_combined_data %>% filter(Type == TF_2),
               aes(x = pos, xend = pos, 
                   y = as.numeric(chr)-0.5+0.33, yend = as.numeric(chr)-0.5 +0.67),
               color = "pink", alpha = 0.5, size = 2) +

  geom_segment(data = filtered_combined_data %>% filter(Type == TF_3),
               aes(x = pos, xend = pos, 
                   y = as.numeric(chr)-0.5+0.67, yend = as.numeric(chr)+0.5),
               color = "orange", alpha = 0.5, size = 2) +

  geom_segment(data = filtered_combined_data %>% filter(Type == TF_4),
               aes(x = pos, xend = pos, 
                   y = as.numeric(chr)-0.5, yend = as.numeric(chr)-0.5+0.33),
               color = "yellow", alpha = 0.5, size = 2) +

  geom_segment(data = filtered_combined_data %>% filter(Type == TF_5),
               aes(x = pos, xend = pos, 
                   y = as.numeric(chr)-0.5+0.33, yend = as.numeric(chr)-0.5 +0.67),
               color = "green", alpha = 0.5, size = 2) +

  geom_segment(data = filtered_combined_data %>% filter(Type == TF_6),
               aes(x = pos, xend = pos, 
                   y = as.numeric(chr)-0.5+0.67, yend = as.numeric(chr)+0.5),
               color = "darkgreen", alpha = 0.5, size = 2) +

  geom_segment(data = filtered_combined_data %>% filter(Type == TF_7),
               aes(x = pos, xend = pos, 
                   y = as.numeric(chr)-0.5, yend = as.numeric(chr)-0.5+0.33),
               color = "blue", alpha = 0.5, size = 2) +

  geom_segment(data = filtered_combined_data %>% filter(Type == TF_8),
               aes(x = pos, xend = pos, 
                   y = as.numeric(chr)-0.5+0.33, yend = as.numeric(chr)-0.5 +0.67),
               color = "black", alpha = 0.5, size = 2) +

  geom_segment(data = filtered_combined_data %>% filter(Type == TF_9),
               aes(x = pos, xend = pos, 
                   y = as.numeric(chr)-0.5+0.67, yend = as.numeric(chr)+0.5),
               color = "cyan", alpha = 0.5, size = 2) +

  scale_y_reverse() +
  labs(x = "Genomic Position (Mbp)", y = "Chromosome") +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")
)
dev.off()


# Create plot
pdf("./plots/TF_types_pos_plot.pdf", width = 1200, height = 700)
print(ggplot() +
  geom_rect(data = plot_data,
            aes(xmin = 0, xmax = End+1000, ymin = as.numeric(chr)-0.5, ymax = as.numeric(chr)+0.5),
            colour="black", fill = NA, alpha = 0.2) +

  # Adjust the y aesthetics to use factor levels for chr
  geom_segment(data = filtered_combined_data %>% filter(Type == TF_1),
               aes(x = pos, xend = pos, 
                   y = as.numeric(chr)-0.5, yend = as.numeric(chr)-0.5+0.33),
               color = "red", alpha = 0.5, size = 2) +

  geom_segment(data = filtered_combined_data %>% filter(Type == TF_2),
               aes(x = pos, xend = pos, 
                   y = as.numeric(chr)-0.5+0.33, yend = as.numeric(chr)-0.5 +0.67),
               color = "pink", alpha = 0.5, size = 2) +

  geom_segment(data = filtered_combined_data %>% filter(Type == TF_3),
               aes(x = pos, xend = pos, 
                   y = as.numeric(chr)-0.5+0.67, yend = as.numeric(chr)+0.5),
               color = "orange", alpha = 0.5, size = 2) +

  geom_segment(data = filtered_combined_data %>% filter(Type == TF_4),
               aes(x = pos, xend = pos, 
                   y = as.numeric(chr)-0.5, yend = as.numeric(chr)-0.5+0.33),
               color = "yellow", alpha = 0.5, size = 2) +

  geom_segment(data = filtered_combined_data %>% filter(Type == TF_5),
               aes(x = pos, xend = pos, 
                   y = as.numeric(chr)-0.5+0.33, yend = as.numeric(chr)-0.5 +0.67),
               color = "green", alpha = 0.5, size = 2) +

  geom_segment(data = filtered_combined_data %>% filter(Type == TF_6),
               aes(x = pos, xend = pos, 
                   y = as.numeric(chr)-0.5+0.67, yend = as.numeric(chr)+0.5),
               color = "darkgreen", alpha = 0.5, size = 2) +

  geom_segment(data = filtered_combined_data %>% filter(Type == TF_7),
               aes(x = pos, xend = pos, 
                   y = as.numeric(chr)-0.5, yend = as.numeric(chr)-0.5+0.33),
               color = "blue", alpha = 0.5, size = 2) +

  geom_segment(data = filtered_combined_data %>% filter(Type == TF_8),
               aes(x = pos, xend = pos, 
                   y = as.numeric(chr)-0.5+0.33, yend = as.numeric(chr)-0.5 +0.67),
               color = "black", alpha = 0.5, size = 2) +

  geom_segment(data = filtered_combined_data %>% filter(Type == TF_9),
               aes(x = pos, xend = pos, 
                   y = as.numeric(chr)-0.5+0.67, yend = as.numeric(chr)+0.5),
               color = "cyan", alpha = 0.5, size = 2) +

  scale_y_reverse() +
  labs(x = "Genomic Position (Mbp)", y = "Chromosome") +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")
)
dev.off()

###############################
# Make some similar plots of Interproscan data 

# join with chr and start Position
library(dplyr)
library(stringr)

Spp_genes0 <- Spp_genes0 %>%
  mutate(gene = str_extract(V9, "(?<=ID=)[^;]+"))

gene_ont_Spp0 <- left_join(gene_ont_Spp, Spp_genes0, by="gene")

# test
gene_ont_Spp0$descrip[grep("GO:0006952", gene_ont_Spp0$GOterm)]


# extract data
geneont_data1 <- gene_ont_Spp0 %>%
  filter(grepl("Homeobox domain", descrip))%>%
  select(chr = V1, pos = V4) %>%
  mutate(Type = "Homeobox")

geneont_data2 <- gene_ont_Spp0 %>%
  filter(grepl("MADS", descrip))%>%
  select(chr = V1, pos = V4) %>%
  mutate(Type = "MADS")

geneont_data3 <- gene_ont_Spp0 %>%
  filter(grepl("transcription factor", descrip))%>%
  select(chr = V1, pos = V4) %>%
  mutate(Type = "transcription factor")

geneont_data4 <- gene_ont_Spp0 %>%
  filter(grepl("GO:0006952", GOterm))%>%
  select(chr = V1, pos = V4) %>%
  mutate(Type = "GO:0006952")

# Define the minimum chromosome length
min_length <- 1e7  # Set your desired minimum length here

# Combine gene and repeat data
combined_data <- bind_rows(
  geneont_data1,
  geneont_data2,
  geneont_data3,
  geneont_data4
)

unique_combined_data <- combined_data %>%
  distinct()

# Set up a plotting area for each chromosome and type
plot_data <- unique_combined_data %>%
  group_by(chr) %>%
  summarise(Start = min(pos), End = max(pos), .groups = 'drop') %>%
  ungroup() %>%
  mutate(Length = End - Start) %>%
  filter(Length >= min_length)

# Filter combined data based on calculated chromosome lengths
filtered_combined_data <- unique_combined_data %>%
  inner_join(plot_data, by = "chr")

# Create a data frame with y positions dynamically
n_types <- length(unique_types)  # Count unique types

# Generate y positions based on the number of unique types
y_positions <- data.frame(
  Type = unique_types,
  ymin = seq(0, (n_types - 1) * 0.1, by = 0.1),  # Each type spaced by 0.1
  ymax = seq(0.1, n_types * 0.1, by = 0.1)        # Each type spaced by 0.1
)

# Convert chr to a factor or numeric if it's not already
filtered_combined_data$chr <- as.factor(filtered_combined_data$chr)
plot_data$chr <- as.factor(plot_data$chr)

# Create plot
png("./plots/Homeobox_TF_defense_pos_plot.png", width = 1200, height = 700)
print(ggplot() +
  geom_rect(data = plot_data,
            aes(xmin = 0, xmax = End+1000, ymin = as.numeric(chr)-0.5, ymax = as.numeric(chr)+0.5),
            colour="black", fill = NA, alpha = 0.2) +

  # Adjust the y aesthetics to use factor levels for chr
  geom_segment(data = filtered_combined_data %>% filter(Type == "Homeobox"),
               aes(x = pos, xend = pos, 
                   y = as.numeric(chr)-0.5, yend = as.numeric(chr)-0.5+0.33),
               color = "tomato3", alpha = 0.5, size = 2) +

  geom_segment(data = filtered_combined_data %>% filter(Type == "MADS"),
               aes(x = pos, xend = pos, 
                   y = as.numeric(chr)-0.5, yend = as.numeric(chr)-0.5+0.33),
               color = "blue", alpha = 0.5, size = 2) +

  geom_segment(data = filtered_combined_data %>% filter(Type == "transcription factor"),
               aes(x = pos, xend = pos, 
                   y = as.numeric(chr)-0.5+0.33, yend = as.numeric(chr)-0.5 +0.67),
               color = "blue", alpha = 0.5, size = 2) +

  geom_segment(data = filtered_combined_data %>% filter(Type == "GO:0006952"),
               aes(x = pos, xend = pos, 
                   y = as.numeric(chr)-0.5+0.67, yend = as.numeric(chr)+0.5),
               color = "green", alpha = 0.5, size = 2) +

  scale_y_reverse() +
  labs(x = "Genomic Position (Mbp)", y = "Chromosome") +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")
)
dev.off()

# Create plot
pdf("./plots/Homeobox_defense_pos_plot.pdf", width = 1200, height = 700)
print(ggplot() +
  geom_rect(data = plot_data,
            aes(xmin = 0, xmax = End+1000, ymin = as.numeric(chr)-0.5, ymax = as.numeric(chr)+0.5),
            colour="black", fill = NA, alpha = 0.2) +

  # Adjust the y aesthetics to use factor levels for chr
  geom_segment(data = filtered_combined_data %>% filter(Type == "Homeobox"),
               aes(x = pos, xend = pos, 
                   y = as.numeric(chr)-0.5, yend = as.numeric(chr)-0.5+0.33),
               color = "tomato3", alpha = 0.5, size = 2) +

  geom_segment(data = filtered_combined_data %>% filter(Type == "MADS"),
               aes(x = pos, xend = pos, 
                   y = as.numeric(chr)-0.5+0.33, yend = as.numeric(chr)-0.5 +0.67),
               color = "blue", alpha = 0.5, size = 2) +

  geom_segment(data = filtered_combined_data %>% filter(Type == "GO:0006952"),
               aes(x = pos, xend = pos, 
                   y = as.numeric(chr)-0.5+0.67, yend = as.numeric(chr)+0.5),
               color = "green", alpha = 0.5, size = 2) +

  scale_y_reverse() +
  labs(x = "Genomic Position (Mbp)", y = "Chromosome") +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")
)
dev.off()


#####################
# plotting microsynteny regions for specific chromosomes












######################
# try joining with DMRs and microsynteny regions










#############################################

# join with chr and start Position
library(dplyr)
library(stringr)

Spp_genes0 <- Spp_genes0 %>%
  mutate(gene = str_extract(V9, "(?<=ID=)[^;]+"))

gene_ont_Spp0 <- left_join(gene_ont_Spp, Spp_genes0, by="gene")

geneont_data1 <- gene_ont_Spp0 %>%
  filter(grepl("methylation", descrip))%>%
  select(chr = V1, pos = V4) %>%
  mutate(Type = "methylation")

geneont_data2 <- gene_ont_Spp0 %>%
  filter(grepl("GO:0048658", GOterm))%>%
  select(chr = V1, pos = V4) %>%
  mutate(Type = "anther")

geneont_data3 <- gene_ont_Spp0 %>%
  filter(grepl("histone", descrip))%>%
  select(chr = V1, pos = V4) %>%
  mutate(Type = "histone")

# Define the minimum chromosome length
min_length <- 1e7  # Set your desired minimum length here

# Combine gene and repeat data
combined_data <- bind_rows(
  geneont_data1,
  geneont_data2,
  geneont_data3
)

unique_combined_data <- combined_data %>%
  distinct()

# Set up a plotting area for each chromosome and type
plot_data <- unique_combined_data %>%
  group_by(chr) %>%
  summarise(Start = 0, End = max(pos), .groups = 'drop') %>%
  ungroup() %>%
  mutate(Length = End - Start) %>%
  filter(Length >= min_length)

# Filter combined data based on calculated chromosome lengths
filtered_combined_data <- unique_combined_data %>%
  inner_join(plot_data, by = "chr")

# Create a data frame with y positions dynamically
n_types <- length(unique_types)  # Count unique types

# Generate y positions based on the number of unique types
y_positions <- data.frame(
  Type = unique_types,
  ymin = seq(0, (n_types - 1) * 0.1, by = 0.1),  # Each type spaced by 0.1
  ymax = seq(0.1, n_types * 0.1, by = 0.1)        # Each type spaced by 0.1
)

# Convert chr to a factor or numeric if it's not already
filtered_combined_data$chr <- as.factor(filtered_combined_data$chr)
plot_data$chr <- as.factor(plot_data$chr)

# Create plot
png("./plots/Epigenetic_plot.png", width = 1500, height = 900)
print(ggplot() +
  geom_rect(data = plot_data,
            aes(xmin = Start, xmax = End, ymin = as.numeric(chr)-0.5, ymax = as.numeric(chr)+0.5),
            colour="black", fill = NA, alpha = 0.2) +

  # Adjust the y aesthetics to use factor levels for chr
  geom_segment(data = filtered_combined_data %>% filter(Type == "methylation"),
               aes(x = pos, xend = pos, 
                   y = as.numeric(chr)-0.5, yend = as.numeric(chr)-0.5+0.33),
               color = "tomato3", alpha = 0.5, size = 2) +

  geom_segment(data = filtered_combined_data %>% filter(Type == "histone"),
               aes(x = pos, xend = pos, 
                   y = as.numeric(chr)-0.5+0.33, yend = as.numeric(chr)-0.5 +0.67),
               color = "blue", alpha = 0.5, size = 2) +

  geom_segment(data = filtered_combined_data %>% filter(Type == "anther"),
               aes(x = pos, xend = pos, 
                   y = as.numeric(chr)-0.5+0.67, yend = as.numeric(chr)+0.5),
               color = "green", alpha = 0.5, size = 2) +

  scale_y_reverse() +
  labs(x = "Genomic Position (Mbp)", y = "Chromosome") +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")
)
dev.off()

# Create plot
pdf("./plots/Epigenetic_plot.pdf", width = 1500, height = 900)
print(ggplot() +
  geom_rect(data = plot_data,
            aes(xmin = Start, xmax = End, ymin = as.numeric(chr)-0.5, ymax = as.numeric(chr)+0.5),
            colour="black", fill = NA, alpha = 0.2) +

  # Adjust the y aesthetics to use factor levels for chr
  geom_segment(data = filtered_combined_data %>% filter(Type == "methylation"),
               aes(x = pos, xend = pos, 
                   y = as.numeric(chr)-0.5, yend = as.numeric(chr)-0.5+0.33),
               color = "tomato3", alpha = 0.5, size = 2) +

  geom_segment(data = filtered_combined_data %>% filter(Type == "histone"),
               aes(x = pos, xend = pos, 
                   y = as.numeric(chr)-0.5+0.33, yend = as.numeric(chr)-0.5 +0.67),
               color = "blue", alpha = 0.5, size = 2) +

  geom_segment(data = filtered_combined_data %>% filter(Type == "anther"),
               aes(x = pos, xend = pos, 
                   y = as.numeric(chr)-0.5+0.67, yend = as.numeric(chr)+0.5),
               color = "green", alpha = 0.5, size = 2) +

  scale_y_reverse() +
  labs(x = "Genomic Position (Mbp)", y = "Chromosome") +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")
)
dev.off()





##########################################
# plotting Kimura distance for repeats

RepeatMasker -pa 2 -s -a -inv -dir ./RepMask -no_is -norna -xsmall -nolow -div 40 -lib EDTA.TElib.fa -cutoff 225 genome.fasta

calcDivergenceFromAlign.pl -s genome.divsum genome.fasta.align


#--------------------
# in R again

# install.packages("reshape")
# install.packages("ggplot2")
# install.packages("viridis")
# install.packages("hrbrthemes")
# install.packages("tidyverse")
# install.packages("gridExtra")

library(reshape)
library(ggplot2)
library(viridis)
library(hrbrthemes)
library(tidyverse)
library(gridExtra)

KimuraDistance <- read.csv("/lustre04/scratch/celphin/Oxyria/EDTA/",sep=" ")

#add here the genome size in bp
genomes_size=230000000

kd_melt = melt(KimuraDistance,id="Div")
kd_melt$norm = kd_melt$value/genomes_size * 100

ggplot(kd_melt, aes(fill=variable, y=norm, x=Div)) + 
  geom_bar(position="stack", stat="identity",color="black") +
  scale_fill_viridis(discrete = T) +
  theme_classic() +
  xlab("Kimura substitution level") +
  ylab("Percent of the genome") + 
  labs(fill = "") +
  coord_cartesian(xlim = c(0, 55)) +
  theme(axis.text=element_text(size=11),axis.title =element_text(size=12))
  