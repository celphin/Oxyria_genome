############################
# Plot gene or repeat densities
# https://genviz.org/module-02-r/0002/03/03/ggplot2_exercises/
# Jan 2025
##############################

# Fir
cd ~/scratch/
mkdir Arctic_comparative_genomes

cd /home/celphin/scratch/Arctic_comparative_genomes

DEST="/home/celphin/scratch/Arctic_comparative_genomes"
mkdir -p -- "$DEST"

files=(
"/home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/genomes/Dryas_octopetala/Dryas_octopetala.gff3"
"/home/celphin/scratch/Oxyria/EDTA/DryOcto_chr_sizes.txt"
"/home/celphin/scratch/Oxyria/EDTA/DoctH0_Main.fasta.mod.EDTA.TEanno.gff3"
"/home/celphin/scratch/Oxyria/DupGen_finder/output/7.wgd.pairs"
"/home/celphin/scratch/Oxyria/DupGen_finder/output/7.tandem.pairs"
"/home/celphin/scratch/Oxyria/DupGen_finder/output/7.proximal.pairs"
"/home/celphin/scratch/Oxyria/DupGen_finder/output/7.transposed.pairs"
"/home/celphin/scratch/Oxyria/DupGen_finder/output/7.dispersed.pairs"

"/home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/genomes/Oxyria_digyna_H1/Oxyria_digyna_H1.gff3"
"/home/celphin/scratch/Oxyria/EDTA/Oxyria_digyna_chr_sizes.txt"
"/home/celphin/scratch/Oxyria/EDTA/Oxyria_digyna.fasta.mod.EDTA.TEanno.gff3"
"/home/celphin/scratch/Oxyria/EDTA/Oxyria_Main.fasta.mod.EDTA.intact.gff3"

"/home/celphin/scratch/Oxyria/DupGen_finder/output/12.wgd.pairs"
"/home/celphin/scratch/Oxyria/DupGen_finder/output/12.tandem.pairs"
"/home/celphin/scratch/Oxyria/DupGen_finder/output/12.proximal.pairs"
"/home/celphin/scratch/Oxyria/DupGen_finder/output/12.transposed.pairs"
"/home/celphin/scratch/Oxyria/DupGen_finder/output/12.dispersed.pairs"

"/home/celphin/scratch/Oxyria/DupGen_finder/data/SequenceIDs.txt"

"/home/celphin/scratch/Oxyria/synteny_quantity/Total_interproscan_output_edited3.tsv"
)

for f in "${files[@]}"; do
  if [ -e "$f" ]; then
    cp -a -- "$f" "$DEST/" || { echo "Failed to copy $f" >&2; }
  else
    echo "Missing: $f" >&2
  fi
done

# Make a copy on Narval via globus

###################################
# Also copy over positive selection data, centromere locations, microsynteny data

cd /home/celphin/scratch/Arctic_comparative_genomes/

# Positive selection data
# Microsynteny
# Centromere locations
# Core and accessory

############################
# plotting

# Narval1
tmux new-session -s gene_plots
tmux attach-session -t gene_plots

# OLD: /lustre04/scratch/celphin/Oxyria/gene_plots
cd /home/celphin/scratch/Arctic_comparative_genomes

tree

├── Centromeres
│     ├── Draniv2_H0-AT_total_possible_range.txt
│     ├── Draniv2_H0-AT_total_possible.txt
│     ├── Dryoct_H0-AT_total_possible_range.txt
│     ├── Dryoct_H0-AT_total_possible.txt
│     ├── DryOcto_H0-AT_total_possible_range.txt
│     ├── DryOcto_H0-AT_total_possible.txt
│     ├── Oxydig_H1-AT_total_possible_range.txt
│     └── Oxydig_H1-AT_total_possible.txt
├── Core_accessory
│     ├── accessory_by_species
│     │     ├── Arabidopsis_lyrata_accessory.bed
│     │     ├── Arabidopsis_thaliana_accessory.bed
│     │     ├── Arabis_alpina_accessory.bed
│     │     ├── Argentina_anserina_accessory.bed
│     │     ├── Brassica_oleracea_accessory.bed
│     │     ├── Capsella_rubella_accessory.bed
│     │     ├── Cochlearia_groenlandica_accessory.bed
│     │     ├── Dryas_octopetala_accessory.bed
│     │     ├── Fagopyrum_escelentum_H2_accessory.bed
│     │     ├── Fagopyrum_tataricum_H1_accessory.bed
│     │     ├── Fragaria_vesca_accessory.bed
│     │     ├── Malus_sylvestris_accessory.bed
│     │     ├── Oxyria_digyna_H1_accessory.bed
│     │     ├── Polygunum_aviculare_H0_accessory.bed
│     │     ├── Prunus_persica_accessory.bed
│     │     ├── Pyrus_bretschneideri_accessory.bed
│     │     ├── Rheum_nobile_H0_accessory.bed
│     │     ├── Rheum_tangaticum_H0_accessory.bed
│     │     ├── Rosa_rugosa_accessory.bed
│     │     └── Thlaspi_arvense_accessory.bed
│     ├── core_by_species
│     │     ├── Arabidopsis_lyrata_core.bed
│     │     ├── Arabidopsis_thaliana_core.bed
│     │     ├── Arabis_alpina_core.bed
│     │     ├── Argentina_anserina_core.bed
│     │     ├── Brassica_oleracea_core.bed
│     │     ├── Capsella_rubella_core.bed
│     │     ├── Cochlearia_groenlandica_core.bed
│     │     ├── Dryas_octopetala_core.bed
│     │     ├── Fagopyrum_escelentum_H2_core.bed
│     │     ├── Fagopyrum_tataricum_H1_core.bed
│     │     ├── Fragaria_vesca_core.bed
│     │     ├── Malus_sylvestris_core.bed
│     │     ├── Oxyria_digyna_H1_core.bed
│     │     ├── Polygunum_aviculare_H0_core.bed
│     │     ├── Prunus_persica_core.bed
│     │     ├── Pyrus_bretschneideri_core.bed
│     │     ├── Rheum_nobile_H0_core.bed
│     │     ├── Rheum_tangaticum_H0_core.bed
│     │     ├── Rosa_rugosa_core.bed
│     │     └── Thlaspi_arvense_core.bed
│     └── singleton_by_species
│         ├── Arabidopsis_lyrata_singleton.bed
│         ├── Arabidopsis_thaliana_singleton.bed
│         ├── Arabis_alpina_singleton.bed
│         ├── Argentina_anserina_singleton.bed
│         ├── Brassica_oleracea_singleton.bed
│         ├── Capsella_rubella_singleton.bed
│         ├── Cochlearia_groenlandica_singleton.bed
│         ├── Dryas_octopetala_singleton.bed
│         ├── Fagopyrum_escelentum_H2_singleton.bed
│         ├── Fagopyrum_tataricum_H1_singleton.bed
│         ├── Fragaria_vesca_singleton.bed
│         ├── Malus_sylvestris_singleton.bed
│         ├── Oxyria_digyna_H1_singleton.bed
│         ├── Polygunum_aviculare_H0_singleton.bed
│         ├── Prunus_persica_singleton.bed
│         ├── Pyrus_bretschneideri_singleton.bed
│         ├── Rheum_nobile_H0_singleton.bed
│         ├── Rheum_tangaticum_H0_singleton.bed
│         ├── Rosa_rugosa_singleton.bed
│         └── Thlaspi_arvense_singleton.bed
├── Dryas_octopetala.gff3
├── DryOcto_chr_sizes.txt
├── Gene_duplications
│     ├── 12.dispersed.pairs
│     ├── 12.proximal.pairs
│     ├── 12.tandem.pairs
│     ├── 12.transposed.pairs
│     ├── 12.wgd.pairs
│     ├── 7.dispersed.pairs
│     ├── 7.proximal.pairs
│     ├── 7.tandem.pairs
│     ├── 7.transposed.pairs
│     ├── 7.wgd.pairs
│     └── SequenceIDs.txt
├── Microsynteny.csv
├── Oxyria_digyna_chr_sizes.txt
├── Oxyria_digyna_H1.gff3
├── Positive_Selection
│     ├── all_significant_genes_guidance_filtered.tsv
│     ├── all_significant_genes_guidance.tsv
│     ├── baseline_guidance_all_genes.tsv
│     ├── Dryasoct_genescores_guidance
│     ├── Dryasoct_interproscan_guidance.tsv
│     ├── Dryasoct_ORA.ermine.results
│     ├── genes_in_4species_filtered_p.txt
│     ├── genes_in_4species_OGs.txt
│     ├── genes_in_4species_with_p.txt
│     ├── Oxydig_genescores_guidance
│     ├── Oxydig_interproscan_guidance.tsv
│     └── Oxydig_ORA_guidance.ermine.results
├── Repeats
│     ├── DoctH0_Main.fasta.mod.EDTA.TEanno.gff3
│     ├── Oxyria_digyna.fasta.mod.EDTA.TEanno.gff3
│     └── Oxyria_Main.fasta.mod.EDTA.intact.gff3
└── Total_interproscan_output_edited3.tsv


#------------------------------
# get chromosome lengths

# remove sequences that are too short
module load StdEnv/2023 seqkit/2.5.1
seqkit seq -m 10000000 ./Dryas/Dry-octo-H2_DoctH0_Main.fasta > DryOcto_chr.fasta
seqkit seq -m 10000000 Oxyria_Main.fasta > Oxyria_digyna_chr.fasta

module load  StdEnv/2020 bioawk/1.0
bioawk -c fastx '{print $name "\t" length($seq)}' DryOcto_chr.fasta  > DryOcto_chr_sizes.txt
bioawk -c fastx '{print $name "\t" length($seq)}' Oxyria_digyna_chr.fasta > Oxyria_digyna_chr_sizes.txt

#-----------------------
# Joining data and plotting in R

module load StdEnv/2023
module load r/4.4.0

R

library(dplyr)
library(ggplot2) 
library(tidyverse)
library(rlang)
library(tidyr)
library(stringr)
library(tibble)

#library(statebins)

# import a text file with gene positions
# Dryas
Dry_genes0 <- read.table("Dryas_octopetala.gff3",sep="\t",header=F)
Dry_chr_sizes0 <- read.table("DryOcto_chr_sizes.txt",sep="\t",header=F)
Dry_TE_repeats0 <- read.table("./Repeats/DoctH0_Main.fasta.mod.EDTA.TEanno.gff3",sep="\t",header=F)

# Oxyria
Oxy_genes0 <- read.table("Oxyria_digyna_H1.gff3",sep="\t",header=F)
Oxy_chr_sizes0 <- read.table("Oxyria_digyna_chr_sizes.txt",sep="\t",header=F)
Oxy_TE_repeats1 <- read.table("./Repeats/Oxyria_digyna.fasta.mod.EDTA.TEanno.gff3",sep="\t",header=F)
Oxy_TE_repeats0 <- read.table("./Repeats/Oxyria_Main.fasta.mod.EDTA.intact.gff3",sep="\t",header=F)

# Gene duplications
Dry_wgd0 <- read.table("./Gene_duplications/7.wgd.pairs",sep="\t",header=T)
Dry_tandem0 <- read.table("./Gene_duplications/7.tandem.pairs",sep="\t",header=T)
Dry_proximal0 <- read.table("./Gene_duplications/7.proximal.pairs",sep="\t",header=T)
Dry_transposed0 <- read.table("./Gene_duplications/7.transposed.pairs",sep="\t",header=T)
Dry_dispersed0 <- read.table("./Gene_duplications/7.dispersed.pairs",sep="\t",header=T)

Oxy_wgd0 <- read.table("./Gene_duplications/12.wgd.pairs",sep="\t",header=T)
Oxy_tandem0 <- read.table("./Gene_duplications/12.tandem.pairs",sep="\t",header=T)
Oxy_proximal0 <- read.table("./Gene_duplications/12.proximal.pairs",sep="\t",header=T)
Oxy_transposed0 <- read.table("./Gene_duplications/12.transposed.pairs",sep="\t",header=T)
Oxy_dispersed0 <- read.table("./Gene_duplications/12.dispersed.pairs",sep="\t",header=T)

# DupGen seq IDs
SequenceIDs <- read.table("./Gene_duplications/SequenceIDs.txt",sep=":",header=F)

# Interproscan data for all species
gene_ont <- read.delim("Total_interproscan_output_edited3.tsv", header = TRUE, sep = "\t", na.strings = "-", colClasses = c("character", "character", "character", "character"))

# Positive selection data
pos_sel <- read.table("./Positive_Selection/all_significant_genes_guidance_filtered.tsv",sep="\t",header=T)
pos_sel_4spp <- read.table("./Positive_Selection/genes_in_4species_filtered_p.txt",sep="\t",header=T)

# Microsynteny
microsyn <- read.table("Microsynteny.csv",sep=",",header=T)

# Centromere locations
Dryas_cent_pos <- read.table("./Centromeres/DryOcto_H0-AT_total_possible_range.txt",sep="\t",header=T)
Oxyria_cent_pos <- read.table("./Centromeres/Oxydig_H1-AT_total_possible_range.txt",sep="\t",header=T)

# Core and accessory
Dryas_core <- read.table("./Core_accessory/core_by_species/Dryas_octopetala_core.bed",sep="\t",header=T)
Oxyria_core <- read.table("./Core_accessory/core_by_species/Oxyria_digyna_H1_core.bed",sep="\t",header=T)

Dryas_accessory <- read.table("./Core_accessory/accessory_by_species/Dryas_octopetala_accessory.bed",sep="\t",header=T)
Oxyria_accessory <- read.table("./Core_accessory/accessory_by_species/Oxyria_digyna_H1_accessory.bed",sep="\t",header=T)

Dryas_singleton <- read.table("./Core_accessory/singleton_by_species/Dryas_octopetala_singleton.bed",sep="\t",header=T)
Oxyria_singleton <- read.table("./Core_accessory/singleton_by_species/Oxyria_digyna_H1_singleton.bed",sep="\t",header=T)

#############################
# Assign species to run plots for

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
# Plot normalized genes

make_gene_density_norm_1Mb <- function(gff_tbl, species_name, binwidth = 1e6,
                                       feature_col="V3", feature_value="gene",
                                       chrcol="V1", poscol="V4",
                                       threshold = 1e7) {
  df <- as.data.frame(gff_tbl)

  df <- df %>% filter(.data[[feature_col]] == feature_value)

  scaffold_lengths <- df %>%
    group_by(.data[[chrcol]]) %>%
    summarise(length = max(.data[["V5"]]), .groups="drop")

  long_scaffolds <- scaffold_lengths %>%
    filter(length > threshold) %>%
    pull(.data[[chrcol]])

  df <- df %>% filter(.data[[chrcol]] %in% long_scaffolds)

  binned <- df %>%
    transmute(chr = .data[[chrcol]],
              pos = as.numeric(.data[[poscol]])) %>%
    mutate(bin_start = floor(pos / binwidth) * binwidth) %>%
    group_by(chr, bin_start) %>%
    summarise(count = n(), .groups="drop") %>%
    group_by(chr) %>%
    mutate(norm = if (max(count) == min(count)) 0
           else (count - min(count)) / (max(count) - min(count))) %>%
    ungroup() %>%
    mutate(dataset = species_name)

  binned
}

plot_one_species <- function(gff_tbl, species_name) {
  dat <- make_gene_density_norm_1Mb(gff_tbl, species_name)

  p <- ggplot(dat, aes(x = bin_start / 1e6, y = norm)) +
    geom_line(linewidth = 0.6, color = "steelblue") +
    facet_wrap(~ chr, ncol = 1, scales = "free_y") +
    xlab("Genomic position (Mb, 1 Mb bins)") +
    ylab("Normalized gene density (0–1)") +
    theme_classic()

  out_base <- paste0("./plots/", species_name, "_gene_density_normalized_0_1_1Mb")
  ggsave(paste0(out_base, ".png"), p, width = 7, height = 12, dpi = 300)
  ggsave(paste0(out_base, ".pdf"), p, width = 7, height = 12)
}

plot_one_species(Dry_genes0, "Dryas_octopetala_H0")
plot_one_species(Oxy_genes0, "Oxyria_digyna")

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
# Normalize from 0-1 again

make_repeat_bp_density_norm_1Mb_by_type <- function(Spp_TE_repeats0, Spp, binwidth = 1e6, threshold = 1e7,
                                                        normalize_bins_to_0_1 = TRUE) {

  scaffold_lengths <- Spp_TE_repeats0 %>%
    as.data.frame() %>%
    group_by(V1) %>%
    summarise(length = max(V5), .groups="drop")

  long_scaffolds <- scaffold_lengths %>%
    filter(length > threshold) %>%
    pull(V1)

  unique_repeat_types <- unique(Spp_TE_repeats0$V3)

  for (repeat_type in unique_repeat_types) {

    filtered_repeats <- Spp_TE_repeats0 %>%
      filter(V1 %in% long_scaffolds, V3 == repeat_type)

    if (nrow(filtered_repeats) == 0) {
      message(paste("No data on more than 1 chromosome for repeat type:", repeat_type))
      next
    }

    df <- as.data.frame(filtered_repeats)

    # start/end
    df <- df %>%
      transmute(
        chr = V1,
        start = as.numeric(V4),
        end   = as.numeric(V5),
        rep_bp = (as.numeric(V5) - as.numeric(V4) + 1),
        pos = as.numeric(V4)
      ) %>%
      filter(is.finite(rep_bp), rep_bp >= 0)

    # Put each repeat into the bin of its START (common/simple approach)
    df <- df %>%
      mutate(bin_start = floor(pos / binwidth) * binwidth)

    binned <- df %>%
      group_by(chr, bin_start) %>%
      summarise(total_rep_bp = sum(rep_bp, na.rm = TRUE), .groups="drop") %>%
      mutate(rep_bp_per_1Mb = total_rep_bp / binwidth)

    if (normalize_bins_to_0_1) {
      mn <- min(binned$rep_bp_per_1Mb, na.rm = TRUE)
      mx <- max(binned$rep_bp_per_1Mb, na.rm = TRUE)
      binned <- binned %>%
        mutate(norm = if (mx == mn) 0 else (rep_bp_per_1Mb - mn) / (mx - mn))
    }

    # Plot
    p <- ggplot(binned, aes(x = bin_start / 1e6)) +
      geom_line(aes(y = if (normalize_bins_to_0_1) norm else rep_bp_per_1Mb),
                linewidth = 0.6) +
      facet_wrap(~ chr, ncol = 1) +
      theme_classic() +
      xlab("Genomic position (Mb, 1 Mb windows)") +
      ylab(if (normalize_bins_to_0_1) "Normalized repeat abundance (0–1)" else "Repeat abundance (bp/1Mb)") +
      ggtitle(paste(Spp, "-", repeat_type))

    filename_png <- paste0("./plots/", Spp, "_repeats_", gsub(" ", "_", repeat_type), "_bp_per_1Mb.png")
    filename_pdf <- paste0("./plots/", Spp, "_repeats_", gsub(" ", "_", repeat_type), "_bp_per_1Mb.pdf")

    ggsave(filename_png, p, width = 7, height = 12, dpi = 300)
    ggsave(filename_pdf, p, width = 7, height = 12)
  }
}

# Example: run for each species separately (Oxyria then Dryas)
make_repeat_bp_density_norm_1Mb_by_type(Spp_TE_repeats0 = Oxy_TE_repeats0, Spp = "Oxyria_digyna", normalize_bins_to_0_1 = TRUE)
make_repeat_bp_density_norm_1Mb_by_type(Spp_TE_repeats0 = Dry_TE_repeats0, Spp = "Dryas_octopetala_H0", normalize_bins_to_0_1 = TRUE)



#----------------------------
# Gene duplicates
#--------------------------
# Normalize gene duplicate counts from 0-1

dup_bins_norm_1Mb <- function(dup_pairs_tbl, dataset_name, binwidth = 1e6, threshold = 1e7) {
  # dup_pairs_tbl has columns Location and Location.1 (as in your code)
  location_data <- dup_pairs_tbl %>%
    select(Location, Location.1) %>%
    pivot_longer(cols = everything(), names_to = "Type", values_to = "Location") %>%
    separate(Location, into = c("chr", "pos"), sep = ":", convert = TRUE)

  scaffold_lengths <- location_data %>%
    as.data.frame() %>%
    group_by(chr) %>%
    summarise(length = max(pos), .groups = "drop")

  long_scaffolds <- scaffold_lengths %>%
    filter(length > threshold) %>%
    pull(chr)

  location_data <- location_data %>% filter(chr %in% long_scaffolds)
  location_data$pos <- as.numeric(location_data$pos)

  binned <- location_data %>%
    mutate(bin_start = floor(pos / binwidth) * binwidth) %>%
    group_by(chr, bin_start) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(chr) %>%
    mutate(norm = if (max(count) == min(count)) 0
                   else (count - min(count)) / (max(count) - min(count))) %>%
    ungroup() %>%
    mutate(dataset = dataset_name)

  binned
}

# ---- Run for each species / duplication class, collect for later combined plotting ----
dup_collect <- function(Spp, wgddata, tanddata, proxdata, transdata, dispdata) {
  bind_rows(
    dup_bins_norm_1Mb(wgddata,  paste0(Spp), 1e6) %>% mutate(track = "WGD"),
    dup_bins_norm_1Mb(tanddata, paste0(Spp), 1e6) %>% mutate(track = "Tandem"),
    dup_bins_norm_1Mb(proxdata, paste0(Spp), 1e6) %>% mutate(track = "Proximal"),
    dup_bins_norm_1Mb(transdata,paste0(Spp), 1e6) %>% mutate(track = "Transposed"),
    dup_bins_norm_1Mb(dispdata, paste0(Spp), 1e6) %>% mutate(track = "Dispersed")
  )
}

dup_oxy <- dup_collect("Oxyria_digyna", Oxy_wgd0, Oxy_tandem0, Oxy_proximal0, Oxy_transposed0, Oxy_dispersed0)
dup_dry <- dup_collect("Dryas_octopetala_H0", Dry_wgd0, Dry_tandem0, Dry_proximal0, Dry_transposed0, Dry_dispersed0)

dup_all <- bind_rows(dup_oxy, dup_dry)

#------------
# Plot gene duplications
plot_dup_species <- function(dup_all, species_name, out_png, out_pdf) {
  dat <- dup_all %>% filter(dataset == species_name)

  p <- ggplot(dat, aes(x = bin_start/1e6, y = norm, color = track, linetype = track)) +
    geom_line(linewidth = 0.6) +
    facet_wrap(~ chr, ncol = 1) +
    xlab("Genomic position (Mb, 1 Mb bins)") +
    ylab("Normalized duplication location density (0–1)") +
    theme_classic() +
    ggtitle(species_name)

  ggsave(out_png, p, width = 7, height = 12, dpi = 300)
  ggsave(out_pdf, p, width = 7, height = 12)
}

# dup_all from the previous step
plot_dup_species(
  dup_all, "Oxyria_digyna",
  "./plots/Oxyria_digyna_duplications_normalized_0_1.png",
  "./plots/Oxyria_digyna_duplications_normalized_0_1.pdf"
)

plot_dup_species(
  dup_all, "Dryas_octopetala_H0",
  "./plots/Dryas_octopetala_H0_duplications_normalized_0_1.png",
  "./plots/Dryas_octopetala_H0_duplications_normalized_0_1.pdf"
)

#-------------------
# Add plot of total gene duplications

dup_total_bins_norm_1Mb <- function(dup_pairs_tbl, dataset_name, binwidth = 1e6, threshold = 1e7) {
  location_data <- dup_pairs_tbl %>%
    select(Location, Location.1) %>%
    pivot_longer(cols = everything(), names_to = "Type", values_to = "Location") %>%
    separate(Location, into = c("chr", "pos"), sep = ":", convert = TRUE)

  scaffold_lengths <- location_data %>%
    as.data.frame() %>%
    group_by(chr) %>%
    summarise(length = max(pos), .groups = "drop")

  long_scaffolds <- scaffold_lengths %>%
    filter(length > threshold) %>%
    pull(chr)

  location_data <- location_data %>% filter(chr %in% long_scaffolds)
  location_data$pos <- as.numeric(location_data$pos)

  binned <- location_data %>%
    mutate(bin_start = floor(pos / binwidth) * binwidth) %>%
    group_by(chr, bin_start) %>%
    summarise(count = n(), .groups="drop") %>%
    group_by(chr) %>%
    mutate(norm = if (max(count) == min(count)) 0
                   else (count - min(count)) / (max(count) - min(count))) %>%
    ungroup() %>%
    mutate(dataset = dataset_name)

  binned
}

# total duplication counts per 1Mb (sum of all duplication types)
# If you want “regardless of type” across WGD/Tandem/Proximal/Transposed/Dispersed,
# bind them first, then bin once.

dup_total_alltypes <- function(Spp, wgddata, tanddata, proxdata, transdata, dispdata) {
  bind_rows(
    wgddata   %>% mutate(track = "WGD"),
    tanddata  %>% mutate(track = "Tandem"),
    proxdata  %>% mutate(track = "Proximal"),
    transdata %>% mutate(track = "Transposed"),
    dispdata  %>% mutate(track = "Dispersed")
  ) %>%
    select(-track) %>%
    dup_total_bins_norm_1Mb(., dataset_name = Spp)
}

dup_total_oxy <- dup_total_alltypes("Oxyria_digyna", Oxy_wgd0, Oxy_tandem0, Oxy_proximal0, Oxy_transposed0, Oxy_dispersed0)
dup_total_dry <- dup_total_alltypes("Dryas_octopetala_H0", Dry_wgd0, Dry_tandem0, Dry_proximal0, Dry_transposed0, Dry_dispersed0)

plot_dup_total_species <- function(dup_total_df, out_png, out_pdf) {
  p <- ggplot(dup_total_df, aes(x = bin_start/1e6, y = norm)) +
    geom_line(linewidth = 0.7, color = "black") +
    facet_wrap(~ chr, ncol = 1) +
    xlab("Genomic position (Mb, 1 Mb bins)") +
    ylab("Normalized total duplication locations (0–1)") +
    theme_classic() +
    ggtitle(unique(dup_total_df$dataset))

  ggsave(out_png, p, width = 7, height = 12, dpi = 300)
  ggsave(out_pdf, p, width = 7, height = 12)
}

plot_dup_total_species(
  dup_total_oxy,
  "./plots/Oxyria_digyna_total_gene_duplications_normalized_0_1.png",
  "./plots/Oxyria_digyna_total_gene_duplications_normalized_0_1.pdf"
)

plot_dup_total_species(
  dup_total_dry,
  "./plots/Dryas_octopetala_H0_total_gene_duplications_normalized_0_1.png",
  "./plots/Dryas_octopetala_H0_total_gene_duplications_normalized_0_1.pdf"
)


#----------------------------
# Positively selected genes
#----------------------------
pos_sel <- read.table("./Positive_Selection/all_significant_genes_guidance_filtered.tsv",
                      sep="\t", header=F)
pos_sel_4spp <- read.table("./Positive_Selection/genes_in_4species_filtered_p.txt",
                           sep=" ", header=F)

colnames(pos_sel) <- c("OG", "Spp", "p-value")
colnames(pos_sel_4spp) <- c("OG", "Spp", "p-value")

#---- helper: bin gene positions to 1Mb and count selected genes per bin ----
count_selected_genes_1Mb <- function(gff_tbl, gene_names, dataset_name,
                                       binwidth=1e6, threshold=1e7) {
  gff <- as.data.frame(gff_tbl)

  gff_genes <- gff %>%
    dplyr::filter(V3 == "gene")

  # case-insensitive matching
  gene_names_up <- toupper(gene_names)

  gene_keep <- gff_genes %>%
    dplyr::mutate(
      gene_id = sub(".*ID=([^;]+).*", "\\1", V9),
      gene_id_up = toupper(gene_id)
    ) %>%
    dplyr::filter(gene_id_up %in% gene_names_up)

  scaffold_lengths <- gene_keep %>%
    dplyr::group_by(V1) %>%
    dplyr::summarise(length = max(V5), .groups="drop")

  long_scaffolds <- scaffold_lengths %>%
    dplyr::filter(length > threshold) %>%
    dplyr::pull(V1)

  gene_keep <- gene_keep %>%
    dplyr::filter(V1 %in% long_scaffolds)

  # ---- total gene counts per bin (for the same long scaffolds) ----
  total_keep <- gff_genes %>%
    dplyr::filter(V1 %in% long_scaffolds) %>%
    dplyr::mutate(
      pos = as.numeric(V4),
      bin_start = floor(pos / binwidth) * binwidth
    ) %>%
    dplyr::group_by(V1, bin_start) %>%
    dplyr::summarise(total_count = dplyr::n(), .groups="drop") %>%
    dplyr::rename(chr = V1)

  # ---- selected gene counts per bin ----
  selected_binned <- gene_keep %>%
    dplyr::mutate(
      pos = as.numeric(V4),
      bin_start = floor(pos / binwidth) * binwidth
    ) %>%
    dplyr::group_by(V1, bin_start) %>%
    dplyr::summarise(selected_count = dplyr::n(), .groups="drop") %>%
    dplyr::rename(chr = V1)

  binned <- selected_binned %>%
    dplyr::full_join(total_keep, by = c("chr", "bin_start")) %>%
    dplyr::mutate(
      selected_count = dplyr::coalesce(selected_count, 0L),
      total_count = dplyr::coalesce(total_count, 0L),
      ratio = dplyr::if_else(total_count > 0, selected_count / total_count, NA_real_)
    ) %>%
    dplyr::group_by(chr) %>%
    dplyr::mutate(
      # keep your original selected-count normalization
      norm = if (max(selected_count) == min(selected_count)) 0
             else (selected_count - min(selected_count)) / (max(selected_count) - min(selected_count)),

      # normalize the ratio to 0–1 for plotting (optional but matches your y scale)
      ratio_norm = if (all(is.na(ratio)) || max(ratio, na.rm=TRUE) == min(ratio, na.rm=TRUE)) 0
                    else (ratio - min(ratio, na.rm=TRUE)) / (max(ratio, na.rm=TRUE) - min(ratio, na.rm=TRUE))
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(dataset = dataset_name)

  binned
}

#---- run: make selected sets from pos_sel ----

dry_genes <- pos_sel %>%
  dplyr::filter(stringr::str_detect(Spp, "^DOCTH0_")) %>%
  dplyr::pull(Spp) %>%
  unique()

oxy_genes <- pos_sel %>%
  dplyr::filter(stringr::str_detect(Spp, "^OXYRIA")) %>%
  dplyr::pull(Spp) %>%
  unique()

# NOTE: these objects must exist in your environment:
#   Dry_genes0, Oxy_genes0
# (They should be your GFF tables for each species, with V1=chr, V3="gene", V4=start, V5=end, V9=attributes.)

Dry_sel_bins <- count_selected_genes_1Mb(
  gff_tbl = Dry_genes0,
  gene_names = dry_genes,
  dataset_name = "Dryas_octopetala_H0"
)

Oxy_sel_bins <- count_selected_genes_1Mb(
  gff_tbl = Oxy_genes0,
  gene_names = oxy_genes,
  dataset_name = "Oxyria_digyna"
)

#---- plotting ----
plot_sel_species <- function(dat, out_png, out_pdf) {
  dat_long <- tidyr::pivot_longer(
    dat,
    cols = c(norm, ratio_norm),
    names_to = "metric",
    values_to = "value"
  ) %>%
    dplyr::mutate(metric = dplyr::recode(
      metric,
      norm = "Selected count (min–max)",
      ratio_norm = "Selected fraction (min–max)"
    ))

  p <- ggplot(dat_long, aes(x = bin_start/1e6, y = value, color = metric)) +
    geom_line(linewidth = 0.7) +
    facet_wrap(~ chr, scales = "fixed", ncol = 1) +
    xlab("Genomic position (Mb, 1 Mb bins)") +
    ylab("Normalized value (0–1)") +
    theme_classic() +
    ggtitle(unique(dat$dataset)) +
    theme(legend.position = "bottom")

  ggsave(out_png, p, width = 7, height = 12, dpi = 300)
  ggsave(out_pdf, p, width = 7, height = 12)
}


plot_sel_species(
  Dry_sel_bins,
  "./plots/Dryas_octopetala_H0_positive_selection_genes_count_1Mb_norm.png",
  "./plots/Dryas_octopetala_H0_positive_selection_genes_count_1Mb_norm.pdf"
)

plot_sel_species(
  Oxy_sel_bins,
  "./plots/Oxyria_digyna_positive_selection_genes_count_1Mb_norm.png",
  "./plots/Oxyria_digyna_positive_selection_genes_count_1Mb_norm.pdf"
)



#----------------------------
# Core and accessory genes
#--------------------------

count_bed_1Mb_with_fraction_from_gff <- function(bed_tbl,
                                                 gff_tbl_all_genes,
                                                 dataset_name,
                                                 binwidth = 1e6,
                                                 min_len = 5e6,
                                                 use = c("midpoint","start","end")) {
  use <- match.arg(use)

  b <- as.data.frame(bed_tbl)
  gff <- as.data.frame(gff_tbl_all_genes)

  gff_genes <- gff %>%
    dplyr::filter(V3 == "gene") %>%
    dplyr::mutate(
      chr = as.character(V1),
      start = as.numeric(V4),
      end = as.numeric(V5)
    )

  # BED coords
  b_start <- as.numeric(b$V2)
  b_end   <- as.numeric(b$V3)
  b_chr   <- as.character(b$V1)

  # keep chromosomes based on BED (same logic as your current BED function)
  chr_len <- data.frame(chr = b_chr, start = b_start, end = b_end) %>%
    dplyr::group_by(chr) %>%
    dplyr::summarise(length = max(end), .groups="drop")

  keep_chr <- chr_len %>% dplyr::filter(length >= min_len) %>% dplyr::pull(chr)

  keep_cat <- b_chr %in% keep_chr
  b_start <- b_start[keep_cat]; b_end <- b_end[keep_cat]; b_chr <- b_chr[keep_cat]

  # denominator gene coords restricted to same long scaffolds
  gff_genes <- gff_genes %>% dplyr::filter(chr %in% keep_chr)

  # positions
  cat_pos <- switch(use,
    midpoint = (b_start + b_end) / 2,
    start    = b_start,
    end      = b_end
  )

  all_pos <- switch(use,
    midpoint = (gff_genes$start + gff_genes$end) / 2,
    start    = gff_genes$start,
    end      = gff_genes$end
  )

  # binning
  cat_df <- tibble(chr = b_chr, pos = cat_pos) %>%
    mutate(bin_start = floor(pos / binwidth) * binwidth) %>%
    group_by(chr, bin_start) %>%
    summarise(selected_count = dplyr::n(), .groups="drop")

  all_df <- tibble(chr = gff_genes$chr, pos = all_pos) %>%
    mutate(bin_start = floor(pos / binwidth) * binwidth) %>%
    group_by(chr, bin_start) %>%
    summarise(total_count = dplyr::n(), .groups="drop")

  out <- cat_df %>%
    dplyr::full_join(all_df, by = c("chr","bin_start")) %>%
    dplyr::mutate(
      selected_count = dplyr::coalesce(selected_count, 0L),
      total_count = dplyr::coalesce(total_count, 0L),
      ratio = dplyr::if_else(total_count > 0, selected_count / total_count, NA_real_)
    ) %>%
    dplyr::group_by(chr) %>%
    dplyr::mutate(
      norm = if (max(selected_count) == min(selected_count)) 0
              else (selected_count - min(selected_count)) / (max(selected_count) - min(selected_count)),

      ratio_norm = if (all(is.na(ratio)) || max(ratio, na.rm=TRUE) == min(ratio, na.rm=TRUE)) 0
                    else (ratio - min(ratio, na.rm=TRUE)) / (max(ratio, na.rm=TRUE) - min(ratio, na.rm=TRUE))
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(dataset = dataset_name)

  out
}


plot_bed_species <- function(dat, out_png, out_pdf) {
  dat_long <- tidyr::pivot_longer(
    dat,
    cols = c(norm, ratio_norm),
    names_to = "metric",
    values_to = "value"
  ) %>%
    dplyr::mutate(metric = dplyr::recode(
      metric,
      norm = "Category count (min–max)",
      ratio_norm = "Category fraction (min–max)"
    ))

  p <- ggplot(dat_long, aes(x = bin_start/1e6, y = value, color = metric)) +
    geom_line(linewidth = 0.7) +
    facet_wrap(~ chr, scales = "fixed", ncol = 1) +
    xlab("Genomic position (Mb, 1 Mb bins)") +
    ylab("Normalized value (0–1)") +
    theme_classic() +
    ggtitle(unique(dat$dataset)) +
    theme(legend.position = "bottom")

  ggsave(out_png, p, width = 7, height = 12, dpi = 300)
  ggsave(out_pdf, p, width = 7, height = 12)
}



#---- read your BEDs ----
Dryas_core <- read.table("./Core_accessory/core_by_species/Dryas_octopetala_core.bed",sep="\t",header=F)
Oxyria_core <- read.table("./Core_accessory/core_by_species/Oxyria_digyna_H1_core.bed",sep="\t",header=F)

Dryas_accessory <- read.table("./Core_accessory/accessory_by_species/Dryas_octopetala_accessory.bed",sep="\t",header=F)
Oxyria_accessory <- read.table("./Core_accessory/accessory_by_species/Oxyria_digyna_H1_accessory.bed",sep="\t",header=F)

Dryas_singleton <- read.table("./Core_accessory/singleton_by_species/Dryas_octopetala_singleton.bed",sep="\t",header=F)
Oxyria_singleton <- read.table("./Core_accessory/singleton_by_species/Oxyria_digyna_H1_singleton.bed",sep="\t",header=F)

#---- bin each category (midpoint) ----
# Core
Dry_core_bins <- count_bed_1Mb_with_fraction_from_gff(
  Dryas_core, Dry_genes0, "Dryas_octopetala_core",
  binwidth = 1e6, min_len = 5e6, use = "midpoint"
)

Oxy_core_bins <- count_bed_1Mb_with_fraction_from_gff(
  Oxyria_core, Oxy_genes0, "Oxyria_digyna_H1_core",
  binwidth = 1e6, min_len = 5e6, use = "midpoint"
)

# Accessory
Dry_acc_bins <- count_bed_1Mb_with_fraction_from_gff(
  Dryas_accessory, Dry_genes0, "Dryas_octopetala_accessory",
  binwidth = 1e6, min_len = 5e6, use = "midpoint"
)

Oxy_acc_bins <- count_bed_1Mb_with_fraction_from_gff(
  Oxyria_accessory, Oxy_genes0, "Oxyria_digyna_H1_accessory",
  binwidth = 1e6, min_len = 5e6, use = "midpoint"
)

# Singleton
Dry_single_bins <- count_bed_1Mb_with_fraction_from_gff(
  Dryas_singleton, Dry_genes0, "Dryas_octopetala_singleton",
  binwidth = 1e6, min_len = 5e6, use = "midpoint"
)

Oxy_single_bins <- count_bed_1Mb_with_fraction_from_gff(
  Oxyria_singleton, Oxy_genes0, "Oxyria_digyna_H1_singleton",
  binwidth = 1e6, min_len = 5e6, use = "midpoint"
)

# Plot
plot_bed_species(Dry_core_bins,
  "./plots/Dryas_octopetala_core_count_and_fraction_1Mb.png",
  "./plots/Dryas_octopetala_core_count_and_fraction_1Mb.pdf"
)

plot_bed_species(Oxy_core_bins,
  "./plots/Oxyria_digyna_H1_core_count_and_fraction_1Mb.png",
  "./plots/Oxyria_digyna_H1_core_count_and_fraction_1Mb.pdf"
)

plot_bed_species(Dry_acc_bins,
  "./plots/Dryas_octopetala_accessory_count_and_fraction_1Mb.png",
  "./plots/Dryas_octopetala_accessory_count_and_fraction_1Mb.pdf"
)

plot_bed_species(Oxy_acc_bins,
  "./plots/Oxyria_digyna_H1_accessory_count_and_fraction_1Mb.png",
  "./plots/Oxyria_digyna_H1_accessory_count_and_fraction_1Mb.pdf"
)

plot_bed_species(Dry_single_bins,
  "./plots/Dryas_octopetala_singleton_count_and_fraction_1Mb.png",
  "./plots/Dryas_octopetala_singleton_count_and_fraction_1Mb.pdf"
)

plot_bed_species(Oxy_single_bins,
  "./plots/Oxyria_digyna_H1_singleton_count_and_fraction_1Mb.png",
  "./plots/Oxyria_digyna_H1_singleton_count_and_fraction_1Mb.pdf"
)



#----------------------------
# Microsynteny genes (count + fraction)
#----------------------------

#---- parse GFF: gene id -> chr + start/end (case-insensitive) ----
gff_gene_coords <- function(gff_tbl) {
  gff <- as.data.frame(gff_tbl)
  gff_genes <- gff %>% dplyr::filter(V3 == "gene")

  tibble(
    gene_id = sub(".*ID=([^;]+).*", "\\1", gff_genes$V9),
    chr = as.character(gff_genes$V1),
    start = as.numeric(gff_genes$V4),
    end = as.numeric(gff_genes$V5)
  ) %>%
    mutate(gene_id_up = toupper(gene_id))
}

Dry_gff <- gff_gene_coords(Dry_genes0)
Oxy_gff <- gff_gene_coords(Oxy_genes0)

#---- bin selected vs total per chr, and compute two metrics ----
count_microsyn_1Mb_with_fraction <- function(sel_chr, sel_pos,
                                              all_chr, all_pos,
                                              dataset_name,
                                              binwidth = 1e6,
                                              min_len = 5e6) {

  sel_chr <- as.character(sel_chr); sel_pos <- as.numeric(sel_pos)
  all_chr <- as.character(all_chr); all_pos <- as.numeric(all_pos)

  chr_len <- tibble(chr = all_chr, pos = all_pos) %>%
    dplyr::group_by(chr) %>%
    dplyr::summarise(length = max(pos, na.rm = TRUE), .groups = "drop")

  keep_chr <- chr_len %>% dplyr::filter(length >= min_len) %>% dplyr::pull(chr)

  sel_keep <- sel_chr %in% keep_chr
  sel_chr <- sel_chr[sel_keep]
  sel_pos <- sel_pos[sel_keep]

  all_keep <- all_chr %in% keep_chr
  all_chr <- all_chr[all_keep]
  all_pos <- all_pos[all_keep]

  cat_df <- tibble(chr = sel_chr, pos = sel_pos) %>%
    dplyr::filter(!is.na(pos)) %>%
    mutate(bin_start = floor(pos / binwidth) * binwidth) %>%
    dplyr::group_by(chr, bin_start) %>%
    dplyr::summarise(selected_count = dplyr::n(), .groups = "drop")

  all_df <- tibble(chr = all_chr, pos = all_pos) %>%
    dplyr::filter(!is.na(pos)) %>%
    mutate(bin_start = floor(pos / binwidth) * binwidth) %>%
    dplyr::group_by(chr, bin_start) %>%
    dplyr::summarise(total_count = dplyr::n(), .groups = "drop")

  out <- cat_df %>%
    full_join(all_df, by = c("chr", "bin_start")) %>%
    mutate(
      selected_count = dplyr::coalesce(selected_count, 0L),
      total_count = dplyr::coalesce(total_count, 0L),
      ratio = dplyr::if_else(total_count > 0, selected_count / total_count, NA_real_)
    ) %>%
    dplyr::group_by(chr) %>%
    dplyr::mutate(
      norm = if (max(selected_count, na.rm = TRUE) == min(selected_count, na.rm = TRUE)) 0
              else (selected_count - min(selected_count, na.rm = TRUE)) /
                   (max(selected_count, na.rm = TRUE) - min(selected_count, na.rm = TRUE)),

      ratio_norm = if (all(is.na(ratio)) ||
                         max(ratio, na.rm = TRUE) == min(ratio, na.rm = TRUE)) 0
                    else (ratio - min(ratio, na.rm = TRUE)) /
                         (max(ratio, na.rm = TRUE) - min(ratio, na.rm = TRUE))
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(dataset = dataset_name)

  out
}


#---- plotting: overlay both metrics within each chromosome ----
plot_microsyn_species <- function(dat, out_png, out_pdf) {
  dat_long <- tidyr::pivot_longer(
    dat,
    cols = c(norm, ratio_norm),
    names_to = "metric",
    values_to = "value"
  ) %>%
    dplyr::mutate(metric = dplyr::recode(
      metric,
      norm = "Microsynteny count (min–max)",
      ratio_norm = "Microsynteny fraction (min–max)"
    ))

  p <- ggplot(dat_long, aes(x = bin_start/1e6, y = value, color = metric)) +
    geom_line(linewidth = 0.7) +
    facet_wrap(~ chr, scales = "fixed", ncol = 1) +
    xlab("Genomic position (Mb, 1 Mb bins)") +
    ylab("Normalized value (0–1)") +
    theme_classic() +
    ggtitle(dat$dataset[1]) +
    theme(legend.position = "bottom")

  ggsave(out_png, p, width = 7, height = 12, dpi = 300)
  ggsave(out_pdf, p, width = 7, height = 12)
}

#----------------------------
# Microsynteny table -> map to GFF -> positions
#----------------------------

# Example columns from your snippet
Dry_gene_col <- "Dryas.octopetala.Gene.1"
Oxy_gene_col <- "Oxyria.digyna.Gene.1"

# microsynteny gene midpoints
Dry_mic_pos <- microsyn %>%
  dplyr::transmute(gene_id = .data[[Dry_gene_col]],
                   gene_id_up = toupper(gene_id)) %>%
  dplyr::inner_join(Dry_gff %>% dplyr::select(gene_id_up, chr, start, end),
                    by = "gene_id_up") %>%
  dplyr::mutate(mid = (start + end)/2) %>%
  dplyr::select(chr, mid)

Oxy_mic_pos <- microsyn %>%
  dplyr::transmute(gene_id = .data[[Oxy_gene_col]],
                   gene_id_up = toupper(gene_id)) %>%
  dplyr::inner_join(Oxy_gff %>% dplyr::select(gene_id_up, chr, start, end),
                    by = "gene_id_up") %>%
  dplyr::mutate(mid = (start + end)/2) %>%
  dplyr::select(chr, mid)

# all gene midpoints (denominator)
Dry_all_pos <- Dry_gff %>%
  dplyr::mutate(mid = (start + end)/2) %>%
  dplyr::select(chr, mid)

Oxy_all_pos <- Oxy_gff %>%
  dplyr::mutate(mid = (start + end)/2) %>%
  dplyr::select(chr, mid)

# bin + metrics
Dry_mic_bins <- count_microsyn_1Mb_with_fraction(
  sel_chr = Dry_mic_pos$chr,
  sel_pos = Dry_mic_pos$mid,
  all_chr = Dry_all_pos$chr,
  all_pos = Dry_all_pos$mid,
  dataset_name = "Dryas_octopetala_microsynteny (Gene.1)",
  binwidth = 1e6,
  min_len = 5e6
)

Oxy_mic_bins <- count_microsyn_1Mb_with_fraction(
  sel_chr = Oxy_mic_pos$chr,
  sel_pos = Oxy_mic_pos$mid,
  all_chr = Oxy_all_pos$chr,
  all_pos = Oxy_all_pos$mid,
  dataset_name = "Oxyria_digyna_microsynteny (Gene.1)",
  binwidth = 1e6,
  min_len = 5e6
)

# plot
plot_microsyn_species(
  Dry_mic_bins,
  "./plots/Dryas_octopetala_microsynteny_gene1_count_and_fraction_1Mb.png",
  "./plots/Dryas_octopetala_microsynteny_gene1_count_and_fraction_1Mb.pdf"
)

plot_microsyn_species(
  Oxy_mic_bins,
  "./plots/Oxyria_digyna_microsynteny_gene1_count_and_fraction_1Mb.png",
  "./plots/Oxyria_digyna_microsynteny_gene1_count_and_fraction_1Mb.pdf"
)

#----------------------------
# Centromere locations - another time
#--------------------------
# # Centromere locations
# Dryas_cent_pos <- read.table("./Centromeres/DryOcto_H0-AT_total_possible_range.txt",sep="\t",header=T)
# Oxyria_cent_pos <- read.table("./Centromeres/Oxydig_H1-AT_total_possible_range.txt",sep="\t",header=T)


#---------------------
# Plot all
# LTRs, total gene duplications, core/accessory 
#---------------------

#-----------------------------
# 1) Standardize outputs
#    Expected columns:
#    chr, bin_start (bp), norm (0-1), track (dataset name), species
#-----------------------------

# Dry/Oxy gene density: your make_gene_density_norm_1Mb already outputs chr, bin_start, norm
gene_tracks <- bind_rows(
  make_gene_density_norm_1Mb(Dry_genes0, "Dryas_octopetala_H0") %>%
    mutate(track = "Gene density", species = "Dryas_octopetala_H0"),
  make_gene_density_norm_1Mb(Oxy_genes0, "Oxyria_digyna") %>%
    mutate(track = "Gene density", species = "Oxyria_digyna")
)

# Duplication (you already have) -> dup_all includes: chr, bin_start, norm, track, dataset
dup_tracks <- dup_all %>%
  transmute(
    chr,
    bin_start,
    count,
    norm,
    track = track,               # WGD/Tandem/...
    species = dataset
  )

# Total duplications -> you need to create a combined df in the same schema
dup_total_tracks <- bind_rows(
  dup_total_oxy %>% transmute(chr, bin_start, count, norm, track = "Total duplications", species = dataset),
  dup_total_dry %>% transmute(chr, bin_start, count, norm, track = "Total duplications", species = dataset)
)

# If you later want repeats-by-type, you can bind them too once you have them as:
# chr, bin_start, norm, track (repeat type), species


#-----------------------------
# 2) Combine everything into one long df
#-----------------------------
all_tracks_raw <- bind_rows(
  gene_tracks,
  dup_tracks,
  dup_total_tracks
)

# uses Oxyrt-1-86582034


# ---- repeats filter subset (edit these) ----
unique(Spp_TE_repeats0$V3)

 # [1] "Mutator_TIR_transposon"       "helitron"
 # [3] "Copia_LTR_retrotransposon"    "CACTA_TIR_transposon"
 # [5] "PIF_Harbinger_TIR_transposon" "LTR_retrotransposon"
 # [7] "hAT_TIR_transposon"           "Tc1_Mariner_TIR_transposon"
 # [9] "Gypsy_LTR_retrotransposon"    "repeat_region"
# [11] "target_site_duplication"      "long_terminal_repeat"

repeat_types_keep <- c(
  "Copia_LTR_retrotransposon", "Gypsy_LTR_retrotransposon", "CACTA_TIR_transposon", "hAT_TIR_transposon"
)

# ---- Helper: re-bin any “track-like” table already in chr/bin_start/norm ----
renorm_per_track_chr01 <- function(dat) {
  dat %>%
    group_by(species, track, chr) %>%
    mutate(norm01 = if (max(norm, na.rm = TRUE) == min(norm, na.rm = TRUE)) 0
                     else (norm - min(norm, na.rm = TRUE)) /
                          (max(norm, na.rm = TRUE) - min(norm, na.rm = TRUE))) %>%
    ungroup()
}

#-----------------------------
# 1) Repeats (subset by repeat type) into the common schema
#-----------------------------
make_repeats_bins_subset <- function(Spp_TE_repeats0, species_name,
                                      repeat_types_keep,
                                      binwidth = 1e6, threshold = 1e7,
                                      normalize_bins_to_0_1 = FALSE) {

  scaffold_lengths <- Spp_TE_repeats0 %>%
    group_by(V1) %>%
    summarise(length = max(V5), .groups = "drop")

  long_scaffolds <- scaffold_lengths %>%
    filter(length > threshold) %>%
    pull(V1)

  df <- Spp_TE_repeats0 %>%
    filter(V1 %in% long_scaffolds, V3 %in% repeat_types_keep)

  # bin on repeat start (V4)
  df2 <- df %>%
    transmute(
      chr = V1,
      start = as.numeric(V4),
      end   = as.numeric(V5),
      pos = as.numeric(V4),
      rep_bp = (as.numeric(V5) - as.numeric(V4) + 1),
      repeat_type = V3
    ) %>%
    filter(is.finite(rep_bp), rep_bp >= 0, is.finite(pos))

  binned <- df2 %>%
    mutate(bin_start = floor(pos / binwidth) * binwidth) %>%
    group_by(chr, bin_start, repeat_type) %>%
    summarise(total_rep_bp = sum(rep_bp, na.rm = TRUE), .groups = "drop") %>%
    mutate(rep_bp_per_1Mb = total_rep_bp / binwidth) %>%
    ungroup()

  # put into common schema: norm is either already 0-1 or raw
  if (normalize_bins_to_0_1) {
    binned <- binned %>%
      group_by(repeat_type, chr) %>%
      mutate(norm = if (max(rep_bp_per_1Mb, na.rm = TRUE) == min(rep_bp_per_1Mb, na.rm = TRUE)) 0
                     else (rep_bp_per_1Mb - min(rep_bp_per_1Mb, na.rm = TRUE)) /
                          (max(rep_bp_per_1Mb, na.rm = TRUE) - min(rep_bp_per_1Mb, na.rm = TRUE))) %>%
      ungroup()
  } else {
    binned <- binned %>% mutate(norm = rep_bp_per_1Mb)
  }

  binned %>%
    transmute(
      chr,
      bin_start,
      total_rep_bp,
      norm,
      track = paste0("Repeat: ", repeat_type),
      species = species_name
    )
}

# ---- build repeat track dataframes (subset) ----
rep_tracks <- bind_rows(
  make_repeats_bins_subset(Oxy_TE_repeats0, "Oxyria_digyna",
                            repeat_types_keep = repeat_types_keep,
                            normalize_bins_to_0_1 = FALSE),
  make_repeats_bins_subset(Dry_TE_repeats0, "Dryas_octopetala_H0",
                            repeat_types_keep = repeat_types_keep,
                            normalize_bins_to_0_1 = FALSE)
)

# uses Oxy-1-8814624


#-----------------------------
# 2) Positive selection tracks -> common schema
#-----------------------------
count_selected_genes_1Mb_to_tracks <- function(gff_tbl, gene_names, dataset_name,
                                                 binwidth = 1e6, threshold = 1e7) {
  gff_genes <- as.data.frame(gff_tbl) %>%
    filter(V3 == "gene") %>%
    mutate(gene_id = sub(".*ID=([^;]+).*", "\\1", V9),
           gene_id_up = toupper(gene_id)) %>%
    filter(toupper(gene_id_up) %in% toupper(gene_names))

  scaffold_lengths <- gff_genes %>%
    group_by(V1) %>%
    summarise(length = max(V5), .groups="drop")

  long_scaffolds <- scaffold_lengths %>%
    filter(length > threshold) %>%
    pull(V1)

  gff_genes <- gff_genes %>% filter(V1 %in% long_scaffolds)

  binned <- gff_genes %>%
    mutate(pos = as.numeric(V4),
           bin_start = floor(pos / binwidth) * binwidth) %>%
    group_by(V1, bin_start) %>%
    summarise(count = dplyr::n(), .groups="drop") %>%
    rename(chr = V1) %>%
    group_by(chr) %>%
    mutate(norm = if (max(count) == min(count)) 0
                   else (count - min(count)) / (max(count) - min(count))) %>%
    ungroup() %>%
    transmute(chr, bin_start, count, norm, 
              track = "Positively selected genes",
              species = dataset_name)

  binned
}

# ---- make dry/oxy selected gene lists exactly like your code ----
pos_sel <- read.table("./Positive_Selection/all_significant_genes_guidance_filtered.tsv",
                      sep="\t", header=F)
colnames(pos_sel) <- c("OG", "Spp", "p-value")

dry_genes <- pos_sel %>%
  filter(stringr::str_detect(Spp, "^DOCTH0_")) %>%
  pull(Spp) %>% unique()

oxy_genes <- pos_sel %>%
  filter(stringr::str_detect(Spp, "^OXYRIA")) %>%
  pull(Spp) %>% unique()

sel_tracks <- bind_rows(
  count_selected_genes_1Mb_to_tracks(Dry_genes0, dry_genes, "Dryas_octopetala_H0"),
  count_selected_genes_1Mb_to_tracks(Oxy_genes0, oxy_genes, "Oxyria_digyna")
)

#-----------------------------
# 3) Core/accessory/singleton tracks -> common schema
#-----------------------------
count_bed_1Mb_to_track <- function(bed_tbl, dataset_name,
                                     track_name,
                                     binwidth = 1e6, min_len = 5e6, use="midpoint") {
  b <- as.data.frame(bed_tbl)

  start <- as.numeric(b$V2)
  end   <- as.numeric(b$V3)
  chr   <- as.character(b$V1)

  chr_len <- data.frame(chr=chr, start=start, end=end) %>%
    group_by(chr) %>%
    summarise(length = max(end), .groups="drop")

  keep_chr <- chr_len %>% filter(length >= min_len) %>% pull(chr)

  keep <- chr %in% keep_chr
  start <- start[keep]; end <- end[keep]; chr <- chr[keep]

  pos <- switch(use,
    midpoint = (start + end)/2,
    start = start,
    end = end
  )

  tibble(chr = chr, pos = pos) %>%
    mutate(bin_start = floor(pos / binwidth) * binwidth) %>%
    group_by(chr, bin_start) %>%
    summarise(count = dplyr::n(), .groups="drop") %>%
    group_by(chr) %>%
    mutate(norm = if (max(count) == min(count)) 0
                   else (count - min(count)) / (max(count) - min(count))) %>%
    ungroup() %>%
    transmute(chr, bin_start, count, norm,
              track = track_name,
              species = dataset_name)
}

Dryas_core <- read.table("./Core_accessory/core_by_species/Dryas_octopetala_core.bed",sep="\t",header=F)
Oxyria_core <- read.table("./Core_accessory/core_by_species/Oxyria_digyna_H1_core.bed",sep="\t",header=F)

Dryas_accessory <- read.table("./Core_accessory/accessory_by_species/Dryas_octopetala_accessory.bed",sep="\t",header=F)
Oxyria_accessory <- read.table("./Core_accessory/accessory_by_species/Oxyria_digyna_H1_accessory.bed",sep="\t",header=F)

Dryas_singleton <- read.table("./Core_accessory/singleton_by_species/Dryas_octopetala_singleton.bed",sep="\t",header=F)
Oxyria_singleton <- read.table("./Core_accessory/singleton_by_species/Oxyria_digyna_H1_singleton.bed",sep="\t",header=F)

core_tracks <- bind_rows(
  count_bed_1Mb_to_track(Dryas_core, "Dryas_octopetala_H0", "Core genes", use="midpoint"),
  count_bed_1Mb_to_track(Oxyria_core, "Oxyria_digyna", "Core genes", use="midpoint"),
  count_bed_1Mb_to_track(Dryas_accessory, "Dryas_octopetala_H0", "Accessory genes", use="midpoint"),
  count_bed_1Mb_to_track(Oxyria_accessory, "Oxyria_digyna", "Accessory genes", use="midpoint"),
  count_bed_1Mb_to_track(Dryas_singleton, "Dryas_octopetala_H0", "Singleton genes", use="midpoint"),
  count_bed_1Mb_to_track(Oxyria_singleton, "Oxyria_digyna", "Singleton genes", use="midpoint")
)

#-----------------------------
# 4) Microsynteny tracks -> common schema
#-----------------------------
gff_gene_coords <- function(gff_tbl) {
  gff_genes <- as.data.frame(gff_tbl) %>% filter(V3 == "gene")
  tibble(
    gene_id = sub(".*ID=([^;]+).*", "\\1", gff_genes$V9),
    chr = as.character(gff_genes$V1),
    start = as.numeric(gff_genes$V4),
    end = as.numeric(gff_genes$V5)
  ) %>%
    mutate(gene_id_up = toupper(gene_id))
}

count_positions_1Mb_to_track <- function(chr, pos, dataset_name, track_name,
                                         binwidth=1e6, min_len=5e6) {
  chr <- as.character(chr)
  pos <- as.numeric(pos)

  chr_len <- tibble(chr = chr, pos = pos) %>%
    group_by(chr) %>%
    summarise(length = max(pos, na.rm = TRUE), .groups="drop")

  keep_chr <- chr_len %>% filter(length >= min_len) %>% pull(chr)

  df <- tibble(chr = chr, pos = pos) %>%
    filter(!is.na(pos), chr %in% keep_chr)

  df %>%
    mutate(bin_start = floor(pos / binwidth) * binwidth) %>%
    group_by(chr, bin_start) %>%
    summarise(count = dplyr::n(), .groups="drop") %>%
    group_by(chr) %>%
    mutate(norm = if (max(count) == min(count)) 0
                   else (count - min(count)) / (max(count) - min(count))) %>%
    ungroup() %>%
    transmute(chr, bin_start, count, norm,
              track = track_name,
              species = dataset_name)
}

microsyn <- read.table("Microsynteny.csv", sep=",", header=TRUE)

Dry_gene_col <- "Dryas.octopetala.Gene.1"
Oxy_gene_col <- "Oxyria.digyna.Gene.1"

Dry_gff <- gff_gene_coords(Dry_genes0)
Oxy_gff <- gff_gene_coords(Oxy_genes0)

Dry_mic_pos <- microsyn %>%
  transmute(gene_id = .data[[Dry_gene_col]],
            gene_id_up = toupper(gene_id)) %>%
  inner_join(Dry_gff %>% select(gene_id_up, chr, start, end), by="gene_id_up") %>%
  mutate(mid = (start + end)/2) %>%
  select(chr, mid)

Oxy_mic_pos <- microsyn %>%
  transmute(gene_id = .data[[Oxy_gene_col]],
            gene_id_up = toupper(gene_id)) %>%
  inner_join(Oxy_gff %>% select(gene_id_up, chr, start, end), by="gene_id_up") %>%
  mutate(mid = (start + end)/2) %>%
  select(chr, mid)

microsynt_tracks <- bind_rows(
  count_positions_1Mb_to_track(Dry_mic_pos$chr, Dry_mic_pos$mid,
                                "Dryas_octopetala_H0",
                                "Microsynteny (Gene.1)",
                                binwidth=1e6, min_len=5e6),
  count_positions_1Mb_to_track(Oxy_mic_pos$chr, Oxy_mic_pos$mid,
                                "Oxyria_digyna",
                                "Microsynteny (Gene.1)",
                                binwidth=1e6, min_len=5e6)
)

#-----------------------------
# 5) Combine with your existing gene+duplication tracks (assumes:
#    all_tracks_raw exists with columns: chr, bin_start, norm, track, species)
#-----------------------------
all_tracks_extended_raw <- bind_rows(
  all_tracks_raw,     # gene density + WGD/Tandem/... + total duplications (from your earlier combined step)
  rep_tracks,
  sel_tracks,
  core_tracks,
  microsynt_tracks
)

 unique(all_tracks_extended_raw$chr)
 # [1] "DoctH0-1"         "DoctH0-2"         "DoctH0-3"         "DoctH0-4"
 # [5] "DoctH0-5"         "DoctH0-6"         "DoctH0-7"         "DoctH0-8"
 # [9] "DoctH0-9"         "Oxyrt-1-86582034" "Oxyrt-2-79714091" "Oxyrt-3-79472951"
# [13] "Oxyrt-4-78410798" "Oxyrt-5-76064323" "Oxyrt-6-73303751" "Oxyrt-7-72361354"
# [17] "Oxy-1-8814624"    "Oxy-2-8078772"    "Oxy-3-7717500"    "Oxy-4-7603626"
# [21] "Oxy-5-7384279"    "Oxy-6-7328924"    "Oxy-7-7013624"


all_tracks_extended_raw <- all_tracks_extended_raw %>%
  mutate(
    chr = case_when(
      grepl("^Oxy\\-\\d+\\-\\d+$", chr) ~ sub("^Oxy\\-(\\d+)\\-\\d+$", "Oxy-\\1", chr),
      grepl("^Oxyrt\\-\\d+\\-\\d+$", chr) ~ sub("^Oxyrt\\-(\\d+)\\-\\d+$", "Oxy-\\1", chr),
      TRUE ~ chr
    )
  )

unique(all_tracks_extended_raw$chr)
 # [1] "DoctH0-1" "DoctH0-2" "DoctH0-3" "DoctH0-4" "DoctH0-5" "DoctH0-6"
 # [7] "DoctH0-7" "DoctH0-8" "DoctH0-9" "Oxy-1"    "Oxy-2"    "Oxy-3"
# [13] "Oxy-4"    "Oxy-5"    "Oxy-6"    "Oxy-7"


all_tracks_extended <- renorm_per_track_chr01(all_tracks_extended_raw) %>%
  rename(norm01 = norm01)

#-----------------------------
# 6) Separate plot per species
#-----------------------------
plot_species_tracks <- function(dat, species_name, out_png, out_pdf) {
  dat <- dat %>% filter(species == species_name)

  p <- ggplot(dat, aes(x = bin_start / 1e6, y = norm01,
                        color = track)) +
    geom_line(linewidth = 0.6, alpha = 0.95) +
    facet_wrap(~ chr, ncol = 1, scales = "fixed") +
    labs(x = "Genomic position (Mb, 1 Mb bins)",
         y = "Normalized value (0–1)",
         title = species_name) +
    theme_classic() +
    theme(legend.position = "right")

  ggsave(out_png, p, width = 7, height = 12, dpi = 300)
  ggsave(out_pdf, p, width = 7, height = 12)
}

plot_species_tracks(
  all_tracks_extended, "Dryas_octopetala_H0",
  "./plots/Dryas_combined_plus_repeats_sel_core_microsynt_norm01.png",
  "./plots/Dryas_combined_plus_repeats_sel_core_microsynt_norm01.pdf"
)

plot_species_tracks(
  all_tracks_extended, "Oxyria_digyna",
  "./plots/Oxyria_combined_plus_repeats_sel_core_microsynt_norm01.png",
  "./plots/Oxyria_combined_plus_repeats_sel_core_microsynt_norm01.pdf"
)

###################################
# Add ratio of genes out of total genes for each category

gene_class_tracks <- c(
  "Core genes", "Accessory genes",
  "Positively selected genes", "Microsynteny (Gene.1)"
)

wide <- all_tracks_extended_raw %>%
  filter(track %in% c("Gene density", gene_class_tracks)) %>%
  mutate(bin_start = as.numeric(bin_start)) %>%
  select(species, chr, bin_start, track, count) %>%
  pivot_wider(names_from = track, values_from = count, values_fill = 0) %>%
  mutate(total_gene_density = .data[["Gene density"]])

counts_long <- wide %>%
  pivot_longer(
    cols = all_of(gene_class_tracks),
    names_to = "track",
    values_to = "count"
  )

ratios_long <- wide %>%
  mutate(across(
    all_of(gene_class_tracks),
    ~ .x / .data[["Gene density"]],
    .names = "ratio_{.col}"
  )) %>%
  pivot_longer(
    cols = starts_with("ratio_"),
    names_to = "ratio_name",
    values_to = "ratio"
  ) %>%
  mutate(track = sub("^ratio_", "", ratio_name)) %>%
  select(-ratio_name)

ratios_bin_chr <- ratios_long %>%
  left_join(
    counts_long %>%
      select(species, chr, bin_start, track, count),
    by = c("species", "chr", "bin_start", "track")
  ) %>%
  select(species, chr, bin_start, track, ratio, count, total_gene_density)


plot_species_tracks_ratio <- function(dat, species_name, out_png, out_pdf) {
  dat <- dat %>% filter(species == species_name)

  p <- ggplot(dat, aes(x = bin_start / 1e6, y = ratio, color = track)) +
    geom_line(linewidth = 0.6, alpha = 0.95) +
    facet_wrap(~ chr, ncol = 1, scales = "fixed") +
    labs(x = "Genomic position (Mb, 1 Mb bins)",
         y = "Ratio",
         title = species_name) +
    theme_classic() +
    theme(legend.position = "right")

  ggsave(out_png, p, width = 7, height = 12, dpi = 300)
  ggsave(out_pdf, p, width = 7, height = 12)
}

plot_species_tracks_ratio(
  ratios_bin_chr, "Dryas_octopetala_H0",
  "./plots/Dryas_with_gene_class_ratios.png",
  "./plots/Dryas_with_gene_class_ratios.pdf"
)

plot_species_tracks_ratio(
  ratios_bin_chr, "Oxyria_digyna",
  "./plots/Oxyria_with_gene_class_ratios.png",
  "./plots/Oxyria_with_gene_class_ratios.pdf"
)


ratios_bin_chr <- ratios_bin_chr %>%
  group_by(chr, track) %>%
  mutate(
    ratio_plot = (ratio - min(ratio, na.rm = TRUE)) /
                  (max(ratio, na.rm = TRUE) - min(ratio, na.rm = TRUE))
  ) %>%
  ungroup() %>%
  mutate(
    ratio_plot = ifelse(is.nan(ratio_plot), 0, ratio_plot)  # handles max==min
  )


plot_species_tracks_rationorm <- function(dat, species_name, out_png, out_pdf) {
  dat <- dat %>% filter(species == species_name)

  p <- ggplot(dat, aes(x = bin_start / 1e6, y = ratio_plot, color = track)) +
    geom_line(linewidth = 0.6, alpha = 0.95) +
    facet_wrap(~ chr, ncol = 1, scales = "fixed") +
    labs(x = "Genomic position (Mb, 1 Mb bins)",
         y = "Ratio Normalized 0-1",
         title = species_name) +
    theme_classic() +
    theme(legend.position = "right")

  ggsave(out_png, p, width = 7, height = 12, dpi = 300)
  ggsave(out_pdf, p, width = 7, height = 12)
}

plot_species_tracks_rationorm(
  ratios_bin_chr, "Dryas_octopetala_H0",
  "./plots/Dryas_with_gene_class_ratiosnorm.png",
  "./plots/Dryas_with_gene_class_ratiosnorm.pdf"
)

plot_species_tracks_rationorm(
  ratios_bin_chr, "Oxyria_digyna",
  "./plots/Oxyria_with_gene_class_ratiosnorm.png",
  "./plots/Oxyria_with_gene_class_ratiosnorm.pdf"
)



# Join
all_tracks_with_ratios <- left_join(all_tracks_extended, ratios_bin_chr)
# Joining with `by = join_by(chr, bin_start, count, track, species)`


#-----------------------------
# Subset to reasonable number of factors

unique(all_tracks_with_ratios$track)
 # [1] "Gene density"                      "WGD"
 # [3] "Tandem"                            "Proximal"
 # [5] "Transposed"                        "Dispersed"
 # [7] "Total duplications"                "Repeat: CACTA_TIR_transposon"
 # [9] "Repeat: Copia_LTR_retrotransposon" "Repeat: Gypsy_LTR_retrotransposon"
# [11] "Repeat: hAT_TIR_transposon"        "Positively selected genes"
# [13] "Core genes"                        "Accessory genes"
# [15] "Singleton genes"                   "Microsynteny (Gene.1)"


# Keep only the tracks you want to show
tracks_keep <- c(
  "Total duplications",
  "Positively selected genes",
  "Core genes","Accessory genes",
  "Microsynteny (Gene.1)",
  "Repeat: Copia_LTR_retrotransposon"
)

all_tracks_subset <- all_tracks_with_ratios %>%
  filter(track %in% tracks_keep)

# then plot as before
plot_species_tracks(
  all_tracks_subset, "Dryas_octopetala_H0",
  "./plots/Dryas_subset_tracks_norm01.png",
  "./plots/Dryas_subset_tracks_norm01.pdf"
)

plot_species_tracks(
  all_tracks_subset, "Oxyria_digyna",
  "./plots/Oxyria_subset_tracks_norm01.png",
  "./plots/Oxyria_subset_tracks_norm01.pdf"
)

# then plot with ratios
plot_species_tracks_ratio(
  all_tracks_subset, "Dryas_octopetala_H0",
  "./plots/Dryas_subset_ratio_norm01.png",
  "./plots/Dryas_subset_ratio_norm01.pdf"
)

plot_species_tracks_ratio(
  all_tracks_subset, "Oxyria_digyna",
  "./plots/Oxyria_subset_ratio_norm01.png",
  "./plots/Oxyria_subset_ratio_norm01.pdf"
)


####################################
# Add centromere locations in and count 
oxychr <- c("Oxy-1", "Oxy-2", "Oxy-3", "Oxy-4",  
 "Oxy-5", "Oxy-6", "Oxy-7")
oxycent <- c("46e+06", "31e+06", "45e+06", "47e+06",  
 "45e+06", "30e+06", "31e+06")
oxychrlength <- c("79e+06", "80e+06", "78e+06", "76e+06",  
 "73e+06", "72e+06", "87e+06")
Oxy_cent <- cbind(oxychr, oxycent, oxychrlength)
# Chr 6: 29Mbp has same repeat!

drychr <- c("DoctH0-1", "DoctH0-2", "DoctH0-3", "DoctH0-4", 
"DoctH0-5", "DoctH0-6", "DoctH0-7", "DoctH0-8", "DoctH0-9")
drycent <- c("27e+06", "8.5e+06", "23e+06", "5e+06", 
"2e+06", "19e+06", "17.5e+06", "17e+06", "2.5e+06")
drychrlength <- c("34e+06", "33e+06", "28e+06", "27e+06", 
"23e+06", "22e+06", "22e+06", "21e+06", "20e+06")
Dry_cent <- cbind(drychr, drycent, drychrlength)

Oxy_cent
     # oxychr  oxycent  oxychrlength
# [1,] "Oxy-1" "79e+06" "79e+06"
# [2,] "Oxy-2" "80e+06" "80e+06"
# [3,] "Oxy-3" "78e+06" "78e+06"
# [4,] "Oxy-4" "76e+06" "76e+06"
# [5,] "Oxy-5" "73e+06" "73e+06"
# [6,] "Oxy-6" "72e+06" "72e+06"
# [7,] "Oxy-7" "87e+06" "87e+06"
Dry_cent
      # drychr     drycent    drychrlength
 # [1,] "DoctH0-1" "27e+06"   "34e+06"
 # [2,] "DoctH0-2" "8.5e+06"  "33e+06"
 # [3,] "DoctH0-3" "23e+06"   "28e+06"
 # [4,] "DoctH0-4" "5e+06"    "27e+06"
 # [5,] "DoctH0-5" "2e+06"    "23e+06"
 # [6,] "DoctH0-6" "19e+06"   "22e+06"
 # [7,] "DoctH0-7" "17.5e+06" "22e+06"
 # [8,] "DoctH0-8" "17e+06"   "21e+06"
 # [9,] "DoctH0-9" "2.5e+06"  "20e+06"

################################
# Compare pericentromeric and telomeric regions

win_size <- 1e6

oxy_cent_tbl <- as_tibble(Oxy_cent, .name_repair = "minimal") %>%
  transmute(chr = oxychr,
            centromere = as.numeric(oxycent))

dry_cent_tbl <- as_tibble(Dry_cent, .name_repair = "minimal") %>%
  transmute(chr = drychr,
            centromere = as.numeric(drycent))

cent_tbl <- bind_rows(
  oxy_cent_tbl %>% mutate(species = "Oxyria_digyna"),
  dry_cent_tbl %>% mutate(species = "Dryas_octopetala_H0")
)

tmp <- ratios_bin_chr %>%
  mutate(bin_start = as.numeric(bin_start),
         bin_center = bin_start + win_size/2) %>%
  left_join(cent_tbl, by = c("species", "chr")) %>%
  mutate(dist_to_cent = abs(bin_center - centromere))

# Label pericentric = closest half of windows within each species
ratios_bin_chr_labeled <- tmp %>%
  group_by(species) %>%
  mutate(
    k = n() %/% 2,
    rank_by_dist = rank(dist_to_cent, ties.method = "first")
  ) %>%
  mutate(pericentric = if_else(rank_by_dist <= k, 1L, 0L)) %>%
  ungroup() %>%
  select(-centromere, -dist_to_cent, -k, -rank_by_dist)


# Plot
ratios_bin_chr_labeled <- ratios_bin_chr_labeled %>%
  mutate(pericentric = factor(pericentric,
                               levels = c(0, 1),
                               labels = c("telomeric", "pericentric")))


p <- ggplot(ratios_bin_chr_labeled,
       aes(x = pericentric, y = ratio, fill = pericentric)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
  facet_grid(track ~ species, scales = "free_y") +
  theme_bw() +
  labs(x = NULL, y = "Ratio") +
  guides(fill = "none")

ggsave(
  filename = "ratio_boxplot_pericentric_telomeric.png",
  plot = p,
  width = 6, height = 12, dpi = 300
)

#####################################
# Look at 2 and 4 Mbp windows

library(dplyr)

make_window_ratios <- function(df, window_bp = 2e6) {
  # df: ratios_bin_chr-like with columns:
  # species, chr, bin_start, track, count, total_gene_density

  window_start_bp <- floor(df$bin_start / window_bp) * window_bp

  df <- df %>%
    mutate(bin_start_win = window_start_bp)

  # Denominator: gene density per window (same for all tracks, so take one copy per window)
  denom <- df %>%
    select(species, chr, bin_start_win, total_gene_density) %>%
    distinct() %>%
    group_by(species, chr, bin_start_win) %>%
    summarise(gene_density_win = sum(total_gene_density), .groups = "drop")

  # Numerator: sum each track's counts within the window
  num <- df %>%
    group_by(species, chr, bin_start_win, track) %>%
    summarise(count_win = sum(count), .groups = "drop")

  out <- num %>%
    left_join(denom, by = c("species", "chr", "bin_start_win")) %>%
    mutate(ratio = count_win / gene_density_win,
           ratio_plot = ratio) %>%   # placeholder; normalize later if you want
    select(species, chr, bin_start = bin_start_win, track, ratio, count = count_win,
           total_gene_density = gene_density_win, ratio_plot)

  out
}

ratios_2Mb <- make_window_ratios(ratios_bin_chr, 2e6)
ratios_4Mb <- make_window_ratios(ratios_bin_chr, 4e6)
ratios_6Mb <- make_window_ratios(ratios_bin_chr, 6e6)


scale_0_1_per_chr_track <- function(df) {
  df %>%
    group_by(chr, track) %>%
    mutate(ratio_plot = (ratio - min(ratio, na.rm = TRUE)) /
                         (max(ratio, na.rm = TRUE) - min(ratio, na.rm = TRUE))) %>%
    ungroup() %>%
    mutate(ratio_plot = ifelse(is.nan(ratio_plot), 0, ratio_plot))
}

ratios_2Mb <- scale_0_1_per_chr_track(ratios_2Mb)
ratios_4Mb <- scale_0_1_per_chr_track(ratios_4Mb)
ratios_6Mb <- scale_0_1_per_chr_track(ratios_6Mb)


#-------------------------------------
# Look at boxplots of different window sizes

# one third closest to centromere = pericentric
label_pericentric <- function(ratios_df, cent_tbl, win_size) {
  ratios_df %>%
    mutate(
      bin_start = as.numeric(bin_start),
      bin_center = bin_start + win_size / 2
    ) %>%
    left_join(cent_tbl, by = c("species", "chr")) %>%
    mutate(dist_to_cent = abs(bin_center - centromere)) %>%
    group_by(species) %>%
    mutate(
      k = n() %/% 2,
      #k = floor(4 * n() / 5),  # closest 3/4 are pericentric
      rank_by_dist = rank(dist_to_cent, ties.method = "first"),
      pericentric = if_else(rank_by_dist <= k, 1L, 0L)
    ) %>%
    ungroup() %>%
    select(-centromere, -dist_to_cent, -k, -rank_by_dist) %>%
    mutate(
      pericentric = factor(pericentric,
                            levels = c(0, 1),
                            labels = c("telomeric", "pericentric"))
    )
}


oxy_cent_tbl <- as_tibble(Oxy_cent, .name_repair = "minimal") %>%
  transmute(chr = oxychr,
            centromere = as.numeric(oxycent))

dry_cent_tbl <- as_tibble(Dry_cent, .name_repair = "minimal") %>%
  transmute(chr = drychr,
            centromere = as.numeric(drycent))

cent_tbl <- bind_rows(
  oxy_cent_tbl %>% mutate(species = "Oxyria_digyna"),
  dry_cent_tbl %>% mutate(species = "Dryas_octopetala_H0")
)


#----------------------
win_size <- 2e6

ratios_2Mb_labeled <- label_pericentric(ratios_2Mb, cent_tbl, win_size)

p2 <- ggplot(ratios_2Mb_labeled,
             aes(x = pericentric, y = ratio, fill = pericentric)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
  facet_grid(track ~ species, scales = "free_y") +
  theme_bw() +
  labs(x = NULL, y = "Ratio") +
  guides(fill = "none")

ggsave(
  filename = "ratio_boxplot_pericentric_telomeric_2Mb.png",
  plot = p2,
  width = 6, height = 12, dpi = 300
)

#--------------------

win_size <- 4e6
ratios_4Mb_labeled <- label_pericentric(ratios_4Mb, cent_tbl, win_size)
# then make p4 and ggsave with a different filename

p3 <- ggplot(ratios_4Mb_labeled,
             aes(x = pericentric, y = ratio, fill = pericentric)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
  facet_grid(track ~ species, scales = "free_y") +
  theme_bw() +
  labs(x = NULL, y = "Ratio") +
  guides(fill = "none")

ggsave(
  filename = "ratio_boxplot_pericentric_telomeric_4Mb.png",
  plot = p3,
  width = 6, height = 12, dpi = 300
)

#--------------------------
win_size <- 6e6
ratios_6Mb_labeled <- label_pericentric(ratios_6Mb, cent_tbl, win_size)
# then make p6 and ggsave with a different filename

p4 <- ggplot(ratios_6Mb_labeled,
             aes(x = pericentric, y = ratio, fill = pericentric)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
  facet_grid(track ~ species, scales = "free_y") +
  theme_bw() +
  labs(x = NULL, y = "Ratio") +
  guides(fill = "none")

ggsave(
  filename = "ratio_boxplot_pericentric_telomeric_6Mb.png",
  plot = p4,
  width = 6, height = 12, dpi = 300
)

############################################
# Check significance of gene type differences

df <- ratios_4Mb_labeled %>%
  filter(!is.na(pericentric)) %>%
  mutate(pericentric = factor(pericentric,
                               levels = c("telomeric","pericentric")))

yvar <- "ratio"   # change to "ratio_plot" if you want to test normalized values

tests <- df %>%
  group_by(species, track) %>%
  group_modify(~{
    tel <- .x[[yvar]][.x$pericentric == "telomeric"]
    per <- .x[[yvar]][.x$pericentric == "pericentric"]

    # need both groups
    if(length(tel) == 0 || length(per) == 0) {
      return(tibble(p_value = NA_real_, n_tel = length(tel), n_per = length(per)))
    }

    wt <- wilcox.test(tel, per, alternative = "two.sided")
    tibble(p_value = wt$p.value, n_tel = length(tel), n_per = length(per))
  }) %>%
  ungroup() %>%
  mutate(q_value = p.adjust(p_value, method = "BH"))


tests %>% arrange(q_value)

# A tibble: 8 × 6
  # species             track                      p_value n_tel n_per  q_value
  # <chr>               <chr>                        <dbl> <int> <int>    <dbl>
# 1 Dryas_octopetala_H0 Positively selected genes 2.00e-15    32    31 1.60e-14
# 2 Dryas_octopetala_H0 Accessory genes           1.77e-13    31    32 7.08e-13
# 3 Dryas_octopetala_H0 Core genes                3.88e-13    31    32 1.03e-12
# 4 Dryas_octopetala_H0 Microsynteny (Gene.1)     1.03e- 9    32    31 2.06e- 9
# 5 Oxyria_digyna       Accessory genes           1.07e- 5    69    69 1.71e- 5
# 6 Oxyria_digyna       Core genes                1.02e- 3    69    69 1.36e- 3
# 7 Oxyria_digyna       Positively selected genes 2.77e- 2    69    69 3.17e- 2
# 8 Oxyria_digyna       Microsynteny (Gene.1)     7.49e- 1    69    69 7.49e- 1


###################################

library(dplyr)
library(coin)

perm_tests <- ratios_4Mb_labeled %>%
  filter(!is.na(pericentric)) %>%
  group_by(species, track) %>%
  group_modify(~{
    # Need both groups present
    if(length(unique(.x$pericentric)) < 2) {
      return(tibble(p_value = NA_real_))
    }

    test <- coin::oneway_test(
      ratio ~ pericentric,
      data = .x,
      distribution = coin::approximate(nresample = 10000)
    )

    tibble(p_value = as.numeric(coin::pvalue(test)))
  }) %>%
  ungroup() %>%
  mutate(q_value = p.adjust(p_value, method = "BH"))


perm_tests %>% arrange(q_value)

  # species             track                     p_value q_value
  # <chr>               <chr>                       <dbl>   <dbl>
# 1 Dryas_octopetala_H0 Accessory genes            0       0
# 2 Dryas_octopetala_H0 Core genes                 0       0
# 3 Dryas_octopetala_H0 Microsynteny (Gene.1)      0       0
# 4 Dryas_octopetala_H0 Positively selected genes  0       0
# 5 Oxyria_digyna       Accessory genes            0       0
# 6 Oxyria_digyna       Core genes                 0.0009  0.0012
# 7 Oxyria_digyna       Positively selected genes  0.0185  0.0211
# 8 Oxyria_digyna       Microsynteny (Gene.1)      0.251   0.251


######################################
# Plot tracks for various window sizes

normalize_chr_track_0_1 <- function(dat) {
  dat %>%
    group_by(chr, track) %>%
    mutate(
      ratio_plot = (ratio - min(ratio, na.rm = TRUE)) /
                    (max(ratio, na.rm = TRUE) - min(ratio, na.rm = TRUE))
    ) %>%
    ungroup() %>%
    mutate(ratio_plot = ifelse(is.nan(ratio_plot), 0, ratio_plot))
}


plot_species_tracks_rationorm <- function(dat, species_name, win_size_bp,
                                            out_png, out_pdf,
                                            cent_tbl) {
  dat <- dat %>% filter(species == species_name)
  cent_sub <- cent_tbl %>% filter(species == species_name)

  win_label <- paste0(win_size_bp/1e6, " Mb")

  p <- ggplot(dat, aes(x = bin_start / 1e6, y = ratio_plot, color = track)) +
    geom_vline(data = cent_sub,
               aes(xintercept = centromere / 1e6),
               linewidth = 0.5, alpha = 0.7, inherit.aes = FALSE) +
    geom_line(linewidth = 0.6, alpha = 0.95) +
    facet_wrap(~ chr, ncol = 1, scales = "fixed") +
    labs(
      x = paste0("Genomic position (Mb, ", win_label, " windows)"),
      y = "Ratio Normalized 0-1",
      title = species_name
    ) +
    theme_classic() +
    theme(legend.position = "right")

  ggsave(out_png, p, width = 7, height = 12, dpi = 300)
  ggsave(out_pdf, p, width = 7, height = 12)
}


ratios_2Mb_plot <- normalize_chr_track_0_1(ratios_2Mb)
ratios_4Mb_plot <- normalize_chr_track_0_1(ratios_4Mb)
ratios_6Mb_plot <- normalize_chr_track_0_1(ratios_6Mb)


plot_species_tracks_rationorm(
  ratios_4Mb_plot, "Dryas_octopetala_H0", 4e6,
  "./plots/Dryas_4Mb_tracks_ratio_norm.png",
  "./plots/Dryas_4Mb_tracks_ratio_norm.pdf",
  cent_tbl
)

plot_species_tracks_rationorm(
  ratios_4Mb_plot, "Oxyria_digyna", 4e6,
  "./plots/Oxyria_4Mb_tracks_ratio_norm.png",
  "./plots/Oxyria_4Mb_tracks_ratio_norm.pdf",
  cent_tbl
)





















###############################
# Check overlap of positively selected genes and microsynteny genes
# Not significantly different than random

head(dry_genes)
# [1] "DOCTH0_CHR100001637" "DOCTH0_CHR300002655" "DOCTH0_CHR800000320"
# [4] "DOCTH0_CHR100008102" "DOCTH0_CHR400005265" "DOCTH0_CHR500006106"


head(microsyn$Dryas.octopetala.Gene.1)
# [1] "DoctH0_Chr600003875" "DoctH0_Chr600003875" "DoctH0_Chr600003912"
# [4] "DoctH0_Chr600003912" "DoctH0_Chr600003916" "DoctH0_Chr600003916"


norm <- function(x) {
  x |>
    toupper()
}

dry_n <- norm(dry_genes)
mic_n <- norm(microsyn$Dryas.octopetala.Gene.1)

ov <- unique(intersect(dry_n, mic_n))


ov
length(ov)
# Overlap 788
length(dry_n)
# Positively selected 3668
length(mic_n)
# Microsynteny related 8910
length(unique(Dry_gff$gene_id_up))
# Total genes 39696

M <- 39696
K <- 3668
N <- 8910
x <- 788

p_upper <- phyper(x - 1, m = K, n = M - K, k = N, lower.tail = FALSE)
p_upper
#  0.9320445
# randomly distributed and overlapping







































#############################
# OLDER



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
  