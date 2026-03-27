############################
# Look at positions of core and accessory genes
# Feb 2025
############################

# NARVAL1
cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Orthogroups

cp Orthogroups.GeneCount.tsv /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/Gene_ontology/core_accessory_genes/
cp Orthogroups.tsv /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/Gene_ontology/core_accessory_genes/

#----------------------------
cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/Gene_ontology/core_accessory_genes/

# Core genes - in all
awk 'NR==1 {print; next} 
{
  present=1;
  for(i=2;i<=NF;i++){
    if($i==0){present=0}
  }
  if(present==1){print}
}' Orthogroups.GeneCount.tsv > core_orthogroups.tsv

wc -l core_orthogroups.tsv
# 3352 core_orthogroups.tsv

#----------------------------
# Accessory genes - in only one
awk 'NR==1 {print; next}
{
  count=0;
  for(i=2;i<=NF;i++){
    if($i>0){count++}
  }
  if(count>0 && count<NF-20){print}
}' Orthogroups.GeneCount.tsv > accessory_orthogroups.tsv

wc -l accessory_orthogroups.tsv
# 8169 accessory_orthogroups.tsv

#----------------------------
# Mid genes - less than 3
awk 'NR==1 {print; next}
{
  count=0;
  for(i=2;i<=NF;i++){
    if($i>0){count++}
  }
  if(count>0 && count<NF-19){print}
}' Orthogroups.GeneCount.tsv > mid-accessory_orthogroups.tsv

wc -l mid-accessory_orthogroups.tsv
# 13477 mid-accessory_orthogroups.tsv


#----------------------------
# extract core gene IDs
cut -f1 core_orthogroups.tsv | tail -n +2 > core_ids.txt
grep -Ff core_ids.txt Orthogroups.tsv > core_genes.tsv
cut -f2- core_genes.tsv | tr ',' '\n' | tr -d ' ' > core_gene_ids.txt

#---
# Extract core gene positions
cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/bed

cat *.bed > /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/Gene_ontology/core_accessory_genes/all_genes.bed

cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/Gene_ontology/core_accessory_genes/

grep -Ff core_gene_ids.txt all_genes.bed > core_genes.bed


# split by species

awk '
$4 ~ /^DoctH0_/                {print > "Dryas_octopetala_core.bed"}
$4 ~ /^Oxyria_/                {print > "Oxyria_digyna_core.bed"}
$4 ~ /^g[0-9]+\.t1$/           {print > "Cochlearia_groenlandica_core.bed"}
$4 ~ /(snap|maker)-/           {print > "Draba_nivalis_core.bed"}
' core_genes.bed

 wc -l *core.bed
  # 2349 Cochlearia_groenlandica_core.bed
  # 2305 Draba_nivalis_core.bed
  # 1873 Dryas_octopetala_core.bed
  # 2450 Oxyria_digyna_core.bed
  # 8977 total

#---------------------------
# extract accessory gene IDs
cut -f1 accessory_orthogroups.tsv | tail -n +2 > accessory_ids.txt
grep -Ff accessory_ids.txt Orthogroups.tsv > accessory_genes.tsv
cut -f2- accessory_genes.tsv | tr ',' '\n' | tr -d ' ' > accessory_gene_ids.txt

#---
# Extract core gene positions
cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/Gene_ontology/core_accessory_genes/
grep -Ff accessory_gene_ids.txt all_genes.bed > accessory_genes.bed


# split by species

awk '
$4 ~ /^DoctH0_/                {print > "Dryas_octopetala_accessory.bed"}
$4 ~ /^Oxyria_/                {print > "Oxyria_digyna_accessory.bed"}
$4 ~ /^g[0-9]+\.t1$/           {print > "Cochlearia_groenlandica_accessory.bed"}
$4 ~ /(snap|maker)-/           {print > "Draba_nivalis_accessory.bed"}
' accessory_genes.bed

 wc -l *_accessory.bed
  # 1604 Cochlearia_groenlandica_accessory.bed
  # 1164 Draba_nivalis_accessory.bed
  # 6905 Dryas_octopetala_accessory.bed
  # 2242 Oxyria_digyna_accessory.bed
 # 11915 total

#-------------------
# Split total genes by species

awk '
$4 ~ /^DoctH0_/                {print > "Dryas_octopetala_all_genes.bed"}
$4 ~ /^Oxyria_/                {print > "Oxyria_digyna_all_genes.bed"}
$4 ~ /^g[0-9]+\.t1$/           {print > "Cochlearia_groenlandica_all_genes.bed"}
$4 ~ /(snap|maker)-/           {print > "Draba_nivalis_all_genes.bed"}
' all_genes.bed


##################################
# Plot in R 

cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/Gene_ontology/core_accessory_genes/

tmux new-session -s R
tmux attach-session -t R


module load StdEnv/2023 r/4.3.1
export R_LIBS_USER=/home/celphin/R/x86_64-pc-linux-gnu-library/4.3

R

install.packages(c("dplyr", "ggplot2", "readr", "tidyr"))

##############################
# Plot raw densities of core and accessory genes

library(ggplot2)
library(readr)
library(dplyr)
library(IRanges)

#---------------------------
# Function to read BED files
read_bed <- function(file, species, type) {
  read_tsv(file,
           col_names = c("chr", "start", "end", "gene"),
           show_col_types = FALSE) %>%
    mutate(species = species,
           type = type)
}

#---------------------------
# Read all files
df <- bind_rows(
  read_bed("Dryas_octopetala_core.bed", "Dryas_octopetala", "Core"),
  read_bed("Dryas_octopetala_accessory.bed", "Dryas_octopetala", "Accessory"),
  
  read_bed("Oxyria_digyna_core.bed", "Oxyria_digyna", "Core"),
  read_bed("Oxyria_digyna_accessory.bed", "Oxyria_digyna", "Accessory"),
  
  read_bed("Draba_nivalis_core.bed", "Draba_nivalis", "Core"),
  read_bed("Draba_nivalis_accessory.bed", "Draba_nivalis", "Accessory")
)

df$type <- factor(df$type, levels = c("Core", "Accessory"))

#---------------------------
# Filter scaffolds < 1 Mbp
chr_lengths <- df %>%
  group_by(species, chr) %>%
  summarise(chr_len = max(end), .groups = "drop") %>%
  filter(chr_len >= 1e6)

df <- df %>%
  semi_join(chr_lengths, by = c("species", "chr"))

#---------------------------
# Parameters for density
window_size <- 1e6  # 1 Mbp windows

#---------------------------
# Compute density per window
density_df <- list()

for (sp in unique(df$species)) {
  for (ch in unique(df$chr[df$species == sp])) {
    for (tp in unique(df$type)) {
      sub_df <- df %>% filter(species == sp, chr == ch, type == tp)
      if (nrow(sub_df) == 0) next
      
      chr_len <- max(sub_df$end)
      bins <- IRanges(start = seq(1, chr_len, by = window_size),
                      width = window_size)
      genes <- IRanges(start = sub_df$start, end = sub_df$end)
      counts <- countOverlaps(bins, genes)
      
      density_df[[length(density_df)+1]] <- data.frame(
        species = sp,
        chr = ch,
        type = tp,
        window_start = start(bins),
        window_end = end(bins),
        count = counts
      )
    }
  }
}

density_df <- bind_rows(density_df)

#---------------------------
# Assign y-axis positions for chromosomes and direction
density_df <- density_df %>%
  group_by(species) %>%
  mutate(
    chr = factor(chr, levels = unique(chr)),
    y_base = as.numeric(chr),
    # Core points down, Accessory points up
    y_min = ifelse(type == "Core", y_base - count / max(count), y_base),
    y_max = ifelse(type == "Core", y_base, y_base + count / max(count))
  ) %>%
  ungroup()

#---------------------------
# Plot densities per species
species_list <- unique(density_df$species)

for (sp in species_list) {
  df_sp <- density_df %>% filter(species == sp)
  
  # Backbone for chromosomes
  backbone_sp <- df_sp %>%
    group_by(chr, y_base) %>%
    summarise(start = min(window_start),
              end   = max(window_end),
              .groups = "drop")
  
  # Create plot
  p <- ggplot() +
    # Backbone lines
    geom_segment(data = backbone_sp,
                 aes(x = start, xend = end,
                     y = y_base, yend = y_base),
                 colour = "grey40", linewidth = 1) +
    
    # Density tracks as filled rectangles
    geom_rect(data = df_sp,
              aes(xmin = window_start, xmax = window_end,
                  ymin = y_min, ymax = y_max,
                  fill = type),
              colour = NA, alpha = 0.9) +
    
    scale_fill_manual(values = c("Core" = "#0B3D91",
                                 "Accessory" = "#8B0000")) +
    
    scale_y_continuous(
      breaks = unique(df_sp$y_base),
      labels = unique(df_sp$chr)
    ) +
    
    theme_bw() +
    labs(
      x = "Genomic position (bp)",
      y = "Chromosome",
      fill = "Gene type",
      title = sp
    )
  
  # Save plot
  ggsave(filename = paste0(sp, "_genes_density_updown.png"),
         plot = p,
         width = 10, height = 10, dpi = 800)
}

#########################
# One plot per species of core and accessory gene positions

library(ggplot2)
library(dplyr)


# List of species
species_list <- unique(df$species)

# Loop over species
for (sp in species_list) {
  
  # Subset data for this species
  df_sp <- df %>% filter(species == sp)
  backbone_sp <- backbone %>% filter(species == sp)
  
  # Create plot
  p <- ggplot() +
    # Backbone lines
    geom_segment(data = backbone_sp,
                 aes(x = start, xend = end,
                     y = y_base, yend = y_base),
                 colour = "grey40", linewidth = 1) +
    
    # Core/accessory genes
    geom_segment(data = df_sp,
             aes(x = start, xend = end,
                 y = y_pos, yend = y_pos,
                 colour = type),
             linewidth = 7,       # thicker lines
             alpha = 1) +         # fully opaque
     scale_colour_manual(values = c("Core" = "#0B3D91",      # very dark blue
                               "Accessory" = "#8B0000")) # dark red
    
    theme_bw() +
    labs(x = "Genomic position",
         y = "Chromosome",
         colour = "Gene type",
         title = sp)
  
  # Save figure for this species
  ggsave(filename = paste0(sp, "_genes_plot.png"), plot = p,
         width = 8, height = 5, dpi = 800)
}



########################
# Relative gene density plot: Core down, Accessory up
# To see releative to total number of genes
#-----------------------------------

library(ggplot2)
library(readr)
library(dplyr)
library(IRanges)

#---------------------------
# Function to read BED files
read_bed <- function(file, species, type = NA) {
  read_tsv(file,
           col_names = c("chr", "start", "end", "gene"),
           show_col_types = FALSE) %>%
    mutate(species = species,
           type = type)
}

#---------------------------
# Read Core/Accessory BEDs
df_genes <- bind_rows(
  read_bed("Dryas_octopetala_core.bed", "Dryas_octopetala", "Core"),
  read_bed("Dryas_octopetala_accessory.bed", "Dryas_octopetala", "Accessory"),
  
  read_bed("Oxyria_digyna_core.bed", "Oxyria_digyna", "Core"),
  read_bed("Oxyria_digyna_accessory.bed", "Oxyria_digyna", "Accessory"),
  
  read_bed("Draba_nivalis_core.bed", "Draba_nivalis", "Core"),
  read_bed("Draba_nivalis_accessory.bed", "Draba_nivalis", "Accessory")
)

df_genes$type <- factor(df_genes$type, levels = c("Core", "Accessory"))

#---------------------------
# Read all genes BED (for total counts)
all_genes <- bind_rows(
  read_bed("Dryas_octopetala_all_genes.bed", "Dryas_octopetala"),
  read_bed("Oxyria_digyna_all_genes.bed", "Oxyria_digyna"),
  read_bed("Draba_nivalis_all_genes.bed", "Draba_nivalis")
)

#---------------------------
# Filter scaffolds < 5 Mbp
chr_lengths <- all_genes %>%
  group_by(species, chr) %>%
  summarise(chr_len = max(end), .groups = "drop") %>%
  filter(chr_len >= 5e6)

all_genes <- all_genes %>%
  semi_join(chr_lengths, by = c("species", "chr"))

df_genes <- df_genes %>%
  semi_join(chr_lengths, by = c("species", "chr"))

#---------------------------
# Parameters
window_size <- 3e6  # 0.5 Mbp

#---------------------------
# Function to compute gene counts per window
compute_window_counts <- function(df, window_size) {
  res_list <- list()
  for (sp in unique(df$species)) {
    for (ch in unique(df$chr[df$species == sp])) {
      for (tp in unique(df$type[df$species == sp & df$chr == ch])) {
        sub_df <- df %>% filter(species == sp, chr == ch, type == tp)
        if (nrow(sub_df) == 0) next
        
        chr_len <- max(sub_df$end)
        bins <- IRanges(start = seq(1, chr_len, by = window_size),
                        width = window_size)
        genes <- IRanges(start = sub_df$start, end = sub_df$end)
        counts <- countOverlaps(bins, genes)
        
        res_list[[length(res_list)+1]] <- data.frame(
          species = sp,
          chr = ch,
          type = tp,
          window_start = start(bins),
          window_end = end(bins),
          count = counts
        )
      }
    }
  }
  bind_rows(res_list)
}

# Core/Accessory counts per window
coreacc_window_counts <- compute_window_counts(df_genes, window_size)

# Total genes per window (no type)
compute_window_counts_total <- function(df, window_size) {
  res_list <- list()
  for (sp in unique(df$species)) {
    for (ch in unique(df$chr[df$species == sp])) {
      sub_df <- df %>% filter(species == sp, chr == ch)
      if (nrow(sub_df) == 0) next
      
      chr_len <- max(sub_df$end)
      bins <- IRanges(start = seq(1, chr_len, by = window_size),
                      width = window_size)
      genes <- IRanges(start = sub_df$start, end = sub_df$end)
      counts <- countOverlaps(bins, genes)
      
      res_list[[length(res_list)+1]] <- data.frame(
        species = sp,
        chr = ch,
        window_start = start(bins),
        window_end = end(bins),
        count_total = counts
      )
    }
  }
  bind_rows(res_list)
}

total_window_counts <- compute_window_counts_total(all_genes, window_size)

#---------------------------
# Merge Core/Accessory with total and compute relative density
density_df <- coreacc_window_counts %>%
  left_join(total_window_counts, by = c("species", "chr", "window_start", "window_end")) %>%
  mutate(
    rel_gene = count / count_total
  )

#---------------------------
# Assign y-axis positions (Core down, Accessory up)
density_df <- density_df %>%
  group_by(species) %>%
  mutate(
    y_base = as.numeric(factor(chr, levels = unique(chr))),
    y_min = ifelse(type == "Core", y_base - rel_gene, y_base),
    y_max = ifelse(type == "Core", y_base, y_base + rel_gene)
  ) %>%
  ungroup()

#---------------------------
# Backbone for chromosomes
backbone <- density_df %>%
  group_by(species, chr, y_base) %>%
  summarise(start = min(window_start),
            end   = max(window_end),
            .groups = "drop")

#---------------------------
# Plot per species
species_list <- unique(density_df$species)

for (sp in species_list) {
  df_sp <- density_df %>% filter(species == sp)
  backbone_sp <- backbone %>% filter(species == sp)
  
  p <- ggplot() +
    geom_segment(data = backbone_sp,
                 aes(x = start, xend = end, y = y_base, yend = y_base),
                 colour = "grey40", linewidth = 1) +
    geom_rect(data = df_sp,
              aes(xmin = window_start, xmax = window_end,
                  ymin = y_min, ymax = y_max,
                  fill = type),
              colour = NA, alpha = 0.9) +
    scale_fill_manual(values = c("Core" = "#0B3D91",
                                 "Accessory" = "#8B0000")) +
    scale_y_continuous(
      breaks = unique(df_sp$y_base),
      labels = unique(df_sp$chr)
    ) +
    theme_bw() +
    labs(
      x = "Genomic position (bp)",
      y = "Chromosome",
      fill = "Gene type",
      title = sp
    )
  
  ggsave(filename = paste0(sp, "_relative_gene_density_updown.png"),
         plot = p,
         width = 10, height = 10, dpi = 800)
}

###################################
# Manually scale
#---------------------------
# Define separate scale factors for visual balance
# Scaling factors per species
scale_table <- tibble::tibble(
  species = c("Dryas_octopetala", "Oxyria_digyna", "Draba_nivalis"),
  core_scale = c(5, 2, 3),       # manually tuned per species
  accessory_scale = c(1, 2, 3)   # usually 1, but you can adjust
)

# Apply scaling
density_df <- density_df %>%
  left_join(scale_table, by = "species") %>%
  mutate(
    y_min_scaled = case_when(
      type == "Core"      ~ y_base - rel_gene * core_scale,
      type == "Accessory" ~ y_base
    ),
    y_max_scaled = case_when(
      type == "Core"      ~ y_base,
      type == "Accessory" ~ y_base + rel_gene * accessory_scale
    )
  )

for (sp in species_list) {
  df_sp <- density_df %>% filter(species == sp)
  backbone_sp <- backbone %>% filter(species == sp)
  
  p <- ggplot() +
    geom_segment(data = backbone_sp,
                 aes(x = start, xend = end, y = y_base, yend = y_base),
                 colour = "grey40", linewidth = 1) +
    geom_rect(data = df_sp,
              aes(xmin = window_start, xmax = window_end,
                  ymin = y_min_scaled, ymax = y_max_scaled,
                  fill = type),
              colour = NA, alpha = 0.9) +
    scale_fill_manual(values = c("Core" = "#0B3D91",
                                 "Accessory" = "#8B0000")) +
    scale_y_continuous(
      breaks = unique(df_sp$y_base),
      labels = unique(df_sp$chr)
    ) +
    theme_bw() +
    labs(
      x = "Genomic position (bp)",
      y = "Chromosome",
      fill = "Gene type",
      title = sp
    )
  
  ggsave(filename = paste0(sp, "_relative_gene_density_updown_scaled.png"),
         plot = p,
         width = 5, height = 10, dpi = 300)
}

############################


# Auto Scale accessory and core genes

#---------------------------
# Compute species-specific scaling factors
scale_table <- density_df %>%
  group_by(species, type) %>%
  summarise(max_rel = max(rel_gene, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = type, values_from = max_rel, values_fill = 0) %>%
  mutate(
    # Core scaled to match Accessory max
    core_scale = ifelse(Core > 0, Accessory / Core, 1),
    accessory_scale = 1  # keep accessory scale = 1
  ) %>%
  select(species, core_scale, accessory_scale)

#---------------------------
# Merge scale factors into density dataframe
density_df <- density_df %>%
  left_join(scale_table, by = "species") %>%
  mutate(
    y_min_scaled = case_when(
      type == "Core"      ~ y_base - rel_gene * core_scale,
      type == "Accessory" ~ y_base
    ),
    y_max_scaled = case_when(
      type == "Core"      ~ y_base,
      type == "Accessory" ~ y_base + rel_gene * accessory_scale
    )
  )

#---------------------------
# Plot per species
for (sp in species_list) {
  df_sp <- density_df %>% filter(species == sp)
  backbone_sp <- backbone %>% filter(species == sp)
  
  p <- ggplot() +
    geom_segment(data = backbone_sp,
                 aes(x = start, xend = end, y = y_base, yend = y_base),
                 colour = "grey40", linewidth = 1) +
    geom_rect(data = df_sp,
              aes(xmin = window_start, xmax = window_end,
                  ymin = y_min_scaled, ymax = y_max_scaled,
                  fill = type),
              colour = NA, alpha = 0.9) +
    scale_fill_manual(values = c("Core" = "#0B3D91",
                                 "Accessory" = "#8B0000")) +
    scale_y_continuous(
      breaks = unique(df_sp$y_base),
      labels = unique(df_sp$chr)
    ) +
    theme_bw() +
    labs(
      x = "Genomic position (bp)",
      y = "Chromosome",
      fill = "Gene type",
      title = sp
    )
  
  ggsave(filename = paste0(sp, "_relative_gene_density_updown_auto_scaled.png"),
         plot = p,
         width = 10, height = 5, dpi = 800)
}
