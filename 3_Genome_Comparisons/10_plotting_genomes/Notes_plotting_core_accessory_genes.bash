############################
# Look at positions of core and accessory genes
# Feb 2025
############################

# NARVAL3
tmux attach-session -t Oxyria1

cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Orthogroups

cp Orthogroups.GeneCount.tsv /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/Gene_ontology/core_accessory_genes/
cp Orthogroups.tsv /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/Gene_ontology/core_accessory_genes/

#----------------------------
cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/Gene_ontology/core_accessory_genes/

# Core genes - in all
awk 'NR==1 {print; next}
{
  count=0;
  for(i=2;i<=NF;i++){
    if($i>0){count++}
  }
  if(count >= NF-1){print}
}' Orthogroups.GeneCount.tsv > core_orthogroups.tsv

wc -l core_orthogroups.tsv
# 3352 core_orthogroups.tsv

#----------------------------
# singleton genes - in only one
awk 'NR==1 {print; next}
{
  count=0;
  for(i=2;i<=NF;i++){
    if($i>0){count++}
  }
  if(count>0 && count<NF-20){print}
}' Orthogroups.GeneCount.tsv > singleton_orthogroups.tsv

wc -l singleton_orthogroups.tsv
# 8169 singleton_orthogroups.tsv

#-------------------------------
# accessory genes - in only one or two
awk 'NR==1 {print; next}
{
  count=0;
  for(i=2;i<=NF;i++){
    if($i>0){count++}
  }
  if(count>0 && count<NF-19){print}
}' Orthogroups.GeneCount.tsv > accessory_orthogroups.tsv

wc -l accessory_orthogroups.tsv
# 13477 accessory_orthogroups.tsv

#----------------------------
# Two spp - in 2 spp
awk 'NR==1 {print; next}
{
  count=0;
  for(i=2;i<=NF;i++){
    if($i>0){count++}
  }
  if(count>2 && count<4){print}
}' Orthogroups.GeneCount.tsv > twospp_orthogroups.tsv

wc -l twospp_orthogroups.tsv
# 5309 twospp_orthogroups.tsv


#----------------------------
# extract core gene IDs
cut -f1 core_orthogroups.tsv | tail -n +2 > core_ids.txt
grep -Ff core_ids.txt Orthogroups.tsv > core_genes.tsv
cut -f2- core_genes.tsv | tr ',' '\n' | tr -d ' ' > core_gene_ids.txt

# extract accessory gene IDs
cut -f1 accessory_orthogroups.tsv | tail -n +2 > accessory_ids.txt
grep -Ff accessory_ids.txt Orthogroups.tsv > accessory_genes.tsv
cut -f2- accessory_genes.tsv | tr ',' '\n' | tr -d ' ' > accessory_gene_ids.txt

# extract singleton gene IDs
cut -f1 singleton_orthogroups.tsv | tail -n +2 > singleton_ids.txt
grep -Ff singleton_ids.txt Orthogroups.tsv > singleton_genes.tsv
cut -f2- singleton_genes.tsv | tr ',' '\n' | tr -d ' ' > singleton_gene_ids.txt

# extract twospp gene IDs
cut -f1 twospp_orthogroups.tsv | tail -n +2 > twospp_ids.txt
grep -Ff twospp_ids.txt Orthogroups.tsv > twospp_genes.tsv
cut -f2- twospp_genes.tsv | tr ',' '\n' | tr -d ' ' > twospp_gene_ids.txt

#--------------------------------
# formatting 
# space to new line
tr -s '[:space:]' '\n' < core_gene_ids.txt > core_flat.txt
tr -s '[:space:]' '\n' < accessory_gene_ids.txt > accessory_flat.txt
tr -s '[:space:]' '\n' < singleton_gene_ids.txt > singleton_flat.txt
tr -s '[:space:]' '\n' < twospp_gene_ids.txt > twospp_flat.txt

# remove lines with ... or back
grep -v -E '^\.\.\.|^back' core_flat.txt > core_clean.txt
grep -v -E '^\.\.\.|^back' accessory_flat.txt > accessory_clean.txt
grep -v -E '^\.\.\.|^back' singleton_flat.txt > singleton_clean.txt
grep -v -E '^\.\.\.|^back' twospp_flat.txt > twospp_clean.txt

# remove blank lines and sort
awk '{$1=$1} NF' core_clean.txt | sort -u > core_gene_ids_sorted.txt
awk '{$1=$1} NF' accessory_clean.txt | sort -u > accessory_gene_ids_sorted.txt
awk '{$1=$1} NF' singleton_clean.txt | sort -u > singleton_gene_ids_sorted.txt
awk '{$1=$1} NF' twospp_clean.txt | sort -u > twospp_gene_ids_sorted.txt

#--------------------------
# Create all genes bed

# working dirs
BEDDIR="/home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/bed"
OUTDIR="/home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/Gene_ontology/core_accessory_genes"
mkdir -p "$OUTDIR"
cd "$BEDDIR"

# 1) Create combined BED with source as first column (species name from filename without .bed)
echo "Creating all_genes.with_src.bed..."
rm -f "$OUTDIR/all_genes.with_src.bed"
for f in *.bed; do
  species="${f%.bed}"
  # preserve original 6-column BED (or whatever), prepend species
  awk -v src="$species" 'BEGIN{OFS="\t"} {print src,$0}' "$f" >> "$OUTDIR/all_genes.with_src.bed"
done

# 2) Remove duplicates: keep longest isoform per gene ID (assumes gene ID is in column 5 of original -> now column 6)
# Adjust column index if your original .bed had different column counts.
echo "Deduplicating to keep longest isoform..."
cd "$OUTDIR"
awk -v OFS="\t" '
{
  # fields: 1=src, 2=chr, 3=start, 4=end, 5=name (if original had 4 cols); adjust if needed
  src=$1
  chr=$2
  start=$3
  end=$4
  name=$5

  # Remove trailing isoform suffix like -... if desired (but keep whole name for uniqueness)
  gene=name
  sub(/-.*$/, "", gene)

  len=end-start
  key = src"\t"gene

  if(!(key in best) || len > best_len[key]){
    best[key]=$0
    best_len[key]=len
  }
}
END {
  for(k in best) print best[k]
}' all_genes.with_src.bed > all_genes.with_src.best_isoform.bed

# 3) Optionally generate a plain all_genes.bed (without src) if you need it
cut -f2- all_genes.with_src.best_isoform.bed > all_genes.best_isoform.bed
# also keep the original all_genes (with src)
cp all_genes.with_src.bed all_genes.bed.with_src.orig

# 4) Extract gene sets (core/accessory/singleton/twospp) by ID lists.
# Assumes ID lists contain the gene IDs in the same form as the NAME column in the bed (before you stripped -R... earlier).
# If your ID lists have no species prefix and your NAME column includes species-specific suffixes, adjust the sub() used below.

extract_set() {
  idsfile=$1
  outprefix=$2
  awk -v idsfile="$idsfile" 'BEGIN{
    FS="\t"; OFS="\t"
    while((getline < idsfile) > 0) ids[$1]=1
    close(idsfile)
  }
  {
    name=$5
    # normalize name to match IDs file by removing trailing -R[A-Z] if your IDs are without that; remove or change if not needed
    n=name
    sub(/-R[A-Z]$/, "", n)
    if(n in ids) print
  }' all_genes.with_src.best_isoform.bed > "${outprefix}.with_src.bed"
}
extract_set core_gene_ids_sorted.txt core_genes
extract_set accessory_gene_ids_sorted.txt accessory_genes
extract_set singleton_gene_ids_sorted.txt singleton_genes
extract_set twospp_gene_ids_sorted.txt twospp_genes


# Check 
wc -l *.bed
    # 66187 accessory_genes.with_src.bed
   # 617927 all_genes.best_isoform.bed
   # 651482 all_genes.with_src.bed
   # 617927 all_genes.with_src.best_isoform.bed
   # 153928 core_genes.with_src.bed
    # 41575 singleton_genes.with_src.bed
    # 24612 twospp_genes.with_src.bed


# 5) Split each set into per-species BED files (remove the source column in output)
split_by_species() {
  infile=$1   # e.g., core_genes.with_src.bed
  prefix=$2   # e.g., core
  mkdir -p "${prefix}_by_species"
  awk -v outdir="${prefix}_by_species" 'BEGIN{FS="\t"; OFS="\t"}
  {
    src=$1
    # rebuild output filename safe (no slashes)
    gsub(/[\/ ]/, "_", src)
    out = outdir"/"src"_" "'"${prefix}"'.bed"
    # print columns 2..end (drop src column)
    $1=""; sub(/^\t/,"")
    print $0 >> out
  }
  END {
    # flush all files (awk typically does this on exit)
  }' "$infile"
  echo "Wrote files to ${prefix}_by_species/"
}

split_by_species core_genes.with_src.bed core
split_by_species accessory_genes.with_src.bed accessory
split_by_species singleton_genes.with_src.bed singleton
split_by_species twospp_genes.with_src.bed twospp

# Accessory
   # 1374 Arabidopsis_lyrata_accessory.bed
    # 718 Arabidopsis_thaliana_accessory.bed
   # 1967 Arabis_alpina_accessory.bed
    # 298 Argentina_anserina_accessory.bed
   # 2429 Brassica_oleracea_accessory.bed
    # 377 Capsella_rubella_accessory.bed
   # 3866 Cochlearia_groenlandica_accessory.bed
  # 10995 Dryas_octopetala_accessory.bed
  # 13379 Fagopyrum_escelentum_H2_accessory.bed
   # 5308 Fagopyrum_tataricum_H1_accessory.bed
    # 635 Fragaria_vesca_accessory.bed
   # 2091 Malus_sylvestris_accessory.bed
   # 5439 Oxyria_digyna_H1_accessory.bed
   # 1398 Polygunum_aviculare_H0_accessory.bed
    # 748 Prunus_persica_accessory.bed
   # 1033 Pyrus_bretschneideri_accessory.bed
   # 6113 Rheum_nobile_H0_accessory.bed
   # 4691 Rheum_tangaticum_H0_accessory.bed
   # 2182 Rosa_rugosa_accessory.bed
   # 1146 Thlaspi_arvense_accessory.bed

# Core

   # 7931 Arabidopsis_lyrata_core.bed
   # 6008 Arabidopsis_thaliana_core.bed
   # 6288 Arabis_alpina_core.bed
   # 5675 Argentina_anserina_core.bed
  # 11343 Brassica_oleracea_core.bed
   # 7582 Capsella_rubella_core.bed
   # 7803 Cochlearia_groenlandica_core.bed
   # 6522 Dryas_octopetala_core.bed
   # 8936 Fagopyrum_escelentum_H2_core.bed
   # 8738 Fagopyrum_tataricum_H1_core.bed
   # 6061 Fragaria_vesca_core.bed
   # 9899 Malus_sylvestris_core.bed
   # 7457 Oxyria_digyna_H1_core.bed
   # 7435 Polygunum_aviculare_H0_core.bed
   # 6302 Prunus_persica_core.bed
   # 9974 Pyrus_bretschneideri_core.bed
   # 7805 Rheum_nobile_H0_core.bed
   # 7678 Rheum_tangaticum_H0_core.bed
   # 7215 Rosa_rugosa_core.bed
   # 7276 Thlaspi_arvense_core.bed

# singleton
    # 546 Arabidopsis_lyrata_singleton.bed
    # 163 Arabidopsis_thaliana_singleton.bed
   # 1467 Arabis_alpina_singleton.bed
     # 95 Argentina_anserina_singleton.bed
   # 1650 Brassica_oleracea_singleton.bed
    # 149 Capsella_rubella_singleton.bed
   # 3011 Cochlearia_groenlandica_singleton.bed
   # 8927 Dryas_octopetala_singleton.bed
   # 9726 Fagopyrum_escelentum_H2_singleton.bed
   # 2782 Fagopyrum_tataricum_H1_singleton.bed
    # 181 Fragaria_vesca_singleton.bed
   # 1018 Malus_sylvestris_singleton.bed
   # 3266 Oxyria_digyna_H1_singleton.bed
    # 733 Polygunum_aviculare_H0_singleton.bed
    # 368 Prunus_persica_singleton.bed
    # 228 Pyrus_bretschneideri_singleton.bed
   # 3218 Rheum_nobile_H0_singleton.bed
   # 2293 Rheum_tangaticum_H0_singleton.bed
   # 1057 Rosa_rugosa_singleton.bed
    # 697 Thlaspi_arvense_singleton.bed


##################################
# Plot in R 

cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/Gene_ontology/core_accessory_genes/

tmux new-session -s R
tmux attach-session -t R


module load StdEnv/2023 r/4.3.1
export R_LIBS_USER=/home/celphin/R/x86_64-pc-linux-gnu-library/4.3

R

install.packages(c("dplyr", "ggplot2", "readr", "tidyr"))


#!/usr/bin/env Rscript
library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(IRanges)
library(fs)
library(purrr)

# Parameters
window_size <- 3e6
min_chr_len <- 1e6

# helper to read all files in a folder into a dataframe with species inferred from filename
read_folder <- function(dir, type_label) {
  files <- dir_ls(dir, glob = "*.bed")
  if(length(files)==0) return(tibble())
  df_list <- lapply(files, function(f) {
    species <- path_file(f) %>%
      path_ext_remove() %>%
      sub(paste0("(?i)_", type_label, "$"), "", ., perl = TRUE)
    dat <- read_tsv(f, col_names = FALSE, show_col_types = FALSE,
                    col_types = cols(.default = "c"))
    if(ncol(dat) < 4) stop("BED has <4 columns: ", f)
    tibble(chr = dat[[1]],
           start = as.integer(dat[[2]]),
           end = as.integer(dat[[3]]),
           gene = as.character(dat[[4]]),
           species = species,
           type = type_label)
  })
  bind_rows(df_list)
}


# Read per-type by-species directories
core_df <- read_folder("core_by_species", "Core")
accessory_df <- read_folder("accessory_by_species", "Accessory")
singleton_df <- read_folder("singleton_by_species", "Singleton")
twospp_df <- read_folder("twospp_by_species", "TwoSpp")

# Combine
df_types <- bind_rows(core_df, accessory_df, singleton_df, twospp_df) %>%
  mutate(type = factor(type, levels = c("Core","Accessory","Singleton","TwoSpp")))

all_genes <- read_tsv("all_genes.best_isoform.bed",
                      col_names = FALSE,
                      show_col_types = FALSE,
                      col_types = cols(.default = "c")) %>%
  as_tibble()

all_genes <- all_genes %>%
  mutate(
    chr = X1,
    start = as.integer(X2),
    end = as.integer(X3),
    gene = pmap_chr(across(everything()), function(...) {
      vals <- list(...)
      vv <- unlist(vals, use.names = FALSE)
      vv <- vv[!is.na(vv) & vv != ""]
      if(length(vv) >= 4) vv[length(vv)] else NA_character_
    })
  ) %>%
  select(chr, start, end, gene) %>%
  filter(!is.na(chr) & !is.na(start) & !is.na(end) & !is.na(gene))

# Build total genes per species by combining unique genes from all type dfs (safer)
# ensure no list-columns
df_types <- df_types %>% mutate(
  chr = as.character(chr),
  start = as.integer(start),
  end = as.integer(end),
  gene = as.character(gene),
  species = as.character(species),
  type = as.character(type)
)

# drop any remaining list-columns (if present)
list_cols <- names(df_types)[vapply(df_types, is.list, logical(1))]
if(length(list_cols) > 0) df_types <- df_types %>% select(-all_of(list_cols))

# now build totals
all_genes_by_species <- df_types %>%
  distinct(species, chr, start, end, gene) %>%
  arrange(species, chr, start)

# Filter chromosomes by length >= min_chr_len (computed from union of all genes for each species)
chr_lengths <- all_genes_by_species %>%
  group_by(species, chr) %>%
  summarise(chr_len = max(end), .groups = "drop") %>%
  filter(chr_len >= min_chr_len)

# keep only entries on large scaffolds
df_types <- df_types %>% semi_join(chr_lengths, by = c("species","chr"))
all_genes_by_species <- all_genes_by_species %>% semi_join(chr_lengths, by = c("species","chr"))

# window counting functions (works with df having species/chr/start/end and optionally type)
compute_window_counts <- function(df, window_size) {
  res_list <- list()
  for(sp in unique(df$species)) {
    chrs <- unique(df$chr[df$species==sp])
    for(ch in chrs) {
      sub_df <- df %>% filter(species==sp, chr==ch)
      chr_len <- max(sub_df$end)
      bins <- IRanges(start = seq(1, chr_len, by = window_size), width = window_size)
      if("type" %in% names(sub_df)) {
        for(tp in unique(sub_df$type)) {
          genes <- IRanges(start = sub_df$start[sub_df$type==tp], end = sub_df$end[sub_df$type==tp])
          counts <- countOverlaps(bins, genes)
          res_list[[length(res_list)+1]] <- data.frame(species=sp, chr=ch, type=tp,
                                                      window_start=start(bins), window_end=end(bins),
                                                      count=counts)
        }
      } else {
        genes <- IRanges(start = sub_df$start, end = sub_df$end)
        counts <- countOverlaps(bins, genes)
        res_list[[length(res_list)+1]] <- data.frame(species=sp, chr=ch,
                                                    window_start=start(bins), window_end=end(bins),
                                                    count_total=counts)
      }
    }
  }
  bind_rows(res_list)
}

# compute per-type window counts
coreacc_window_counts <- compute_window_counts(df_types, window_size)

# compute total genes per window
total_window_counts <- compute_window_counts(all_genes_by_species, window_size)

# For each species, build merged table and plot
species_list <- unique(df_types$species)

for(sp in species_list) {
  message("Processing ", sp)
  cc <- coreacc_window_counts %>% filter(species==sp)
  tot <- total_window_counts %>% filter(species==sp)
  if(nrow(cc)==0 || nrow(tot)==0) {
    warning("No data for ", sp); next
  }
  density_wide <- cc %>%
    pivot_wider(names_from = type, values_from = count, values_fill = 0) %>%
    left_join(tot, by = c("species","chr","window_start","window_end"))
  # ensure Core and Accessory columns exist
  if(!"Core" %in% names(density_wide)) density_wide$Core <- 0
  if(!"Accessory" %in% names(density_wide)) density_wide$Accessory <- 0
  # proportions
  p_core <- sum(density_wide$Core) / sum(density_wide$count_total)
  p_accessory <- sum(density_wide$Accessory) / sum(density_wide$count_total)
  density_wide <- density_wide %>%
    mutate(expected_core = count_total * p_core,
           expected_accessory = count_total * p_accessory,
           log2_core_enrichment = log2((Core + 1e-6) / expected_core),
           log2_accessory_enrichment = log2((Accessory + 1e-6) / expected_accessory))
  # plots
  p1 <- ggplot(density_wide, aes(x = window_start, y = log2_core_enrichment)) +
    geom_line() + geom_hline(yintercept=0, linetype="dashed") +
    facet_wrap(~chr, ncol=1, scales="free") + theme_bw() + labs(title=paste0(sp," Core Enrichment"))
  ggsave(paste0(sp,"_core.png"), p1, width=5, height=10, dpi=800)
  p2 <- ggplot(density_wide, aes(x = window_start, y = log2_accessory_enrichment)) +
    geom_line() + geom_hline(yintercept=0, linetype="dashed") +
    facet_wrap(~chr, ncol=1, scales="free") + theme_bw() + labs(title=paste0(sp," Accessory Enrichment"))
  ggsave(paste0(sp,"_accessory.png"), p2, width=5, height=10, dpi=800)
  p3 <- ggplot(density_wide, aes(x = window_start)) +
    geom_line(aes(y = log2_core_enrichment, color = "Core")) +
    geom_line(aes(y = log2_accessory_enrichment, color = "Accessory")) +
    facet_wrap(~ chr, ncol = 1, scales = "free") +
    geom_hline(yintercept = 0, linetype = "dashed") + theme_bw() +
    labs(y = "log2 enrichment", color = "Type", title = paste0(sp," Core vs Accessory"))
  ggsave(paste0(sp,"_total.png"), p3, width=5, height=10, dpi=800)

  # Singleton plot if present
  if("Singleton" %in% unique(df_types$type)) {
    if("Singleton" %in% colnames(density_wide)) {
      density_wide <- density_wide %>% mutate(expected_singleton = count_total * (sum(Singleton, na.rm=TRUE) / sum(count_total)),
                                              log2_singleton_enrichment = log2((Singleton + 1e-6)/ expected_singleton))
      p4 <- ggplot(density_wide, aes(x = window_start, y = log2_singleton_enrichment)) +
        geom_line() + geom_hline(yintercept=0, linetype="dashed") +
        facet_wrap(~chr, ncol=1, scales="free") + theme_bw() + labs(title=paste0(sp," Singleton Enrichment"))
      ggsave(paste0(sp,"_singleton.png"), p4, width=5, height=10, dpi=800)
    }
  }
}

# Contig only:
# Processing Arabidopsis_lyrata
# Processing Cochlearia_groenlandica
# Processing Fagopyrum_escelentum_H2
# Processing Pyrus_bretschneideri

# Worked!





















#############################

library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(IRanges)

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
# Filter out scaffolds < 1 Mbp
chr_lengths <- df %>%
  group_by(species, chr) %>%
  summarise(chr_len = max(end), .groups = "drop") %>%
  filter(chr_len >= 1e6)

df <- df %>%
  semi_join(chr_lengths, by = c("species", "chr"))

#---------------------------
# Read all genes BED (for total counts)
all_genes <- bind_rows(
  read_bed("Dryas_octopetala_all_genes.bed", "Dryas_octopetala", "ALL"),
  read_bed("Oxyria_digyna_all_genes.bed", "Oxyria_digyna", "ALL"),
  read_bed("Draba_nivalis_all_genes.bed", "Draba_nivalis", "ALL")
)

all_genes <- all_genes %>%
  semi_join(chr_lengths, by = c("species", "chr"))

#---------------------------
# Parameters
window_size <- 3e6  # 0.5 Mbp

#------------------------------
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

#---------------------------------
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

#-------------------------------------
run_species_pipeline <- function(species_name, df, window_size, coreacc_window_counts) {

  df_sp <- df %>% filter(species == species_name)
  all_genes_sp <- all_genes %>% filter(species == species_name)

  # ----------------------------
  # Core/Accessory counts per window
  coreacc_window_counts <- compute_window_counts(df_sp, window_size)

  # Total counts per window
  total_window_counts <- compute_window_counts_total(all_genes_sp, window_size)

  # ----------------------------
  # merge core/accessory + total
  density_counts_wide <- coreacc_window_counts %>%
    filter(species == species_name) %>%
    pivot_wider(
      names_from = type,
      values_from = count,
      values_fill = 0
    ) %>%
    left_join(total_window_counts,
              by = c("species", "chr", "window_start", "window_end"))

  # ----------------------------
  # proportions
  p_core <- sum(density_counts_wide$Core) / sum(density_counts_wide$count_total)
  p_accessory <- sum(density_counts_wide$Accessory) / sum(density_counts_wide$count_total)

  density_counts_wide <- density_counts_wide %>%
    mutate(
      expected_core = count_total * p_core,
      expected_accessory = count_total * p_accessory,
      log2_core_enrichment = log2((Core + 1e-6) / expected_core),
      log2_accessory_enrichment = log2((Accessory + 1e-6) / expected_accessory)
    )

  # ----------------------------
  # plots
  p1 <- ggplot(density_counts_wide, aes(x = window_start, y = log2_core_enrichment)) +
    geom_line() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_wrap(~ chr, ncol=1, scales = "free") +
    theme_bw() +
    labs(title = paste0(species_name, " Core Enrichment"))

  ggsave(paste0(species_name, "_core.png"), p1, width = 5, height = 10, dpi = 800)

  p2 <- ggplot(density_counts_wide, aes(x = window_start, y = log2_accessory_enrichment)) +
    geom_line() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_wrap(~ chr, ncol=1, scales = "free") +
    theme_bw() +
    labs(title = paste0(species_name, " Accessory Enrichment"))

  ggsave(paste0(species_name, "_accessory.png"), p2, width = 5, height = 10, dpi = 800)

p3 <- ggplot(density_counts_wide, aes(x = window_start)) +
  geom_line(aes(y = log2_core_enrichment, color = "Core")) +
  geom_line(aes(y = log2_accessory_enrichment, color = "Accessory")) +
  facet_wrap(~ chr, ncol = 1, scales = "free") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_bw() +
  labs(
    y = "log2 enrichment",
    color = "Type"
  )
  
  ggsave(paste0(species_name, "_total.png"), p3, width = 5, height = 10, dpi = 800)

  return(density_counts_wide)
}


species_list <- unique(df$species)

results <- lapply(
  species_list,
  function(sp) run_species_pipeline(sp, df, window_size, coreacc_window_counts)
)

names(results) <- species_list


######################################































####################################
# Older
#---------------------------

# clean up duplicates - take longest gene for each ID
awk '
{
  id=$4
  gene=id
  sub(/-.*$/, "", gene)

  len=$3-$2

  if(!(gene in best) || len > best_len[gene]){
    best[gene]=$0
    best_len[gene]=len
  }
}
END {
  for(g in best) print best[g]
}' all_genes.bed > all_genes.best_isoform.bed

wc -l all_genes*
  # 651482 all_genes.bed
  # 617694 all_genes.best_isoform.bed

#---------------------------
# Extract gene positions

awk '
NR==FNR {ids[$1]=1; next}
{
  id=$4
  sub(/-R[A-Z]$/, "", id)

  if(id in ids) print
}' core_gene_ids_sorted.txt all_genes.bed > core_genes.bed

awk '
NR==FNR {ids[$1]=1; next}
{
  id=$4

  sub(/-R[A-Z]$/, "", id)

  if(id in ids) print
}'  accessory_gene_ids_sorted.txt all_genes.bed > accessory_genes.bed

awk '
NR==FNR {ids[$1]=1; next}
{
  id=$4

  sub(/-R[A-Z]$/, "", id)

  if(id in ids) print
}' singleton_gene_ids_sorted.txt all_genes.bed > singleton_genes.bed

awk '
NR==FNR {ids[$1]=1; next}
{
  id=$4

  sub(/-R[A-Z]$/, "", id)

  if(id in ids) print
}' twospp_gene_ids_sorted.txt all_genes.bed > twospp_genes.bed

wc -l *_genes.bed
  # 69351 accessory_genes.bed
  # 651482 all_genes.bed
  # 161788 core_genes.bed
  # 25424 twospp_genes.bed
  # 43927 singleton_genes.bed

#-------------------------
# split by species

# Core
awk '
$4 ~ /DoctH0/       {print > "Dryas_octopetala_core.bed"}
$4 ~ /Oxyria/       {print > "Oxyria_digyna_core.bed"}
$4 ~ /g[0-9]+\.t1$/ {print > "Cochlearia_groenlandica_core.bed"}
$4 ~ /(snap|maker)/  {print > "Draba_nivalis_core.bed"}
' core_genes.bed

wc -l *_core.bed
   # 7401 Cochlearia_groenlandica_core.bed
   # 7860 Draba_nivalis_core.bed
   # 6522 Dryas_octopetala_core.bed
   # 7457 Oxyria_digyna_core.bed


# old
  # 2349 Cochlearia_groenlandica_core.bed
  # 2305 Draba_nivalis_core.bed
  # 1873 Dryas_octopetala_core.bed
  # 2450 Oxyria_digyna_core.bed
  # 8977 total

#----
# Accessory

awk '
$4 ~ /DoctH0/                {print > "Dryas_octopetala_accessory.bed"}
$4 ~ /Oxyria/                {print > "Oxyria_digyna_accessory.bed"}
$4 ~ /g[0-9]+\.t1$/           {print > "Cochlearia_groenlandica_accessory.bed"}
$4 ~ /(snap|maker)/            {print > "Draba_nivalis_accessory.bed"}
' accessory_genes.bed

 wc -l *_accessory.bed

   # 3665 Cochlearia_groenlandica_accessory.bed
   # 3164 Draba_nivalis_accessory.bed
  # 10995 Dryas_octopetala_accessory.bed
   # 5439 Oxyria_digyna_accessory.bed

 # old
  # 1604 Cochlearia_groenlandica_accessory.bed
  # 1164 Draba_nivalis_accessory.bed
  # 6905 Dryas_octopetala_accessory.bed
  # 2242 Oxyria_digyna_accessory.bed

#----
# Singleton

awk '
$4 ~ /DoctH0/                {print > "Dryas_octopetala_singleton.bed"}
$4 ~ /Oxyria/                {print > "Oxyria_digyna_singleton.bed"}
$4 ~ /g[0-9]+\.t1$/           {print > "Cochlearia_groenlandica_singleton.bed"}
$4 ~ /(snap|maker)/           {print > "Draba_nivalis_singleton.bed"}
' singleton_genes.bed

wc -l *_singleton.bed

  # 2848 Cochlearia_groenlandica_singleton.bed
  # 2352 Draba_nivalis_singleton.bed
  # 8927 Dryas_octopetala_singleton.bed
  # 3266 Oxyria_digyna_singleton.bed

#----
# Twospp

awk '
$4 ~ /DoctH0/                {print > "Dryas_octopetala_twospp.bed"}
$4 ~ /Oxyria/                {print > "Oxyria_digyna_twospp.bed"}
$4 ~ /g[0-9]+\.t1$/           {print > "Cochlearia_groenlandica_twospp.bed"}
$4 ~ /(snap|maker)/           {print > "Draba_nivalis_twospp.bed"}
' twospp_genes.bed

wc -l *_twospp.bed
   # 817 Cochlearia_groenlandica_twospp.bed
   # 812 Draba_nivalis_twospp.bed
  # 2068 Dryas_octopetala_twospp.bed
  # 2173 Oxyria_digyna_twospp.bed

#-------------------
# Split total genes by species

awk '
$4 ~ /DoctH0/                {print > "Dryas_octopetala_all_genes.bed"}
$4 ~ /Oxyria/                {print > "Oxyria_digyna_all_genes.bed"}
$4 ~ /g[0-9]+\.t1$/           {print > "Cochlearia_groenlandica_all_genes.bed"}
$4 ~ /(snap|maker)/           {print > "Draba_nivalis_all_genes.bed"}
' all_genes.bed

wc -l *all_genes.bed
   # 29675 Cochlearia_groenlandica_all_genes.bed
   # 33557 Draba_nivalis_all_genes.bed
   # 39696 Dryas_octopetala_all_genes.bed
   # 33799 Oxyria_digyna_all_genes.bed


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
window_size <- 2e6  # 1 Mbp windows


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
window_size <- 2e6  # 0.5 Mbp

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
         width = 5, height = 10, dpi = 800)
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
# Newer


library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(IRanges)

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
# Filter out scaffolds < 1 Mbp
chr_lengths <- df %>%
  group_by(species, chr) %>%
  summarise(chr_len = max(end), .groups = "drop") %>%
  filter(chr_len >= 1e6)

df <- df %>%
  semi_join(chr_lengths, by = c("species", "chr"))

#---------------------------
# Read all genes BED (for total counts)
all_genes <- bind_rows(
  read_bed("Dryas_octopetala_all_genes.bed", "Dryas_octopetala", "ALL"),
  read_bed("Oxyria_digyna_all_genes.bed", "Oxyria_digyna", "ALL"),
  read_bed("Draba_nivalis_all_genes.bed", "Draba_nivalis", "ALL")
)

all_genes <- all_genes %>%
  semi_join(chr_lengths, by = c("species", "chr"))

#---------------------------
# Parameters
window_size <- 3e6  # 0.5 Mbp

#------------------------------
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

#---------------------------------
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

#-------------------------------------
run_species_pipeline <- function(species_name, df, window_size, coreacc_window_counts) {

  df_sp <- df %>% filter(species == species_name)
  all_genes_sp <- all_genes %>% filter(species == species_name)

  # ----------------------------
  # Core/Accessory counts per window
  coreacc_window_counts <- compute_window_counts(df_sp, window_size)

  # Total counts per window
  total_window_counts <- compute_window_counts_total(all_genes_sp, window_size)

  # ----------------------------
  # merge core/accessory + total
  density_counts_wide <- coreacc_window_counts %>%
    filter(species == species_name) %>%
    pivot_wider(
      names_from = type,
      values_from = count,
      values_fill = 0
    ) %>%
    left_join(total_window_counts,
              by = c("species", "chr", "window_start", "window_end"))

  # ----------------------------
  # proportions
  p_core <- sum(density_counts_wide$Core) / sum(density_counts_wide$count_total)
  p_accessory <- sum(density_counts_wide$Accessory) / sum(density_counts_wide$count_total)

  density_counts_wide <- density_counts_wide %>%
    mutate(
      expected_core = count_total * p_core,
      expected_accessory = count_total * p_accessory,
      log2_core_enrichment = log2((Core + 1e-6) / expected_core),
      log2_accessory_enrichment = log2((Accessory + 1e-6) / expected_accessory)
    )

  # ----------------------------
  # plots
  p1 <- ggplot(density_counts_wide, aes(x = window_start, y = log2_core_enrichment)) +
    geom_line() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_wrap(~ chr, ncol=1, scales = "free") +
    theme_bw() +
    labs(title = paste0(species_name, " Core Enrichment"))

  ggsave(paste0(species_name, "_core.png"), p1, width = 5, height = 10, dpi = 800)

  p2 <- ggplot(density_counts_wide, aes(x = window_start, y = log2_accessory_enrichment)) +
    geom_line() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_wrap(~ chr, ncol=1, scales = "free") +
    theme_bw() +
    labs(title = paste0(species_name, " Accessory Enrichment"))

  ggsave(paste0(species_name, "_accessory.png"), p2, width = 5, height = 10, dpi = 800)

p3 <- ggplot(density_counts_wide, aes(x = window_start)) +
  geom_line(aes(y = log2_core_enrichment, color = "Core")) +
  geom_line(aes(y = log2_accessory_enrichment, color = "Accessory")) +
  facet_wrap(~ chr, ncol = 1, scales = "free") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_bw() +
  labs(
    y = "log2 enrichment",
    color = "Type"
  )
  
  ggsave(paste0(species_name, "_total.png"), p3, width = 5, height = 10, dpi = 800)

  return(density_counts_wide)
}


species_list <- unique(df$species)

results <- lapply(
  species_list,
  function(sp) run_species_pipeline(sp, df, window_size, coreacc_window_counts)
)

names(results) <- species_list