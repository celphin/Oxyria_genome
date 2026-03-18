##################################
# From the Arctic and non-Arctic plant RepeatOBserver output
# Predict one or two positions per chromosome based on the longest region
# Compare length of repeat regions
# March 2026
###############################

# First fixing initial ranges to be more complete for length

# Copy all RepeatOBserver results to Narval
cd /home/celphin/scratch/Arctic_centromere_lengths/raw_data/Summary_output

ls *_total_possible_range.txt
# Argentinaanserina_H0-AT_total_possible_range.txt
# Brassica_H0-AT_total_possible_range.txt
# COLCEN_H0-AT_total_possible_range.txt
# Draniv2_H0-AT_total_possible_range.txt
# DryOcto_H0-AT_total_possible_range.txt
# Fagoesc_H1-AT_total_possible_range.txt
# Fagoesc_H2-AT_total_possible_range.txt
# Fagotat_H1-AT_total_possible_range.txt
# Fagotat_H2-AT_total_possible_range.txt
# Fragariavesca_H0-AT_total_possible_range.txt
# Malussylvestris_H0-AT_total_possible_range.txt
# MN47_H0-AT_total_possible_range.txt
# Oxydig_H1-AT_total_possible_range.txt
# Polavi_H0-AT_total_possible_range.txt
# Prunuspersica_H0-AT_total_possible_range.txt
# Rhunob_H0-AT_total_possible_range.txt
# Rhutan_H0-AT_total_possible_range.txt
# Rubusidaeus_H0-AT_total_possible_range.txt
# Strawberry_H0-AT_total_possible_range.txt


#################################
# How to get best positions and best length estimates?

cd /home/celphin/scratch/Arctic_centromere_lengths/raw_data

# narval2
tmux new-session -s R
tmux attach-session -t R


module load StdEnv/2023 r/4.5.0

R

############################################################
# Libraries
############################################################

library(RepeatOBserverV1)
library(zoo)
library(dplyr)
library(stringr)
library(ChemoSpecUtils)

############################################################
# Helper functions
############################################################

read_two_col <- function(files, value_name) {

  bind_rows(lapply(files, function(f) {

    df <- read.table(f, header = FALSE)

    if (ncol(df) != 2) {
      stop(paste("File does not have 2 columns:", f))
    }

    colnames(df) <- c("Genome_position", value_name)

    chr <- str_extract(basename(f), "Chr[0-9]+")

    if (is.na(chr)) {
      stop(paste("Cannot extract chromosome from filename:", f))
    }

    df$Chr <- chr
    df
  }))
}


select_peaks <- function(values, n = 3, decreasing = FALSE, min_dist = 1e6) {

  ord <- order(values, decreasing = decreasing)
  selected <- c()

  for (idx in ord) {

    if (length(selected) == 0) {
      selected <- c(selected, idx)
    } else {

      if (all(abs(idx - selected) > min_dist)) {
        selected <- c(selected, idx)
      }
    }

    if (length(selected) == n) break
  }

  selected
}

############################################################
# Centromere detection
############################################################

calculate_centromeres <- function(fname,
                                  outpath,
                                  rep_bin,
                                  shannon_bin) {

  message("====================================")
  message("Starting genome: ", fname)
  message("====================================")

  summary_path <- file.path(outpath, "Summary_output", "output_data")
  file_list <- list.files(summary_path, full.names = TRUE)

  message("Searching files...")

  rep_files <- file_list[
    grepl("Repeat_abundance_sum", file_list) &
    grepl(fname, file_list)
  ]

  sh_files <- file_list[
    grepl("Shannon_div", file_list) &
    grepl(fname, file_list)
  ]

  message("Repeat files found: ", length(rep_files))
  message("Shannon files found: ", length(sh_files))

  if (length(rep_files) == 0 || length(sh_files) == 0) {
    warning("Missing files for: ", fname)
    return(NULL)
  }

  message("Loading data...")

  Repeat_total  <- read_two_col(rep_files, "Repeat")
  Shannon_total <- read_two_col(sh_files, "Shannon")

  Repeat_total$Genome_position  <- as.numeric(Repeat_total$Genome_position)
  Repeat_total$Repeat           <- as.numeric(Repeat_total$Repeat)

  Shannon_total$Genome_position <- as.numeric(Shannon_total$Genome_position)
  Shannon_total$Shannon         <- as.numeric(Shannon_total$Shannon)

  chromosomes <- intersect(unique(Repeat_total$Chr),
                           unique(Shannon_total$Chr))

  message("Chromosomes to process: ", length(chromosomes))

  final_ranges <- list()

############################################################
# Chromosome loop
############################################################

  for (i in seq_along(chromosomes)) {

    chr <- chromosomes[i]

    message("Processing ", chr,
            " (", i, "/", length(chromosomes), ")")

	edge_buffer <- 2e6  # 2 Mbp

	rep_chr <- Repeat_total %>%
	  filter(Chr == chr,
			 Genome_position > edge_buffer,
			 Genome_position < max(Genome_position) - edge_buffer) %>%
	  arrange(Genome_position)

	sh_chr <- Shannon_total %>%
	  filter(Chr == chr,
			 Genome_position > edge_buffer,
			 Genome_position < max(Genome_position) - edge_buffer) %>%
	  arrange(Genome_position)

	if (nrow(rep_chr) < rep_bin || nrow(sh_chr) < shannon_bin) {
	  message("Skipping ", chr, " (insufficient data after edge trimming)")
	  next
	}

############################################################
# Rolling means
############################################################

    rep_chr$roll <- zoo::rollmean(rep_chr$Repeat,
                                  k = rep_bin,
                                  fill = NA,
                                  align = "center")

    sh_chr$roll <- zoo::rollmean(sh_chr$Shannon,
                                 k = shannon_bin,
                                 fill = NA,
                                 align = "center")

    rep_chr <- rep_chr %>% filter(!is.na(roll))
    sh_chr  <- sh_chr  %>% filter(!is.na(roll))

############################################################
# Find top 3 min values - Shannon rolling mean
############################################################

    sh_mean <- mean(sh_chr$roll, trim = 0.1, na.rm = TRUE)
    sh_extreme_min <- min(sh_chr$roll, na.rm = TRUE)
    sh_thresh <- sh_mean - 0.5 * (sh_mean - sh_extreme_min)

    sh_min_idx <- select_peaks(sh_chr$roll, n = 3, decreasing = FALSE)

    sh_min_pos <- sh_chr$Genome_position[sh_min_idx]

############################################################
# Find top 3 min values - Repeat Abundance rolling mean
############################################################

    rep_mean <- mean(rep_chr$roll, trim = 0.1, na.rm = TRUE)
    rep_extreme_min <- min(rep_chr$roll, na.rm = TRUE)
    rep_thresh_min <- rep_mean - 0.5 * (rep_mean - rep_extreme_min)

    rep_min_idx <- select_peaks(rep_chr$roll, n = 3, decreasing = FALSE)

    rep_min_pos <- rep_chr$Genome_position[rep_min_idx]

############################################################
# Find top 3 max values  - Repeat Abundance rolling mean
############################################################

    rep_extreme_max <- max(rep_chr$roll, na.rm = TRUE)
    rep_thresh_max <- rep_mean + 0.5 * (rep_extreme_max - rep_mean)

    rep_max_idx <- select_peaks(rep_chr$roll, n = 3, decreasing = TRUE)

    rep_max_pos <- rep_chr$Genome_position[rep_max_idx]

############################################################
# Shannon diversity ranges
# Find the range of genome positions from each min that stay below the mean 
# Once a position has the same value as the mean no longer include in the range
############################################################

	sh_ranges <- lapply(sh_min_idx, function(idx) {

	  left <- idx
	  right <- idx

	while (left > 1 && sh_chr$roll[left] < sh_thresh) {
	  left <- left - 1
	}
	while (right < nrow(sh_chr) && sh_chr$roll[right] < sh_thresh) {
	  right <- right + 1
	}

	  start <- sh_chr$Genome_position[left + 1]
	  end   <- sh_chr$Genome_position[right - 1]

	  data.frame(
		Species = fname,
		Chr = chr,
		Start = start,
		End = end,
		Length = end - start,
		Type = "Shannon",
		Dist = sh_mean - sh_chr$roll[idx]
	  )
	})


############################################################
# Repeat abundance ranges 
# Find the range of genome positions from each min and max that stay below/above the mean 
# Once a position has the same value as the mean no longer include in the range
############################################################

	rep_ranges_min <- lapply(rep_min_idx, function(idx) {

	  left <- idx
	  right <- idx

      while (left > 1 && rep_chr$roll[left] < rep_thresh_min) { left <- left - 1 }
      while (right < nrow(rep_chr) && rep_chr$roll[right] < rep_thresh_min) { right <- right + 1 }

	  start <- rep_chr$Genome_position[left + 1]
	  end   <- rep_chr$Genome_position[right - 1]

	  data.frame(
		Species = fname,
		Chr = chr,
		Start = start,
		End = end,
		Length = end - start,
		Type = "RepMin",
		Dist = rep_mean - rep_chr$roll[idx]
	  )
	})


	rep_ranges_max <- lapply(rep_max_idx, function(idx) {

	  left <- idx
	  right <- idx

      while (left > 1 && rep_chr$roll[left] > rep_thresh_max) { left <- left - 1 }
      while (right < nrow(rep_chr) && rep_chr$roll[right] > rep_thresh_max) { right <- right + 1 }

	  start <- rep_chr$Genome_position[left + 1]
	  end   <- rep_chr$Genome_position[right - 1]

	  data.frame(
		Species = fname,
		Chr = chr,
		Start = start,
		End = end,
		Length = end - start,
		Type = "RepMax",
		Dist = rep_chr$roll[idx] - rep_mean
	  )
	})


############################################################
# Find the longest range for each chromosome
# Filter for the type furthest from the mean and then filter by length within that type
# Record the Spp, Chr, Start, End,  Length, Type: Shannon or RepMin or RepMax
############################################################

	chr_ranges <- do.call(
	  rbind,
	  c(sh_ranges, rep_ranges_min, rep_ranges_max)
	)

	if (nrow(chr_ranges) > 0) {

	  best_type <- chr_ranges %>%
		group_by(Type) %>%
		summarise(Dist = max(Dist)) %>%
		arrange(desc(Dist)) %>%
		slice(1) %>%
		pull(Type)

	type_candidates <- chr_ranges %>%
	  filter(Type == best_type)

	# normalize within the chromosome/type
	if(nrow(type_candidates) > 1){
	  type_candidates <- type_candidates %>%
		mutate(
		  Dist_norm   = (Dist - min(Dist)) / (max(Dist) - min(Dist)),
		  Length_norm = (Length - min(Length)) / (max(Length) - min(Length)),
		  Score = 0.7*Dist_norm + 0.3*Length_norm
		)
	} else {
	  type_candidates <- type_candidates %>%
		mutate(
		  Dist_norm = 1,
		  Length_norm = 1,
		  Score = 1
		)
	}

	longest <- type_candidates %>%
	  arrange(desc(Score)) %>%
	  slice(1)
	  
	  final_ranges[[length(final_ranges) + 1]] <- longest
	}


}

############################################################
# Output results
############################################################

  if (length(final_ranges) > 0) {

    final_ranges <- do.call(rbind, final_ranges)

  } else {

    final_ranges <- data.frame(
      Species = character(),
      Chr = character(),
      Start = numeric(),
      End = numeric(),
      Length = numeric(),
      Type = character()
    )
  }

  out_file <- file.path(outpath,
                        paste0(fname,
                               "_Centromere_candidates.txt"))

  message("Writing results to: ", out_file)
  message("Number of candidate regions: ", nrow(final_ranges))

  write.table(final_ranges,
              file = out_file,
              row.names = FALSE,
              col.names = TRUE,
              sep = "\t")

  message("Finished genome: ", fname)

  return(final_ranges)
}

############################################################
# Run all genomes
############################################################

outpath <- "/home/celphin/scratch/Arctic_centromere_lengths/raw_data"

fnames <- c(
"Argentinaanserina_H0-AT",
"Brassica_H0-AT",
"COLCEN_H0-AT",
"Draniv2_H0-AT",
"DryOcto_H0-AT",
"Fagoesc_H1-AT",
"Fagoesc_H2-AT",
"Fagotat_H1-AT",
"Fagotat_H2-AT",
"Malussylvestris_H0-AT",
"MN47_H0-AT",
"Oxydig_H1-AT",
"Polavi_H0-AT",
"Prunuspersica_H0-AT",
"Rhunob_H0-AT",
"Rhutan_H0-AT",
"Rubusidaeus_H0-AT",
"Strawberry_H0-AT"
)

results_all <- lapply(fnames, function(x)
  calculate_centromeres(x, outpath, rep_bin = 100, shannon_bin = 100))

results_all <- do.call(rbind, results_all)


#########################
# Get chromosome lengths

extract_chr_lengths <- function(fname, outpath) {

  summary_path <- file.path(outpath, "Summary_output", "output_data")
  file_list <- list.files(summary_path, full.names = TRUE)

  rep_files <- file_list[
    grepl("Repeat_abundance_sum", file_list) &
    grepl(fname, file_list)
  ]

  sh_files <- file_list[
    grepl("Shannon_div", file_list) &
    grepl(fname, file_list)
  ]

  if (length(rep_files) == 0 & length(sh_files) == 0) {
    warning("No files found for: ", fname)
    return(NULL)
  }

  # Prefer Shannon but fall back to repeat if needed
  if (length(sh_files) > 0) {
    df <- read_two_col(sh_files, "Shannon")
  } else {
    df <- read_two_col(rep_files, "Repeat")
  }

  df$Genome_position <- as.numeric(df$Genome_position)

  chr_lengths <- df %>%
    dplyr::group_by(Chr) %>%
    dplyr::summarise(Chr_length = max(Genome_position, na.rm = TRUE)) %>%
    dplyr::mutate(Genome = fname) %>%
    dplyr::ungroup()

  return(chr_lengths)
}

chr_lengths_all <- lapply(fnames, function(x)
  extract_chr_lengths(x, outpath)
)

chr_lengths_all <- dplyr::bind_rows(chr_lengths_all)

head(chr_lengths_all)

colnames(chr_lengths_all) <- c("Chr", "ChrLength", "Spp")

##########################
# Look at data with plots for all spp

##########################
# copy *_total_possible_range.txt to one folder

# Narval
cd /home/celphin/scratch/Arctic_centromere_lengths/raw_data
cat *_Centromere_candidates.txt > Arctic_genomes_total_possible_range.txt

wc -l Arctic_genomes_total_possible_range.txt
# 2624 Arctic_genomes_total_possible_range.txt

# Copy to local machine
# ~/Github/Oxyria_genome/3_Genome_Comparisons/6_RepeatOBserver

#############################
# read in cent_pred

tmux attach-session -t R

fname="Arctic_genomes"
outpath <- "/home/celphin/scratch/Arctic_centromere_lengths/raw_data"

# Wrap in function
#automated_centromere_detection <- function(fname=fname,  outpath=outpath) {

# read in the data
df0 <- utils::read.table(paste0(outpath,"/",fname,"_total_possible_range.txt"), sep = "\t", header=TRUE, check.names = FALSE)

df <- df0[-which(df0$Species=="Species"),]

df_sorted <- df

colnames(df_sorted) <-  c("Spp", "Chr","Start", "End","Length",
"Type","Dist", "Dist_norm", "Length_norm", "Score")


unique(df_sorted$Spp)

arctic_spp <- c("Draniv2_H0-AT", "DryOcto_H0-AT", "Dryoct_H0-AT", "Oxydig_H1-AT")   
alpine_spp <- c("Rhunob_H0-AT")  

df2 <- df_sorted %>%
  mutate(
    Habitat = case_when(
      Spp %in% arctic_spp ~ "Arctic",
      Spp %in% alpine_spp ~ "Alpine",
      TRUE ~ "Other"
    )
  )

#-----------------------
# add in plant family

unique(df2$Spp)
#[1] Arabidopsis_H0-AT     Brassica_H0-AT        COLCEN_H0-AT          Draniv2_H0-AT         DryOcto_H0-AT       
#[7] Fagoesc_H1-AT         Fagoesc_H2-AT         Fagotat_H1-AT         MN47_H0-AT            Malussylvestris_H0-AT Oxydig_H1-AT
#[13] Polavi_H0-AT          Prunuspersica_H0-AT   Rhunob_H0-AT          Rhutan_H0-AT          Rubusidaeus_H0-AT     Strawberry_H0-AT

family_lookup <- data.frame(
  Spp = c("Argentinaanserina_H0-AT",
    "Brassica_H0-AT",
    "COLCEN_H0-AT",
    "Draniv2_H0-AT",
    "DryOcto_H0-AT",
    "Fagoesc_H1-AT",
    "Fagoesc_H2-AT",
    "Fagotat_H1-AT",
    "Fagotat_H2-AT",
    "MN47_H0-AT",
    "Malussylvestris_H0-AT",
    "Oxydig_H1-AT",
    "Polavi_H0-AT",
    "Prunuspersica_H0-AT",
    "Rhunob_H0-AT",
    "Rhutan_H0-AT",
    "Rubusidaeus_H0-AT",
    "Strawberry_H0-AT"
  ),
  Family = c(
    "Brassicaceae",
    "Brassicaceae",
    "Brassicaceae",
    "Brassicaceae",
    "Rosaceae",
    "Polygonaceae",
    "Polygonaceae",
    "Polygonaceae",
    "Polygonaceae",
    "Brassicaceae",
    "Rosaceae",
    "Polygonaceae",
    "Polygonaceae",
    "Rosaceae",
    "Polygonaceae",
    "Polygonaceae",
    "Rosaceae",
    "Rosaceae"
  )
)

df2 <- df2 %>%
  left_join(family_lookup, by = "Spp")

df2$Family <- as.factor(df2$Family)


#################################
# calculate the average centromere length per species

df2$Length <- as.numeric(df2$Length)

avg_length <- df2 %>%
  group_by(Spp) %>%
  summarise(mean_length = mean(Length))

avg_length_habitat <- df2 %>%
  group_by(Habitat) %>%
  summarise(mean_length = mean(Length))

df2 <- df2 %>%
  mutate(
    HabitatGroup = ifelse(Habitat %in% c("Arctic","Alpine"),
                          "Arctic/Alpine",
                          "Other")
  )


write.csv(df2, paste0(outpath,"/",fname,"_final_cent_pred_reptypes_habitat.csv"), row.names = FALSE)


#############################
# add in ChrLengths

df3 <- merge(
  df2,
  chr_lengths_all,
  by = c("Spp", "Chr"),
  all.x = TRUE
)

# relative centromere sizes

df3 <- df3 %>%
  dplyr::mutate(RelLength = Length / ChrLength)

df3$Spp <- with(df3, reorder(Spp, RelLength, FUN = mean))

df2$Spp <- with(df2, reorder(Spp, Length, FUN = mean))


df2_filtered <- df2 %>%
  filter(!Spp %in% c("Fagoesc_H2-AT", "Fagotat_H2-AT"))

df3_filtered <- df3 %>%
  filter(!Spp %in% c("Fagoesc_H2-AT", "Fagotat_H2-AT"))


########################
# Plot the data
library(ggplot2)

# Looking at Repeat lengths

p1 <- ggplot(df2_filtered, aes(x = HabitatGroup, y = Length, fill = HabitatGroup)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.6) +
  labs(
    x = "Habitat group",
    y = "Repeat region length",
    title = "Chromosome repeat lengths: Arctic/Alpine vs Other"
  ) +
  theme_classic()

ggsave("habitat_boxplot.png", plot = p1, width = 6, height = 4, dpi = 300)

#----------------------
# spp patterns

# Order spp by family
df2 <- df2 %>%
  mutate(Spp = reorder(Spp, Length, FUN = mean))

p2 <- ggplot(df2_filtered, aes(x = Spp, y = Length, fill = HabitatGroup)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  labs(
    x = "Species",
    y = "Repeat region length",
    fill = "Group",
    title = "Repeat lengths per species"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("species_repeat_lengths.png", plot = p2, width = 8, height = 5, dpi = 300)

#--------------------------------
# spp centromere sizes split by plant family

p3 <- ggplot(df2_filtered, aes(x = Spp, y = Length, fill = HabitatGroup)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  facet_wrap(~Family, scales = "free_x") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("family_repeat_lengths.png", plot = p3, width = 10, height = 6, dpi = 300)


#------------------------------

p4 <- ggplot(df3_filtered, aes(x = HabitatGroup, y = RelLength, fill = HabitatGroup)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.6) +
  labs(
    x = "Habitat group",
    y = "Relative Repeat region length",
    title = "Chromosome repeat lengths: Arctic/Alpine vs Other"
  ) +
  theme_classic()

ggsave("rel_habitat_boxplot.png", plot = p4, width = 6, height = 4, dpi = 300)


#------------------------------
# relative spp centromere sizes split by plant family
p5 <- ggplot(df3_filtered, aes(x = Spp, y = RelLength, fill = HabitatGroup)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  facet_wrap(~Family, scales = "free_x") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("rel_family_repeat_lengths.png", plot = p5, width = 10, height = 6, dpi = 300)


#------------------------------
# relative spp centromere sizes split by plant family
p6 <- ggplot(df3_filtered, aes(x = Spp, y = RelLength, fill = HabitatGroup)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  facet_wrap(~HabitatGroup, scales = "free_x") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("rel_Habitat_SPP_repeat_lengths.png", plot = p6, width = 10, height = 6, dpi = 300)


##########################
# Test

table(df3$HabitatGroup)

# Arctic/Alpine         Other
           # 35           143


kruskal.test(RelLength ~ HabitatGroup, data = df3_filtered)

# Kruskal-Wallis rank sum test
# data:  RelLength by HabitatGroup
# Kruskal-Wallis chi-squared = 13.55, df = 1, p-value = 0.0002323

#########################
# Correct for spp groupings
library(lme4)

# Linear mixed model
mod <- lmer(RelLength ~ HabitatGroup + (1 | Spp), data = df3_filtered)

# Test significance of HabitatGroup
library(lmerTest)
anova(mod)

# Type III Analysis of Variance Table with Satterthwaite's method
                # Sum Sq   Mean Sq NumDF  DenDF F value Pr(>F)
# HabitatGroup 0.0016764 0.0016764     1 13.657  2.1858  0.162





###################################
# Old version for raw output

# remove header rows
df <- df0[-which(df0$Label=="Label"),]

list(unique(df$Spp))


#-----------------
# merge close regions within 1Mbp
library(dplyr)

# Convert Start and End to numeric for comparison
df$Start <- as.numeric(df$Start)
df$End <- as.numeric(df$End)

# Function to merge overlapping/close regions
merge_regions <- function(df) {
  df$Length <- df$End - df$Start

  # Remove rows with NA values in Start or End
  df <- df[!is.na(df$Start) & !is.na(df$End), ]

  # Check if there's any data left after removing NAs
  if (nrow(df) == 0) {
    return(df)  # Return empty dataframe if there's no data
  }

  # If there's only one region, return it as is
  if (nrow(df) == 1) {
    return(df)
  }

  df <- df[order(df$Start), ]  # Sort by Start position
  merged <- list()
  current <- df[1, ]  # Start with the first row

  for (i in 2:nrow(df)) {
    if (df$Start[i] - current$End <= 2e6) {
      # If the next region is within the max length, merge it
      current$End <- max(current$End, df$End[i])
    } else {
      merged <- append(merged, list(current))
      current <- df[i, ]  # Start a new region
    }
  }
  merged <- append(merged, list(current))  # Append the last region
  do.call(rbind, merged)
}

# Apply merging function by chromosome and repeat type
df_merged0 <- df %>%
  group_by(Spp, Chr, Label) %>%
  do(merge_regions(.)) %>%
  ungroup()

# View the merged data
head(df_merged0)

df_merged <- df_merged0 %>%
  arrange(Spp, Chr)  # Sort by the 'Chr' column

#-----------------------
# Find two longest regions

# Calculate the length of each region
df_merged$Length <- df_merged$End - df_merged$Start

# Find the two longest regions for each chromosome
df_longest <- df_merged %>%
  group_by(Spp, Chr, Label) %>%
  arrange(desc(Length)) %>%
  slice_head(n = 1) %>%
  ungroup()

# View the  longest regions for each chromosome and repeat type
head(df_longest)

#------------------------
# Compare lengths

# Step 1: Sort by Length in descending order, then calculate the PercentLength
df_final <- df_longest %>%
  group_by(Spp, Chr) %>%
  arrange(desc(Length)) %>%  # Sort by Length within each chromosome (largest first)
  mutate(
    # If there are at least two regions, calculate PercentLength as the ratio between second and first largest
    PercentLength = ifelse(n() > 1, Length[2] / Length[1], 1)  # Default to 1 if there's only one region
  ) %>%
  ungroup()

# Step 1: Sort by Length in descending order, then filter based on PercentLength and Length
df_filtered <- df_final %>%
  group_by(Spp, Chr) %>%
  arrange(desc(Length)) %>%  # Sort by Length in descending order within each chromosome
  mutate(keep_rows = ifelse(PercentLength > 0.85, row_number() <= 1, row_number() == 1)) %>%
  ungroup() %>%
  filter(keep_rows) %>%
  select(-keep_rows)  # Remove the temporary column used for filtering

# View the final filtered data
head(df_filtered)


#----------------------
# add repeat type column

df_transformed <- df_filtered %>%
  mutate(RepeatType = case_when(
    Label == "Shannon" ~ "TandemRepeat",
    Label == "MinRepeatAbund" ~ "TandemRepeat",
    Label == "MaxRepeatAbund" ~ "Retrotransposon",
    TRUE ~ Label  # In case there are any other labels you want to keep as is
  ))

# View the transformed data
head(df_transformed)

#-------------------------
# Write out file

# Sort the data by 'Chr' column
df_sorted <- df_transformed %>%
  arrange(Spp, Chr)  # Sort by the 'Chr' column


# Write the final dataframe to a CSV file
write.csv(df_sorted, paste0(outpath,"/",fname,"_final_cent_pred_reptypes.csv"), row.names = FALSE)

#}

df_sorted

########################
# Identify habitat types in dataset

unique(df_sorted$Spp)

arctic_spp <- c("Draniv2_H0-AT", "DryOcto_H0-AT", "Dryoct_H0-AT", "Oxydig_H1-AT")   # replace with your Arctic species
alpine_spp <- c("Rhunob_H0-AT")   # replace with Alpine species

df2 <- df_sorted %>%
  mutate(
    Habitat = case_when(
      Spp %in% arctic_spp ~ "Arctic",
      Spp %in% alpine_spp ~ "Alpine",
      TRUE ~ "Other"
    )
  )

#-----------------------
# add in plant family

unique(df2$Spp)
#[1] Arabidopsis_H0-AT     Brassica_H0-AT        COLCEN_H0-AT          Draniv2_H0-AT         DryOcto_H0-AT         Dryoct_H0-AT
#[7] Fagoesc_H1-AT         Fagoesc_H2-AT         Fagotat_H1-AT         MN47_H0-AT            Malussylvestris_H0-AT Oxydig_H1-AT
#[13] Polavi_H0-AT          Prunuspersica_H0-AT   Rhunob_H0-AT          Rhutan_H0-AT          Rubusidaeus_H0-AT     Strawberry_H0-AT

family_lookup <- data.frame(
  Spp = c(
    "Arabidopsis_H0-AT",
    "Brassica_H0-AT",
    "COLCEN_H0-AT",
    "Draniv2_H0-AT",
    "DryOcto_H0-AT",
    "Dryoct_H0-AT",
    "Fagoesc_H1-AT",
    "Fagoesc_H2-AT",
    "Fagotat_H1-AT",
    "MN47_H0-AT",
    "Malussylvestris_H0-AT",
    "Oxydig_H1-AT",
    "Polavi_H0-AT",
    "Prunuspersica_H0-AT",
    "Rhunob_H0-AT",
    "Rhutan_H0-AT",
    "Rubusidaeus_H0-AT",
    "Strawberry_H0-AT"
  ),
  Family = c(
    "Brassicaceae",
    "Brassicaceae",
    "Brassicaceae",
    "Brassicaceae",
    "Rosaceae",
    "Rosaceae",
    "Polygonaceae",
    "Polygonaceae",
    "Polygonaceae",
    "Brassicaceae",
    "Rosaceae",
    "Polygonaceae",
    "Polygonaceae",
    "Rosaceae",
    "Polygonaceae",
    "Polygonaceae",
    "Rosaceae",
    "Rosaceae"
  )
)

df2 <- df2 %>%
  left_join(family_lookup, by = "Spp")

df2$Family <- as.factor(df2$Family)

#################################
# calculate the average centromere length per species

avg_length <- df2 %>%
  group_by(Spp) %>%
  summarise(mean_length = mean(Length))

avg_length_habitat <- df2 %>%
  group_by(Habitat) %>%
  summarise(mean_length = mean(Length))

df2 <- df2 %>%
  mutate(
    HabitatGroup = ifelse(Habitat %in% c("Arctic","Alpine"),
                          "Arctic/Alpine",
                          "Other")
  )


write.csv(df2, paste0(outpath,"/",fname,"_final_cent_pred_reptypes_habitat.csv"), row.names = FALSE)

########################
# Plot the data
library(ggplot2)

# Looking at Repeat lengths

ggplot(df2, aes(x = HabitatGroup, y = Length, fill = HabitatGroup)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.6) +
  labs(
    x = "Habitat group",
    y = "Repeat region length",
    title = "Chromosome repeat lengths: Arctic/Alpine vs Other"
  ) +
  theme_classic()

#----------------------
# spp patterns

# Order spp by family
df2 <- df2 %>%
  mutate(Spp = reorder(Spp, Length, FUN = mean))

ggplot(df2, aes(x = Spp, y = Length, fill = HabitatGroup)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  labs(
    x = "Species",
    y = "Repeat region length",
    fill = "Group",
    title = "Repeat lengths per species"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#--------------------------------
# spp centromere sizes split by plant family

ggplot(df2, aes(x = Spp, y = Length, fill = HabitatGroup)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  facet_wrap(~Family, scales = "free_x") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))






























##########################
# Older version

calculate_ranges <- function(fname=fname, outpath=outpath,
                             Shannon_bin_size=800, Shannon_SD=1.5,
                             RepAbundMax_bin_size=200, RepAbundMax_SD=1.5,
                             RepAbundMin_bin_size=800, RepAbundMin_SD=1.5){

  #-----------------------
  # Read in data for all chromosomes and plot
  summary_path <- paste0(outpath,"/", fname, "/", "Summary_output/output_data")
  file_list <- list.files(summary_path, full.names=TRUE)
  rollsumhist_list <- file_list[grep("Repeat_abundance_sum", file_list)]

  # read in Repeat Abundance sums
  lsd <- lapply(rollsumhist_list, read.table)
  sd_chr_list0 <- basename(rollsumhist_list)
  sd_chr_list1 <- stringr::str_split(sd_chr_list0, "_", simplify =TRUE)
  sd_chr_list2 <- sd_chr_list1[,3]
  names(lsd) <- sd_chr_list2
  RepeatAbundance_total <- dplyr::bind_rows(lsd, .id = 'chromosome')
  colnames(RepeatAbundance_total) <- c("Chromosome", "Genome_position", "RepeatAbundance")

  RepeatAbundance_total_parts <- RepeatAbundance_total[grep("-", RepeatAbundance_total$Chromosome),]

  #-----------------------------
  # join chromosome parts if >400Mbp
  if (nrow(RepeatAbundance_total_parts)>0){
    RepeatAbundance_total$Chrnum0 <- as.factor(stringr::str_split(RepeatAbundance_total$Chromosome, "r", simplify =TRUE)[,2])
    RepeatAbundance_total$Chrnum <- as.factor(stringr::str_split(RepeatAbundance_total$Chrnum0, "-", simplify =TRUE)[,1])
    RepeatAbundance_total$Chrpart <- as.numeric(stringr::str_split(RepeatAbundance_total$Chrnum0, "-", simplify =TRUE)[,2])
    RepeatAbundance_total$Chrpart[which(is.na(RepeatAbundance_total$Chrpart))] <-1
    RepeatAbundance_total$Genome_position <- RepeatAbundance_total$Genome_position + (RepeatAbundance_total$Chrpart -1)*4e8
  }else{
    RepeatAbundance_total$Chrnum <- as.factor(stringr::str_split(RepeatAbundance_total$Chromosome, "r", simplify =TRUE)[,2])
  }

  #--------------------------
  # formatting

  RepeatAbundance_total$Genome_position <- as.numeric(as.character(RepeatAbundance_total$Genome_position))
  RepeatAbundance_total$RepeatAbundance <- as.numeric(as.character(RepeatAbundance_total$RepeatAbundance))
  RepeatAbundance_total$Chrnum <- as.numeric(as.character(RepeatAbundance_total$Chrnum))

  #-------------------------------------
  # Calculate max, min and ranges from roll sum data

  # find start and end of highly repeating regions based 1.5 SD from mean
  RepeatAbund_cent_min <- data.frame()
  RepeatAbund_cent_max <- data.frame()
  RepeatAbund_min <- data.frame()
  RepeatAbund_max <- data.frame()
  RepeatAbund_length <- data.frame()
  RepeatAbundTotal<- data.frame()

  for (chromosome in unique(RepeatAbundance_total$Chrnum)){
    RepeatAbundance_chr <- RepeatAbundance_total[which(RepeatAbundance_total$Chrnum == chromosome),]
    subsetRepeatAbundance <- RepeatAbundance_chr$RepeatAbundance
    subsetRepeatAbundance[c((length(subsetRepeatAbundance)-400):length(subsetRepeatAbundance))] <- NA

    # run rolling sum
    RepeatAbund100 <-  zoo::rollsum(subsetRepeatAbundance, k=RepAbundMax_bin_size,  fill = NA, align = "center")

    # run rolling sum
    RepeatAbund500 <-  zoo::rollsum(subsetRepeatAbundance, k=RepAbundMin_bin_size, fill = NA, align = "center")

    # join sums with genome positions
    Repeat_abund_chr <- cbind(RepeatAbundance_chr, RepeatAbund100, RepeatAbund500)

    #plot
    grDevices::png(filename = paste0(outpath,"/", fname, "/Summary_output/histograms/", fname,"_Chr", chromosome, "_", RepAbundMax_bin_size,"_Repeat_abundance_rollsum_max.png"), width = 1400, height = 700)
    plot(RepeatAbundance_chr$Genome_position, RepeatAbund100)
    grDevices::dev.off()

    #plot
    grDevices::png(filename = paste0(outpath,"/", fname, "/Summary_output/histograms/", fname,"_Chr", chromosome, "_", RepAbundMin_bin_size,"_Repeat_abundance_rollsum_min.png"), width = 1400, height = 700)
    plot(RepeatAbundance_chr$Genome_position, RepeatAbund500)
    grDevices::dev.off()

    # write out total file
    utils::write.table(x=Repeat_abund_chr, file=paste0(outpath,"/", fname, "/Summary_output/output_data/", fname,"_Chr", chromosome, "_Repeat_abundance_rollsum.txt"), sep = "\t", dec = ".",row.names = FALSE, col.names = TRUE)

    Repeat_abund_chr <- as.data.frame(Repeat_abund_chr)

    Repeat_abund_chr$Genome_position <- as.numeric(Repeat_abund_chr$Genome_position)
    Repeat_abund_chr$RepeatAbundance <- as.numeric(Repeat_abund_chr$RepeatAbundance)
    Repeat_abund_chr$RepeatAbund100 <- as.numeric(Repeat_abund_chr$RepeatAbund100)
    Repeat_abund_chr$RepeatAbund500 <- as.numeric(Repeat_abund_chr$RepeatAbund500)

    RepeatAbundTotal <- rbind(RepeatAbundTotal, Repeat_abund_chr)

    #---------------------------------

    # find min position
    cent_min <- Repeat_abund_chr$Genome_position[which(Repeat_abund_chr$RepeatAbund500 == min(Repeat_abund_chr$RepeatAbund500, na.rm=TRUE))]
    cent_max <- Repeat_abund_chr$Genome_position[which(Repeat_abund_chr$RepeatAbund100 == max(Repeat_abund_chr$RepeatAbund100, na.rm=TRUE))]

    # find SD of data
    SD_repeatAbund_up <- sd(Repeat_abund_chr$RepeatAbund100, na.rm=TRUE)
    SD_repeatAbund_low <- sd(Repeat_abund_chr$RepeatAbund500, na.rm=TRUE)

    thres_upper = mean(Repeat_abund_chr$RepeatAbund100, na.rm=TRUE) + (RepAbundMax_SD* SD_repeatAbund_up)
    thres_lower = mean(Repeat_abund_chr$RepeatAbund500, na.rm=TRUE) - (RepAbundMin_SD* SD_repeatAbund_low)

    #-------------------------
    # find positions of - two SD less than mean
    cent_range_wind <- Repeat_abund_chr$Genome_position[which(Repeat_abund_chr$RepeatAbund500 <= thres_lower )]/50000
    cent_range_wind <- Genome_position / 100000
    cent_range_wind_max <- Genome_position / 10000

    # find range of these values
    # ChemoSpecUtils
    # https://rdrr.io/cran/ChemoSpecUtils/man/check4Gaps.html

    #library(ChemoSpecUtils)
    if (length(cent_range_wind)>0) {
      cent_range <- ChemoSpecUtils::check4Gaps(cent_range_wind)
      cent_range[nrow(cent_range)+1,] <- c(0,0,0,0,0)

      cent_range_pos_start <- cent_range[,1]*50000
      cent_range_pos_end <- cent_range[,2]*50000

      SPP_l <- rep(fname, length(cent_range_pos_start))
      Chr_l <- rep(chromosome, length(cent_range_pos_start))

      RepeatAbund_cent_chr_min <- cbind(SPP_l, Chr_l, cent_range_pos_start, cent_range_pos_end)
      RepeatAbund_cent_min <- rbind(RepeatAbund_cent_min, RepeatAbund_cent_chr_min)
    }
    #-------------------------------
    # find positions of + two SD more than mean
    cent_range_wind_max <- Repeat_abund_chr$Genome_position[which(Repeat_abund_chr$RepeatAbund100 >= thres_upper)]/5000

    if (length(cent_range_wind_max)>0) {
      cent_range_max <- ChemoSpecUtils::check4Gaps(cent_range_wind_max)
      cent_range_max[nrow(cent_range_max)+1,] <- c(0,0,0,0,0)

      cent_range_pos_start_max <- cent_range_max[,1]*5000
      cent_range_pos_end_max <- cent_range_max[,2]*5000

      SPP_l_max <- rep(fname, length(cent_range_pos_start_max))
      Chr_l_max <- rep(chromosome, length(cent_range_pos_start_max))

      RepeatAbund_cent_chr_max <- cbind(SPP_l_max, Chr_l_max, cent_range_pos_start_max, cent_range_pos_end_max)
      RepeatAbund_cent_max <- rbind(RepeatAbund_cent_max, RepeatAbund_cent_chr_max)
    }

    #------------------------------

    RepeatAbund_min_chr <- c("MinRepeatAbund", fname, chromosome, cent_min)
    RepeatAbund_min <- rbind(RepeatAbund_min, RepeatAbund_min_chr)

    RepeatAbund_max_chr <- c("MaxRepeatAbund", fname, chromosome, cent_max)
    RepeatAbund_max <- rbind(RepeatAbund_max, RepeatAbund_max_chr)

    chr_length <- max(Repeat_abund_chr$Genome_position)
    RepeatAbund_length_chr <- c("Length", fname, chromosome, chr_length)
    RepeatAbund_length <- rbind(RepeatAbund_length, RepeatAbund_length_chr)

  }

  # remove zeros if files are not empty
  if (length(RepeatAbund_cent_min)>0){
    RepeatAbund_cent_min <- RepeatAbund_cent_min[-which(RepeatAbund_cent_min[,3] == 0),]
    RepeatAbund_cent_min <- as.data.frame(RepeatAbund_cent_min)
    RepeatAbund_cent_min$Label <- rep("MinRepeatAbund", nrow(RepeatAbund_cent_min))
  }
  if(length(RepeatAbund_cent_max)>0){
    RepeatAbund_cent_max <- RepeatAbund_cent_max[-which(RepeatAbund_cent_max[,3] == 0),]
    RepeatAbund_cent_max <- as.data.frame(RepeatAbund_cent_max)
    RepeatAbund_cent_max$Label <- rep("MaxRepeatAbund", nrow(RepeatAbund_cent_max))
  }
  if (length(RepeatAbund_cent_max)>0 && length(RepeatAbund_cent_min)>0){
    colnames(RepeatAbund_cent_max) <-  colnames(RepeatAbund_cent_min)
  }

  print(RepeatAbund_max)
  print(RepeatAbund_min)

  RepeatAbund_cent <- rbind(RepeatAbund_cent_min, RepeatAbund_cent_max)

  # output final files
  utils::write.table(x=RepeatAbund_min, file=paste0(outpath,"/", fname,"/Summary_output/histograms/", fname,  "_RepeatAbund_centromere_prediction_min.txt"), sep = "\t", dec = ".",row.names = FALSE, col.names = FALSE)
  utils::write.table(x=RepeatAbund_max, file=paste0(outpath,"/", fname,"/Summary_output/histograms/", fname,  "_RepeatAbund_centromere_prediction_max.txt"), sep = "\t", dec = ".",row.names = FALSE, col.names = FALSE)
  utils::write.table(x=RepeatAbund_length, file=paste0(outpath,"/", fname,"/Summary_output/histograms/", fname, "_RepeatAbund_centromere_prediction_length.txt"), sep = "\t", dec = ".",row.names = FALSE, col.names = FALSE)
  utils::write.table(x=RepeatAbund_cent_max, file=paste0(outpath,"/", fname,"/Summary_output/histograms/", fname,  "_RepeatAbund_centromere_maxrange.txt"), sep = "\t", dec = ".",row.names = FALSE, col.names = FALSE)
  utils::write.table(x=RepeatAbund_cent_min, file=paste0(outpath,"/", fname,"/Summary_output/histograms/", fname,  "_RepeatAbund_centromere_minrange.txt"), sep = "\t", dec = ".",row.names = FALSE, col.names = FALSE)
  utils::write.table(x=RepeatAbund_cent, file=paste0(outpath,"/", fname,"/Summary_output/histograms/", fname,  "_RepeatAbund_centromere_range.txt"), sep = "\t", dec = ".",row.names = FALSE, col.names = FALSE)

  #-----------------------------------------------
  # get ranges for Shannon diversity

  Shannon_list <- file_list[grep("Shannon", file_list)]

  # get Shannon div data
  lsd <- lapply(Shannon_list, read.table)
  sd_chr_list0 <- basename(Shannon_list)
  sd_chr_list1 <- stringr::str_split(sd_chr_list0, "_", simplify =TRUE)
  sd_chr_list2 <- sd_chr_list1[,3]
  names(lsd) <- sd_chr_list2
  Shannon_div_total <- dplyr::bind_rows(lsd, .id = 'chromosome')
  colnames(Shannon_div_total) <- c("Chromosome", "Genome_position", "Shannon_div")

  #--------------------------------------
  # join parts
  Shannon_div_total_parts <- Shannon_div_total[grep("-", Shannon_div_total$Chromosome),]

  if (nrow(Shannon_div_total_parts)>0){
    Shannon_div_total$Chrnum0 <- as.factor(stringr::str_split(Shannon_div_total$Chromosome, "r", simplify =TRUE)[,2])
    Shannon_div_total$Chrnum <- as.factor(stringr::str_split(Shannon_div_total$Chrnum0, "-", simplify =TRUE)[,1])
    Shannon_div_total$Chrpart <- as.numeric(stringr::str_split(Shannon_div_total$Chrnum0, "-", simplify =TRUE)[,2])
    Shannon_div_total$Chrpart[which(is.na(Shannon_div_total$Chrpart))] <-1
    Shannon_div_total$Genome_position <- Shannon_div_total$Genome_position + (Shannon_div_total$Chrpart -1)*4e8
  }else{
    Shannon_div_total$Chrnum <- as.factor(stringr::str_split(Shannon_div_total$Chromosome, "r", simplify =TRUE)[,2])
  }

  Shannon_div_total$Genome_position <- as.numeric(as.character(Shannon_div_total$Genome_position))
  Shannon_div_total$Shannon_div <- as.numeric(as.character(Shannon_div_total$Shannon_div))

  #---------------------------------------------

  # find start and end of highly repeating regions based 1SD from min
  Shannon_cent <- data.frame()
  Shannon_min <- data.frame()
  Shannon_length <- data.frame()
  Shannon <- data.frame()

  for (chromosome in unique(Shannon_div_total$Chrnum)){

    Shannon_div_chr0 <- Shannon_div_total[which(Shannon_div_total$Chrnum == chromosome),]

    Shannon_div_chr <- Shannon_div_chr0[which(!is.na(Shannon_div_chr0$Shannon_div)),]

    grDevices::png(filename = paste0(outpath,"/", fname, "/Summary_output/Shannon_div/", fname,"_Chr", chromosome, "_Shannon_div.png"), width = 1400, height = 700)
    plot(Shannon_div_chr$Genome_position, Shannon_div_chr$Shannon_div)
    grDevices::dev.off()

    Shannon_div_chr$roll_mean_Shannon <- zoo::rollmean(Shannon_div_chr$Shannon_div, k=Shannon_bin_size, fill = NA, align = "center")

    # find min position
    cent_min <- Shannon_div_chr$Genome_position[which(Shannon_div_chr$roll_mean_Shannon == min(Shannon_div_chr$roll_mean_Shannon, na.rm=TRUE))]
    print(min(Shannon_div_chr$roll_mean_Shannon, na.rm=TRUE))

    grDevices::png(filename = paste0(outpath,"/", fname, "/Summary_output/Shannon_div/", fname,"_Chr", chromosome, "_Shannon_div_rollmean_",Shannon_bin_size,".png"), width = 1400, height = 700)
    plot(Shannon_div_chr$Genome_position, Shannon_div_chr$roll_mean_Shannon)
    grDevices::dev.off()

    # find SD of data
    SD_repeatAbund <- sd(Shannon_div_chr$roll_mean_Shannon, na.rm=TRUE)

    thres = mean(Shannon_div_chr$roll_mean_Shannon, na.rm=TRUE) - (Shannon_SD * SD_repeatAbund)

    # find positions of + 1 SD from mean
    cent_range_wind <- Shannon_div_chr$Genome_position[which(Shannon_div_chr$roll_mean_Shannon <= thres)]/5000

    # find range of these values
    # ChemoSpecUtils
    # https://rdrr.io/cran/ChemoSpecUtils/man/check4Gaps.html

    #library(ChemoSpecUtils)

    #cent_range <- ChemoSpecUtils::check4Gaps(cent_range_wind)
    #cent_range[nrow(cent_range)+1,] <- c(0,0,0,0,0)

    if (length(cent_range_wind)>0) {
        cent_range <- ChemoSpecUtils::check4Gaps(cent_range_wind)
    } else {cent_range <- data.frame(0,0,0,0)}
    cent_range[nrow(cent_range)+1,] <- c(0,0,0,0,0)

    cent_range_pos_start <- cent_range[,1]*5000
    cent_range_pos_end <- cent_range[,2]*5000

    SPP_l <- rep(fname, length(cent_range_pos_start))
    Chr_l <- rep(chromosome, length(cent_range_pos_start))
    typel <- rep("Shannon", length(cent_range_pos_start))

    Shannon_cent_chr <- cbind(typel, SPP_l, Chr_l, cent_range_pos_start, cent_range_pos_end)
    Shannon_cent <- rbind(Shannon_cent, Shannon_cent_chr)

    Shannon_min_chr <- c("Shannon", fname, chromosome, cent_min)
    Shannon_min <- rbind(Shannon_min, Shannon_min_chr)

    chr_length <- max(Shannon_div_chr$Genome_position)
    Shannon_length_chr <- c("Length", fname, chromosome, chr_length)
    Shannon_length <- rbind(Shannon_length, Shannon_length_chr)
    Shannon <- rbind(Shannon, Shannon_div_chr)

  }

  # remove zeros
  Shannon_cent <- Shannon_cent[-which(Shannon_cent[,4] == 0),]
  print(Shannon_min)

  # output final file
  utils::write.table(x=Shannon_cent, file=paste0(outpath,"/", fname,"/Summary_output/Shannon_div/", fname, "_Shannon_centromere_range.txt"), sep = "\t", dec = ".",row.names = FALSE, col.names = FALSE)
  utils::write.table(x=Shannon_min, file=paste0(outpath,"/", fname,"/Summary_output/Shannon_div/", fname, "_Shannon_centromere_prediction_min.txt"), sep = "\t", dec = ".",row.names = FALSE, col.names = FALSE)
  utils::write.table(x=Shannon_length, file=paste0(outpath,"/", fname,"/Summary_output/Shannon_div/", fname, "_Shannon_centromere_prediction_length.txt"), sep = "\t", dec = ".",row.names = FALSE, col.names = FALSE)

}

cent_finalize <- function(fname=fname, outpath=outpath){

  # read in RepAbund data
  summary_path <- paste0(outpath,"/", fname, "/", "Summary_output/output_data")
  file_list <- list.files(summary_path, full.names=TRUE)

  rollsumhist_list <- file_list[grep("Repeat_abundance_sum", file_list)]
  lsd <- lapply(rollsumhist_list, read.table)
  sd_chr_list0 <- basename(rollsumhist_list)
  sd_chr_list1 <- stringr::str_split(sd_chr_list0, "_", simplify =TRUE)
  sd_chr_list2 <- sd_chr_list1[,3]
  names(lsd) <- sd_chr_list2
  RepeatAbundance_total <- dplyr::bind_rows(lsd, .id = 'chromosome')
  colnames(RepeatAbundance_total) <- c("Chromosome", "Genome_position", "RepeatAbundance")

  #--------------------------------------
  # Read in centromere predictions
  RepAbund_maxcent <- base::as.data.frame(base::as.matrix(utils::read.table(paste0(outpath,"/", fname, "/", "Summary_output/histograms/", fname,"_RepeatAbund_centromere_prediction_max.txt"), header = FALSE, check.names = FALSE)))
  RepAbund_mincent <- base::as.data.frame(base::as.matrix(utils::read.table(paste0(outpath,"/", fname, "/", "Summary_output/histograms/", fname, "_RepeatAbund_centromere_prediction_min.txt"), header = FALSE, check.names = FALSE)))
  Shannon_cent <- base::as.data.frame(base::as.matrix(utils::read.table(paste0(outpath,"/", fname, "/", "Summary_output/Shannon_div/", fname, "_Shannon_centromere_prediction_min.txt"), header = FALSE, check.names = FALSE)))
  Hist_cent <- base::as.data.frame(base::as.matrix(utils::read.table(paste0(outpath,"/", fname, "/", "Summary_output/histograms/", fname, "_Centromere_histograms_summary.txt"), sep=" ", fill=T, header = FALSE, check.names = FALSE)))

  # add column names
  colnames(RepAbund_maxcent) <- c("Label", "Spp", "Chr", "Centromere")
  colnames(RepAbund_mincent) <- c("Label", "Spp", "Chr", "Centromere")
  colnames(Shannon_cent) <- c("Label", "Spp", "Chr", "Centromere")
  colnames(Hist_cent) <- c("Spp", "Chr", "Centromere", "Length")

  # add/change label for some

  RepAbund_maxcent$Chr <- as.numeric(RepAbund_maxcent$Chr)
  RepAbund_mincent$Chr <- as.numeric(RepAbund_mincent$Chr)
  Shannon_cent$Chr <- as.numeric(Shannon_cent$Chr)

  RepAbund_maxcent$Chr <- paste0("Chr", RepAbund_maxcent$Chr)
  RepAbund_mincent$Chr <- paste0("Chr", RepAbund_mincent$Chr)
  Shannon_cent$Chr <- paste0("Chr", Shannon_cent$Chr)

  Hist_cent$Label <- rep("Histogram", nrow(Hist_cent))

  # subset columns
  Hist_cent <- cbind(Hist_cent$Label, Hist_cent$Spp, Hist_cent$Chr, Hist_cent$Centromere)
  colnames(Hist_cent) <- c("Label", "Spp", "Chr", "Centromere")

  # make a total possible cent positions - with all methods
  Total_cent <- rbind(RepAbund_maxcent, RepAbund_mincent, Shannon_cent, Hist_cent)

  utils::write.table(x=Total_cent, file=paste0(outpath,"/", fname, "/Summary_output/", fname,"_total_possible.txt"), sep = "\t", dec = ".",row.names = FALSE, col.names = TRUE)

  #-------------------------------------------
  # read in ranges
  RepAbund_cent_range <- base::as.data.frame(base::as.matrix(utils::read.table(paste0(outpath,"/", fname, "/", "Summary_output/histograms/", fname, "_RepeatAbund_centromere_range.txt"), header = FALSE, check.names = FALSE)))
  Shannon_cent_range <- base::as.data.frame(base::as.matrix(utils::read.table(paste0(outpath,"/", fname, "/", "Summary_output/Shannon_div/", fname, "_Shannon_centromere_range.txt"), header = FALSE, check.names = FALSE)))

  # rearrange columns for joining
  RepAbund_cent_range_ed <- cbind(RepAbund_cent_range$V5, RepAbund_cent_range$V1,RepAbund_cent_range$V2,RepAbund_cent_range$V3,RepAbund_cent_range$V4)

  colnames(RepAbund_cent_range_ed) <-  c("Label", "Spp", "Chr", "Start", "End")
  colnames(Shannon_cent_range) <-  c("Label", "Spp", "Chr", "Start", "End")

  # make a total possible cent ranges - with all methods
  Total_cent_range <- rbind(RepAbund_cent_range_ed, Shannon_cent_range)

  Total_cent_range$Chr <- as.numeric(Total_cent_range$Chr)
  Total_cent_range$Chr <- paste0("Chr", Total_cent_range$Chr)

  utils::write.table(x=Total_cent_range, file=paste0(outpath,"/", fname, "/Summary_output/", fname,"_total_possible_range.txt"), sep = "\t", dec = ".",row.names = FALSE, col.names = TRUE)

}

