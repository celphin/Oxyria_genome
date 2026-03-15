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

library(RepeatOBserverV1) 
library(zoo)
library(dplyr)

# Function

calculate_centromeres <- function(fname,
                                  outpath,
                                  Shannon_bin_size = 800,
                                  Shannon_SD = 1.5,
                                  Rep_bin_size = 200,
                                  Rep_SD = 1.5,
                                  min_region_length = 1e6) {


  message("Processing: ", fname)

  summary_path <- file.path(outpath,"Summary_output/output_data")
  file_list <- list.files(summary_path, full.names = TRUE)

  # ---------------------------
  # READ REPEAT DATA
  # ---------------------------

  repeat_files <- file_list[grepl("Repeat_abundance_sum", file_list)]

  repeat_list <- lapply(repeat_files, read.table)

  names(repeat_list) <- gsub(".*_(Chr[0-9]+).*", "\\1",
                             basename(repeat_files))

  Repeat_total <- bind_rows(repeat_list, .id = "Chr")

  colnames(Repeat_total) <- c("Chromosome",
                              "Genome_position",
                              "RepeatAbundance")

  Repeat_total$Genome_position <- as.numeric(Repeat_total$Genome_position)
  Repeat_total$RepeatAbundance <- as.numeric(Repeat_total$RepeatAbundance)

  results <- list()

  # ===========================
  # PER CHROMOSOME ANALYSIS
  # ===========================

  for (chr in unique(Repeat_total$Chromosome)) {

    df <- Repeat_total %>%
      filter(Chromosome == chr) %>%
      arrange(Genome_position)

    # Rolling mean
    df$roll <- rollmean(df$RepeatAbundance,
                        k = Rep_bin_size,
                        fill = NA,
                        align = "center")

    # Z-score (more stable than mean±SD alone)
    df$z <- (df$roll - mean(df$roll, na.rm=TRUE)) /
            sd(df$roll, na.rm=TRUE)

    # Detect dips
    dip <- df$z <= -Rep_SD

    # Run length encoding (fast edge detection)
    r <- rle(dip)
    ends <- cumsum(r$lengths)
    starts <- ends - r$lengths + 1

    for (i in which(r$values)) {

      start_pos <- df$Genome_position[starts[i]]
      end_pos   <- df$Genome_position[ends[i]]

      region_length <- end_pos - start_pos

      if (!is.na(region_length) &&
          region_length >= min_region_length) {

        results[[length(results)+1]] <-
          data.frame(Spp = fname,
                     Chr = chr,
                     Method = "Repeat",
                     Start = start_pos,
                     End = end_pos,
                     Length = region_length)
      }
    }
  }

  # ---------------------------
  # SHANNON ANALYSIS
  # ---------------------------

  shannon_files <- file_list[grepl("Shannon", file_list)]

  if (length(shannon_files) > 0) {

    shannon_list <- lapply(shannon_files, read.table)
    names(shannon_list) <- gsub(".*_(Chr[0-9]+).*",
                                "\\1",
                                basename(shannon_files))

    Shannon_total <- bind_rows(shannon_list, .id="Chr")
    colnames(Shannon_total) <- c("Chromosome",
                                "Genome_position",
                                "Shannon")

    Shannon_total$Genome_position <- as.numeric(Shannon_total$Genome_position)
    Shannon_total$Shannon <- as.numeric(Shannon_total$Shannon)

    for (chr in unique(Shannon_total$Chromosome)) {

      df <- Shannon_total %>%
        filter(Chromosome == chr) %>%
        arrange(Genome_position)

      df$roll <- rollmean(df$Shannon,
                          k = Shannon_bin_size,
                          fill=NA,
                          align="center")

      df$z <- (df$roll - mean(df$roll, na.rm=TRUE)) /
              sd(df$roll, na.rm=TRUE)

      valley <- df$z <= -Shannon_SD

      r <- rle(valley)
      ends <- cumsum(r$lengths)
      starts <- ends - r$lengths + 1

      for (i in which(r$values)) {

        start_pos <- df$Genome_position[starts[i]]
        end_pos   <- df$Genome_position[ends[i]]

        region_length <- end_pos - start_pos

        if (!is.na(region_length) &&
            region_length >= min_region_length) {

          results[[length(results)+1]] <-
            data.frame(Spp = fname,
                       Chr = chr,
                       Method = "Shannon",
                       Start = start_pos,
                       End = end_pos,
                       Length = region_length)
        }
      }
    }
  }

  # ---------------------------
  # COMBINE + PICK BEST PER CHR
  # ---------------------------

  if (length(results) == 0) {
    message("No regions detected.")
    return(NULL)
  }

  final_df <- bind_rows(results)

  # Pick strongest region per chromosome
  final_best <- final_df %>%
    group_by(Chr) %>%
    slice_max(order_by = Length, n = 1) %>%
    ungroup()

  # ---------------------------
  # WRITE OUTPUT
  # ---------------------------

  out_file <- file.path(paste0(outpath, fname,
                               "_centromere_best.txt"))

  write.table(final_best,
              file = out_file,
              sep = "\t",
              row.names = FALSE,
              quote = FALSE)

  return(final_best)
}


# -----------------------------
# Run for all genomes
outpath <- "/home/celphin/scratch/Arctic_centromere_lengths/raw_data/"

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
  calculate_centromeres(x, outpath))

results_all <- do.call(rbind, results_all)


# Combine all species
final_table <- bind_rows(all_results)

# Save master file
write.table(final_table,
            file=paste0(outpath,
                        "ALL_species_centromere_ranges.txt"),
            sep="\t",
            row.names=FALSE,
            quote=FALSE)


##########################
# Look at data with plots for all spp

##########################
# copy *_total_possible_range.txt to one folder

# Narval
# cd /home/celphin/scratch/Arctic_centromere_lengths/
# cat *_total_possible_range.txt > Arctic_genomes_total_possible_range.txt

# wc -l Arctic_genomes_total_possible_range.txt
# 2624 Arctic_genomes_total_possible_range.txt

# Copy to local machine
# ~/Github/Oxyria_genome/3_Genome_Comparisons/6_RepeatOBserver

#############################
# read in cent_pred

fname="Arctic_genomes"
outpath="~/Github/Oxyria_genome/3_Genome_Comparisons/6_RepeatOBserver"

# Wrap in function
#automated_centromere_detection <- function(fname=fname,  outpath=outpath) {

# read in the data
df0 <- utils::read.table(paste0(outpath,"/",fname,"_total_possible_range.txt"), sep = "\t", header=TRUE, check.names = FALSE)

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

