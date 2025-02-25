############################
# Check Genespace data for matrix of space between gene positions and other genes
# Oct 2024
############################

# Combined BLASTS between genomes
# For each species
# Find genes that blast to all other genomes
# Append go term 
# For each pair of orthogroups for each species
# Are they on the same chromosome?
# If yes, how far apart are the starts?
# If no, make diff 200Mbp or large value
# Also how far from centromere - last column
# Sum matrices between species of a type
# What gene pairs have low value for all spp 
# Sort by gene ontology term and ortholog group

#---------------------------
tmux new-session -s synteny
tmux attach-session -t synteny

cd /home/celphin/scratch/Oxyria/synteny_quantity

# find OG in all
grep "Oxy" /home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/results/Orthogroups.tsv \
| grep "DoctH0" | grep "\-lg" | grep "AALP" | grep Polavi | grep FT | grep FEH | grep Rno | grep Rta >> Total_Orthogroups_in_all.tsv

wc -l Total_Orthogroups_in_all.tsv
# 6909 

# copy over BLAST combined data
cp /home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/results/combBed.txt ./Total_combBed.txt

wc -l Total_combBed.txt
# 651560 Total_combBed.txt

#------------------
# Join Interpro files but keep file name
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


############################################
# Write code to run in 500x500 parts

# each 500 x500 part takes about 2-3 hours
# Oxyria - Time difference of 2.810299 hours
# Dryas - Time difference of 1.896123 hours

mkdir /home/celphin/scratch/Oxyria/synteny_quantity
cd /home/celphin/scratch/Oxyria/synteny_quantity
# mkdir tmp

tmux new-session -s synteny
tmux attach-session -t synteny

salloc -c40 --time 7:00:00 --mem 191000m --account def-rieseber

module load  StdEnv/2020 r/4.2.2
export R_LIBS_USER=/home/celphin/R/x86_64-pc-linux-gnu-library/4.1.0/

R

path="/home/celphin/scratch/Oxyria/synteny_quantity"
og_file="Total_Orthogroups_in_all.tsv"
comb_file="Total_combBed.txt"
x_cpu=40
spp_hap="Oxyria_digyna_H1"

orthogroup_distances(spp_hap="Oxyria_digyna_H1", path=path, og_file=og_file, comb_file=comb_file, x_cpu=x_cpu)
orthogroup_distances(spp_hap="Dryas_octopetala", path=path, og_file=og_file, comb_file=comb_file, x_cpu=x_cpu)
orthogroup_distances(spp_hap="Rheum_nobile_H0", path=path, og_file=og_file, comb_file=comb_file, x_cpu=x_cpu)
orthogroup_distances(spp_hap="Arabis_alpina", path=path, og_file=og_file, comb_file=comb_file, x_cpu=x_cpu)
orthogroup_distances(spp_hap="Polygunum_aviculare_H0", path=path, og_file=og_file, comb_file=comb_file, x_cpu=x_cpu)
orthogroup_distances(spp_hap="Draba_nivalis", path=path, og_file=og_file, comb_file=comb_file, x_cpu=x_cpu)


##############################
# Functions

orthogroup_distances <- function(spp_hap, path, og_file, comb_file, x_cpu) {
  
  # Load files
  orthogroups_in_all <- read.table(file.path(path, og_file), sep = "\t", header = FALSE, check.names = FALSE)
  
  combined_genes <- read.table(file.path(path, comb_file), sep = "\t", header = TRUE, check.names = FALSE)
  
  colnames(combined_genes) <- c("chr", "start", "end", "gene", "ofID", "pepLen", "ord", "genome", 
                                 "arrayID", "isArrayRep", "orthogroup", "globHOG", "noAnchor", "og")
  
  colnames(orthogroups_in_all)[1] <- "orthogroup"
  
  # Join files
  combined_genes_orthogroups <- dplyr::left_join(orthogroups_in_all, combined_genes, by = "orthogroup")
  
  # Filter species of interest
  orthogroups_in_spp <- combined_genes_orthogroups[combined_genes_orthogroups$genome == spp_hap, ]
  orthogroup_list_total <- unique(orthogroups_in_spp$orthogroup)
  
  n_groups <- length(orthogroup_list_total)
  og_groups_list <- seq(0, ((n_groups %/% 500)-1))
  X_seq_along <- outer(og_groups_list, og_groups_list, function(i, j) paste0(i, "_", j))
  
  cl <- parallel::makeCluster(x_cpu)
  on.exit(parallel::stopCluster(cl))  # Ensure cluster stops after execution
  
  results <- parallel::parSapply(cl, seq_along(X_seq_along), 
                                  orthogroup_loop, X_seq_along, spp_hap = spp_hap, path = path, 
                                  orthogroups_in_spp = orthogroups_in_spp,
                                  orthogroup_list_total = orthogroup_list_total)
  
# Initialize lists to hold matrices with dimensions
pos_diff_matrices <- lapply(1:(n_groups %/% 500), function(x) matrix(NA, nrow = 500, ncol = 0))
ord_diff_matrices <- lapply(1:(n_groups %/% 500), function(x) matrix(NA, nrow = 500, ncol = 0))
pos_diff_max_matrices <- lapply(1:(n_groups %/% 500), function(x) matrix(NA, nrow = 500, ncol = 0))
gene_names1_matrices <- lapply(1:(n_groups %/% 500), function(x) matrix(NA, nrow = 500, ncol = 0))
gene_names2_matrices <- lapply(1:(n_groups %/% 500), function(x) matrix(NA, nrow = 500, ncol = 0))

# Fill matrices with results
for (index in seq_along(X_seq_along)) {
  # Extract the relevant orthogroup indices from X_seq_along
  split_values <- as.numeric(unlist(strsplit(X_seq_along[index], "_")))
  i <- split_values[1] +1
  j <- split_values[2] +1
  
  # Define the row and column indices for the matrix
  row_indices <- ((i - 1) * 500 + 1):((i - 1) * 500 + 9)
  col_indices <- ((j - 1) * 500 + 1):((j - 1) * 500 + 9)
  
  # Fill the matrices for the corresponding index
  pos_diff_matrices[[i]] <- cbind(pos_diff_matrices[[i]], as.data.frame(results[(index - 1) * 5 + 1, drop = FALSE]))
  ord_diff_matrices[[i]] <- cbind(ord_diff_matrices[[i]], as.data.frame(results[(index - 1) * 5 + 2, drop = FALSE]))
  pos_diff_max_matrices[[i]] <- cbind(pos_diff_max_matrices[[i]], as.data.frame(results[(index - 1) * 5 + 3, drop = FALSE]))
  gene_names1_matrices[[i]] <- cbind(gene_names1_matrices[[i]], as.data.frame(results[(index - 1) * 5 + 4, drop = FALSE]))
  gene_names2_matrices[[i]] <- cbind(gene_names2_matrices[[i]], as.data.frame(results[(index - 1) * 5 + 5, drop = FALSE]))
}

# Combine matrices from the lists into final matrices
pos_diff_matrix <- do.call(rbind, pos_diff_matrices)
ord_diff_matrix <- do.call(rbind, ord_diff_matrices)
pos_diff_max_matrix <- do.call(rbind, pos_diff_max_matrices)
gene_names1_matrix <- do.call(rbind, gene_names1_matrices)
gene_names2_matrix <- do.call(rbind, gene_names2_matrices)

  # Set row and column names
  # rownames(pos_diff_matrix) <- orthogroup_list_total
  # colnames(pos_diff_matrix) <- orthogroup_list_total
  # rownames(ord_diff_matrix) <- orthogroup_list_total
  # colnames(ord_diff_matrix) <- orthogroup_list_total
  # rownames(pos_diff_max_matrix) <- orthogroup_list_total
  # colnames(pos_diff_max_matrix) <- orthogroup_list_total
  # rownames(gene_names1_matrix) <- orthogroup_list_total
  # colnames(gene_names1_matrix) <- orthogroup_list_total
  # rownames(gene_names2_matrix) <- orthogroup_list_total
  # colnames(gene_names2_matrix) <- orthogroup_list_total

  # Write output matrices to files
  utils::write.table(pos_diff_matrix, file = file.path(path, "tmp", sprintf("%s_synteny_matrix.txt", spp_hap)), sep = "\t", row.names = TRUE, col.names = TRUE)
  utils::write.table(ord_diff_matrix, file = file.path(path, "tmp", sprintf("%s_synteny_matrix_ord.txt", spp_hap)), sep = "\t", row.names = TRUE, col.names = TRUE)
  utils::write.table(pos_diff_max_matrix, file = file.path(path, "tmp", sprintf("%s_synteny_matrix_max.txt", spp_hap)), sep = "\t", row.names = TRUE, col.names = TRUE)
  utils::write.table(gene_names1_matrix, file = file.path(path, "tmp", sprintf("%s_synteny_matrix_genes1.txt", spp_hap)), sep = "\t", row.names = TRUE, col.names = TRUE)
  utils::write.table(gene_names2_matrix, file = file.path(path, "tmp", sprintf("%s_synteny_matrix_genes2.txt", spp_hap)), sep = "\t", row.names = TRUE, col.names = TRUE)

}


orthogroup_loop <- function(X, X_seq_along, spp_hap, path, orthogroups_in_spp, orthogroup_list_total) {

  XX <- X_seq_along[X]
  split_values <- strsplit(XX, "_")[[1]]
  i <- as.numeric(split_values[1])
  j <- as.numeric(split_values[2])
  
  start_orthog1 = (i * 500) + 1
  numb_orthogroups1 = start_orthog1 + 499
  
  start_orthog2 = (j * 500) + 1
  numb_orthogroups2 = start_orthog2 + 499
  
  # Ensure indices do not exceed the length of the list
  start_orthog1 <- min(start_orthog1, (length(orthogroup_list_total)-9))
  numb_orthogroups1 <- min(numb_orthogroups1, length(orthogroup_list_total))
  start_orthog2 <- min(start_orthog2, (length(orthogroup_list_total)-9))
  numb_orthogroups2 <- min(numb_orthogroups2, length(orthogroup_list_total))
  
    multiResultClass <- function(pos_diff=NULL, ord_diff=NULL, pos_diff_max=NULL, gene_names1=NULL, gene_names2=NULL)
  {
    me <- list(
      pos_diff = pos_diff,
      ord_diff = ord_diff,
      pos_diff_max = pos_diff_max,
      gene_names1 = gene_names1,
      gene_names2 = gene_names2
    )

    ## Set the name for the class
    class(me) <- append(class(me),"multiResultClass")
    return(me)
  }
  
  # Get this set of orthogroups
  orthogroup_list1 <- orthogroup_list_total[start_orthog1:numb_orthogroups1]
  orthogroup_list2 <- orthogroup_list_total[start_orthog2:numb_orthogroups2]
  
  if (file.exists(paste0(path,"/tmp/", spp_hap, "_", start_orthog1, "_", start_orthog2,"_synteny_matrix_ord.txt"))) {
    print("File exists already")
  }else{
  # Initialize matrices
  n1 <- length(orthogroup_list1)
  n2 <- length(orthogroup_list2)
  pos_diff <- matrix(NA, nrow = n1, ncol = n2)
  ord_diff <- matrix(NA, nrow = n1, ncol = n2)
  pos_diff_max <- matrix(NA, nrow = n1, ncol = n2)
  gene_names1 <- matrix(NA, nrow = n1, ncol = n2)
  gene_names2 <- matrix(NA, nrow = n1, ncol = n2)

  # Use outer for pairwise operations
  for (og1 in 1:n1) {
    for (og2 in 1:n2) {

      gene_list1 <- orthogroups_in_spp$gene[which(orthogroups_in_spp$orthogroup==orthogroup_list1[og1])]
      gene_list2 <- orthogroups_in_spp$gene[which(orthogroups_in_spp$orthogroup==orthogroup_list2[og2])]
      
      chr_list1 <- orthogroups_in_spp$chr[which(orthogroups_in_spp$orthogroup==orthogroup_list1[og1])]
      chr_list2 <- orthogroups_in_spp$chr[which(orthogroups_in_spp$orthogroup==orthogroup_list2[og2])]
      
      gene_diff <- outer(gene_list1, gene_list2, Vectorize(function(g1, g2) {
        chr1 <- chr_list1[which(gene_list1 == g1)]
        chr2 <- chr_list2[which(gene_list2 == g2)]
        if (chr1 == chr2) {
          start1 <- orthogroups_in_spp$start[orthogroups_in_spp$gene == g1]
          start2 <- orthogroups_in_spp$start[orthogroups_in_spp$gene == g2]
          abs(start1 - start2)
        } else {
          400e6  # Large distance if chromosomes are different
        }
      }))
      
      min_diff_idx <- which(gene_diff == min(as.numeric(unlist(gene_diff)), na.rm = TRUE), arr.ind = TRUE)
      max_diff_idx <- which(gene_diff == max(gene_diff, na.rm = TRUE), arr.ind = TRUE)
	  min_diff_idx <- rbind(min_diff_idx, newrow = c(0, 0))
	  max_diff_idx <- rbind(max_diff_idx, newrow = c(0, 0))

      pos_diff[og1, og2] <- gene_diff[min_diff_idx[1,1], min_diff_idx[1,2]]
      ord_diff[og1, og2] <- abs(orthogroups_in_spp$ord[which(orthogroups_in_spp$gene == gene_list1[min_diff_idx[1,1]])] - 
                                orthogroups_in_spp$ord[which(orthogroups_in_spp$gene == gene_list2[min_diff_idx[1,2]])])
      pos_diff_max[og1, og2] <- gene_diff[max_diff_idx[1,1], max_diff_idx[1,2]]
      
      gene_names1[og1, og2] <- gene_list1[min_diff_idx[1,1]]
      gene_names2[og1, og2] <- gene_list2[min_diff_idx[1,2]]
    }
  }
  # Set row and column names
  rownames(pos_diff) <- orthogroup_list1
  colnames(pos_diff) <- orthogroup_list2
  rownames(ord_diff) <- orthogroup_list1
  colnames(ord_diff) <- orthogroup_list2
  rownames(pos_diff_max) <- orthogroup_list1
  colnames(pos_diff_max) <- orthogroup_list2
  rownames(gene_names1) <- orthogroup_list1
  colnames(gene_names1) <- orthogroup_list2
  rownames(gene_names2) <- orthogroup_list1
  colnames(gene_names2) <- orthogroup_list2
  
  # Return result
  result <- multiResultClass()
  result$pos_diff <- pos_diff
  result$ord_diff <- ord_diff
  result$pos_diff_max <- pos_diff_max
  result$gene_names1 <- gene_names1
  result$gene_names2 <- gene_names2
  
  utils::write.table(x=pos_diff, file=paste0(path,"/tmp/", spp_hap, "_", start_orthog1, "_", start_orthog2,"_synteny_matrix.txt"), append = FALSE, sep = "\t", dec = ".", row.names = TRUE, col.names = TRUE)
  utils::write.table(x=ord_diff, file=paste0(path,"/tmp/", spp_hap, "_", start_orthog1, "_", start_orthog2,"_synteny_matrix_ord.txt"), append = FALSE, sep = "\t", dec = ".", row.names = TRUE, col.names = TRUE)
  utils::write.table(x=pos_diff_max, file=paste0(path,"/tmp/", spp_hap, "_", start_orthog1, "_", start_orthog2,"_synteny_matrix_max.txt"), append = FALSE, sep = "\t", dec = ".", row.names = TRUE, col.names = TRUE)
  utils::write.table(x=gene_names1, file=paste0(path,"/tmp/", spp_hap, "_", start_orthog1, "_", start_orthog2,"_synteny_matrix_genes1.txt"), append = FALSE, sep = "\t", dec = ".", row.names = TRUE, col.names = TRUE)
  utils::write.table(x=gene_names2, file=paste0(path,"/tmp/", spp_hap, "_", start_orthog1, "_", start_orthog2,"_synteny_matrix_genes2.txt"), append = FALSE, sep = "\t", dec = ".", row.names = TRUE, col.names = TRUE)
  return(result)
}
}

#-----------------------
# once all have run move important files up
cp Rheum_nobile_H0_synteny_matrix_ord.txt ..
cp Oxyria_digyna_H1_synteny_matrix_ord.txt ..
cp Dryas_octopetala_synteny_matrix_ord.txt ..
cp Arabis_alpina_synteny_matrix_ord.txt ..
cp Polygunum_aviculare_H0_synteny_matrix_ord.txt ..
cp Draba_nivalis_synteny_matrix_ord.txt ..

cp Rheum_nobile_H0_synteny_matrix_genes1.txt ..
cp Oxyria_digyna_H1_synteny_matrix_genes1.txt ..
cp Dryas_octopetala_synteny_matrix_genes1.txt ..
cp Arabis_alpina_synteny_matrix_genes1.txt ..
cp Polygunum_aviculare_H0_synteny_matrix_genes1.txt ..
cp Draba_nivalis_synteny_matrix_genes1.txt ..

cp Rheum_nobile_H0_synteny_matrix_genes2.txt ..
cp Oxyria_digyna_H1_synteny_matrix_genes2.txt ..
cp Dryas_octopetala_synteny_matrix_genes2.txt ..
cp Arabis_alpina_synteny_matrix_genes2.txt ..
cp Polygunum_aviculare_H0_synteny_matrix_genes2.txt ..
cp Draba_nivalis_synteny_matrix_genes2.txt ..


#################################
# compare matrices for different species

tmux new-session -s synteny1
tmux attach-session -t synteny1

cd /home/celphin/scratch/Oxyria/synteny_quantity

salloc -c1 --time 2:00:00 --mem 191000m --account def-rieseber

module load  StdEnv/2020 r/4.2.2
export R_LIBS_USER=/home/celphin/R/x86_64-pc-linux-gnu-library/4.1.0/

R

path="/home/celphin/scratch/Oxyria/synteny_quantity"
numb_orthogroups=500
comb_file="Total_combBed.txt"

spp_hap1="Oxyria_digyna_H1"
spp_hap2="Dryas_octopetala"
spp_hap3="Rheum_nobile_H0"
spp_hap4="Arabis_alpina"
spp_hap5="Polygunum_aviculare_H0"
spp_hap6="Draba_nivalis"

#--------------------------
# get gene positions
combined_genes <- base::as.data.frame(utils::read.table(paste0(path,"/", comb_file), sep="\t", header = TRUE, check.names = FALSE))
colnames(combined_genes) <- c( "chr", "start", "end", "gene", "ofID", "pepLen", "ord", "genome", 
"arrayID","isArrayRep", "orthogroup", "globHOG", "noAnchor", "og")
gene_pos <- cbind(combined_genes$gene, combined_genes$chr, combined_genes$start, combined_genes$ord, combined_genes$orthogroup, combined_genes$genome)

#load files
gene1_names_spp1 <- base::as.matrix(utils::read.table(paste0(path,"/", spp_hap1, "_synteny_matrix_genes1.txt"), row.names = 1, header = TRUE, check.names = FALSE))
gene1_names_spp2 <- base::as.matrix(utils::read.table(paste0(path,"/", spp_hap2, "_synteny_matrix_genes1.txt"), row.names = 1, header = TRUE, check.names = FALSE))
gene1_names_spp3 <- base::as.matrix(utils::read.table(paste0(path,"/", spp_hap3, "_synteny_matrix_genes1.txt"), row.names = 1, header = TRUE, check.names = FALSE))
gene1_names_spp4 <- base::as.matrix(utils::read.table(paste0(path,"/", spp_hap4, "_synteny_matrix_genes1.txt"), row.names = 1, header = TRUE, check.names = FALSE))
gene1_names_spp5 <- base::as.matrix(utils::read.table(paste0(path,"/", spp_hap5, "_synteny_matrix_genes1.txt"), row.names = 1, header = TRUE, check.names = FALSE))
gene1_names_spp6 <- base::as.matrix(utils::read.table(paste0(path,"/", spp_hap6, "_synteny_matrix_genes1.txt"), row.names = 1, header = TRUE, check.names = FALSE))

gene2_names_spp1 <- base::as.matrix(utils::read.table(paste0(path,"/", spp_hap1, "_synteny_matrix_genes2.txt"), row.names = 1, header = TRUE, check.names = FALSE))
gene2_names_spp2 <- base::as.matrix(utils::read.table(paste0(path,"/", spp_hap2, "_synteny_matrix_genes2.txt"), row.names = 1, header = TRUE, check.names = FALSE))
gene2_names_spp3 <- base::as.matrix(utils::read.table(paste0(path,"/", spp_hap3, "_synteny_matrix_genes2.txt"), row.names = 1, header = TRUE, check.names = FALSE))
gene2_names_spp4 <- base::as.matrix(utils::read.table(paste0(path,"/", spp_hap4, "_synteny_matrix_genes2.txt"), row.names = 1, header = TRUE, check.names = FALSE))
gene2_names_spp5 <- base::as.matrix(utils::read.table(paste0(path,"/", spp_hap5, "_synteny_matrix_genes2.txt"), row.names = 1, header = TRUE, check.names = FALSE))
gene2_names_spp6 <- base::as.matrix(utils::read.table(paste0(path,"/", spp_hap6, "_synteny_matrix_genes2.txt"), row.names = 1, header = TRUE, check.names = FALSE))

ord_diff1 <- base::as.matrix(utils::read.table(paste0(path,"/", spp_hap1, "_synteny_matrix_ord.txt"), row.names = 1, header = TRUE, check.names = FALSE))
ord_diff2 <- base::as.matrix(utils::read.table(paste0(path,"/", spp_hap2, "_synteny_matrix_ord.txt"), row.names = 1, header = TRUE, check.names = FALSE))
ord_diff3 <- base::as.matrix(utils::read.table(paste0(path,"/", spp_hap3, "_synteny_matrix_ord.txt"), row.names = 1, header = TRUE, check.names = FALSE))
ord_diff4 <- base::as.matrix(utils::read.table(paste0(path,"/", spp_hap4, "_synteny_matrix_ord.txt"), row.names = 1, header = TRUE, check.names = FALSE))
ord_diff5 <- base::as.matrix(utils::read.table(paste0(path,"/", spp_hap5, "_synteny_matrix_ord.txt"), row.names = 1, header = TRUE, check.names = FALSE))
ord_diff6 <- base::as.matrix(utils::read.table(paste0(path,"/", spp_hap6, "_synteny_matrix_ord.txt"), row.names = 1, header = TRUE, check.names = FALSE))

#---------------------------
# find min diff and min sum
synteny_sum = ord_diff2 + ord_diff3 + ord_diff4 + ord_diff5 + ord_diff6 + ord_diff1 

#---------------------------
# determine which pos_diff are less than X Mbp - maybe 5Mbp
ord_binary1_small <- ord_diff1
ord_binary2_small <- ord_diff2
ord_binary3_small <- ord_diff3
ord_binary4_small <- ord_diff4
ord_binary5_small <- ord_diff5
ord_binary6_small <- ord_diff6

ord_binary1_small[which(ord_diff1==0, arr.ind=TRUE)] <- NA
ord_binary2_small[which(ord_diff2==0, arr.ind=TRUE)] <- NA
ord_binary3_small[which(ord_diff3==0, arr.ind=TRUE)] <- NA
ord_binary4_small[which(ord_diff4==0, arr.ind=TRUE)] <- NA
ord_binary5_small[which(ord_diff5==0, arr.ind=TRUE)] <- NA
ord_binary6_small[which(ord_diff6==0, arr.ind=TRUE)] <- NA

gene_count=40

ord_binary1_small[which(ord_diff1 < gene_count, arr.ind=TRUE)] <- 1
ord_binary1_small[which(ord_diff1 > gene_count, arr.ind=TRUE)] <- 0

ord_binary2_small[which(ord_diff2 < gene_count, arr.ind=TRUE)] <- 1
ord_binary2_small[which(ord_diff2 > gene_count, arr.ind=TRUE)] <- 0

ord_binary3_small[which(ord_diff3 < gene_count, arr.ind=TRUE)] <- 1
ord_binary3_small[which(ord_diff3 > gene_count, arr.ind=TRUE)] <- 0

ord_binary4_small[which(ord_diff4 < gene_count, arr.ind=TRUE)] <- 1
ord_binary4_small[which(ord_diff4 > gene_count, arr.ind=TRUE)] <- 0

ord_binary5_small[which(ord_diff5 < gene_count, arr.ind=TRUE)] <- 1
ord_binary5_small[which(ord_diff5 > gene_count, arr.ind=TRUE)] <- 0

ord_binary6_small[which(ord_diff6 < gene_count, arr.ind=TRUE)] <- 1
ord_binary6_small[which(ord_diff6 > gene_count, arr.ind=TRUE)] <- 0

# add binary
binary_sum_ord <- ord_binary2_small + ord_binary3_small + ord_binary4_small + ord_binary5_small + ord_binary6_small + ord_binary1_small


###############################
# explore data

grDevices::png(paste0(path,"/Multi-spp_all_ortho_", gene_count, "_heatmap_binarysum_ord.png"), height=5000, width=5000)
  print(
  RepeatOBserverV1::imagenan(binary_sum_ord,lnumr=numb_orthogroups,lnumc=numb_orthogroups)
  )
  grDevices::dev.off()
  
grDevices::png(paste0(path,"/Multi-spp_all_ortho_", gene_count, "_heatmap_syntenysum.png"), height=5000, width=5000)
  print(
  RepeatOBserverV1::imagenan(synteny_sum,lnumr=numb_orthogroups,lnumc=numb_orthogroups)
  )
  grDevices::dev.off()
  
# which ones found in all
binary_sum_ord6 <- binary_sum_ord
binary_sum_ord6[which(binary_sum_ord6 <6, arr.ind=TRUE)]<-NA

# too large?? - Error: cannot allocate vector of size 161.2 Mb - need allocation
grDevices::png(paste0(path,"/Multi-spp_all_ortho_", gene_count, "_heatmap_binarysum6_ord.png"), height=5000, width=5000)
  print(
  RepeatOBserverV1::imagenan(binary_sum_ord6,lnumr=numb_orthogroups,lnumc=numb_orthogroups)
  )
  grDevices::dev.off()
  
#----------------------------
inds_bin <- as.data.frame(which(binary_sum_ord >0, arr.ind=TRUE))
nrow(inds_bin)
# 1 453 934 with all orthogroups for within 40 genes

# 33 946/2 # on same chromosome in some - total 50k/2

inds_bin <- as.data.frame(which(synteny_sum < 50, arr.ind=TRUE))
nrow(inds_bin)
# 10 980 - 40 with all

# 730 have less than 50 summed
# 1342 have less than 150 summed

#----------------------
# how many ==2 in 70 gene distance
inds_bin <- as.data.frame(which(binary_sum_ord == 1, arr.ind=TRUE))
nrow(inds_bin)
# 17218/2
inds_bin <- as.data.frame(which(binary_sum_ord == 2, arr.ind=TRUE))
nrow(inds_bin)
#6990 /2
inds_bin <- as.data.frame(which(binary_sum_ord == 3, arr.ind=TRUE))
nrow(inds_bin)
#5044/2
inds_bin <- as.data.frame(which(binary_sum_ord == 4, arr.ind=TRUE))
nrow(inds_bin)
#2578/2
inds_bin <- as.data.frame(which(binary_sum_ord == 5, arr.ind=TRUE))
nrow(inds_bin)
# 1372/2 for first 500
# 19 446/2 for all, 19446/(6505*6505)= 0.00045955307 = 0.046%
inds_bin <- as.data.frame(which(binary_sum_ord == 6, arr.ind=TRUE))
nrow(inds_bin)
# 17144 with all

#1206
# with 150 gene distance

#--------------------------
# # how many ==2 in 150 gene distance
# inds_bin <- as.data.frame(which(binary_sum_ord == 1, arr.ind=TRUE))
# nrow(inds_bin)
# # 20 594/2
# inds_bin <- as.data.frame(which(binary_sum_ord == 2, arr.ind=TRUE))
# nrow(inds_bin)
# #12922/2
# inds_bin <- as.data.frame(which(binary_sum_ord == 3, arr.ind=TRUE))
# nrow(inds_bin)
# #10120/2
# inds_bin <- as.data.frame(which(binary_sum_ord == 4, arr.ind=TRUE))
# nrow(inds_bin)
# #5766/2
# inds_bin <- as.data.frame(which(binary_sum_ord == 5, arr.ind=TRUE))
# nrow(inds_bin)
# #2902/2
# inds_bin <- as.data.frame(which(binary_sum_ord == 6, arr.ind=TRUE))
# nrow(inds_bin)

############################
# what do these genes found near each other in all do

spp1genes_inall <- gene1_names_spp1[which(binary_sum_ord == 6, arr.ind=TRUE)]
spp2genes_inall <- gene1_names_spp2[which(binary_sum_ord == 6, arr.ind=TRUE)]

genes_inall<- cbind(spp1genes_inall, spp2genes_inall)

genes_inall <- as.data.frame(genes_inall)
colnames(genes_inall) <- c("spp1genes", "spp2genes")

#-------------------------
# load GO ont data

Gene_ont_file <- "Total_interproscan_output_edited2.tsv"
gene_ont <- read.delim(paste0(path,"/", Gene_ont_file), header = FALSE, sep = "\t", na.strings = "-", colClasses = c("character", "character", "character", "character"))

#wc -l Total_interproscan_output_edited2.tsv
#507 846 Total_interproscan_output_edited2.tsv

colnames(gene_ont) <- c("spp", "gene", "Pfam", "INTPRO", "descrip", "GOterm")
length(unique(gene_ont$INTPRO))
# 10201 # all 6spp

nrow(gene_ont)
#[1]  724 425 # getting read in properly, half of the rows missing before

#-----------------------------
# formatting Interproscan to have no duplicates of genes - one row per gene
library(dplyr)
library(tidyr)

# collapse GO terms
collapsed_go_terms1 <- gene_ont %>%
  group_by(spp, gene) %>%
  summarize(
    descrip = paste(unique(descrip), collapse = ","),  # Collapse unique descriptions
    GOterm = paste(unique(GOterm), collapse = ","),    # Collapse unique GO terms
    INTPRO = paste(unique(INTPRO), collapse = ","), .groups="keep"      # Collapse unique IPR terms
  )
  
# remove duplicate values in a list
cleaned_tibble <- collapsed_go_terms1 %>%
  separate_rows(GOterm, sep = ",") %>%  # Split the GOterm string into multiple rows
  separate_rows(INTPRO, sep = ",") %>%  # Split the INTPRO string into multiple rows
  separate_rows(descrip, sep = ",") %>%  # Split the descrip string into multiple rows
  filter(GOterm != "NA") %>%              # Remove rows with '-'
  distinct(gene, GOterm, INTPRO,descrip,.keep_all = TRUE ) #%>% # Keep unique terms with gene info

# recollapse 
collapsed_go_terms2 <- cleaned_tibble %>%
  group_by(spp, gene) %>%
  summarize(
    INTPRO = paste(sort(unique(INTPRO)), collapse = ","),           # Collapse IPR terms with unique values
    descrip = paste(sort(unique(descrip)), collapse = ","),        # Collapse descriptions with unique values
    GOterm = paste(sort(unique(GOterm)), collapse = ",") , .groups="keep"         # Collapse GO terms with unique values
  )

# format for ermineJ
collapsed_go_terms <- collapsed_go_terms2 %>%
  #mutate(gene2 = gene) %>%
  select(spp, gene, everything())

collapsed_go_terms_df <- as.data.frame(collapsed_go_terms)

gene_ont <- collapsed_go_terms_df

# write out file
utils::write.table(x=gene_ont , file=paste0(path,"/Total_interproscan_output_edited3.tsv"), append = FALSE, sep = "\t", dec = ".", row.names = FALSE, col.names = TRUE)

#----------------------
# subset for spp 1 and join to data
unique(gene_ont$spp)
gene_ont_spp1 <- gene_ont[which(gene_ont$spp=="Oxyria_digyna_H1_interproscan_output.tsv"),]
colnames(gene_ont_spp1) <- c("spp", "spp1genes", "INTPRO", "descrip", "GOterm")
genes_inall_go1 <- dplyr::left_join(genes_inall, gene_ont_spp1, by="spp1genes")

#-----------------------------------
# subset for spp 1 and join to data
unique(gene_ont$spp)
gene_ont_spp2 <- gene_ont[which(gene_ont$spp=="Dryas_octopetala_interproscan_output.tsv"),]
colnames(gene_ont_spp2) <- c("spp", "spp2genes", "INTPRO", "descrip", "GOterm")
genes_inall_go12 <- dplyr::left_join(genes_inall_go1, gene_ont_spp2, by="spp2genes")

#---------------------------------
unique(genes_inall_go12$INTPRO.x)
unique(genes_inall_go12$INTPRO.y)

length(unique(genes_inall_go12$INTPRO.x))
length(unique(genes_inall_go12$INTPRO.y))
# 3332 Oxyria and Dryas for all
# 3290

unique(genes_inall_go12$descrip.x)
unique(genes_inall_go12$descrip.y)

# write out genes_inall_go12
utils::write.table(x=genes_inall_go12, file=paste0(path,"/All_spp_gene_ontology_spp1_spp2.txt"), append = FALSE, sep = "\t", dec = ".", row.names = TRUE, col.names = TRUE)


#############################
# what is the probability distribution when randomly arranging 30000 objects in a line that
# two objects fall within 70 objects of each other 5 times in a row

#The final probability of arranging 30,000 objects in a line such that two specific objects
# fall within 70 objects of each other five times in a row is very low, roughly 5.3×10−125.3×10−12.

##############################
# Add Arctic spp and add non-Arctic species

# any comb only in one or other?

Oxy_Dryas <- ord_binary1_small + ord_binary2_small + ord_binary6_small
Rhe_Pol <- ord_binary3_small + ord_binary4_small + ord_binary5_small

Oxy_Dryas[which(Oxy_Dryas != 3, arr.ind=TRUE)] <- 0
Rhe_Pol[which(Rhe_Pol != 3, arr.ind=TRUE)] <- 0

Arctic_NonArctic <- Oxy_Dryas - Rhe_Pol

grDevices::png(paste0(path,"/Multi-spp_", numb_orthogroups,"_", gene_count, "_heatmap_synteny_Arctic-non.png"), height=5000, width=5000)
  print(
  RepeatOBserverV1::imagenan(Arctic_NonArctic,lnumr=numb_orthogroups,lnumc=numb_orthogroups)
  )
  grDevices::dev.off()


Arctic_NonArctic[which(Arctic_NonArctic ==0, arr.ind=TRUE)] <- NA

grDevices::png(paste0(path,"/Multi-spp_", numb_orthogroups,"_", gene_count, "_heatmap_synteny_Arctic-non_rm0.png"), height=5000, width=5000)
  print(
  RepeatOBserverV1::imagenan(Arctic_NonArctic,lnumr=numb_orthogroups,lnumc=numb_orthogroups)
  )
  grDevices::dev.off()

# count
inds_bin <- as.data.frame(which(Arctic_NonArctic >0, arr.ind=TRUE))
nrow(inds_bin)

# 8910 for all orthogroups Arctic specific

# what is the probability distribution when randomly arranging 30000 
# objects in a line that two objects fall within 70 objects of each other 5 times in a row

inds_bin <- as.data.frame(which(Arctic_NonArctic <0, arr.ind=TRUE))
nrow(inds_bin)
# 17974 - Arabis, Polygonum, Rheum but not Arctic spp

######################################
# Plot the most extreme values on the chromosomes of one species
# get the genes that are extreme

# matrix to list
ord_diff1_ls <- base::as.list(ord_diff1)
gene1_names_spp1_ls <- base::as.list(gene1_names_spp1)
gene2_names_spp1_ls <- base::as.list(gene2_names_spp1)
Arctic_NonArctic_ls <- base::as.list(Arctic_NonArctic)

# unlist and cbind
spp1_data <- cbind(as.numeric(unlist(ord_diff1_ls)), unlist(gene1_names_spp1_ls),  unlist(gene2_names_spp1_ls), as.numeric(unlist(Arctic_NonArctic_ls)))

# match colnames
colnames(spp1_data) <- c("ord_diff", "gene",  "geney", "arctic_score")

spp1_data <- as.data.frame(spp1_data)
gene_pos <- as.data.frame(gene_pos)

# remove NA rows
spp1_data<- spp1_data[which(spp1_data$arctic_score ==3, arr.ind=TRUE),]

# join with gene start pos for gene 1 spp1 and gene 1 spp2
colnames(gene_pos) <- c("gene", "chr", "start", "ord", "og", "spp")
spp1_data_starts0 <- dplyr::left_join(spp1_data, gene_pos, by="gene")

colnames(gene_pos) <- c("geney", "chry", "starty", "ordy", "ogy", "sppy")
spp1_data_starts <- dplyr::left_join(spp1_data_starts0, gene_pos, by="geney")

spp1_data_starts$start <- as.numeric(spp1_data_starts$start)
spp1_data_starts$ord_diff <- as.numeric(spp1_data_starts$ord_diff)
spp1_data_starts$chr <- as.factor(spp1_data_starts$chr)

# colour by go term here
grDevices::png(paste0(path,"/Arctic_spp_", numb_orthogroups, "_", gene_count,"_synteny_on_genome.png"), height=2500, width=5000)
  print(
    ggplot2::ggplot(data=spp1_data_starts, ggplot2::aes(x=start, y=ord_diff))+
      ggplot2::geom_point(ggplot2::aes(x=start, y=ord_diff))+
      ggplot2::facet_wrap(~chr, scales = "free")+
      ggplot2::theme_classic()
  )
  grDevices::dev.off()



grDevices::png(paste0(path,"/Arctic_spp_", numb_orthogroups, "_", gene_count,"_synteny_on_genome_histogram.png"), height=2500, width=5000)
  print(
    ggplot2::ggplot(data=spp1_data_starts, ggplot2::aes(x=start))+
      ggplot2::geom_histogram(binwidth=1e5)+
      ggplot2::facet_wrap(~chr, scales = "free")+
      ggplot2::theme_classic()
  )
  grDevices::dev.off()

######################################
# Plot the most extreme values on the chromosomes of the other species
# get the genes that are extreme

# matrix to list
ord_diff2_ls <- base::as.list(ord_diff2)
gene1_names_spp2_ls <- base::as.list(gene1_names_spp2)
gene2_names_spp2_ls <- base::as.list(gene2_names_spp2)
Arctic_NonArctic_ls <- base::as.list(Arctic_NonArctic)

# unlist and cbind
spp2_data <- cbind(as.numeric(unlist(ord_diff2_ls)), unlist(gene1_names_spp2_ls), unlist(gene2_names_spp2_ls), as.numeric(unlist(Arctic_NonArctic_ls)))

# match colnames
colnames(spp2_data) <- c("ord_diff", "gene", "geney", "arctic_score")
colnames(gene_pos) <- c("gene", "chr", "start", "ord", "og", "spp")

spp2_data <- as.data.frame(spp2_data)
gene_pos <- as.data.frame(gene_pos)

# remove NA rows
spp2_data<- spp2_data[which(spp2_data$arctic_score ==3, arr.ind=TRUE),]

# join with gene start pos for gene 1 spp1 
colnames(gene_pos) <- c("gene", "chr", "start", "ord", "og", "spp")
spp2_data_starts0 <- dplyr::left_join(spp2_data, gene_pos, by="gene")

colnames(gene_pos) <- c("geney", "chry", "starty", "ordy", "ogy", "sppy")
spp2_data_starts <- dplyr::left_join(spp2_data_starts0, gene_pos, by="geney")

# formatting 
spp2_data_starts$start <- as.numeric(spp2_data_starts$start )
spp2_data_starts$ord_diff <- as.numeric(spp2_data_starts$ord_diff)
spp2_data_starts$chr <- as.factor(spp2_data_starts$chr)


grDevices::png(paste0(path,"/Arctic_spp2_", numb_orthogroups, "_", gene_count,"_synteny_on_genome.png"), height=2500, width=5000)
  print(
    ggplot2::ggplot(data=spp2_data_starts, ggplot2::aes(x=start, y=ord_diff))+
      ggplot2::geom_point(ggplot2::aes(x=start, y=ord_diff))+
      ggplot2::facet_wrap(~chr, scales = "free")+
      ggplot2::theme_classic()
  )
  grDevices::dev.off()



grDevices::png(paste0(path,"/Arctic_spp2_", numb_orthogroups, "_", gene_count,"_synteny_on_genome_histogram.png"), height=2500, width=5000)
  print(
    ggplot2::ggplot(data=spp2_data_starts, ggplot2::aes(x=start))+
      ggplot2::geom_histogram(binwidth=1e5)+
      ggplot2::facet_wrap(~chr, scales = "free")+
      ggplot2::theme_classic()
  )
  grDevices::dev.off()




##########################################
# what do these genes found near each other in clusters do

# Gene ontology data - gene_ont file

# subset for spp 1 and join to data
unique(gene_ont$spp)
gene_ont_spp1 <- gene_ont[which(gene_ont$spp=="Oxyria_digyna_H1_interproscan_output.tsv"),]
colnames(gene_ont_spp1) <- c("spp", "gene", "INTPRO", "descrip", "GOterm")
genes_Arctic_go1 <- dplyr::left_join(spp1_data, gene_ont_spp1, by="gene")

#unique(genes_Arctic_go1$INTPRO)
#unique(genes_Arctic_go1$descrip)
# 643 go terms
colnames(gene_pos) <- c("gene", "chr", "start", "ord", "og", "spp")
genes_Arctic_go1_starts0 <- dplyr::left_join(genes_Arctic_go1, gene_pos, by="gene")
colnames(gene_pos) <- c("geney", "chry", "starty", "ordy", "ogy", "sppy")
genes_Arctic_go1_starts <- dplyr::left_join(genes_Arctic_go1_starts0, gene_pos, by="geney")

# write out genes_Arctic_go1
utils::write.table(x=genes_Arctic_go1_starts, file=paste0(path,"/Arctic_gene_ontology_starts_spp1.txt"), append = FALSE, sep = "\t", dec = ".", row.names = FALSE, col.names = TRUE)


#-----------------------------
# plot GO terms as colours on points plot

genes_Arctic_go1_starts$start <- as.numeric(genes_Arctic_go1_starts$start )
genes_Arctic_go1_starts$ord_diff <- as.numeric(genes_Arctic_go1_starts$ord_diff)
genes_Arctic_go1_starts$chr <- as.factor(genes_Arctic_go1_starts$chr)
genes_Arctic_go1_starts$GOterm <- as.factor(genes_Arctic_go1_starts$GOterm)
genes_Arctic_go1_starts$INTPRO <- as.factor(genes_Arctic_go1_starts$INTPRO)

library(ggplot2)
library(RColorBrewer)

# Create a color palette based on the number of unique GO terms
unique_terms <- unique(genes_Arctic_go1_starts$INTPRO)
num_terms <- length(unique_terms)
palette <- brewer.pal(n = min(num_terms, 12), name = "Set3")  # Adjust if more than 12 terms
if (num_terms > 12) {
  palette <- colorRampPalette(brewer.pal(12, "Set3"))(num_terms)  # Generate more colors
}

# Create a named vector for color mapping
color_mapping <- setNames(palette, unique_terms)


grDevices::png(paste0(path,"/Arctic_spp_", numb_orthogroups, "_", gene_count,"_synteny_on_genome_colourINTERPRO.png"), height=2500, width=5000)
print(
  ggplot2::ggplot(data = genes_Arctic_go1_starts, ggplot2::aes(x = start, y = ord_diff)) +
    geom_point(ggplot2::aes(color = as.factor(INTPRO))) +  # Map colors to GO terms
    ggplot2::facet_wrap(~chr, scales = "free") +
    ggplot2::scale_color_manual(values = color_mapping) +  # Use the color mapping
    ggplot2::theme_classic()+
	theme(legend.position = "none") 

)
  grDevices::dev.off()


##########################
# Run same for spp2

# subset for spp 1 and join to data
unique(gene_ont$spp)
gene_ont_spp2 <- gene_ont[which(gene_ont$spp=="Dryas_octopetala_interproscan_output.tsv"),]
colnames(gene_ont_spp2) <- c("spp", "gene", "INTPRO", "descrip", "GOterm")
genes_Arctic_go2 <- dplyr::left_join(spp2_data, gene_ont_spp2, by="gene")

#unique(genes_Arctic_go2$INTPRO)
#unique(genes_Arctic_go1$descrip)
# 643 go terms

colnames(gene_pos) <- c("gene", "chr", "start", "ord", "og", "spp")
genes_Arctic_go2_starts0 <- dplyr::left_join(genes_Arctic_go2, gene_pos, by="gene")
colnames(gene_pos) <- c("geney", "chry", "starty", "ordy", "ogy", "sppy")
genes_Arctic_go2_starts <- dplyr::left_join(genes_Arctic_go2_starts0, gene_pos, by="geney")

# write out genes_Arctic_go1
utils::write.table(x=genes_Arctic_go2_starts, file=paste0(path,"/Arctic_gene_ontology_starts_spp2.txt"), append = FALSE, sep = "\t", dec = ".", row.names = FALSE, col.names = TRUE)


#-----------------------------
# plot GO terms as colours on points plot

genes_Arctic_go2_starts$start <- as.numeric(genes_Arctic_go2_starts$start )
genes_Arctic_go2_starts$ord_diff <- as.numeric(genes_Arctic_go2_starts$ord_diff)
genes_Arctic_go2_starts$chr <- as.factor(genes_Arctic_go2_starts$chr)
genes_Arctic_go2_starts$GOterm <- as.factor(genes_Arctic_go2_starts$GOterm)
genes_Arctic_go2_starts$INTPRO <- as.factor(genes_Arctic_go2_starts$INTPRO)

library(ggplot2)
library(RColorBrewer)

# Create a color palette based on the number of unique GO terms
unique_terms <- unique(genes_Arctic_go2_starts$INTPRO)
num_terms <- length(unique_terms)
palette <- brewer.pal(n = min(num_terms, 12), name = "Set3")  # Adjust if more than 12 terms
if (num_terms > 12) {
  palette <- colorRampPalette(brewer.pal(12, "Set3"))(num_terms)  # Generate more colors
}

# Create a named vector for color mapping
color_mapping <- setNames(palette, unique_terms)


grDevices::png(paste0(path,"/Arctic_spp2_", numb_orthogroups, "_", gene_count,"_synteny_on_genome_colourINTERPRO.png"), height=2500, width=5000)
print(
  ggplot2::ggplot(data = genes_Arctic_go2_starts, ggplot2::aes(x = start, y = ord_diff)) +
    geom_point(ggplot2::aes(color = as.factor(INTPRO))) +  # Map colors to GO terms
    ggplot2::facet_wrap(~chr, scales = "free") +
    ggplot2::scale_color_manual(values = color_mapping) +  # Use the color mapping
    ggplot2::theme_classic()+
	theme(legend.position = "none") 

)
  grDevices::dev.off()

##############################################
# Plot one genome against the other genome and see where they match
# genes_Arctic_go1_starts vs genes_Arctic_go1_starts 

colnames(genes_Arctic_go1_starts) <- c("ord_diff1","gene1", "geney1" , "arctic_score1", "spp.x1",
"INTPRO1","descrip1" ,"GOterm1" ,"chr1" ,"start1","ord1","og1" ,"spp.y1" , "chry1","starty1","ordy1", "ogy1" ,"sppy1")

colnames(genes_Arctic_go2_starts) <- c("ord_diff2","gene2", "geney2" , "arctic_score2", "spp.x2",
"INTPRO2","descrip2" ,"GOterm2" ,"chr2" ,"start2","ord2","og2" ,"spp.y2" , "chry2","starty2","ordy2", "ogy2" ,"sppy2")

genes_Arctic_go12_starts <- cbind(genes_Arctic_go1_starts, genes_Arctic_go2_starts)

genes_Arctic_go12_starts$start1 <- as.numeric(genes_Arctic_go12_starts$start1)
genes_Arctic_go12_starts$start2 <- as.numeric(genes_Arctic_go12_starts$start2)
genes_Arctic_go12_starts$GOterm2 <- as.factor(genes_Arctic_go12_starts$GOterm2)
genes_Arctic_go12_starts$INTPRO2 <- as.factor(genes_Arctic_go12_starts$INTPRO2)

# write out genes_Arctic_go12
utils::write.table(x=genes_Arctic_go12_starts, file=paste0(path,"/Arctic_gene_ontology_starts_spp1and2.txt"), append = FALSE, sep = "\t", dec = ".", row.names = FALSE, col.names = TRUE)

#------------------------------
# alternatively read in file
tmux new-session -s synteny1
tmux attach-session -t synteny1

cd /home/celphin/scratch/Oxyria/synteny_quantity

salloc -c1 --time 2:00:00 --mem 191000m --account def-rieseber

module load  StdEnv/2020 r/4.2.2
export R_LIBS_USER=/home/celphin/R/x86_64-pc-linux-gnu-library/4.1.0/

R

library(ggplot2)
library(RColorBrewer)

path="/home/celphin/scratch/Oxyria/synteny_quantity"
numb_orthogroups=500
comb_file="Total_combBed.txt"

spp_hap1="Oxyria_digyna_H1"
spp_hap2="Dryas_octopetala"
spp_hap3="Rheum_nobile_H0"
spp_hap4="Arabis_alpina"
spp_hap5="Polygunum_aviculare_H0"
spp_hap6="Draba_nivalis"

genes_Arctic_go12_starts <- base::as.data.frame(utils::read.table(paste0(path,"/Arctic_gene_ontology_starts_spp1and2.txt"),sep="\t", header = TRUE, check.names = FALSE))
gene_ont <- base::as.data.frame(utils::read.table(paste0(path,"/Total_interproscan_output_edited3.tsv"),  sep="\t", header = TRUE, check.names = FALSE))


#-----------------------------
# join GOtermsy for second gene in every pair
colnames(genes_Arctic_go12_starts)
colnames(gene_ont) <- c("sppy1",  "geney1",    "INTPROy1",  "descripy1", "GOtermy1")
genes_Arctic_second_go0 <- dplyr::left_join(genes_Arctic_go12_starts, gene_ont, by="geney1")
colnames(gene_ont) <- c("sppy2",  "geney2",    "INTPROy2",  "descripy2", "GOtermy2")
genes_Arctic_second_go <- dplyr::left_join(genes_Arctic_second_go0, gene_ont, by="geney2")


# write out genes_Arctic_go12
utils::write.table(x=genes_Arctic_second_go, file=paste0(path,"/Arctic_gene_ontology_starts_spp1and2_bothgo.txt"), append = FALSE, sep = "\t", dec = ".", row.names = FALSE, col.names = TRUE)

#----------------------------------
# Set up the plot to plot with text
grDevices::png(paste0(path, "/Arctic_spp1_spp2_genomevsgenome_allGOTermstext.png"), height=4000, width=7000)
print(
  ggplot2::ggplot(data = genes_Arctic_go12_starts, ggplot2::aes(x = start1, y = start2)) +
    geom_point(size = 2) +  # Points can be sized for visibility
    geom_text(ggplot2::aes(label = GOterm2), hjust = 0.5, vjust = 0.5, size = 2.5, check_overlap = TRUE) +  # Smaller, centered text
    ggplot2::facet_wrap(~chr1 + chr2, scales = "free") +
    ggplot2::theme_classic() +
    theme(legend.position = "none")
)
grDevices::dev.off()

# Set up the plot to plot with text
grDevices::png(paste0(path, "/Arctic_spp1_spp2_genomevsgenome_alldescriptext.png"), height=4000, width=7000)
print(
  ggplot2::ggplot(data = genes_Arctic_go12_starts, ggplot2::aes(x = start1, y = start2)) +
    geom_point(size = 2) +  # Points can be sized for visibility
    geom_text(ggplot2::aes(label = descrip2), hjust = 0.5, vjust = 0.5, size = 2.5, check_overlap = TRUE) +  # Smaller, centered text
    ggplot2::facet_wrap(~chr1 + chr2, scales = "free") +
    ggplot2::theme_classic() +
    theme(legend.position = "none")
)
grDevices::dev.off()

#########################################
# Define the GO terms you want to highlight
temp_terms <- c("GO:1900034", "GO:0009409", "GO:0042752")
light_terms <- c("GO:0071486", "GO:0009409", "GO:0042752", "GO:0015995", "GO:0033014", "GO:0015979", "GO:0016628", "GO:0045550", "GO:0071949", "GO:0010207", "GO:0009785", "GO:0009767", "GO:0009765", "GO:0009639")
water_terms <- c("GO:0048364", "GO:0010374", "GO:0010274", "GO:0010090", "GO:0010088", "GO:0009415", "GO:0009269")
growth_terms <- c("GO:0045927", "GO:0009734", "GO:0009733")
pollen_terms <- c("GO:0007142", "GO:0010183")
biotic_defense <- c("GO:0050832", "GO:0006952", "GO:1900150", "GO:0098542")
metal_terms <- c("GO:0006826", "GO:0010038")

highlighted_terms <- c(water_terms, temp_terms, light_terms, growth_terms, pollen_terms, biotic_defense, metal_terms)

# Create a color palette based on the number of unique highlighted GO terms
num_highlighted_terms <- length(unique(highlighted_terms))
palette <- RColorBrewer::brewer.pal(n = min(num_highlighted_terms, 12), name = "Set3")
if (num_highlighted_terms > 12) {
  palette <- colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(num_highlighted_terms)
}

# Create a named vector for color mapping
color_mapping <- setNames(rep("grey", length(unique_terms)), unique_terms)  # Default all to grey

# Assign colors to the highlighted terms
color_mapping[highlighted_terms] <- palette[seq_along(highlighted_terms)]

# Create a logical column for highlighted terms
genes_Arctic_go12_starts$highlighted <- sapply(genes_Arctic_go12_starts$GOterm2, function(go) {
  any(highlighted_terms %in% go)
})

# Filter the data to include only highlighted terms
genes_filtered <- genes_Arctic_go12_starts[genes_Arctic_go12_starts$highlighted, ]

# Set up the plot
grDevices::png(paste0(path, "/Arctic_spp1_spp2_genomevsgenome_colourspecificGO.png"), height=2500, width=5000)
print(
  ggplot2::ggplot(data = genes_filtered, ggplot2::aes(x = start1, y = start2)) +
    ggplot2::geom_point(ggplot2::aes(color = as.factor(GOterm2))) +  # Map colors to GO terms
    ggplot2::facet_wrap(~chr1 + chr2, scales = "free") +
    ggplot2::scale_color_manual(values = color_mapping) +  # Use the color mapping
    ggplot2::theme_classic() +
    theme(legend.position = "none")
)
grDevices::dev.off()

# Print the color mapping
print(color_mapping)

# Set up the plot to plot with text
grDevices::png(paste0(path, "/Arctic_spp1_spp2_genomevsgenome_GOTermstext.png"), height=4000, width=7000)
print(
  ggplot2::ggplot(data = genes_filtered, ggplot2::aes(x = start1, y = start2)) +
    geom_point(size = 2) +  # Points can be sized for visibility
    geom_text(ggplot2::aes(label = GOterm2), hjust = 0.5, vjust = 0.5, size = 2.5, check_overlap = TRUE) +  # Smaller, centered text
    ggplot2::facet_wrap(~chr1 + chr2, scales = "free") +
    ggplot2::theme_classic() +
    theme(legend.position = "none")
)
grDevices::dev.off()


#-------------------------------
# Define highlighted descriptions
highlighted_descrip <- c("Heat", "HEAT", "Aquaporin", "light", "cold", "Cold", "temp", "Homeobox")

# Create a logical column for highlighted terms using grepl for partial matching
genes_Arctic_go12_starts$highlighted <- sapply(genes_Arctic_go12_starts$descrip1, function(descrip) {
  any(grepl(paste(highlighted_descrip, collapse = "|"), descrip, ignore.case = TRUE))
})

# Filter the data to include only highlighted terms
genes_filtered <- genes_Arctic_go12_starts[genes_Arctic_go12_starts$highlighted, ]

# Set up the plot to plot with text
grDevices::png(paste0(path, "/Arctic_spp1_spp2_genomevsgenome_descrip_text.png"), height=4000, width=7000)
print(
  ggplot2::ggplot(data = genes_filtered, ggplot2::aes(x = start1, y = start2)) +
    geom_point(size = 2) +  # Points can be sized for visibility
    geom_text(ggplot2::aes(label = descrip1), hjust = 0.5, vjust = 0.5, size = 2.5, check_overlap = TRUE) +  # Smaller, centered text
    ggplot2::facet_wrap(~chr1 + chr2, scales = "free") +
    ggplot2::theme_classic() +
    theme(legend.position = "none")
)
grDevices::dev.off()


#------------------------------
# Clusters - share the same chromosome

go_terms_summary <- genes_Arctic_go12_starts %>%
  group_by(chr1, chr2) %>%
  summarize(start1 = list(unique(start1)),  go_terms2 = list(unique(GOterm2)),.groups = 'drop')  # Use list to store unique GO terms

# View the result
print(go_terms_summary)

# Convert the list of GO terms into a single string for each row
go_terms_summary <- go_terms_summary %>%
  mutate(go_terms2 = sapply(go_terms2, function(x) paste(x, collapse = ","))) %>%# Join GO terms into a single string
  mutate(start1 = sapply(start1, function(x) paste(x, collapse = ",")))

# Define the output file path
output_file_path <- paste0(path, "/GOterms_per_chr_combo.txt")  # Update your path

# Write the modified data frame to a text file
utils::write.table(x = go_terms_summary, file = output_file_path, append = FALSE, sep = "\t", dec = ".", row.names = FALSE, col.names = TRUE)

#################################################
# What are the gene clusters that are found in many Arctic species?

#Function to find peaks in histogram
find_peaks_in_hist <- function(data) {
  hist_data <- hist(data, plot = FALSE)
  peaks <- pracma::findpeaks(hist_data$counts)
  peak_locations <- hist_data$mids[peaks[, 2]]
  return(peak_locations)
}

# Apply the peak finding function to each group
peak_results <- lapply(split(genes_Arctic_go1_starts$start1, genes_Arctic_go1_starts$chr1), find_peaks_in_hist)

peak_results
names(peak_results)

spp1_chr6_50_60Mbp <- genes_Arctic_go1_starts[which((genes_Arctic_go1_starts$chr1 == "Oxyrt-6-73303751") &  (genes_Arctic_go1_starts$start1 > 5.0e+07) &  (genes_Arctic_go1_starts$start1 < 6.0e+07)),]
utils::write.table(x=spp1_chr6_50_60Mbp, file=paste0(path,"/spp1_chr6_50_60Mbp.txt"), append = FALSE, sep = "\t", dec = ".", row.names = FALSE, col.names = TRUE)

#----------------------------

# Apply the peak finding function to each group
peak_results <- lapply(split(genes_Arctic_go2_starts$start2, genes_Arctic_go2_starts$chr2), find_peaks_in_hist)

peak_results
names(peak_results)

spp2_chr1_17_18Mbp <- genes_Arctic_go2_starts[which((genes_Arctic_go2_starts$chr2 == "DoctH0-1") &  (genes_Arctic_go2_starts$start2 >  15000000) &  (genes_Arctic_go2_starts$start2 <  20000000)),]
utils::write.table(x=spp2_chr1_17_18Mbp, file=paste0(path,"/spp2_chr1_17_18Mbp.txt"), append = FALSE, sep = "\t", dec = ".", row.names = FALSE, col.names = TRUE)

spp2_chr6_6_8Mbp <- genes_Arctic_go2_starts[which((genes_Arctic_go2_starts$chr2 == "DoctH0-6") &  (genes_Arctic_go2_starts$start2 >  6000000) &  (genes_Arctic_go2_starts$start2 <  8000000)),]
utils::write.table(x=spp2_chr6_6_8Mbp, file=paste0(path,"/spp2_chr6_6_8Mbp.txt"), append = FALSE, sep = "\t", dec = ".", row.names = FALSE, col.names = TRUE)

spp2_chr8_4_6Mbp <- genes_Arctic_go2_starts[which((genes_Arctic_go2_starts$chr2 == "DoctH0-8") &  (genes_Arctic_go2_starts$start2 >  4000000) &  (genes_Arctic_go2_starts$start2 <  6000000)),]
utils::write.table(x=spp2_chr8_4_6Mbp, file=paste0(path,"/spp2_chr8_4_6Mbp.txt"), append = FALSE, sep = "\t", dec = ".", row.names = FALSE, col.names = TRUE)

spp2_chr3_9_11Mbp <- genes_Arctic_go2_starts[which((genes_Arctic_go2_starts$chr2 == "DoctH0-3") &  (genes_Arctic_go2_starts$start2 >  9000000) &  (genes_Arctic_go2_starts$start2 <  11000000)),]
utils::write.table(x=spp2_chr3_9_11Mbp, file=paste0(path,"/spp2_chr3_9_11Mbp.txt"), append = FALSE, sep = "\t", dec = ".", row.names = FALSE, col.names = TRUE)

########################
# Look at GO terms in Revigo
# http://revigo.irb.hr/Results?jobid=424806203

tmux new-session -s Enrichment
tmux attach-session -t Enrichment

#---------------------------------
# Prepare gene score file for Arctic synteny pairs
# Arctic_gene_ontology_starts_spp1and2_bothgo.txt

cd /home/celphin/scratch/Oxyria/synteny_quantity
mkdir enrichment; cd enrichment

cp ../Arctic_gene_ontology_starts_spp1and2_bothgo.txt .
cp ~/scratch/Oxyria/CAFE/enrichment_analysis/Oxydig_interproscan_edited.tsv .
cp ~/scratch/Oxyria/CAFE/enrichment_analysis/Dryasoct_interproscan_edited.tsv .
cp ~/scratch/Oxyria/CAFE/enrichment_analysis/Oxydig_GO_mappings.ermineJ.txt .
cp ~/scratch/Oxyria/CAFE/enrichment_analysis/Dryasoct_GO_mappings.ermineJ.txt .

#---------------------------------------
cd /home/celphin/scratch/Oxyria/synteny_quantity/enrichment/

# first make geneID list
awk 'BEGIN{FS="\t"}{print $2}' Arctic_gene_ontology_starts_spp1and2_bothgo.txt > Arctic_synteny_genelist1.txt
sed "s/\"//g" Arctic_synteny_genelist1.txt |sort -u >  Oxydig_synteny_genelist.txt
awk 'BEGIN{FS="\t"}{print $20}' Arctic_gene_ontology_starts_spp1and2_bothgo.txt > Arctic_synteny_genelist2.txt
sed "s/\"//g" Arctic_synteny_genelist2.txt  |sort -u > Dryasoct_synteny_genelist.txt

for taxon in Oxydig Dryasoct ; \
do cat "$taxon"_interproscan_edited.tsv |sed 's/ //g' | awk 'BEGIN{FS="\t"}{print $2,"0"}' |sort -u > "$taxon"_total_genesets ; done

wc -l *_total_genesets
# 19364 Dryasoct_total_genesets
# 20406 Oxydig_total_genesets

# change 0 to 1 for given genes
for taxon in Oxydig Dryasoct ; \
do cat "$taxon"_synteny_genelist.txt | \
 while read gene ; do sed -i "s/$gene 0/$gene 1/g" "$taxon"_total_genesets ; done ; done

for taxon in Oxydig Dryasoct ; \
do sed -i 's/ /\t/g' "$taxon"_total_genesets; done

#-----------------------------
salloc -c1 --time 3:00:00 --mem 120000m --account def-rieseber

cd /home/celphin/scratch/Oxyria/synteny_quantity/enrichment

# https://erminej.msl.ubc.ca/help/tutorials/erminej-cli/
ERMINEJ_HOME=/home/celphin/ermineJ-3.2
export JAVA_HOME=/cvmfs/soft.computecanada.ca/easybuild/software/2020/Core/java/13.0.2/

#-------------------
# find enriched GO terms

for taxon in Oxydig Dryasoct  ; do $ERMINEJ_HOME/bin/ermineJ.sh \
-a "$taxon"_GO_mappings.ermineJ.txt \
-s "$taxon"_total_genesets \
-c /home/celphin/ermineJ.data/go.obo \
--genesOut -aspects BCM  \
-o "$taxon"_Arctic_synteny_genesets.ermine.results -y 5 -b ; done

# simplify list of enriched terms:  http://revigo.irb.hr/

grep "!" *_Arctic_synteny_genesets.ermine.results  | awk -F '\t' '$7 <0.05 { print }' | awk  -F $'\t' '{print $3}' \
| sort | uniq -c | awk '$1 > 1' | awk  -F " " '{print $2, $1}' > Revigio_Arctic_synteny.txt

cat *_Arctic_synteny_genesets.ermine.results | awk -F '\t' '$7 <1  { print $3 , $7}'  > Revigio_Arctic_synteny_pval.txt


##########################################################
# Plot GO term data

# try plotting on genome 

module load  StdEnv/2020 r/4.2.2
export R_LIBS_USER=/home/celphin/R/x86_64-pc-linux-gnu-library/4.1.0/

R

file_path <- "/home/celphin/scratch/Oxyria/synteny_quantity/enrichment/Oxydig_Arctic_synteny_genesets.ermine.results"

# Read the file into R
mydata <- read.table(file_path, 
                   header = TRUE,       # Use TRUE if the first line contains headers
                   sep = "\t",         # Use tab as a separator (adjust if different)
                   skip = 25,           # Skip the first 25lines 
                   fill = TRUE,        # Fill in missing values
                   stringsAsFactors = FALSE)  # Prevent automatic conversion of strings to factors

# View the imported data
head(mydata)

colnames(mydata) <- c("X", "descrip", "GOterm", "probescore", "genescore", "score", "pval", "cor-pval", "MFP", "cor-MFP", "X1", "genes", "X2")
GO_genes <- cbind(mydata$GOterm, mydata$descrip, mydata$pval, mydata$genes) 


#----------------
# Read in gene start positions
comb_file = "/home/celphin/scratch/Oxyria/synteny_quantity/Total_combBed.txt"
combined_genes <- base::as.data.frame(utils::read.table(comb_file, sep="\t", header = TRUE, check.names = FALSE))
colnames(combined_genes) <- c( "chr", "start", "end", "gene", "ofID", "pepLen", "ord", "genome", 
"arrayID","isArrayRep", "orthogroup", "globHOG", "noAnchor", "og")
gene_pos <- cbind(combined_genes$gene, combined_genes$chr, combined_genes$start, combined_genes$ord, combined_genes$orthogroup, combined_genes$genome)

#---------------------
# subset for enriched GO terms


#--------------------
# join enriched genes with start values and chr




#------------------------
# plot genes coloured by GO terms, y value = p-val,  on genome



###############################################
module load  StdEnv/2020 r/4.2.2
export R_LIBS_USER=/home/celphin/R/x86_64-pc-linux-gnu-library/4.1.0/

R

library(GOplot)

chord_dat(data, genes, process)

# Generate the plotting object
circ <- circle_dat(EC$david, EC$genelist)





##################################################
# Exploring

# Dryas1/Oxy2 N-terminal,GBF-interacting protein 1,UBA-like superfamily	Cold-regulated 413 protein
# 6mbp and 75Mbp

#Cold-regulated 413 protein	CCG-binding protein 1



# Oxy6 30-40
# Dry2 25-30

# large region
# GO:0010090 - Zinc finger C2H2 superfamily,Zinc finger C2H2-type,Zinc finger protein GIS3/ZFP5/ZFP6
# GO:0005975 -  all-beta, C-terminal beta-sheet, catalytic domain, family 13,Alpha amylase,Alpha-amylase,Glycoside hydrolase superfamily,Glycosyl hydrolase
# GO:0016192 - coiled-coil homology domain,Longin domain,Longin-like domain superfamily,Synaptobrevin-like,v-SNARE
# GO:0004672 - active site, ATP binding site,Protein kinase,Protein kinase domain,Protein kinase-like domain superfamily,Serine/threonine-protein kinase
# GO:0008235 - M28 peptidase domain,Endoplasmic reticulum metallopeptidase 1-like,Peptidase M28,Peptidase M28 family




# Oxy4 0-20
# Dry9 10-20

# Oxy1 38-50
# Dry2 25-30

# Oxy1 15-17
# Dry2 20-25

# Oxy2 59-65
# Dry6 6-10

# Oxy2 75-80
# Dry6 4-6

# Oxy2 0-15
# Dry9 5-15

# heat
Oxy1 and Dry8
GO:0003700,GO:0006355,GO:0043565 - DNA-binding,Heat shock factor (HSF)-type,Heat shock transcription factor family,Winged helix DNA-binding domain superfamily,Winged helix-like DNA-binding domain superfamily
GO:0015631 -  TOG domain, type 2,Armadillo-like helical,Armadillo-type fold,CLASP N-terminal domain,HEAT,TOG domain,XMAP215/Dis1/CLASP

# Oxy chr1 Dry chr2
GO:0005975 -  HEAT repeat region,Armadillo-like helical,Armadillo-type fold,Laa1/Sip1/HEATR5-like,Protein SWEETIE

# Oxy ch1 Dry chr4
GO:0003684,GO:0005515,GO:0006289,GO:0043161- Heat shock chaperonin-binding,UBA-like superfamily,Ubiquitin-associated domain,Ubiquitin-like domain,Ubiquitin-like domain superfamily,UV excision repair protein Rad23,XPC-binding domain,XPC-binding domain superfamily

