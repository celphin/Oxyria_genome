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
# copy over Interproscan files
cp -v ../Interproscan/*_interproscan_output.tsv .

# join 
cat *_interproscan_output.tsv > Total_interproscan_output.tsv

# format 
awk -v FS="\t" '{print $1 "\t" $4 "\t" $12 "\t" $13 "\t" $14}' Total_interproscan_output.tsv | wc -l 
# 1669585

awk -v FS="\t" '{print $1 "\t" $4 "\t" $12 "\t" $13 "\t" $14}' Total_interproscan_output.tsv | sort | uniq
# 1112967

awk -v FS="\t" '{print $1 "\t" $12 "\t" $13 "\t" $14}' Total_interproscan_output.tsv | sort | uniq |wc -l
# 669657

awk -v FS="\t" '{print $1 "\t" $12 "\t" $13 "\t" $14}' Total_interproscan_output.tsv | sort | uniq > Total_interproscan_output_edited.tsv

grep -v $'\t''-'$'\t''-'$'\t' Total_interproscan_output_edited.tsv > Total_interproscan_output_edited1.tsv

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

library(parallel)
library(foreach)
library(iterators)
library(doParallel)

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
  orthogroups_in_all <- read.table(file.path(path, og_file), sep="\t", header = FALSE, check.names = FALSE)
  combined_genes <- read.table(file.path(path, comb_file), sep="\t", header = TRUE, check.names = FALSE)
  
  colnames(combined_genes) <- c("chr", "start", "end", "gene", "ofID", "pepLen", "ord", "genome", 
                                 "arrayID", "isArrayRep", "orthogroup", "globHOG", "noAnchor", "og")
  
  colnames(orthogroups_in_all)[1] <- "orthogroup"
  
  # Join files
  combined_genes_orthogroups <- dplyr::left_join(orthogroups_in_all, combined_genes, by="orthogroup")
  
  # Select only the species of interest
  orthogroups_in_spp <- combined_genes_orthogroups[combined_genes_orthogroups$genome == spp_hap, ]
  orthogroup_list_total <- unique(orthogroups_in_spp$orthogroup)
  
  og_groups_list <- 0:(round(length(orthogroup_list_total)/500)-1)
  
  X_seq_along <- NULL
  for (i in 0:(length(og_groups_list)-1)) {
    for (j in 0:(length(og_groups_list)-1)) {
      X_seq_along <- c(X_seq_along, paste0(i, "_", j))
    }
  }
  
  cl <- parallel::makeCluster(x_cpu)
  results <- parallel::parSapply(cl, seq_along(X_seq_along), 
                                  orthogroup_loop, X_seq_along=X_seq_along,
                                  spp_hap = spp_hap, path = path, 
                                  orthogroups_in_spp = orthogroups_in_spp,
                                  orthogroup_list_total = orthogroup_list_total)
  
  parallel::stopCluster(cl)
}

# test version

orthogroup_loop <- function(X, X_seq_along, spp_hap, path, orthogroups_in_spp, orthogroup_list_total) {
  
  XX <- X_seq_along[X]
  
  # Split each value in X_seq_along by "_"
  split_values <- strsplit(XX, "_")[[1]]
  i <- as.numeric(split_values[1])
  j <- as.numeric(split_values[2])
  
  start_orthog = (i * 500) + 1
  numb_orthogroups = (j * 500) + 500
  
  # Ensure indices do not exceed the length of the list
  if (start_orthog > length(orthogroup_list_total)) {
    start_orthog = length(orthogroup_list_total)
  }
  if (numb_orthogroups > length(orthogroup_list_total)) {
    numb_orthogroups = length(orthogroup_list_total)
  }
  
  # Get this set of orthogroups
  orthogroup_list <- orthogroup_list_total[start_orthog:numb_orthogroups]
  
  # Initialize matrices
  pos_diff <- matrix(NA, nrow = length(orthogroup_list), ncol = length(orthogroup_list))
  ord_diff <- matrix(NA, nrow = length(orthogroup_list), ncol = length(orthogroup_list))
  
  library(parallel)
  library(foreach)
  library(iterators)
  library(doParallel)

  output <- foreach::foreach(
    og1 = 1:length(orthogroup_list),
    og2 = 1:length(orthogroup_list),
    .combine = rbind,
    .packages = "dplyr"
  ) %do% {
    orthogroup1 <- orthogroup_list[og1]
    orthogroup2 <- orthogroup_list[og2]
    gene_list1 <- orthogroups_in_spp$gene[orthogroups_in_spp$orthogroup == orthogroup1]
    gene_list2 <- orthogroups_in_spp$gene[orthogroups_in_spp$orthogroup == orthogroup2]
    
    data.frame(genes1 = length(gene_list1), genes2 = length(gene_list2))
  }
  
  # Writing output
  tryCatch({
    write.table(output, file.path(path, "tmp", paste0(spp_hap, "_", start_orthog, "_", numb_orthogroups, "_gene_counts.txt")),
                sep = "\t", row.names = TRUE, col.names = TRUE)
  }, error = function(e) {
    message("Error writing output: ", e)
  })
}

orthogroup_loop <- function(X, X_seq_along, spp_hap, path, orthogroups_in_spp, orthogroup_list_total) {
  
  XX <- X_seq_along[X]
  
  # Split each value in X_seq_along by "_"
  split_values <- strsplit(XX, "_")[[1]]
  i <- as.numeric(split_values[1])
  j <- as.numeric(split_values[2])
  
  start_orthog = (i * 500) + 1
  numb_orthogroups = (j * 500) + 500
  
  # Ensure indices do not exceed the length of the list
  if (start_orthog > length(orthogroup_list_total)) {
    start_orthog = length(orthogroup_list_total)
  }
  if (numb_orthogroups > length(orthogroup_list_total)) {
    numb_orthogroups = length(orthogroup_list_total)
  }
  
  # Get this set of orthogroups
  orthogroup_list <- orthogroup_list_total[start_orthog:numb_orthogroups]
  
  # Initialize matrices
  pos_diff <- matrix(NA, nrow= length(orthogroup_list), ncol=length(orthogroup_list))
  ord_diff <- matrix(NA, nrow= length(orthogroup_list), ncol=length(orthogroup_list))
  pos_diff_max <- matrix(NA, nrow= length(orthogroup_list), ncol=length(orthogroup_list))
  gene_names1 <- matrix(NA, nrow= length(orthogroup_list), ncol=length(orthogroup_list))
  gene_names2 <- matrix(NA, nrow= length(orthogroup_list), ncol=length(orthogroup_list))

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

# https://cran.r-project.org/web/packages/foreach/vignettes/foreach.html
output <- NULL
old <- Sys.time()
output <- foreach::foreach(
      og1 = og.df$og1x,
      og2 = og.df$og2x,
     .combine = 'c',
     .packages="foreach") %do% {
    orthogroup1 <- orthogroup_list[og1]
    orthogroup2 <- orthogroup_list[og2]
    #print(og1)
    #print(og2)
    gene_list1 <- orthogroups_in_spp$gene[which(orthogroups_in_spp$orthogroup==orthogroup1)]
    gene_list2 <- orthogroups_in_spp$gene[which(orthogroups_in_spp$orthogroup==orthogroup2)]
    chr_list1 <- orthogroups_in_spp$chr[which(orthogroups_in_spp$orthogroup==orthogroup1)]
    chr_list2 <- orthogroups_in_spp$chr[which(orthogroups_in_spp$orthogroup==orthogroup2)]
    gene_diff <- matrix(NA, nrow= length(gene_list1), ncol=length(gene_list2))
    gene_ord <- matrix(NA, nrow= length(gene_list1), ncol=length(gene_list2))
    gene_name1 <- matrix(NA, nrow= length(gene_list1), ncol=length(gene_list2))
    gene_name2 <- matrix(NA, nrow= length(gene_list1), ncol=length(gene_list2))
    #print(length(gene_list1))
    #print(length(gene_list2))
    for (g1 in c(1:length(gene_list1))){
      for (g2 in c(1:length(gene_list2))){
        gene1 = gene_list1[g1]
        gene2 = gene_list2[g2]
        chr1 = chr_list1[g1]
        chr2 = chr_list2[g2]
        #print(gene1)
        #print(gene2)
        start1 <- orthogroups_in_spp$start[which(orthogroups_in_spp$gene==gene1)]
        start2 <- orthogroups_in_spp$start[which(orthogroups_in_spp$gene==gene2)]
        ord1 <- orthogroups_in_spp$ord[which(orthogroups_in_spp$gene==gene1)]
        ord2 <- orthogroups_in_spp$ord[which(orthogroups_in_spp$gene==gene2)]
        #print(start1)
        #print(start2)
        if (chr1==chr2){
          gene_diff[g1, g2] <- abs(start1-start2)
          gene_ord[g1, g2] <- abs(ord1-ord2)
          gene_name1[g1, g2] <- gene1
          gene_name2[g1, g2] <- gene2
        }else{
           gene_diff[g1,g2] <- 400e6
           gene_name1[g1, g2] <- gene1
           gene_name2[g1, g2] <- gene2
        }
    }
  }
  if (length(which(gene_diff == min(gene_diff, na.rm=TRUE)))!=0){
    inds <- data.frame()
    inds <- as.data.frame(which(gene_diff == min(gene_diff, na.rm=TRUE), arr.ind=TRUE))
    inds <- rbind(inds, newrow = c(0, 0))
    inds <- inds[1,]
    pos_diff[og1,og2] <- gene_diff[as.numeric(inds[1]),as.numeric(inds[2])]
    ord_diff[og1,og2] <- gene_ord[as.numeric(inds[1]),as.numeric(inds[2])]
    gene_names1[og1,og2] <- gene_name1[as.numeric(inds[1]),as.numeric(inds[2])]
    gene_names2[og1,og2] <- gene_name2[as.numeric(inds[1]),as.numeric(inds[2])]
  }
  if (length(which(gene_diff == max(gene_diff, na.rm=TRUE)))!=0){
    inds <- data.frame()
    inds <- as.data.frame(which(gene_diff == max(gene_diff, na.rm=TRUE), arr.ind=TRUE))
    inds <- rbind(inds, newrow = c(0, 0))
    inds <- inds[1,]
    pos_diff_max[og1,og2] <- gene_diff[as.numeric(inds[1]),as.numeric(inds[2])]
  }
  #print(pos_diff[og1,og2])
  #print(pos_diff_max[og1,og2]) 
  result <- multiResultClass()
  result$pos_diff <- pos_diff[og1,og2]
  result$ord_diff <- ord_diff[og1,og2]
  result$pos_diff_max <- pos_diff_max[og1,og2]
  result$gene_names1 <- gene_names1[og1,og2]
  result$gene_names2 <- gene_names2[og1,og2]
  # return 
  return(result)
}
new <- Sys.time()- old
print(new)

# for 10 in par Time difference of 35.2029 secs
# for 10 no par Time difference of 52.03195 secs
# for 10 cpu par for 500 orthoglogs Oxyria and Dryas Time difference of 9.361634 mins
# Rheum 10 cpu 500 = 7.461092 mins
 
df.output <- matrix(unlist(output), ncol = 5, byrow = TRUE)

pos_diff <- matrix(df.output[,1], ncol = (numb_orthogroups), byrow = TRUE)
ord_diff <- matrix(df.output[,2], nrow = (numb_orthogroups), byrow = TRUE)
pos_diff_max <- matrix(df.output[,3], nrow = (numb_orthogroups), byrow = TRUE)
gene_names1 <- matrix(df.output[,4], nrow = (numb_orthogroups), byrow = TRUE)
gene_names2 <- matrix(df.output[,5], nrow = (numb_orthogroups), byrow = TRUE)

rownames(pos_diff) <- orthogroup_list
colnames(pos_diff) <- orthogroup_list

rownames(ord_diff) <- orthogroup_list
colnames(ord_diff) <- orthogroup_list

rownames(pos_diff_max) <- orthogroup_list
colnames(pos_diff_max) <- orthogroup_list

rownames(gene_names1) <- orthogroup_list
colnames(gene_names1) <- orthogroup_list

rownames(gene_names2) <- orthogroup_list
colnames(gene_names2) <- orthogroup_list

utils::write.table(x=pos_diff, file=paste0(path,"/tmp/", spp_hap, "_", start_orthog, "_", numb_orthogroups,"_synteny_matrix.txt"), append = FALSE, sep = "\t", dec = ".", row.names = TRUE, col.names = TRUE)
utils::write.table(x=ord_diff, file=paste0(path,"/tmp/", spp_hap, "_", start_orthog, "_", numb_orthogroups,"_synteny_matrix_ord.txt"), append = FALSE, sep = "\t", dec = ".", row.names = TRUE, col.names = TRUE)
utils::write.table(x=pos_diff_max, file=paste0(path,"/tmp/", spp_hap, "_", start_orthog, "_", numb_orthogroups,"_synteny_matrix_max.txt"), append = FALSE, sep = "\t", dec = ".", row.names = TRUE, col.names = TRUE)
utils::write.table(x=gene_names1, file=paste0(path,"/tmp/", spp_hap, "_", start_orthog, "_", numb_orthogroups,"_synteny_matrix_genes1.txt"), append = FALSE, sep = "\t", dec = ".", row.names = TRUE, col.names = TRUE)
utils::write.table(x=gene_names2, file=paste0(path,"/tmp/", spp_hap, "_", start_orthog, "_", numb_orthogroups,"_synteny_matrix_genes2.txt"), append = FALSE, sep = "\t", dec = ".", row.names = TRUE, col.names = TRUE)

}

##################################
# load and join the 14 matrices

tmp_path <- paste0(path,"/tmp/")
file_list <- list.files(summary_path, full.names=TRUE)
ord_list <- file_list[grep("synteny_matrix_ord", file_list)]
gene1_list <- file_list[grep("synteny_matrix_genes1", file_list)]

#-----------------------------------
# join Shannon div data
lsd <- lapply(ord_list, read.table)
start_end_list0 <- basename(ord_list)
start_end_list1 <- stringr::str_split(start_end_list0, "_", simplify =TRUE)

start_list <- start_end_list1[,3]
end_list <- start_end_list1[,3]

names(lsd_start) <- start_list
names(lsd_end) <- end_list

Shannon_div_total <- dplyr::bind_rows(lsd, .id = 'chromosome')
Shannon_div_total <- dplyr::bind_cols(lsd, .id = 'chromosome')

colnames(Shannon_div_total) <- c("Chromosome", "Genome_position", "Shannon_div")

Shannon_div_total_parts <- Shannon_div_total[grep("-", Shannon_div_total$Chromosome),]




###########################################
# run with just 500 to see how it works

mkdir /home/celphin/scratch/Oxyria/synteny_quantity
cd /home/celphin/scratch/Oxyria/synteny_quantity

tmux new-session -s synteny
tmux attach-session -t synteny

# find OG in all
grep "Oxy" /home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/results/Orthogroups.tsv \
| grep "DoctH0" | grep "\-lg" | grep "AALP" | grep Polavi | grep FT | grep FEH | grep Rno | grep Rta >> Total_Orthogroups_in_all.tsv

wc -l Total_Orthogroups_in_all.tsv
# 6909 

# copy over BLAST combined data
cp /home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/results/combBed.txt ./Total_combBed.txt

wc -l Total_combBed.txt
# 651560 Total_combBed.txt

#---------------------------
# make parallel
# https://www.blasbenito.com/post/02_parallelizing_loops_with_r/
# in R

salloc -c40 --time 7:00:00 --mem 191000m --account def-rieseber

module load  StdEnv/2020 r/4.2.2
export R_LIBS_USER=/home/celphin/R/x86_64-pc-linux-gnu-library/4.1.0/

R

path="/home/celphin/scratch/Oxyria/synteny_quantity"
og_file="Total_Orthogroups_in_all.tsv"
comb_file="Total_combBed.txt"

orthogroup_distances(spp_hap="Oxyria_digyna_H1", path=path, og_file=og_file, comb_file=comb_file)
orthogroup_distances(spp_hap="Dryas_octopetala", path=path, og_file=og_file, comb_file=comb_file)
orthogroup_distances(spp_hap="Rheum_nobile_H0", path=path, og_file=og_file, comb_file=comb_file)
orthogroup_distances(spp_hap="Arabis_alpina", path=path, og_file=og_file, comb_file=comb_file)
orthogroup_distances(spp_hap="Polygunum_aviculare_H0", path=path, og_file=og_file, comb_file=comb_file)
orthogroup_distances(spp_hap="Draba_nivalis", path=path, og_file=og_file, comb_file=comb_file)


#----------------
orthogroup_distances <- function(spp_hap, path=path, og_file=og_file, comb_file=comb_file){

# load files above in R
orthogroups_in_all <- base::as.data.frame(utils::read.table(paste0(path,"/", og_file), sep="\t", header = FALSE, check.names = FALSE))
combined_genes <- base::as.data.frame(utils::read.table(paste0(path,"/", comb_file), sep="\t", header = TRUE, check.names = FALSE))

colnames(combined_genes) <- c( "chr", "start", "end", "gene", "ofID", "pepLen", "ord", "genome", 
"arrayID","isArrayRep", "orthogroup", "globHOG", "noAnchor", "og")

colnames(orthogroups_in_all)[1] <- "orthogroup"

#----------------------
# join files
combined_genes_orthogroups <- dplyr::left_join(orthogroups_in_all, combined_genes, by="orthogroup")
unique(combined_genes_orthogroups$genome)

orthogroups_in_spp <-combined_genes_orthogroups[which(combined_genes_orthogroups$genome==spp_hap), ]
orthogroup_list_total <- unique(orthogroups_in_spp$orthogroup)

nrow(orthogroups_in_spp)
#13916 Polygonum
# same Oxyria

length(orthogroup_list_total)
# 6909

#------------------------------
# inputs
#numb_orthogroups <- (length(orthogroup_list))
numb_orthogroups <- 500
start_orthog <- 1

orthogroup_list <- orthogroup_list_total[c(start_orthog:numb_orthogroups)]

# formatting
orthogroups_in_spp$start <- as.numeric(orthogroups_in_spp$start)
orthogroups_in_spp$ord <- as.numeric(orthogroups_in_spp$ord)

og.df <- expand.grid(og1x = c(1:numb_orthogroups),
                     og2x = c(1:numb_orthogroups))

pos_diff <- matrix(NA, nrow= length(orthogroup_list), ncol=length(orthogroup_list))
ord_diff <- matrix(NA, nrow= length(orthogroup_list), ncol=length(orthogroup_list))
pos_diff_max <- matrix(NA, nrow= length(orthogroup_list), ncol=length(orthogroup_list))
gene_names1 <- matrix(NA, nrow= length(orthogroup_list), ncol=length(orthogroup_list))
gene_names2 <- matrix(NA, nrow= length(orthogroup_list), ncol=length(orthogroup_list))

# install.packages("foreach")
# install.packages( "doParallel")
# install.packages("iterators")
library(parallel)
library(foreach)
library(iterators)
library(doParallel)

#create and register cluster
# n.cores=40
# my.cluster <- parallel::makeCluster(n.cores)
# doParallel::registerDoParallel(cl = my.cluster)

# to run
#for (og1 in c(1:(length(orthogroup_list)-1))) {
#  for (og2 in c((og1+1):length(orthogroup_list))){

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

# https://cran.r-project.org/web/packages/foreach/vignettes/foreach.html
output <- NULL
old <- Sys.time()
output <- foreach::foreach(
      og1 = og.df$og1x,
      og2 = og.df$og2x,
     .combine = 'c',
     .packages="foreach") %do% {
    orthogroup1 <- orthogroup_list[og1]
    orthogroup2 <- orthogroup_list[og2]
    #print(og1)
    #print(og2)
    gene_list1 <- orthogroups_in_spp$gene[which(orthogroups_in_spp$orthogroup==orthogroup1)]
    gene_list2 <- orthogroups_in_spp$gene[which(orthogroups_in_spp$orthogroup==orthogroup2)]
    chr_list1 <- orthogroups_in_spp$chr[which(orthogroups_in_spp$orthogroup==orthogroup1)]
    chr_list2 <- orthogroups_in_spp$chr[which(orthogroups_in_spp$orthogroup==orthogroup2)]
    gene_diff <- matrix(NA, nrow= length(gene_list1), ncol=length(gene_list2))
    gene_ord <- matrix(NA, nrow= length(gene_list1), ncol=length(gene_list2))
    gene_name1 <- matrix(NA, nrow= length(gene_list1), ncol=length(gene_list2))
    gene_name2 <- matrix(NA, nrow= length(gene_list1), ncol=length(gene_list2))
    #print(length(gene_list1))
    #print(length(gene_list2))
    for (g1 in c(1:length(gene_list1))){
      for (g2 in c(1:length(gene_list2))){
        gene1 = gene_list1[g1]
        gene2 = gene_list2[g2]
        chr1 = chr_list1[g1]
        chr2 = chr_list2[g2]
        #print(gene1)
        #print(gene2)
        start1 <- orthogroups_in_spp$start[which(orthogroups_in_spp$gene==gene1)]
        start2 <- orthogroups_in_spp$start[which(orthogroups_in_spp$gene==gene2)]
        ord1 <- orthogroups_in_spp$ord[which(orthogroups_in_spp$gene==gene1)]
        ord2 <- orthogroups_in_spp$ord[which(orthogroups_in_spp$gene==gene2)]
        #print(start1)
        #print(start2)
        if (chr1==chr2){
          gene_diff[g1, g2] <- abs(start1-start2)
          gene_ord[g1, g2] <- abs(ord1-ord2)
          gene_name1[g1, g2] <- gene1
          gene_name2[g1, g2] <- gene2
        }else{
           gene_diff[g1,g2] <- 400e6
           gene_name1[g1, g2] <- gene1
           gene_name2[g1, g2] <- gene2
        }
    }
  }
  if (length(which(gene_diff == min(gene_diff, na.rm=TRUE)))!=0){
    inds <- data.frame()
    inds <- as.data.frame(which(gene_diff == min(gene_diff, na.rm=TRUE), arr.ind=TRUE))
    inds <- rbind(inds, newrow = c(0, 0))
    inds <- inds[1,]
    pos_diff[og1,og2] <- gene_diff[as.numeric(inds[1]),as.numeric(inds[2])]
    ord_diff[og1,og2] <- gene_ord[as.numeric(inds[1]),as.numeric(inds[2])]
    gene_names1[og1,og2] <- gene_name1[as.numeric(inds[1]),as.numeric(inds[2])]
    gene_names2[og1,og2] <- gene_name2[as.numeric(inds[1]),as.numeric(inds[2])]
  }
  if (length(which(gene_diff == max(gene_diff, na.rm=TRUE)))!=0){
    inds <- data.frame()
    inds <- as.data.frame(which(gene_diff == max(gene_diff, na.rm=TRUE), arr.ind=TRUE))
    inds <- rbind(inds, newrow = c(0, 0))
    inds <- inds[1,]
    pos_diff_max[og1,og2] <- gene_diff[as.numeric(inds[1]),as.numeric(inds[2])]
  }
  #print(pos_diff[og1,og2])
  #print(pos_diff_max[og1,og2]) 
  result <- multiResultClass()
  result$pos_diff <- pos_diff[og1,og2]
  result$ord_diff <- ord_diff[og1,og2]
  result$pos_diff_max <- pos_diff_max[og1,og2]
  result$gene_names1 <- gene_names1[og1,og2]
  result$gene_names2 <- gene_names2[og1,og2]
  # return 
  return(result)
}
new <- Sys.time()- old
print(new)

# for 10 in par Time difference of 35.2029 secs
# for 10 no par Time difference of 52.03195 secs
# for 10 cpu par for 500 orthoglogs Oxyria and Dryas Time difference of 9.361634 mins
# Rheum 10 cpu 500 = 7.461092 mins
 
df.output <- matrix(unlist(output), ncol = 5, byrow = TRUE)

pos_diff <- matrix(df.output[,1], ncol = (numb_orthogroups), byrow = TRUE)
ord_diff <- matrix(df.output[,2], nrow = (numb_orthogroups), byrow = TRUE)
pos_diff_max <- matrix(df.output[,3], nrow = (numb_orthogroups), byrow = TRUE)
gene_names1 <- matrix(df.output[,4], nrow = (numb_orthogroups), byrow = TRUE)
gene_names2 <- matrix(df.output[,5], nrow = (numb_orthogroups), byrow = TRUE)

rownames(pos_diff) <- orthogroup_list
colnames(pos_diff) <- orthogroup_list

rownames(ord_diff) <- orthogroup_list
colnames(ord_diff) <- orthogroup_list

rownames(pos_diff_max) <- orthogroup_list
colnames(pos_diff_max) <- orthogroup_list

rownames(gene_names1) <- orthogroup_list
colnames(gene_names1) <- orthogroup_list

rownames(gene_names2) <- orthogroup_list
colnames(gene_names2) <- orthogroup_list

# check heatmaps
class(pos_diff) <- "numeric"
class(ord_diff) <- "numeric"
class(pos_diff_max) <- "numeric"

utils::write.table(x=pos_diff, file=paste0(path,"/", spp_hap, "_", start_orthog, "_", numb_orthogroups,"_synteny_matrix.txt"), append = FALSE, sep = "\t", dec = ".", row.names = TRUE, col.names = TRUE)
utils::write.table(x=ord_diff, file=paste0(path,"/", spp_hap, "_", start_orthog, "_", numb_orthogroups,"_synteny_matrix_ord.txt"), append = FALSE, sep = "\t", dec = ".", row.names = TRUE, col.names = TRUE)
utils::write.table(x=pos_diff_max, file=paste0(path,"/", spp_hap, "_", start_orthog, "_", numb_orthogroups,"_synteny_matrix_max.txt"), append = FALSE, sep = "\t", dec = ".", row.names = TRUE, col.names = TRUE)
utils::write.table(x=gene_names1, file=paste0(path,"/", spp_hap, "_", start_orthog, "_", numb_orthogroups,"_synteny_matrix_genes1.txt"), append = FALSE, sep = "\t", dec = ".", row.names = TRUE, col.names = TRUE)
utils::write.table(x=gene_names2, file=paste0(path,"/", spp_hap, "_", start_orthog, "_", numb_orthogroups,"_synteny_matrix_genes2.txt"), append = FALSE, sep = "\t", dec = ".", row.names = TRUE, col.names = TRUE)

}

#################################
# compare matrices for different species

tmux new-session -s synteny1
tmux attach-session -t synteny1

cd /home/celphin/scratch/Oxyria/synteny_quantity

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
gene1_names_spp1 <- base::as.matrix(utils::read.table(paste0(path,"/", spp_hap1, "_", numb_orthogroups,"_synteny_matrix_genes1.txt"), row.names = 1, header = TRUE, check.names = FALSE))
gene1_names_spp2 <- base::as.matrix(utils::read.table(paste0(path,"/", spp_hap2, "_", numb_orthogroups,"_synteny_matrix_genes1.txt"), row.names = 1, header = TRUE, check.names = FALSE))
gene1_names_spp3 <- base::as.matrix(utils::read.table(paste0(path,"/", spp_hap3, "_", numb_orthogroups,"_synteny_matrix_genes1.txt"), row.names = 1, header = TRUE, check.names = FALSE))
gene1_names_spp4 <- base::as.matrix(utils::read.table(paste0(path,"/", spp_hap4, "_", numb_orthogroups,"_synteny_matrix_genes1.txt"), row.names = 1, header = TRUE, check.names = FALSE))
gene1_names_spp5 <- base::as.matrix(utils::read.table(paste0(path,"/", spp_hap5, "_", numb_orthogroups,"_synteny_matrix_genes1.txt"), row.names = 1, header = TRUE, check.names = FALSE))
gene1_names_spp6 <- base::as.matrix(utils::read.table(paste0(path,"/", spp_hap6, "_", numb_orthogroups,"_synteny_matrix_genes1.txt"), row.names = 1, header = TRUE, check.names = FALSE))

ord_diff1 <- base::as.matrix(utils::read.table(paste0(path,"/", spp_hap1, "_", numb_orthogroups,"_synteny_matrix_ord.txt"), row.names = 1, header = TRUE, check.names = FALSE))
ord_diff2 <- base::as.matrix(utils::read.table(paste0(path,"/", spp_hap2, "_", numb_orthogroups,"_synteny_matrix_ord.txt"), row.names = 1, header = TRUE, check.names = FALSE))
ord_diff3 <- base::as.matrix(utils::read.table(paste0(path,"/", spp_hap3, "_", numb_orthogroups,"_synteny_matrix_ord.txt"), row.names = 1, header = TRUE, check.names = FALSE))
ord_diff4 <- base::as.matrix(utils::read.table(paste0(path,"/", spp_hap4, "_", numb_orthogroups,"_synteny_matrix_ord.txt"), row.names = 1, header = TRUE, check.names = FALSE))
ord_diff5 <- base::as.matrix(utils::read.table(paste0(path,"/", spp_hap5, "_", numb_orthogroups,"_synteny_matrix_ord.txt"), row.names = 1, header = TRUE, check.names = FALSE))
ord_diff6 <- base::as.matrix(utils::read.table(paste0(path,"/", spp_hap6, "_", numb_orthogroups,"_synteny_matrix_ord.txt"), row.names = 1, header = TRUE, check.names = FALSE))


#---------------------------
# find min diff and min sum
synteny_sum = ord_diff1 + ord_diff2 + ord_diff3 + ord_diff4 + ord_diff5

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

gene_count=70

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
binary_sum_ord <- ord_binary1_small + ord_binary2_small + ord_binary3_small + ord_binary4_small + ord_binary5_small + ord_binary6_small


###############################
# explore data

grDevices::png(paste0(path,"/Multi-spp_", numb_orthogroups, "_", gene_count, "_heatmap_binarysum_ord.png"), height=5000, width=5000)
  print(
  RepeatOBserverV1::imagenan(binary_sum_ord,lnumr=numb_orthogroups,lnumc=numb_orthogroups)
  )
  grDevices::dev.off()
  
grDevices::png(paste0(path,"/Multi-spp_", numb_orthogroups,"_", gene_count, "_heatmap_syntenysum.png"), height=5000, width=5000)
  print(
  RepeatOBserverV1::imagenan(synteny_sum,lnumr=numb_orthogroups,lnumc=numb_orthogroups)
  )
  grDevices::dev.off()
  
# which ones found in all
binary_sum_ord6 <- binary_sum_ord
binary_sum_ord6[which(binary_sum_ord6 <6, arr.ind=TRUE)]<-NA

grDevices::png(paste0(path,"/Multi-spp_", numb_orthogroups, "_", gene_count, "_heatmap_binarysum6_ord.png"), height=5000, width=5000)
  print(
  RepeatOBserverV1::imagenan(binary_sum_ord6,lnumr=numb_orthogroups,lnumc=numb_orthogroups)
  )
  grDevices::dev.off()
  
#----------------------------
inds_bin <- as.data.frame(which(binary_sum_ord >0, arr.ind=TRUE))
nrow(inds_bin)
# 33946/2 # on same chromosome in some
# total 50k/2

inds_bin <- as.data.frame(which(synteny_sum < 50, arr.ind=TRUE))
nrow(inds_bin)
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
#1372/2 
inds_bin <- as.data.frame(which(binary_sum_ord == 6, arr.ind=TRUE))
nrow(inds_bin)
#1206
# with 150 gene distance

#--------------------------
# how many ==2 in 70 gene distance
inds_bin <- as.data.frame(which(binary_sum_ord == 1, arr.ind=TRUE))
nrow(inds_bin)
# 20 594/2
inds_bin <- as.data.frame(which(binary_sum_ord == 2, arr.ind=TRUE))
nrow(inds_bin)
#12922/2
inds_bin <- as.data.frame(which(binary_sum_ord == 3, arr.ind=TRUE))
nrow(inds_bin)
#10120/2
inds_bin <- as.data.frame(which(binary_sum_ord == 4, arr.ind=TRUE))
nrow(inds_bin)
#5766/2
inds_bin <- as.data.frame(which(binary_sum_ord == 5, arr.ind=TRUE))
nrow(inds_bin)
#2902/2
inds_bin <- as.data.frame(which(binary_sum_ord == 6, arr.ind=TRUE))
nrow(inds_bin)

############################
# what do these genes found near each other in all do

spp1genes_inall <- gene1_names_spp1[which(binary_sum_ord == 5, arr.ind=TRUE)]
spp2genes_inall <- gene1_names_spp2[which(binary_sum_ord == 5, arr.ind=TRUE)]

genes_inall<- cbind(spp1genes_inall, spp2genes_inall)

genes_inall <- as.data.frame(genes_inall)
colnames(genes_inall) <- c("spp1genes", "spp2genes")

#-------------------------
# load GO ont data

Gene_ont_file <- "Total_interproscan_output_edited1.tsv"

gene_ont <- base::as.data.frame(utils::read.table(paste0(path,"/", Gene_ont_file), sep="\t",  header = FALSE, check.names = FALSE))

colnames(gene_ont) <- c("gene", "INTPRO", "descrip", "GOterm")
length(unique(gene_ont$INTPRO))
# 2099

#----------------------
# join with genes 
colnames(gene_ont) <- c("spp1genes", "INTPRO", "descrip", "GOterm")
genes_inall_go1 <- dplyr::distinct(dplyr::left_join(genes_inall, gene_ont, by="spp1genes"))

colnames(gene_ont) <- c("spp2genes", "INTPRO", "descrip", "GOterm")
genes_inall_go12 <- dplyr::distinct(dplyr::left_join(genes_inall_go1, gene_ont, by="spp2genes"))

unique(genes_inall_go12$INTPRO.x)
unique(genes_inall_go12$INTPRO.y)
#~790 Oxyria and Dryas for all

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
# 1128 with 70 gene dist and Draba, Dryas and Oxyria but not between Polygonaceae

# 1542 # in the two Arctic plants but not in the Polygonaceae - 70 gene dist
# 3374 - same but 150 gene dist
# 10186 # 500 gene dist

# what is the probability distribution when randomly arranging 30000 
# objects in a line that two objects fall within 70 objects of each other 5 times in a row


inds_bin <- as.data.frame(which(Arctic_NonArctic <0, arr.ind=TRUE))
nrow(inds_bin)
# 1750 - Arabis, Polygonum, Rheum but not Arctic spp

# 8414 #in the Polygonaceae plants and not Arctic - 70 gene dist
# 16084 # 150 gene dist
# 29526 # 500 gene dist

######################################
# Plot the most extreme values on the chromosomes of one species
# get the genes that are extreme

# matrix to list
ord_diff1_ls <- base::as.list(ord_diff1)
gene1_names_spp1_ls <- base::as.list(gene1_names_spp1)
Arctic_NonArctic_ls <- base::as.list(Arctic_NonArctic)

# unlist and cbind
spp1_data <- cbind(as.numeric(unlist(ord_diff1_ls)), unlist(gene1_names_spp1_ls), as.numeric(unlist(Arctic_NonArctic_ls)))

# match colnames
colnames(spp1_data) <- c("ord_diff", "gene", "arctic_score")
colnames(gene_pos) <- c("gene", "chr", "start", "ord", "og", "spp")

spp1_data <- as.data.frame(spp1_data)
gene_pos <- as.data.frame(gene_pos)

# remove NA rows
spp1_data<- spp1_data[which(spp1_data$arctic_score ==3, arr.ind=TRUE),]

# join with gene start pos for gene 1 spp1 and gene 1 spp2
spp1_data_starts <- dplyr::left_join(spp1_data, gene_pos, by="gene")

spp1_data_starts$start <- as.numeric(spp1_data_starts$start )
spp1_data_starts$ord_diff <- as.numeric(spp1_data_starts$ord_diff)
spp1_data_starts$chr <- as.factor(spp1_data_starts$chr)


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
Arctic_NonArctic_ls <- base::as.list(Arctic_NonArctic)

# unlist and cbind
spp2_data <- cbind(as.numeric(unlist(ord_diff2_ls)), unlist(gene1_names_spp2_ls), as.numeric(unlist(Arctic_NonArctic_ls)))

# match colnames
colnames(spp2_data) <- c("ord_diff", "gene", "arctic_score")
colnames(gene_pos) <- c("gene", "chr", "start", "ord", "og", "spp")

spp2_data <- as.data.frame(spp2_data)
gene_pos <- as.data.frame(gene_pos)

# remove NA rows
spp2_data<- spp2_data[which(spp2_data$arctic_score ==3, arr.ind=TRUE),]

# join with gene start pos for gene 1 spp1 and gene 1 spp2
spp2_data_starts <- dplyr::left_join(spp2_data, gene_pos, by="gene")

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


######################################
# Plot the most extreme values non-Arctic values on a non-Arctic plant
# get the genes that are extreme

# matrix to list
ord_diff4_ls <- base::as.list(ord_diff4)
gene1_names_spp4_ls <- base::as.list(gene1_names_spp4)
Arctic_NonArctic_ls <- base::as.list(Arctic_NonArctic)

# unlist and cbind
spp4_data <- cbind(as.numeric(unlist(ord_diff4_ls)), unlist(gene1_names_spp4_ls), as.numeric(unlist(Arctic_NonArctic_ls)))

# match colnames
colnames(spp4_data) <- c("ord_diff", "gene", "arctic_score")
colnames(gene_pos) <- c("gene", "chr", "start", "ord", "og", "spp")

spp4_data <- as.data.frame(spp4_data)
gene_pos <- as.data.frame(gene_pos)

# remove NA rows
spp4_data<- spp4_data[which(spp4_data$arctic_score ==-3, arr.ind=TRUE),]

# join with gene start pos
spp4_data_starts <- dplyr::left_join(spp4_data, gene_pos, by="gene")

spp4_data_starts$start <- as.numeric(spp4_data_starts$start )
spp4_data_starts$ord_diff <- as.numeric(spp4_data_starts$ord_diff)
spp4_data_starts$chr <- as.factor(spp4_data_starts$chr)


grDevices::png(paste0(path,"/NonArctic_spp4_", numb_orthogroups, "_", gene_count,"_synteny_on_genome.png"), height=2500, width=5000)
  print(
    ggplot2::ggplot(data=spp4_data_starts, ggplot2::aes(x=start, y=ord_diff))+
      ggplot2::geom_point(ggplot2::aes(x=start, y=ord_diff))+
      ggplot2::facet_wrap(~chr, scales = "free")+
      ggplot2::theme_classic()
  )
  grDevices::dev.off()



grDevices::png(paste0(path,"/NonArctic_spp4_", numb_orthogroups, "_", gene_count,"_synteny_on_genome_histogram.png"), height=2500, width=5000)
  print(
    ggplot2::ggplot(data=spp4_data_starts, ggplot2::aes(x=start))+
      ggplot2::geom_histogram(binwidth=1e5)+
      ggplot2::facet_wrap(~chr, scales = "free")+
      ggplot2::theme_classic()
  )
  grDevices::dev.off()

##########################################
# what do these genes found near each other in clusters do

# Gene ontology data - gene_ont file

colnames(gene_ont) <- c("gene", "INTPRO", "descrip", "GOterm")
length(unique(gene_ont$INTPRO))
# 2099

#----------------------
# join with genes of interest
colnames(gene_ont) <- c("gene", "INTPRO", "descrip", "GOterm")
genes_Arctic_go1 <- dplyr::distinct(dplyr::left_join(spp1_data, gene_ont, by="gene"))

unique(genes_Arctic_go1$INTPRO)
#unique(genes_Arctic_go1$descrip)
# 643 go terms

# write out genes_Arctic_go1
utils::write.table(x=genes_Arctic_go1, file=paste0(path,"/Arctic_gene_ontology_spp1.txt"), append = FALSE, sep = "\t", dec = ".", row.names = FALSE, col.names = TRUE)


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
peak_results <- lapply(split(spp1_data_starts$start, spp1_data_starts$chr), find_peaks_in_hist)

peak_results
names(peak_results)

spp1_chr6_50_60Mbp <- dplyr::distinct(spp1_data_starts[which( (spp1_data_starts$chr == "Oxyrt-6-73303751") &  (spp1_data_starts$start > 5.0e+07) &  (spp1_data_starts$start < 6.0e+07)),])
spp1_chr6_50_60Mbp_go <- dplyr::distinct(dplyr::left_join(spp1_chr6_50_60Mbp , gene_ont, by="gene"))
utils::write.table(x=spp1_chr6_50_60Mbp_go, file=paste0(path,"/spp1_chr6_50_60Mbp.txt"), append = FALSE, sep = "\t", dec = ".", row.names = FALSE, col.names = TRUE)

#----------------------------

# Apply the peak finding function to each group
peak_results <- lapply(split(spp2_data_starts$start, spp2_data_starts$chr), find_peaks_in_hist)

peak_results
names(peak_results)

spp2_chr1_17_18Mbp <- dplyr::distinct(spp2_data_starts[which( (spp2_data_starts$chr == "DoctH0-1") &  (spp2_data_starts$start >  15000000) &  (spp2_data_starts$start <  20000000)),])
spp2_chr1_17_18Mbp_go <- dplyr::distinct(dplyr::left_join(spp2_chr1_17_18Mbp , gene_ont, by="gene"))
utils::write.table(x=spp2_chr1_17_18Mbp_go, file=paste0(path,"/spp2_chr1_17_18Mbp.txt"), append = FALSE, sep = "\t", dec = ".", row.names = FALSE, col.names = TRUE)

spp2_chr6_6_8Mbp <- dplyr::distinct(spp2_data_starts[which( (spp2_data_starts$chr == "DoctH0-6") &  (spp2_data_starts$start >  6000000) &  (spp2_data_starts$start <  8000000)),])
spp2_chr6_6_8Mbp_go <- dplyr::distinct(dplyr::left_join(spp2_chr6_6_8Mbp , gene_ont, by="gene"))
utils::write.table(x=spp2_chr6_6_8Mbp_go, file=paste0(path,"/spp2_chr6_6_8Mbp.txt"), append = FALSE, sep = "\t", dec = ".", row.names = FALSE, col.names = TRUE)

spp2_chr8_4_6Mbp <- dplyr::distinct(spp2_data_starts[which( (spp2_data_starts$chr == "DoctH0-8") &  (spp2_data_starts$start >  4000000) &  (spp2_data_starts$start <  6000000)),])
spp2_chr8_4_6Mbp_go <- dplyr::distinct(dplyr::left_join(spp2_chr8_4_6Mbp , gene_ont, by="gene"))
utils::write.table(x=spp2_chr8_4_6Mbp_go, file=paste0(path,"/spp2_chr8_4_6Mbp.txt"), append = FALSE, sep = "\t", dec = ".", row.names = FALSE, col.names = TRUE)

spp2_chr3_9_11Mbp <- dplyr::distinct(spp2_data_starts[which( (spp2_data_starts$chr == "DoctH0-3") &  (spp2_data_starts$start >  9000000) &  (spp2_data_starts$start <  11000000)),])
spp2_chr3_9_11Mbp_go <- dplyr::distinct(dplyr::left_join(spp2_chr3_9_11Mbp , gene_ont, by="gene"))
utils::write.table(x=spp2_chr3_9_11Mbp_go, file=paste0(path,"/spp2_chr3_9_11Mbp.txt"), append = FALSE, sep = "\t", dec = ".", row.names = FALSE, col.names = TRUE)


##########################################################















#######################
#OLDER exploration
# only two species

#--------------------------
# get gene positions
combined_genes <- base::as.data.frame(utils::read.table(paste0(path,"/", comb_file), sep="\t", header = TRUE, check.names = FALSE))
colnames(combined_genes) <- c( "chr", "start", "end", "gene", "ofID", "pepLen", "ord", "genome", 
"arrayID","isArrayRep", "orthogroup", "globHOG", "noAnchor", "og")
gene_pos <- cbind(combined_genes$gene, combined_genes$chr, combined_genes$start, combined_genes$ord, combined_genes$orthogroup, combined_genes$genome)


#load files
gene1_names_spp1 <- base::as.matrix(utils::read.table(paste0(path,"/", spp_hap1, "_", numb_orthogroups,"_synteny_matrix_genes1.txt"), row.names = 1, header = TRUE, check.names = FALSE))
gene1_names_spp2 <- base::as.matrix(utils::read.table(paste0(path,"/", spp_hap2, "_", numb_orthogroups,"_synteny_matrix_genes1.txt"), row.names = 1, header = TRUE, check.names = FALSE))

ord_diff1 <- base::as.matrix(utils::read.table(paste0(path,"/", spp_hap1, "_", numb_orthogroups,"_synteny_matrix_ord.txt"), row.names = 1, header = TRUE, check.names = FALSE))
ord_diff2 <- base::as.matrix(utils::read.table(paste0(path,"/", spp_hap2, "_", numb_orthogroups,"_synteny_matrix_ord.txt"), row.names = 1, header = TRUE, check.names = FALSE))


#---------------------------
# find min diff and min sum
synteny_diff = abs(pos_diff1 - pos_diff2)
synteny_sum = pos_diff1 + pos_diff2

# determine which pos_diff are less than X Mbp - maybe 5Mbp
pos_binary1_small <- pos_diff1
pos_binary2_small <- pos_diff2

pos_binary1_small[which(pos_diff1 < 1e6, arr.ind=TRUE)] <- 1
pos_binary1_small[which(pos_diff1 > 1e6, arr.ind=TRUE)] <- 0

pos_binary2_small[which(pos_diff2 < 1e6, arr.ind=TRUE)] <- 1
pos_binary2_small[which(pos_diff2 > 1e6, arr.ind=TRUE)] <- 0

# add binary
binary_sum_small <- pos_binary1_small+pos_binary2_small

#---------------------------
# determine which pos_diff are on diff chromosomes
pos_binary1_large <- pos_diff1
pos_binary2_large <- pos_diff2

pos_binary1_large[which(pos_diff1 > 4e8, arr.ind=TRUE)] <- 1
pos_binary1_large[which(pos_diff1 < 4e8, arr.ind=TRUE)] <- 0

pos_binary2_large[which(pos_diff2 > 4e8, arr.ind=TRUE)] <- 1
pos_binary2_large[which(pos_diff2 < 4e8, arr.ind=TRUE)] <- 0

# add binary
binary_sum_large <- pos_binary1_large+pos_binary2_large

#---------------------------
# determine which pos_diff are less than X Mbp - maybe 5Mbp
ord_binary1_small <- ord_diff1
ord_binary2_small <- ord_diff2

ord_binary1_small[which(ord_diff1==0, arr.ind=TRUE)] <- NA
ord_binary2_small[which(ord_diff2==0, arr.ind=TRUE)] <- NA

ord_binary1_small[which(ord_diff1 < 50, arr.ind=TRUE)] <- 1
ord_binary1_small[which(ord_diff1 > 50, arr.ind=TRUE)] <- 0

ord_binary2_small[which(ord_diff2 < 50, arr.ind=TRUE)] <- 1
ord_binary2_small[which(ord_diff2 > 50, arr.ind=TRUE)] <- 0

# add binary
binary_sum_ord <- ord_binary1_small + ord_binary2_small


#----------------------
# use heatmaps in imagenan
#install.packages("devtools")
#library(devtools)
#install_github("celphin/RepeatOBserverV1") #to install the package
# Select 1:All to install all the required packages

grDevices::png(paste0(path,"/", spp_hap1, "_", numb_orthogroups,"_heatmap_posdiff1.png"), height=5000, width=5000)
  print(
  RepeatOBserverV1::imagenan(pos_diff1,lnumr=numb_orthogroups,lnumc=numb_orthogroups)
  )
  grDevices::dev.off()

  
grDevices::png(paste0(path,"/", spp_hap2, "_", numb_orthogroups,"_heatmap_posdiff2.png"), height=5000, width=5000)
  print(
 RepeatOBserverV1::imagenan(pos_diff2,lnumr=numb_orthogroups,lnumc=numb_orthogroups)
  )
  grDevices::dev.off()

grDevices::png(paste0(path,"/", spp_hap2, "_",spp_hap1, "_", numb_orthogroups,"_heatmap_syntenydiff.png"), height=5000, width=5000)
  print(
  RepeatOBserverV1::imagenan(synteny_diff,lnumr=numb_orthogroups,lnumc=numb_orthogroups)
  )
  grDevices::dev.off()
  
grDevices::png(paste0(path,"/", spp_hap2, "_",spp_hap1, "_", numb_orthogroups,"_heatmap_syntenysum.png"), height=5000, width=5000)
  print(
  RepeatOBserverV1::imagenan(synteny_sum,lnumr=numb_orthogroups,lnumc=numb_orthogroups)
  )
  grDevices::dev.off()

 
grDevices::png(paste0(path,"/", spp_hap2, "_",spp_hap1, "_", numb_orthogroups,"_heatmap_binarysumsmall.png"), height=5000, width=5000)
  print(
  RepeatOBserverV1::imagenan(binary_sum_small,lnumr=numb_orthogroups,lnumc=numb_orthogroups)
  )
  grDevices::dev.off()


grDevices::png(paste0(path,"/", spp_hap2, "_",spp_hap1, "_", numb_orthogroups,"_heatmap_binarysumlarge.png"), height=5000, width=5000)
  print(
  RepeatOBserverV1::imagenan(binary_sum_large,lnumr=numb_orthogroups,lnumc=numb_orthogroups)
  )
  grDevices::dev.off()
  
grDevices::png(paste0(path,"/", spp_hap2, "_",spp_hap1, "_", numb_orthogroups,"_heatmap_binarysum_ord.png"), height=5000, width=5000)
  print(
  RepeatOBserverV1::imagenan(binary_sum_ord,lnumr=numb_orthogroups,lnumc=numb_orthogroups)
  )
  grDevices::dev.off()
  
##########################
# subset to remove values > 400
pos_diff1_sub <- pos_diff1
pos_diff1_sub[which(pos_diff1 > 3e8, arr.ind=TRUE)] <- NA

pos_diff2_sub <- pos_diff2
pos_diff2_sub[which(pos_diff2 > 3e8, arr.ind=TRUE)] <- NA

synteny_diff_sub <- synteny_diff
synteny_diff_sub[which(synteny_diff > 3e8, arr.ind=TRUE)] <- NA

synteny_sum_sub <- synteny_sum
synteny_sum_sub[which(synteny_sum > 3e8, arr.ind=TRUE)] <- NA

# use heatmaps in imagenan
#install.packages("devtools")
#library(devtools)
#install_github("celphin/RepeatOBserverV1") #to install the package
# Select 1:All to install all the required packages

grDevices::png(paste0(path,"/", spp_hap1, "_", numb_orthogroups,"_heatmap_posdiff1_sub.png"), height=5000, width=5000)
  print(
  RepeatOBserverV1::imagenan(pos_diff1_sub,lnumr=numb_orthogroups,lnumc=numb_orthogroups)
  )
  grDevices::dev.off()

grDevices::png(paste0(path,"/", spp_hap2, "_", numb_orthogroups,"_heatmap_posdiff2_sub.png"), height=5000, width=5000)
  print(
 RepeatOBserverV1::imagenan(pos_diff2_sub,lnumr=numb_orthogroups,lnumc=numb_orthogroups)
  )
  grDevices::dev.off()

grDevices::png(paste0(path,"/", spp_hap2, "_",spp_hap1, "_", numb_orthogroups,"_heatmap_syntenydiff_sub.png"), height=5000, width=5000)
  print(
  RepeatOBserverV1::imagenan(synteny_diff_sub,lnumr=numb_orthogroups,lnumc=numb_orthogroups)
  )
  grDevices::dev.off()
  
grDevices::png(paste0(path,"/", spp_hap2, "_",spp_hap1, "_", numb_orthogroups,"_heatmap_syntenysum_sub.png"), height=5000, width=5000)
  print(
  RepeatOBserverV1::imagenan(synteny_sum_sub,lnumr=numb_orthogroups,lnumc=numb_orthogroups)
  )
  grDevices::dev.off()


#----------------------------
# how many ==2

inds_bin <- as.data.frame(which(binary_sum_small == 2, arr.ind=TRUE))
nrow(inds_bin)
# 7970 - Rheum/Oxyria
# 8052 - Dryas/Oxyria
# 14120 - Polygonum/Oxyria

inds_bin <- as.data.frame(which(binary_sum_ord == 2, arr.ind=TRUE))
nrow(inds_bin)
# 3506 - Dryas and Oxyria within 50
# 9724 - Polygonum and Oxyria


#---------------------------------
# find positions in matrix with zero diff - count
inds_diff <- data.frame()
inds_diff <- as.data.frame(which(synteny_diff == min(synteny_diff, na.rm=TRUE), arr.ind=TRUE))
nrow(inds_diff)
# 19860 Dryas, Oxy
# 19594, Rheum, Oxy

# remove zeros along diagonal
synteny_sum[which(synteny_sum == 0)] <- NA

inds_sum <- data.frame()
inds_sum <- as.data.frame(which(synteny_sum == min(synteny_sum, na.rm=TRUE), arr.ind=TRUE))
inds_sum <- rbind(inds_sum, newrow = c(0, 0))
synteny_sum_min <- synteny_sum[as.numeric(inds_sum[,1]),as.numeric(inds_sum[,2])]

# Dryas, Oxy
          # OG0000181 OG0000586
# OG0000586      6451        NA
# OG0000181        NA      6451

#Rheum/Oxy
          # OG0000034 OG0000208
# OG0000208      3610        NA
# OG0000034        NA      3610

#################################################
# which are the corresponding ogs and genes
synteny_Sppgene1_minsum <- gene1_names_spp1[as.numeric(inds_sum[,1]),as.numeric(inds_sum[,2])]
synteny_Sppgene2_minsum <- gene1_names_spp2[as.numeric(inds_sum[,1]),as.numeric(inds_sum[,2])]

          # OG0000181                  OG0000586
# OG0000586 "Oxyria_NCBI_Chr400000961" "Oxyria_NCBI_Chr200000937"
# OG0000181 "Oxyria_NCBI_Chr100002866" "Oxyria_NCBI_Chr400000962"

             # OG0000181             OG0000586
# OG0000586 "DoctH0_Chr900001699" "DoctH0_Chr200009208"
# OG0000181 "DoctH0_Chr100001737" "DoctH0_Chr900001696"

#################################
# calculate total sum of the differences
totalsum_score_diff <- sum(synteny_diff)
totalsum_score_diff
# 3.75155e+13 , Dryas/Oxy
#  2.528466e+13, Rheum/Oxy

##################################
# Plot sum values on genome

# make each matrix long 
# column for gene1, gene2, OG1, OG2, and diff

pos_diff1_ls <- base::as.list(pos_diff1)
pos_diff2_ls <- base::as.list(pos_diff2)
gene1_names_spp1_ls <- base::as.list(gene1_names_spp1)
gene1_names_spp2_ls <- base::as.list(gene1_names_spp2)
synteny_sum_ls <- base::as.list(synteny_sum)

spp1_data <- cbind(as.numeric(unlist(pos_diff1_ls)), unlist(gene1_names_spp1_ls), as.numeric(unlist(synteny_sum_ls)))
spp2_data <- cbind(as.numeric(unlist(pos_diff2_ls)), unlist(gene1_names_spp2_ls), as.numeric(unlist(synteny_sum_ls)))

# match colnames
colnames(spp1_data) <- c("pos_diff", "gene", "synteny_score")
colnames(spp2_data) <- c("pos_diff", "gene", "synteny_score")
colnames(gene_pos) <- c("gene", "chr", "start", "og", "spp")

spp1_data <- as.data.frame(spp1_data)
spp2_data <- as.data.frame(spp2_data)
gene_pos <- as.data.frame(gene_pos)

# join with gene start pos for gene 1 spp1 and gene 1 spp2
spp1_data_starts <- dplyr::left_join(spp1_data, gene_pos, by="gene")
spp2_data_starts <- dplyr::left_join(spp2_data, gene_pos, by="gene")


#------------------------------------------
# plot start pos and synteny sum for spp1 and spp2

library(ggplot2)

spp2_data_starts$start <- as.numeric(spp2_data_starts$start )
spp2_data_starts$synteny_score <- as.numeric(spp2_data_starts$synteny_score)
spp2_data_starts$chr <- as.factor(spp2_data_starts$chr)

spp1_data_starts$start <- as.numeric(spp1_data_starts$start )
spp1_data_starts$synteny_score <- as.numeric(spp1_data_starts$synteny_score)
spp1_data_starts$chr <- as.factor(spp1_data_starts$chr)

grDevices::png(paste0(path,"/", spp_hap1, "_",spp_hap2, "_", numb_orthogroups,"_synteny_on_genome.png"), height=2500, width=5000)
  print(
    ggplot2::ggplot(data=spp1_data_starts, ggplot2::aes(x=start, y=synteny_score))+
      ggplot2::geom_point(ggplot2::aes(x=start, y=synteny_score))+
      ggplot2::facet_wrap(~chr, scales = "free")+
      ggplot2::theme_classic()
  )
  grDevices::dev.off()

grDevices::pdf(paste0(path,"/", spp_hap2, "_",spp_hap1, "_", numb_orthogroups,"_synteny_on_genome.pdf"))
  print(
    ggplot2::ggplot(data=spp2_data_starts, ggplot2::aes(x=start, y=synteny_score))+
      ggplot2::geom_point(ggplot2::aes(x=start, y=synteny_score))+
      ggplot2::facet_wrap(~chr, scales = "free")+
      ggplot2::theme_classic()
  )
  grDevices::dev.off()
  
#####################################
# try min value for each start position
# https://dplyr.tidyverse.org/reference/summarise.html

# group by the start value/gene name
spp1_data_gb_genes <- dplyr::group_by(spp1_data_starts, gene, chr, start, spp, og)

# summarize the min value
min_syn_scores_spp1 <- dplyr::summarize(spp1_data_gb_genes, min=min(synteny_score, na.rm=TRUE))

# plot again
min_syn_scores_spp1$start <- as.numeric(min_syn_scores_spp1$start )
min_syn_scores_spp1$min<- as.numeric(min_syn_scores_spp1$min)
min_syn_scores_spp1$chr <- as.factor(min_syn_scores_spp1$chr)

# remove values greater than chromosome size  (>350Mbp)
length(which(min_syn_scores_spp1$min>3.5e8))
# 15 Rheum and Oxyria

min_syn_scores_spp1_rmlg <- min_syn_scores_spp1[-which(min_syn_scores_spp1$min>3.5e8),]

grDevices::png(paste0(path,"/", spp_hap1, "_", spp_hap2, "_",numb_orthogroups,"_min_synteny_on_genome.png"), height=2500, width=5000)
  print(
    ggplot2::ggplot(data=min_syn_scores_spp1_rmlg, ggplot2::aes(x=start, y=min))+
      ggplot2::geom_point(ggplot2::aes(x=start, y=min))+
      ggplot2::facet_wrap(~chr, scales = "free")+
      ggplot2::theme_classic()
  )
  grDevices::dev.off()

#------------------------------
# group by the start value/gene name
spp2_data_gb_genes <- dplyr::group_by(spp2_data_starts, gene, chr, start, spp, og)

# summarize the min value
min_syn_scores_spp2 <- dplyr::summarize(spp2_data_gb_genes, min=min(synteny_score, na.rm=TRUE))

# plot again
min_syn_scores_spp2$start <- as.numeric(min_syn_scores_spp2$start )
min_syn_scores_spp2$min <- as.numeric(min_syn_scores_spp2$min)
min_syn_scores_spp2$chr <- as.factor(min_syn_scores_spp2$chr)

# remove values greater than chromosome size  (>350Mbp)
length(which(min_syn_scores_spp2$min>3.5e8))
# 7 Rheum and Oxyria
# 15 Dryas and Oxyria
min_syn_scores_spp2_rmlg <- min_syn_scores_spp2[-which(min_syn_scores_spp2$min>3.5e8),]

grDevices::png(paste0(path,"/", spp_hap2, "_",spp_hap1, "_", numb_orthogroups,"_min_synteny_on_genome.png"), height=2500, width=5000)
  print(
    ggplot2::ggplot(data=min_syn_scores_spp2_rmlg, ggplot2::aes(x=start, y=min))+
      ggplot2::geom_point(ggplot2::aes(x=start, y=min))+
      ggplot2::facet_wrap(~chr, scales = "free")+
      ggplot2::theme_classic()
  )
  grDevices::dev.off()

#############################
# running rolling mean over values

bin_size=10
min_syn_scores_spp1_rmlg$roll_min <- zoo::rollapply(min_syn_scores_spp1_rmlg$min, width = bin_size, FUN=mean, fill = NA, partial=(bin_size/2))
min_syn_scores_spp2_rmlg$roll_min <- zoo::rollapply(min_syn_scores_spp2_rmlg$min, width = bin_size, FUN=mean, fill = NA, partial=(bin_size/2))

grDevices::png(paste0(path,"/", spp_hap1, "_",spp_hap2, "_", numb_orthogroups,"_rollmin_synteny_on_genome.png"), height=2500, width=5000)
  print(
    ggplot2::ggplot(data=min_syn_scores_spp1_rmlg, ggplot2::aes(x=start, y=roll_min))+
      ggplot2::geom_line(ggplot2::aes(x=start, y=roll_min))+
      ggplot2::facet_wrap(~chr, scales = "free")+
      ggplot2::theme_classic()
  )
  grDevices::dev.off()


grDevices::png(paste0(path,"/", spp_hap2, "_", spp_hap1, "_",numb_orthogroups,"_rollmin_synteny_on_genome.png"), height=2500, width=5000)
  print(
    ggplot2::ggplot(data=min_syn_scores_spp2_rmlg, ggplot2::aes(x=start, y=roll_min))+
      ggplot2::geom_line(ggplot2::aes(x=start, y=roll_min))+
      ggplot2::facet_wrap(~chr, scales = "free")+
      ggplot2::theme_classic()
  )
  grDevices::dev.off()


#################################
# run for difference not sum to get a score like value 

# synteny_diff

pos_diff1_ls <- base::as.list(pos_diff1)
pos_diff2_ls <- base::as.list(pos_diff2)
gene1_names_spp1_ls <- base::as.list(gene1_names_spp1)
gene1_names_spp2_ls <- base::as.list(gene1_names_spp2)
synteny_diff_ls <- base::as.list(synteny_diff)

spp1_data <- cbind(as.numeric(unlist(pos_diff1_ls)), unlist(gene1_names_spp1_ls), as.numeric(unlist(synteny_diff_ls)))
spp2_data <- cbind(as.numeric(unlist(pos_diff2_ls)), unlist(gene1_names_spp2_ls), as.numeric(unlist(synteny_diff_ls)))

# match colnames
colnames(spp1_data) <- c("pos_diff", "gene", "synteny_score")
colnames(spp2_data) <- c("pos_diff", "gene", "synteny_score")
colnames(gene_pos) <- c("gene", "chr", "start", "og", "spp")

spp1_data <- as.data.frame(spp1_data)
spp2_data <- as.data.frame(spp2_data)
gene_pos <- as.data.frame(gene_pos)

# join with gene start pos for gene 1 spp1 and gene 1 spp2
spp1_data_starts <- dplyr::left_join(spp1_data, gene_pos, by="gene")
spp2_data_starts <- dplyr::left_join(spp2_data, gene_pos, by="gene")


#------------------------------------------
# plot start pos and synteny sum for spp1 and spp2

library(ggplot2)

spp2_data_starts$start <- as.numeric(spp2_data_starts$start )
spp2_data_starts$synteny_score <- as.numeric(spp2_data_starts$synteny_score)
spp2_data_starts$chr <- as.factor(spp2_data_starts$chr)

spp1_data_starts$start <- as.numeric(spp1_data_starts$start )
spp1_data_starts$synteny_score <- as.numeric(spp1_data_starts$synteny_score)
spp1_data_starts$chr <- as.factor(spp1_data_starts$chr)

grDevices::png(paste0(path,"/", spp_hap1, "_",spp_hap2, "_", numb_orthogroups,"_synteny_on_genome.png"), height=2500, width=5000)
  print(
    ggplot2::ggplot(data=spp1_data_starts, ggplot2::aes(x=start, y=synteny_score))+
      ggplot2::geom_point(ggplot2::aes(x=start, y=synteny_score))+
      ggplot2::facet_wrap(~chr, scales = "free")+
      ggplot2::theme_classic()
  )
  grDevices::dev.off()

grDevices::png(paste0(path,"/", spp_hap2, "_",spp_hap1, "_", numb_orthogroups,"_synteny_on_genome.png"), height=2500, width=5000)
  print(
    ggplot2::ggplot(data=spp2_data_starts, ggplot2::aes(x=start, y=synteny_score))+
      ggplot2::geom_point(ggplot2::aes(x=start, y=synteny_score))+
      ggplot2::facet_wrap(~chr, scales = "free")+
      ggplot2::theme_classic()
  )
  grDevices::dev.off()

#####################################
# try mean value for each start position
# https://dplyr.tidyverse.org/reference/summarise.html

# group by the start value/gene name
spp1_data_gb_genes <- dplyr::group_by(spp1_data_starts, gene, chr, start, spp, og)

# summarize the mean value
mean_syn_scores_spp1 <- dplyr::summarize(spp1_data_gb_genes, mean=mean(synteny_score, na.rm=TRUE))

# plot again
mean_syn_scores_spp1$start <- as.numeric(mean_syn_scores_spp1$start )
mean_syn_scores_spp1$mean<- as.numeric(mean_syn_scores_spp1$mean)
mean_syn_scores_spp1$chr <- as.factor(mean_syn_scores_spp1$chr)


grDevices::png(paste0(path,"/", spp_hap1, "_", spp_hap2, "_",numb_orthogroups,"_mean_syntenydiff_on_genome.png"), height=2500, width=5000)
  print(
    ggplot2::ggplot(data=mean_syn_scores_spp1, ggplot2::aes(x=start, y=mean))+
      ggplot2::geom_point(ggplot2::aes(x=start, y=mean))+
      ggplot2::facet_wrap(~chr, scales = "free")+
      ggplot2::theme_classic()
  )
  grDevices::dev.off()

#-----------------------
# group by the start value/gene name
spp2_data_gb_genes <- dplyr::group_by(spp2_data_starts, gene, chr, start, spp, og)

# summarize the mean value
mean_syn_scores_spp2 <- dplyr::summarize(spp2_data_gb_genes, mean=mean(synteny_score, na.rm=TRUE))

# plot again
mean_syn_scores_spp2$start <- as.numeric(mean_syn_scores_spp2$start )
mean_syn_scores_spp2$mean <- as.numeric(mean_syn_scores_spp2$mean)
mean_syn_scores_spp2$chr <- as.factor(mean_syn_scores_spp2$chr)



grDevices::png(paste0(path,"/", spp_hap2, "_",spp_hap1, "_", numb_orthogroups,"_mean_syntenydiff_on_genome.png"), height=2500, width=5000)
  print(
    ggplot2::ggplot(data=mean_syn_scores_spp2, ggplot2::aes(x=start, y=mean))+
      ggplot2::geom_point(ggplot2::aes(x=start, y=mean))+
      ggplot2::facet_wrap(~chr, scales = "free")+
      ggplot2::theme_classic()
  )
  grDevices::dev.off()

##############################
# remove values with difference more than 5 Mbp on average
length(which(mean_syn_scores_spp1$mean>5e6))
# 2699 Dryas and Oxyria
#  2863 Rheum and Oxyria

mean_syn_scores_spp1_rmlg <- mean_syn_scores_spp1[which(mean_syn_scores_spp1$mean<5e6),]

grDevices::png(paste0(path,"/", spp_hap1, "_", spp_hap2, "_",numb_orthogroups,"_mean_syntenydiff_on_genome_rmlg.png"), height=2500, width=5000)
  print(
    ggplot2::ggplot(data=mean_syn_scores_spp1_rmlg, ggplot2::aes(x=start, y=mean))+
      ggplot2::geom_point(ggplot2::aes(x=start, y=mean))+
      ggplot2::facet_wrap(~chr, scales = "free")+
      ggplot2::theme_classic()
  )
  grDevices::dev.off()

#------------------------------
# remove values with a difference more than 5 Mbp
length(which(mean_syn_scores_spp2$mean>5e6))

# 2184 Dryas and Oxyria
mean_syn_scores_spp2_rmlg <- mean_syn_scores_spp2[which(mean_syn_scores_spp2$mean<5e6),]

grDevices::png(paste0(path,"/", spp_hap2, "_",spp_hap1, "_", numb_orthogroups,"_mean_syntenydiff_on_genome_rmlg.png"), height=2500, width=5000)
  print(
    ggplot2::ggplot(data=mean_syn_scores_spp2_rmlg, ggplot2::aes(x=start, y=mean))+
      ggplot2::geom_point(ggplot2::aes(x=start, y=mean))+
      ggplot2::facet_wrap(~chr, scales = "free")+
      ggplot2::theme_classic()
  )
  grDevices::dev.off()
  



#################################
# try running with GO terms not orthogroups
# seems like some genes are in many orthogroups

# join gff3 and interproscan output 

# same as orthologs but with GO terms