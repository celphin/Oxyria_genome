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


################################################
# # Results_Aug18
# cd /home/celphin/scratch/Oxyria/GeneSpace/Polygonaceae_genomes/results
# more combBed.txt
# chr     start              end              id           ofID   pepLen  ord     genome                     arrayID      isArrayRep  globOG         globHOG               noAnchor   og
# Chr1    37076943        37077498        FEHAP201918.t1  0_2049  184     2050    Fagopyrum_escelentum_H2    Arr000001       TRUE    OG0000006       N0.HOG0000006 OG0000006 FALSE   29631
# Chr1    99260663        99261385        FEHAP203619.t1  0_3805  223     3806    Fagopyrum_escelentum_H2    Arr000002       TRUE    OG0000006       N0.HOG0000006 OG0000006 FALSE   34058
# Chr1    145550788       145551679       FEHAP205055.t1  0_5313  259     5314    Fagopyrum_escelentum_H2    Arr000003       TRUE    OG0000006       N0.HOG0000006 OG0000006 FALSE   37741

# #------------------------
# grep Oxy Orthogroups.tsv | grep Polavi | grep FT | grep FEH | grep Rno | grep Rta |wc -l
# 11683

# # take  OG groups in all from above - extract start positions for each species
# awk -F '\t' '$11 == "OG0015429"' combBed.txt

# Chr8                    139721165       139726687       FEHAP248785.t1                  0_51638 536     51639   Fagopyrum_escelentum_H2    Arr026996       TRUE    OG0015429       N0.HOG0018908 OG0015429 FALSE   18275
# Chr2                    13755697        13762221        FT01Gene05295.t1                1_5654  536     5655    Fagopyrum_tataricum_H1     Arr034548       TRUE    OG0015429       N0.HOG0018908 OG0015429 FALSE   18275
# Oxyrt-7-72361354        36549423        36554224        Oxyria_NCBI_Chr700002107        2_30883 538     31229   Oxyria_digyna_H1           Arr077401       TRUE    OG0015429       N0.HOG0018908 OG0015429 FALSE   18275
# Polavi-6                22400038        22404469        Polavi_Chr600002219             3_18540 538     15706   Polygunum_aviculare_H0     Arr093293       TRUE    OG0015429       N0.HOG0018908 OG0015429 FALSE   18275
# RnoChr04                80507799        80513772        RnoG0031967.1                   4_11646 537     11647   Rheum_nobile_H0            Arr110490       TRUE    OG0015429       N0.HOG0018908 OG0015429 FALSE   18275
# RtaChr04                124186686       124192887       RtaG0015388.1                   5_11649 527     11650   Rheum_tangaticum_H0        Arr137927       TRUE    OG0015429       N0.HOG0018908 OG0015429 FALSE   18275

###########################################
mkdir /home/celphin/scratch/Oxyria/synteny_quantity
cd /home/celphin/scratch/Oxyria/synteny_quantity

tmux new-session -s synteny
tmux attach-session -t synteny

# find OG in all
grep Oxy /home/celphin/scratch/Oxyria/GeneSpace/Polygonaceae_genomes/results/Orthogroups.tsv \
| grep Polavi | grep FT | grep FEH | grep Rno | grep Rta >> Polygonaceae_Orthogroups_in_all.tsv

# copy over BLAST combined data
cp /home/celphin/scratch/Oxyria/GeneSpace/Polygonaceae_genomes/results/combBed.txt .

# in R
module load  StdEnv/2020 r/4.2.2
export R_LIBS_USER=/home/celphin/R/x86_64-pc-linux-gnu-library/4.1.0/

R

path="/home/celphin/scratch/Oxyria/synteny_quantity"
og_file="Polygonaceae_Orthogroups_in_all.tsv"
comb_file="Polygonaceae_combBed.txt"
spp_hap="Oxyria_digyna_H1"
#----------------
#orthogroup_distances <- function(path=path, og_file=og_file, comb_file=comb_file){

# load files above in R
orthogroups_in_all <- base::as.data.frame(utils::read.table(paste0(path,"/", og_file), sep="\t", header = FALSE, check.names = FALSE))
combined_genes <- base::as.data.frame(utils::read.table(paste0(path,"/", comb_file), sep="\t", header = TRUE, check.names = FALSE))

colnames(combined_genes) <- c( "chr", "start", "end", "gene", "ofID", "pepLen", "ord", "genome", 
"arrayID","isArrayRep", "orthogroup", "globHOG", "noAnchor", "og")

colnames(orthogroups_in_all) <- c("orthogroup", "Spp1", "Spp2", "Spp3", "Spp4", "Spp5", "Spp6")

#----------------------
# join files

combined_genes_orthogroups <- dplyr::left_join(orthogroups_in_all, combined_genes, by="orthogroup")
orthogroups_in_spp <-combined_genes_orthogroups[which(combined_genes_orthogroups$genome==spp_hap), ]

nrow(orthogroups_in_spp)
#[1] 16236

length(orthogroup_list)
#[1] 11683

#------------------------------
orthogroup_list <- unique(orthogroups_in_spp$orthogroup)

pos_diff <- data.frame(nrow( orthogroups_in_spp), nrow( orthogroups_in_spp), rep(NA)
gene_names1 <- data.frame(nrow( orthogroups_in_spp), nrow( orthogroups_in_spp), NA)
gene_names2 <- data.frame(nrow( orthogroups_in_spp), nrow( orthogroups_in_spp), NA)

for (og1 in c(1:(length(orthogroup_list)-1))) {
  for (og2 in c((og1+1):length(orthogroup_list))){
    orthogroup1 <- orthogroup_list[og1]
    orthogroup2 <- orthogroup_list[og2]
    print(orthogroup1)
    print(orthogroup2)
    gene_list1 <- orthogroups_in_spp$gene[which(orthogroups_in_spp$orthogroup==orthogroup1)]
    gene_list2 <- orthogroups_in_spp$gene[which(orthogroups_in_spp$orthogroup==orthogroup2)]
    chr_list1 <- orthogroups_in_spp$chr[which(orthogroups_in_spp$orthogroup==orthogroup1)]
    chr_list2 <- orthogroups_in_spp$chr[which(orthogroups_in_spp$orthogroup==orthogroup2)]
    gene_diff <- data.frame()
    gene_name1 <- data.frame()
    gene_name2 <- data.frame()
    print(length(gene_list1))
    print(length(gene_list2))
    for (g1 in c(1:length(gene_list1))){
      for (g2 in c(1:length(gene_list2))){
        gene1 <- gene_list1[g1]
        gene2 <- gene_list2[g2]
        chr1 <- chr_list1[g1]
        chr2 <- chr_list2[g2]
        #print(gene1)
        #print(gene2)
        start1 <- orthogroups_in_spp$start[which(orthogroups_in_spp$gene==gene1)]
        start2 <- orthogroups_in_spp$start[which(orthogroups_in_spp$gene==gene2)]
        #print(start1)
        #print(start2)
        if (chr1==chr2){
          gene_diff[gene1, gene2] <- abs(start1-start2)
          gene_name1[gene1, gene2] <- gene1
          gene_name2[gene1, gene2] <- gene2
        } else{
          gene_diff[gene1, gene2] <- NA
          gene_name1[gene1, gene2] <- NA
          gene_name2[gene1, gene2] <- NA
          }
    }
  }
  if (length(which(gene_diff == min(gene_diff, na.rm=TRUE)))!=0){
    inds <- data.frame()
    inds <- as.data.frame(which(gene_diff == min(gene_diff, na.rm=TRUE), arr.ind=TRUE))
    inds <- rbind(inds, newrow = c(0, 0))
    inds <- inds[1,]
    pos_diff[og1,og2] <- gene_diff[as.numeric(inds[1]),as.numeric(inds[2])]
    gene_names1[og1,og2] <- gene_name1[as.numeric(inds[1]),as.numeric(inds[2])]
    gene_names2[og1,og2] <- gene_name2[as.numeric(inds[1]),as.numeric(inds[2])]
  }else{
    pos_diff[og1,og2] <- 400e6
    gene_names1[og1,og2] <- "X"
    gene_names2[og1,og2] <- "X"
  }
  print(pos_diff[og1,og2])
 }
}

utils::write.table(x=pos_diff, file=paste0(path,"/Oxyria_synteny_matrix.txt"), append = FALSE, sep = "\t", dec = ".", row.names = TRUE, col.names = TRUE)

####################################
# once working  - make parallel
# https://www.blasbenito.com/post/02_parallelizing_loops_with_r/

mkdir /home/celphin/scratch/Oxyria/synteny_quantity
cd /home/celphin/scratch/Oxyria/synteny_quantity

tmux new-session -s synteny1
tmux attach-session -t synteny1

more /home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/results/Orthogroups.tsv
# gene starts: LOC, FEHAP, FT, Polavi
# Fagopyrum_tataricum_H1, Fagopyrum_escelentum_H2, Polygunum_aviculare_H0


# find OG in all
grep "Oxyria"  /home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/results/Orthogroups.tsv \
| grep "Rno" | grep "DoctH0" | grep "\-lg" | grep "AALP" | grep "FT" | grep "FEHAP" | grep "Polavi" >> Total_Orthogroups_in_all.tsv

wc -l Total_Orthogroups_in_all.tsv
# 7553 Total_Orthogroups_in_all.tsv
# 14662 more species

# copy over BLAST combined data
grep -E 'Oxyria_digyna|Rheum_nobile|Dryas_octopetala|Draba_nivalis|Arabis_alpina|Fagopyrum_tataricum_H1|Fagopyrum_escelentum_H2|Polygunum_aviculare_H0' \
/home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/results/combBed.txt \
 >> Total_combBed.txt

wc -l Total_combBed.txt
# 163359 Total_combBed.txt
# 444628 Total_combBed.txt

#---------------------------
# in R

salloc -c10 --time 0:30:00 --mem 120000m --account def-rieseber

module load  StdEnv/2020 r/4.2.2
export R_LIBS_USER=/home/celphin/R/x86_64-pc-linux-gnu-library/4.1.0/

R

path="/home/celphin/scratch/Oxyria/synteny_quantity"
og_file="Total_Orthogroups_in_all.tsv"
comb_file="Total_combBed.txt"
#spp_hap="Oxyria_digyna_H1"
#spp_hap="Dryas_octopetala"
#spp_hap="Rheum_nobile_H0"
#spp_hap="Arabis_alpina"
#spp_hap="Draba_nivalis" # did not run
spp_hap="Fagopyrum_tataricum_H1"
#spp_hap="Fagopyrum_escelentum_H2"
#spp_hap="Polygunum_aviculare_H0"

#----------------
#orthogroup_distances <- function(path=path, og_file=og_file, comb_file=comb_file){

# load files above in R
orthogroups_in_all <- base::as.data.frame(utils::read.table(paste0(path,"/", og_file), sep="\t", header = FALSE, check.names = FALSE))
combined_genes <- base::as.data.frame(utils::read.table(paste0(path,"/", comb_file), sep="\t", header = FALSE, check.names = FALSE))

colnames(combined_genes) <- c( "chr", "start", "end", "gene", "ofID", "pepLen", "ord", "genome", 
"arrayID","isArrayRep", "orthogroup", "globHOG", "noAnchor", "og")

colnames(orthogroups_in_all)[1] <- "orthogroup"

#----------------------
# join files
combined_genes_orthogroups <- dplyr::left_join(orthogroups_in_all, combined_genes, by="orthogroup")
unique(combined_genes_orthogroups$genome)

orthogroups_in_spp <-combined_genes_orthogroups[which(combined_genes_orthogroups$genome==spp_hap), ]
orthogroup_list <- unique(orthogroups_in_spp$orthogroup)

nrow(orthogroups_in_spp)
#[1] 28385

length(orthogroup_list)
#[1] 7553

#------------------------------
# inputs
pos_diff <- matrix(NA, nrow= length(orthogroup_list), ncol=length(orthogroup_list))
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
n.cores=10
my.cluster <- parallel::makeCluster(n.cores)
doParallel::registerDoParallel(cl = my.cluster)

# to run
#for (og1 in c(1:(length(orthogroup_list)-1))) {
#  for (og2 in c((og1+1):length(orthogroup_list))){

multiResultClass <- function(pos_diff=NULL,pos_diff_max=NULL, gene_names1=NULL, gene_names2=NULL)
{
  me <- list(
    pos_diff = pos_diff,
    pos_diff_max = pos_diff_max,
    gene_names1 = gene_names1,
    gene_names2 = gene_names2
  )

  ## Set the name for the class
  class(me) <- append(class(me),"multiResultClass")
  return(me)
}

#numb_orthogroups <- (length(orthogroup_list))
numb_orthogroups <- 500

og.df <- expand.grid(og1x = c(1:numb_orthogroups),
                     og2x = c(1:numb_orthogroups))

output <- NULL
old <- Sys.time()
output <- foreach::foreach(
      og1 = og.df$og1x,
      og2 = og.df$og2x,
     .combine = 'c',
     .packages="foreach") %dopar% {
    orthogroup1 <- orthogroup_list[og1]
    orthogroup2 <- orthogroup_list[og2]
    #print(og1)
    #print(og2)
    gene_list1 <- orthogroups_in_spp$gene[which(orthogroups_in_spp$orthogroup==orthogroup1)]
    gene_list2 <- orthogroups_in_spp$gene[which(orthogroups_in_spp$orthogroup==orthogroup2)]
    chr_list1 <- orthogroups_in_spp$chr[which(orthogroups_in_spp$orthogroup==orthogroup1)]
    chr_list2 <- orthogroups_in_spp$chr[which(orthogroups_in_spp$orthogroup==orthogroup2)]
    gene_diff <- matrix(NA, nrow= length(gene_list1), ncol=length(gene_list2))
    gene_name1 <- matrix(NA, nrow= length(gene_list1), ncol=length(gene_list2))
    gene_name2 <- matrix(NA, nrow= length(gene_list1), ncol=length(gene_list2))
    #print(length(gene_list1))
    #print(length(gene_list2))
    for (g1 in c(1:length(gene_list1))){
      for (g2 in c(1:length(gene_list2))){
        gene1 <- gene_list1[g1]
        gene2 <- gene_list2[g2]
        chr1 <- chr_list1[g1]
        chr2 <- chr_list2[g2]
        #print(gene1)
        #print(gene2)
        start1 <- orthogroups_in_spp$start[which(orthogroups_in_spp$gene==gene1)]
        start2 <- orthogroups_in_spp$start[which(orthogroups_in_spp$gene==gene2)]
        #print(start1)
        #print(start2)
        if (chr1==chr2){
          gene_diff[g1, g2] <- abs(start1-start2)
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
 
df.output <- matrix(unlist(output), ncol = 4, byrow = TRUE)

pos_diff <- matrix(df.output[,1], ncol = numb_orthogroups, byrow = TRUE)
pos_diff_max <- matrix(df.output[,2], ncol = numb_orthogroups, byrow = TRUE)
gene_names1 <- matrix(df.output[,3], ncol = numb_orthogroups, byrow = TRUE)
gene_names2 <- matrix(df.output[,4], ncol = numb_orthogroups, byrow = TRUE)

rownames(pos_diff) <- orthogroup_list[c(1:numb_orthogroups)]
colnames(pos_diff) <- orthogroup_list[c(1:numb_orthogroups)]

rownames(pos_diff_max) <- orthogroup_list[c(1:numb_orthogroups)]
colnames(pos_diff_max) <- orthogroup_list[c(1:numb_orthogroups)]

rownames(gene_names1) <- orthogroup_list[c(1:numb_orthogroups)]
colnames(gene_names1) <- orthogroup_list[c(1:numb_orthogroups)]

rownames(gene_names2) <- orthogroup_list[c(1:numb_orthogroups)]
colnames(gene_names2) <- orthogroup_list[c(1:numb_orthogroups)]


utils::write.table(x=pos_diff, file=paste0(path,"/", spp_hap, "_", numb_orthogroups,"_synteny_matrix.txt"), append = FALSE, sep = "\t", dec = ".", row.names = TRUE, col.names = TRUE)
utils::write.table(x=pos_diff_max, file=paste0(path,"/", spp_hap, "_", numb_orthogroups,"_synteny_matrix_max.txt"), append = FALSE, sep = "\t", dec = ".", row.names = TRUE, col.names = TRUE)
utils::write.table(x=gene_names1, file=paste0(path,"/", spp_hap, "_", numb_orthogroups,"_synteny_matrix_genes1.txt"), append = FALSE, sep = "\t", dec = ".", row.names = TRUE, col.names = TRUE)
utils::write.table(x=gene_names2, file=paste0(path,"/", spp_hap, "_", numb_orthogroups,"_synteny_matrix_genes2.txt"), append = FALSE, sep = "\t", dec = ".", row.names = TRUE, col.names = TRUE)


#################################
# compare matrices for different species

mkdir /home/celphin/scratch/Oxyria/synteny_quantity
cd /home/celphin/scratch/Oxyria/synteny_quantity

tmux new-session -s synteny2
tmux attach-session -t synteny2

module load  StdEnv/2020 r/4.2.2
export R_LIBS_USER=/home/celphin/R/x86_64-pc-linux-gnu-library/4.1.0/

R

path="/home/celphin/scratch/Oxyria/synteny_quantity"
numb_orthogroups=500
spp_hap1="Oxyria_digyna_H1"
#spp_hap2="Dryas_octopetala"
spp_hap2="Rheum_nobile_H0"

comb_file="Total_combBed.txt"

#--------------------------
#load files
pos_diff1 <- base::as.matrix(utils::read.table(paste0(path,"/", spp_hap1, "_", numb_orthogroups,"_synteny_matrix.txt"), row.names = 1, header = TRUE, check.names = FALSE))
pos_diff2 <- base::as.matrix(utils::read.table(paste0(path,"/", spp_hap2, "_", numb_orthogroups,"_synteny_matrix.txt"), row.names = 1, header = TRUE, check.names = FALSE))

gene1_names_spp1 <- base::as.matrix(utils::read.table(paste0(path,"/", spp_hap1, "_", numb_orthogroups,"_synteny_matrix_genes1.txt"), row.names = 1, header = TRUE, check.names = FALSE))
gene1_names_spp2 <- base::as.matrix(utils::read.table(paste0(path,"/", spp_hap2, "_", numb_orthogroups,"_synteny_matrix_genes1.txt"), row.names = 1, header = TRUE, check.names = FALSE))


combined_genes <- base::as.data.frame(utils::read.table(paste0(path,"/", comb_file), sep="\t", header = TRUE, check.names = FALSE))
colnames(combined_genes) <- c( "chr", "start", "end", "gene", "ofID", "pepLen", "ord", "genome", 
"arrayID","isArrayRep", "orthogroup", "globHOG", "noAnchor", "og")

gene_pos <- cbind(combined_genes$gene, combined_genes$chr, combined_genes$start, combined_genes$orthogroup, combined_genes$genome)

#---------------------------
# find min diff and min sum
synteny_diff = abs(pos_diff1 - pos_diff2)
synteny_sum = pos_diff1 + pos_diff2

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

#-----------------------------
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

grDevices::png(paste0(path,"/", spp_hap2, "_",spp_hap1, "_", numb_orthogroups,"_synteny_on_genome.png"), height=2500, width=5000)
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
  

####################################
# make heatmap of synteny diff and synteny sum

grDevices::png(paste0(path,"/", spp_hap2, "_",spp_hap1, "_", numb_orthogroups,"_heatmap_syntenydiff.png"), height=4000, width=5000)
  print(
  heatmap(synteny_diff)
  )
  grDevices::dev.off()
  
grDevices::png(paste0(path,"/", spp_hap2, "_",spp_hap1, "_", numb_orthogroups,"_heatmap_syntenysum.png"), height=4000, width=5000)
  print(
  heatmap(synteny_sum)
  )
  grDevices::dev.off()

#################################
# try running with GO terms not orthogroups
# seems like some genes are in many orthogroups

# join gff3 and interproscan output 

# same as orthologs but with GO terms