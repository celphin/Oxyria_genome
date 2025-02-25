##Interpro_merger
#remotes::install_github("vitkl/ProtDomSeq")
# read InterProScan result, download and add InterPro Entry Types information, extract from relevant columns and add names metadata and sequence length information
#InterProScan_result = readInterProGFF3("./processed_data_files/all_human_viral_protein_domains.gff3.gz")
#InterProScan_result = addInterProEntryTypes(InterProScan_result, "./data_files/entry.list")
# create a subset that contains "Domain", "Active_site", "Binding_site", "Conserved_site", "PTM" signatures
#InterProScan_domains = SubsetByInterProEntryType(InterProScan_result)

#library(ProtDomSeq)
library(dplyr)

setwd("/DATA/home/jmlazaro/Projects/Mimulus/CE10/")

gff.step2_filtered_genes <- as.data.frame(readGFF("FINAL_ANNOTATION/FINAL_CE10.AED_0.8.sorted.gff3"))

setwd("/DATA/home/jmlazaro/Projects/Mimulus/CE10/")


INTERPRO_Data <- read.csv("INTERPRO/CE10.AED_0.8_protein.interpro.fasta_curated", header = FALSE)
INTERPRO_Data$Pfam <- paste0(INTERPRO_Data$V2,":",INTERPRO_Data$V3)
INTERPRO_Data$IPR <- paste0(INTERPRO_Data$V4,":",INTERPRO_Data$V5)
INTERPRO_Data$GO <- INTERPRO_Data$V6
INTERPRO_Data$Name <- INTERPRO_Data$V1

INTERPRO_Data_Curated <- INTERPRO_Data[,c("Name","Pfam","IPR","GO")]
INTERPRO_Data_Curated <-INTERPRO_Data_Curated %>% distinct(Name, .keep_all = TRUE)

###Flag2 Inter Assembly Annotation
gff.with_functional_annotation <- merge(gff.step2_filtered_genes, INTERPRO_Data_Curated,
                         by ="Name", all.x = TRUE)

gff.with_functional_annotation$Parent <- sub("character","",gff.with_functional_annotation$Parent)
gff.with_functional_annotation$Parent <- sub("(0)","",gff.with_functional_annotation$Parent)
gff.with_functional_annotation$Parent <- sub("\\(","",gff.with_functional_annotation$Parent)
gff.with_functional_annotation$Parent <- sub("\\)","",gff.with_functional_annotation$Parent)

gff.with_functional_annotation_sorted <- gff.with_functional_annotation %>% arrange(seqid, start)
gff.with_functional_annotation_sorted$Alias <-gff.with_functional_annotation_sorted$ID

export(gff.with_functional_annotation,"CE10.with_functional_annotation.gff3",format="gff3")

system("/DATA/home/jmlazaro/github/gff3sort/gff3sort.pl --chr_order original CE10.with_functional_annotation.gff3 > CE10.with_functional_annotation.sorted.gff3")
