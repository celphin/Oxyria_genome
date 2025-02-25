#bedtools maskfasta -fi ../Ref/PK_Scaffold_Split/PK_hap2.scaffold_10.fasta -bed PK_hap2.scaffold_10.fasta.mod.EDTA.intact.gff3 -fo PK_hap2.scaffo ld_10.masked.fasta              

library(GenomicFeatures)
library(Biostrings)
library(dplyr)
library(ORFik)
library(optparse)
library(GenomicRanges)
library(rtracklayer)
library(BSgenome)

#usage: R --vanilla < Post_Maker_ORF_Analysis_Terminal.R 
#--args -g MAKER_Filtered_PK_hap2.scaffold_10.all.AED_0.8.gff3 -a PK_hap2.scaffold_10.fasta -m PK_hap2.scaffold_10.masked.fasta -o PK_hap2.scaffold_10

#https://cran.r-project.org/web/packages/optparse/readme/README.html

option_list = list(
  make_option(c("-g", "--gff_file"), type="character", default=NULL, 
              help="Required -GFF3 file", metavar="Fasta_File"),
  make_option(c("-a", "--assembly"), type="character", default=NULL, 
              help="Required unmasked Assembly Name file", metavar="Assembly_File"),
  make_option(c("-m", "--massembly"), type="character", default=NULL, 
              help="Required masked Assembly Name file", metavar="Assembly_File"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="Output Name file", metavar="Assembly_File")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

File_Name = opt$gff_file
Assembly_Name = opt$assembly
Masked_Assembly_Name = opt$massembly
Output_Name = opt$output


my_columns <- c("seqid","start","end","strand","type")
my_filter <- list(type="gene")

Maker_scaffold_gff3_raw <- as.data.frame(readGFF(File_Name))
Maker_scaffold_gff3_genes<- as.data.frame(readGFF(File_Name,columns = my_columns, filter = my_filter))
###Import Fasta file and transform into DNAstring Dataset,names correlated with GFF3
Assembly_DNA <- readDNAStringSet(Assembly_Name)
Masked_Assembly_DNA <- readDNAStringSet(Masked_Assembly_Name)

Maker_Genes_gr <- with(Maker_scaffold_gff3_genes, GRanges(seqid, IRanges(start, end),strand, id = Name))



Maker_All_genes <- getSeq(Masked_Assembly_DNA,Maker_Genes_gr)
Maker_Masked_genes <- as.data.frame(alphabetFrequency(Maker_All_genes))
Maker_Masked_genes$length <- seqlengths(Maker_All_genes)
Maker_Masked_genes$TE_Overlap_Percentage <- (Maker_Masked_genes$N/Maker_Masked_genes$length)*100
Maker_Masked_genes$Name <- Maker_Genes_gr@elementMetadata@listData$id


#Filter Genes less than 100bp and overlap below 75%
Maker_Masked_genes_filtered <- filter(Maker_Masked_genes, TE_Overlap_Percentage <=74.99)
Maker_Masked_genes_filtered <- filter(Maker_Masked_genes_filtered, length>=150)
Maker_Masked_genes_filtered_list <- Maker_Masked_genes_filtered[,c("Name","TE_Overlap_Percentage")]

Maker_genes_Curated_file <- merge(Maker_scaffold_gff3_genes,Maker_Masked_genes_filtered_list, by = "Name", all.y = TRUE)
####END STEP1
MAKER_raw.gr <- makeGRangesFromDataFrame(Maker_scaffold_gff3_raw, ignore.strand = FALSE,
                                          seqinfo = NULL,
                                          seqnames.field = "seqid",
                                          start.field = "start",
                                          end.field = "end",
                                          strand.field = "strand",
                                          keep.extra.columns = TRUE)


step1_filtered_genes.gr <- makeGRangesFromDataFrame(Maker_genes_Curated_file , ignore.strand = FALSE,
                                                    seqinfo = NULL,
                                                    seqnames.field = "seqid",
                                                    start.field = "start",
                                                    end.field = "end",
                                                    strand.field = "strand",
                                                    keep.extra.columns = TRUE)

MAKER_Filtered_GFF <- subsetByOverlaps(MAKER_raw.gr, step1_filtered_genes.gr)

###CDS Filtering
gff.MAKER_CDS <- filter(as.data.frame(MAKER_Filtered_GFF),type=="CDS")
MAKER_CDS_gr <- with(gff.MAKER_CDS, GRanges(seqnames, IRanges(start, end),strand, id = ID))


###OBTAINING DNA SEQUENCE FROM MASKED FASTA IN ALL THE GFF REGIONS
MAKER_GFF_DNA <- getSeq(Masked_Assembly_DNA,MAKER_Filtered_GFF)

#ORF Analysis START
##################
## make TxDb from GTF file 
txdb <- makeTxDbFromGRanges(MAKER_Filtered_GFF)

## get intron information
all.introns <- intronicParts(txdb)
all.exons <- exonicParts(txdb)

ORF_filtered_genes.gr <-GRanges()
Proteins_list <- list()

for (i in 1:length(step1_filtered_genes.gr)) 
{
  Gene_Number <- i
  print(Gene_Number)
  Sequence_to_Analyze <- getSeq(Masked_Assembly_DNA, step1_filtered_genes.gr[Gene_Number])
  Sequence_to_Analyze <- DNAString(as.character(Sequence_to_Analyze))
  
  Sequence_Strand <- as.character(step1_filtered_genes.gr[Gene_Number]@strand)
  Reverse_CDS=FALSE
  
  #Reverse complement if in the negative Strand
  if (Sequence_Strand=="-"){
    Sequence_to_Analyze <- reverseComplement(Sequence_to_Analyze)
    Reverse_CDS=TRUE
  }
  
  #ORF bigger than 150bp (50aa)
  Possible_ORFS <- as.data.frame(findORFs(as.character(Sequence_to_Analyze)))
  
  Possible_ORFS <- filter(Possible_ORFS, width > 150)
  
  Gene_Name <- as.character(step1_filtered_genes.gr[Gene_Number]@seqnames)
  Gene_Pos <- step1_filtered_genes.gr[Gene_Number]@ranges
  #ORFS count
  No_ORFS <- as.numeric(nrow(Possible_ORFS))
  
  #Exons and Introns Count
  Lookup_Gene <- as.data.frame(step1_filtered_genes.gr[Gene_Number])
  
  Lookup_Gene <- with(Lookup_Gene, GRanges(seqnames, IRanges(start, end),strand, id))
  
  #table(!is.na(findOverlaps(all.exons, Lookup_Gene, select="arbitrary")))
  
  CDS_region <- as.data.frame(subsetByOverlaps(MAKER_CDS_gr, Lookup_Gene))
  
  CDS_region <- with(CDS_region, GRanges(seqnames, IRanges(start, end),strand, id))
  ###Order the SEQUENCE
  CDS_region <- sortSeqlevels(CDS_region)
  CDS_region <- sort(CDS_region, decreasing = Reverse_CDS)
  
  if(length(CDS_region)>0)
  {
    CDS_Sequence <- getSeq(Assembly_DNA,CDS_region)
    CDS_Sequence <- do.call(xscat, CDS_Sequence)
    
    Protein_translated <- as.character(translate(CDS_Sequence, if.fuzzy.codon="solve"))
  } else
  {
    Protein_translated <- "NNNNNNNNNNNNNNNNNNNNNNN"
  }
  #Adding to protein list
  Proteins_list <- append(Proteins_list,c(paste0(">",step1_filtered_genes.gr[Gene_Number]$ID),Protein_translated))
  
  Exons_Count <- as.numeric(nrow(as.data.frame(subsetByOverlaps(all.exons, Lookup_Gene))))
  
  #table(!is.na(findOverlaps(all.introns, Lookup_Gene, select="arbitrary")))
  Introns_Count <- as.numeric(nrow(as.data.frame(subsetByOverlaps(all.introns, Lookup_Gene))))
  
  #Granges Object Creations
  Gene_Granges <- GRanges(seqnames = Gene_Name,ranges = Gene_Pos, strand = Sequence_Strand, 
                          ORFS_Number = No_ORFS, No_Exons = Exons_Count, No_Introns = Introns_Count )
  ORF_filtered_genes.gr <- append(ORF_filtered_genes.gr,Gene_Granges)
}

ORF_filtered_genes.DF <- as.data.frame(ORF_filtered_genes.gr)
Pseudogenes <- ORF_filtered_genes.DF

ORF_filtered_genes.DF <- filter(ORF_filtered_genes.DF,ORFS_Number >0)

ORF_filtered_genes.gr <- with(ORF_filtered_genes.DF, GRanges(seqnames, IRanges(start, end),
                                                             strand, ORFS_Number))

###Getting Step2 Sequence filtered GFF3
MAKER_Step2_Filtered_GFF <- as.data.frame(subsetByOverlaps(MAKER_raw.gr, ORF_filtered_genes.gr))
MAKER_Step2_Filtered_GFF <- MAKER_Step2_Filtered_GFF[,c("seqnames","start","end","width","strand","source","type","score","phase", "ID","Name","Alias","Parent")]

export(MAKER_Step2_Filtered_GFF,paste0("MAKER_ORF_Filtered_",Output_Name,".gff3"),format="gff3")

sink(paste0("predicted_proteins_MAKER_ORF_Filtered_",Output_Name,".fasta"))
writeLines(unlist(lapply(Proteins_list, paste, collapse=" ")))
sink()

####Pseudogenes
Pseudogenes$Pseudo <- ifelse(Pseudogenes$ORFS_Number==0,"PseudoGene","Gene")
Pseudogenes$Pseudo <- ifelse(Pseudogenes$No_Introns==0,"Likely_PseudoGene",Pseudogenes$Pseudo)
Pseudogenes$Gene_Name <- step1_filtered_genes.gr$ID

export(MAKER_Step2_Filtered_GFF,paste0("predicted_Pseudogenes_MAKER_ORF_Filtered_",Output_Name,".gff3"),format="gff3")


