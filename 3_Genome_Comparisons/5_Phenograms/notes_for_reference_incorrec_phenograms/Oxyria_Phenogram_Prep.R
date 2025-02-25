#module load  StdEnv/2020 r/4.1.0
#export R_LIBS_USER=/home/msandler/R/x86_64-pc-linux-gnu-library/4.1.0/

#R
library(dplyr)
library(tidyr)

#Make go terms - go name file:
process_go_terms_file <- function(gotermfile){
    goterm_file_lines <- readLines(gotermsfile)
    pattern <- c("GO:")
    filtered_goterm_lines <- goterm_file_lines[grep(pattern, goterm_file_lines)]
    filtered_goterms <- data.frame(filtered_goterm_lines)
    colnames(filtered_goterms) <- c("Lines")
    split_file <- data.frame(separate(filtered_goterms, Lines, into = c("IPRscan", "Name", "GO_Term"), sep = ">|;"))
    goterms <- split_file[,c("Name", "GO_Term")]
    goterms$GO_Terms <- gsub("\\s", "", goterms$GO_Term)
    no_dub_goterms <- distinct(goterms, GO_Term, Name)
    print(head(no_dub_goterms))
    return(no_dub_goterms)

}

gotermsfile <- "interpro2go"
gotermstoname <- process_go_terms_file(gotermfile)


process_annotation_goterms <- function(annotation_file, genetogofile) {
    annotation_file <- read.table(annotation_file, header=FALSE, sep='\t')
    genetogo <- read.table(genetogofile, header=FALSE, sep='\t')
    colnames(annotation_file) <- c("type", "Start", "GeneID")
    colnames(genetogo) <- c("GeneID", "GO_Term")

    annotation_file$GeneID <- gsub("^.*ID=([^;]+);.*$", "\\1", annotation_file$GeneID)
    genes_annotation_file <- subset(annotation_file, grepl("^gene", annotation_file$type))
    genes_annotation_file <- unique(genes_annotation_file)
    genes_annotation_file$CHR <- as.integer(gsub("^.*Chr(\\d).*$", "\\1", genes_annotation_file$GeneID)) #get chromosome number

    
    # Merge df1 with filtered_df2
    merged_df <- merge(genes_annotation_file, genetogo, by = "GeneID")
    print(head(merged_df))
#Uncollapse GO_Terms

    phenogram_df <- merged_df %>% separate_rows(GO_Term, sep = ",")  # trim whitespace if necessary
    phenogram_df <- data.frame(phenogram_df)
    phenogram_table <- phenogram_df[, c("CHR", "Start", "GO_Term")]
    print(head(phenogram_table))
    return(phenogram_table)
}

Oxydig_expanded <- process_annotation_goterms("Oxydig_Phenograms/Oxydig_expanded_genes_annotation", "Oxydig_Phenograms/Oxydig_expanded_goterms")
Oxydig_contracted <- process_annotation_goterms("Oxydig_Phenograms/Oxydig_contracted_genes_annotation", "Oxydig_Phenograms/Oxydig_contracted_goterms") 
OxyRheum_expanded <- process_annotation_goterms("OxyRheum_Phenogram/OxyRheum_expanded_genes_annotation", "OxyRheum_Phenogram/OxyRheum_expanded_goterms")
OxyRheum_contracted <- process_annotation_goterms("OxyRheum_Phenogram/OxyRheum_contracted_genes_annotation", "OxyRheum_Phenogram/OxyRheum_contracted_goterms")







make_phenogram <- function(gotermstoname, gene_table, filename){
    gotermstoname$GO_Term_id <- gsub("\\s+", "", gotermstoname$GO_Term)
    gene_table$GO_Term_id <- gsub("\\s+", "", gene_table$GO_Term)
    phenogram_table_interpretable <- merge(gotermstoname, gene_table, by = "GO_Term_id")
    phenogram_final <- phenogram_table_interpretable[,c("CHR", "Start", "Name")]
    colnames(phenogram_final) <- c("CHR", "POS", "PHENOTYPE")
    write.table(phenogram_final, file = filename, sep = "\t", quote = FALSE, row.names = FALSE)
}

make_phenogram(gotermstoname, Oxydig_expanded, "Oxydig_Phenograms/Oxydig_expanded.txt")
make_phenogram(gotermstoname, Oxydig_contracted, "Oxydig_Phenograms/Oxydig_contracted.txt")
make_phenogram(gotermstoname, OxyRheum_expanded, "OxyRheum_Phenogram/OxyRheum_expanded.txt")
make_phenogram(gotermstoname, OxyRheum_contracted, "OxyRheum_Phenogram/OxyRheum_contracted.txt")



    






