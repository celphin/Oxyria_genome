#################################################################
#May 2024
#Visualization of genes associated with expanded GO terms as identified in ermine
#!!! TO DO : clean up like in Oxyria_Phenogram_Prep.R in reference notes folders
#################################################################
#1 Install Phenogram:
mkdir Phenograms; cd Phenograms
wget https://ritchielab.org/files/RL_software/ruby_install.sh
wget https://ritchielab.org/files/RL_software/pheno_gram.rb
sh ruby_install.sh
ruby pheno_gram.rb

#######################################################################################################################################
#2 Make file of:
    # Go terms and their names (modified from ermine) -  OxyRheumnob_Expanded_Name_GOTerm.txt
    # Genes associated with those GO terms, and their GO terms:  OxyRheumnob_expanded_GOtermmap.txt
    # Start and end of genes Identified above: genes_expanded_annotation_rnooxy.txt

mkdir OxyRheum_Phenogram2; cd OxyRheum_Phenogram2
 
head -n 43 ~/scratch/Oxyria/CAFE/enrichment_analysis/Oxydig_Rheumnob_expanded_genesets.ermine.results | tail -n 21 | cut -f2,3,8 >> OxyRheumnob_Expanded_Name_GOTerm.txt
cut -f2 OxyRheumnob_Expanded_Name_GOTerm.txt  > OxyRheumnob_GOterms.txt
grep -f OxyRheumnob_GOterms.txt  ~/scratch/Oxyria/CAFE/enrichment_analysis/Interproscan_files/Oxyria_AED0.6_interproscan.tsv | cut -f1,14 > OxyRheumnob_expanded_GOtermmap.txt
cut -f1 OxyRheumnob_expanded_GOtermmap.txt | uniq > expanded_rnooxygenelist

grep -f expanded_rnooxygenelists  ~/scratch/Oxyria/CAFE/enrichment_analysis/gff3/Oxyria_digyna_H1.gff3 | cut -f3,4,9 | grep "gene"  > genes_expanded_annotation_rnooxy.txt

#######################################################################################################################################
#Format for phenogram
#!!! TO DO : clean up like in Oxyria_Phenogram_Prep.R
module load  StdEnv/2020 r/4.1.0
export R_LIBS_USER=/home/msandler/R/x86_64-pc-linux-gnu-library/4.1.0/
R

library(dplyr)
library(tidyr)

genesfile <- "OxyRheumnob_expanded_GOtermmap.txt"
annotation_file <- "genes_expanded_annotation_rnooxy.txt"
gotermstoname_file <- "OxyRheumnob_Expanded_Name_GOTerm.txt"


annotation_file <- read.table(annotation_file, header=FALSE, sep='\t')
genetogo <- read.table(genesfile, header=FALSE, sep='\t')


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

    phenogram_df <- merged_df %>% separate_rows(GO_Term, sep = "\\|")  # trim whitespace if necessary
    phenogram_df <- data.frame(phenogram_df)
    phenogram_table <- phenogram_df[, c("CHR", "Start", "GO_Term")]
    print(head(phenogram_table))

  gene_table <- phenogram_table

  gotermstoname <-  read.table(gotermstoname_file, header=FALSE, sep='\t')
   colnames(gotermstoname) <- c("Name", "GO_Term", "P-val")
    gotermstoname$GO_Term_id <- gsub("\\s+", "", gotermstoname$GO_Term)
    gene_table$GO_Term_id <- gsub("\\s+", "", gene_table$GO_Term)
    phenogram_table_interpretable <- merge(gotermstoname, gene_table, by = "GO_Term_id")
    phenogram_final <- phenogram_table_interpretable[,c("CHR", "Start", "Name")]
    colnames(phenogram_final) <- c("CHR", "POS", "PHENOTYPE")
    write.table(phenogram_final, file = "OxyRheumnob_phenogram_input.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#############################################################################################
#4: remove 8th chromosome, remove duplicates
first_line=$(head -n 1 OxyRheum_Phenogram2/OxyRheumnob_phenogram_input.txt)

# Extract all lines except the first, sort them, and remove duplicates
sorted_unique_lines=$(tail -n +2 OxyRheum_Phenogram2/OxyRheumnob_phenogram_input.txt | awk '$1 != 8' | sort | uniq)

# Combine the first line with the sorted unique lines
echo "$first_line" > "OxyRheum_Phenogram2/OxyRheum_expanded_input2.txt"
echo "$sorted_unique_lines" >> "OxyRheum_Phenogram2/OxyRheum_expanded_input2.txt"

##############################################################################################
#5: Plot phenograms 
ruby pheno_gram.rb -i OxyRheum_Phenogram2/OxyRheum_expanded_input2.txt -g Oxyria_Genome.txt -t "GO terms of Expanded Orthogroup for R. nobile and O. digyna" -o OxyRheum_Phenogram2/oxyriarheum_expanded2 -f jpg


##############################################################################################
#Copy to locals
scp -v msandler@graham.computecanada.ca:/home/msandler/scratch/Oxyria/Phenograms/*/*.jpg .