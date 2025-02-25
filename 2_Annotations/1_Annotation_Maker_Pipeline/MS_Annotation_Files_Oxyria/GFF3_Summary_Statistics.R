library(GenomicFeatures)
library(Biostrings)
library(dplyr)
library(optparse)
library(GenomicRanges)
library(rtracklayer)
library(BSgenome)
library("stargazer")
library(reshape2)

option_list = list(
  make_option(c("-g", "--gff"), type="character", default=NULL, 
              help="Required -GFF3 file", metavar="GFF3_file"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="Required Output Name file", metavar="Output_Name"))

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

GFF3_File_Name = opt$gff
Output_Name = opt$output
			  
GFF3_file <- as.data.frame(readGFF(GFF3_File_Name))

GFF3_file_summary <- group_by(GFF3_file,seqid, type)
GFF3_file_summary <- summarize(GFF3_file_summary, Value = n())

GFF3_file_summary_wide <- dcast(GFF3_file_summary,seqid~type, value.var = "Value")
ncol(GFF3_file_summary_wide)


#Table to text file
stargazer(GFF3_file_summary_wide,                 # Export txt
          summary = FALSE,
          type = "text",
          out = paste0(Output_Name,"_GFF3_summary.txt"))


sink(paste0(Output_Name,"_GFF3_summary.txt"), append = TRUE)
colSums(GFF3_file_summary_wide[,2:ncol(GFF3_file_summary_wide)], na.rm = TRUE)
sink()
