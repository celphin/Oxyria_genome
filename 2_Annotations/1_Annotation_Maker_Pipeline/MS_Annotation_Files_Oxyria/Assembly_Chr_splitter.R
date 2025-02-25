library(GenomicRanges)
library(Biostrings)
library(optparse)
library(rtracklayer)

#usage: R --vanilla < Assembly_Chr_splitter.R --args -f PK_hap2.09082021.fasta  -a PK_hap2
#https://cran.r-project.org/web/packages/optparse/readme/README.html

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="Required -GFF3 file", metavar="Fasta_File"))

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

File_Name = opt$file

#Read Assembly (UNMASKED)
Assembly <- readDNAStringSet(File_Name)

for (i in 1:7)
  {
  fasta_name = paste0(names(Assembly)[i],".fasta")
  writeXStringSet(Assembly[i],file=fasta_name)
}


#print("Scaffold 8 contains all minor scaffolds")
#fasta_name = "scaffold_8.fasta"
#writeXStringSet(Assembly[8:length(Assembly)],file=fasta_name)




