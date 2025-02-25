########################################################################
#NCBI Annotation of Oxyria digyna :
#Last updated: Jan 2024
#######################################################################
#January 2024
#Set up:
#Change Chromosome naming conventions
#Change all instances of Oxyria to Oxyria_NCBI
#Change genome name
########################################################################
cp Snakefile_Oxyria_Rheum Snakefile_Oxyria_NCBI

#----------------------------
#Format chromosomes:
module load StdEnv/2020 bioawk/1.0
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}'  Oxyria_digyna.fasta  |\
awk -F '\t' '{printf("%d\t%s\n",length($2),$0);}' |\
sort -k1,1nr | cut -f 2- | tr "\t" "\n" > Oxyria_ragtag.scaffold_sorted.fasta
bioawk -c fastx '{ print ">Oxyrt-" ++i  "-" length($seq) "\n" $seq}'  < Oxyria_ragtag.scaffold_sorted.fasta > Oxyria_ncbi.fasta
awk 'BEGIN{FS="[> ]"} /^>/{val=$2;next}  {print val,length($0)}'  Oxyria_ncbi.fasta | head -n 23

grep ">" Oxyria_ncbi.fasta | head -n 20
>Oxyrt-1-86582034
>Oxyrt-2-79714091
>Oxyrt-3-79472951
>Oxyrt-4-78410798
>Oxyrt-5-76064323
>Oxyrt-6-73303751
>Oxyrt-7-72361354
>Oxyrt-8-755699
>Oxyrt-9-538353
>Oxyrt-10-287654
>Oxyrt-11-284167
>Oxyrt-12-166496
>Oxyrt-13-127425
>Oxyrt-14-124062
>Oxyrt-15-123133

#----------------------------
#In snake file replace:
#Oxyria.fasta -> Oxyria_ncbi.fasta
#project Oxyria: -> Oxyria_NCBI

tmux new-session -s Oxyria_NCBI
module load  StdEnv/2020 r/4.1.0
export R_LIBS_USER=/home/msandler/R/x86_64-pc-linux-gnu-library/4.1.0/
R
library(GenomicRanges)
library(Biostrings)
library(optparse)
library(GenomicFeatures)
library(Biostrings)
library(ORFik)
library(BSgenome)
library(rtracklayer)
#q()


module load StdEnv/2020 intel/2020.1.217 metagenome-atlas/2.5.0 
snakemake --snakefile Snakefile_Oxyria_NCBI --profile cc-slurm 

#--------------------------------
# Make new protein file for current plant family


