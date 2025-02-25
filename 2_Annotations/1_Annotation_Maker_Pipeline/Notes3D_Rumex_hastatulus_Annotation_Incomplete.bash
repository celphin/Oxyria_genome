########################################################################
#NCBI Annotation of Rumex hastatulus :
#Last updated: Jan 2024
#INCOMPLETE
#######################################################################
#January 2024
#Set up:
#Change Chromosome naming conventions
#Change all instances of Oxyria to Rheum
#Change genome name
#Change chromosome splitter

#Note: incomplete, scaffolds + 1 chromosome didn't properly run with EDTA

########################################################################
module load StdEnv/2020 bioawk/1.0
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}'  R_hastatulus_linkage.fasta  |\
awk -F '\t' '{printf("%d\t%s\n",length($2),$0);}' |\
sort -k1,1nr | cut -f 2- | tr "\t" "\n" > R_hastatulus.scaffold_sorted.fasta
bioawk -c fastx '{ print ">Rhast-" ++i  "-" length($seq) "\n" $seq}'  < R_hastatulus.scaffold_sorted.fasta > R_hastatulus.fasta
awk 'BEGIN{FS="[> ]"} /^>/{val=$2;next}  {print val,length($0)}'  R_hastatulus.fasta | head -n 23


cp Snakefile_Oxyria Snakefile_Rumex
cp Assembly_Chr_splitter_Oxyria.R Assembly_Chr_splitter_Rumex.R
#Change Assembly splitter Rheum
for (i in 1:5)
  {
  fasta_name = paste0(names(Assembly)[i],".fasta")
  writeXStringSet(Assembly[i],file=fasta_name)
}


print("Scaffold 6 contains all minor scaffolds")
fasta_name = "scaffolds.fasta"
writeXStringSet(Assembly[6:length(Assembly)],file=fasta_name)

#Adjust first parameters:
CHRS = '1 2 3 4 5'.split()
PROJECT = "R_hastatalus"
REFERENCE = "R_hastatulus.fasta"
NANOPORE_FASTQ = "Nanopore.fastq"
AED_FILTER = '0.6 1'.split()

rule Chr_splitting:
	input:
		rules.Init.output.Main_fasta,
	output:
		expand("{Project}/Ref/scaffold_{Chrs}.fasta",Project=PROJECT,Chrs=CHRS),
	params:
		project=PROJECT,
	threads: 1
	resources:
		mem_mb=4000, #4GB
		time=30 #minutes
	shell:
		"""
		echo Assembly_Chr_splitter_Rumex.R --args -f {input} 
		echo Assembly split into Chromosomes
		module load  StdEnv/2020 r/4.1.0
		R --vanilla < Assembly_Chr_splitter_Rumex.R --args -f {input} &&
		mv Rhast-*-*.fasta {params.project}/Ref/
		mv scaffolds.fasta {params.project}/Ref/
		echo "finished split"
		mv {params.project}/Ref/Rhast-1-*.fasta {params.project}/Ref/scaffold_1.fasta
		mv {params.project}/Ref/Rhast-2-*.fasta {params.project}/Ref/scaffold_2.fasta
		mv {params.project}/Ref/Rhast-3-*.fasta {params.project}/Ref/scaffold_3.fasta
		mv {params.project}/Ref/Rhast-4-*.fasta {params.project}/Ref/scaffold_4.fasta
		mv {params.project}/Ref/Rhast-5-*.fasta {params.project}/Ref/scaffold_5.fasta
		mv {params.project}/Ref/scaffolds.fasta {params.project}/Ref/scaffold_6.fasta
		module unload r/4.1.0
		"""
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
snakemake --snakefile Snakefile_Rumex --profile cc-slurm 
####
module load seqkit/2.3.1
seqkit seq -m 1000 R_hastatulus.fasta > R_hastatulus_cut.fasta 