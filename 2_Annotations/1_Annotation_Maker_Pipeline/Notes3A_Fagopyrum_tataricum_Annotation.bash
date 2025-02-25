########################################################################
#Fagopyrum tataricum Annotation specific modification notes:
#Last updated: February 2024
#Final cleaned annotation notes for fagopyrum tataricum
#######################################################################
#Set up:
cd scratch/Oxyria/Polygonaceae_Genomes/Fagopyrum_tataricum
mkdir Annotation; cd Annotation

#Copy essential files:
cp /home/msandler/scratch/Oxyria/Annotation/CAP_Snakemake/* . 
mv Snakefile_Oxyria Snakefile_Fagopyrum_Rheum

module load StdEnv/2020 bioawk/1.0
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}'  Fagopyrum_genome.fasta  |\
awk -F '\t' '{printf("%d\t%s\n",length($2),$0);}' |\
sort -k1,1nr | cut -f 2- | tr "\t" "\n" > Fagopyrum.scaffold_sorted.fasta
bioawk -c fastx '{ print ">Fagopyrum-" ++i  "-" length($seq) "\n" $seq}'  < Fagopyrum.scaffold_sorted.fasta > Fagopyrum.fasta
awk 'BEGIN{FS="[> ]"} /^>/{val=$2;next}  {print val,length($0)}'  Fagopyrum.fasta 
head -n 23 Fagopyrum.fasta

grep ">" Fagopyrum.fasta


>Fagopyrum-1-68031765
>Fagopyrum-2-61235386
>Fagopyrum-3-57706077
>Fagopyrum-4-56655744
>Fagopyrum-5-53883329
>Fagopyrum-6-52287906
>Fagopyrum-7-51545819
>Fagopyrum-8-49982843



##########################################################

#Adjust first parameters:
CHRS = '1 2 3 4 5 6 7 8 9'.split()
PROJECT = "Fagopyrum"
REFERENCE = "Fagopyrum.fasta"
NANOPORE_FASTQ = "Nanopore.fastq"
AED_FILTER = '0.6 1'.split()

#Adjust Assembly splitter file 
cp Assembly_Chr_splitter.R Assembly_Chr_splitter_Fagopyrum.R
nano Assembly_Chr_splitter_Fagopyrum.R


for (i in 1:8)
  {
  fasta_name = paste0(names(Assembly)[i],".fasta")
  writeXStringSet(Assembly[i],file=fasta_name)
}


print("Scaffold 9 contains all minor scaffolds")
fasta_name = "scaffold_9.fasta"
writeXStringSet(Assembly[9:length(Assembly)],file=fasta_name)






#In snake file, change init function:

#Change init file:
		echo $PWD
		mkdir -p {params.project}
		mkdir -p {params.project}/Ref
		module load StdEnv/2020 intel/2020.1.217 metagenome-atlas/2.5.0
		#snakemake --dag | dot -Tsvg > {params.project}/dag.svg
		unzip -o {input.maker_files} -d {params.project}

		mv {params.project}/Maker_Files_Polygon/* {params.project}/Maker_Files
		cp -v  {input.reference} {output.Main_fasta}
		#Index FASTA file
		module load  StdEnv/2020 samtools/1.17
		samtools faidx {output.Main_fasta} 
		module unload samtools/1.17
		cd {params.project}

  #Replace all instaces of Maker_Files.zip -> Maker_Files_Polygon.zip
#Replace all: csa.trans.Protein.10072011.fasta with Rheum_tanguticum_Proteins.fa


#Replace all instances of with chromosome  names
#modify moving (based on how fasta is built)
mv {params.project}/Ref/Fagopyrum-1-68031765*.fasta {params.project}/Ref/scaffold_1.fasta
mv {params.project}/Ref/Fagopyrum-2-61235386*.fasta {params.project}/Ref/scaffold_2.fasta
mv {params.project}/Ref/Fagopyrum-3-57706077*.fasta {params.project}/Ref/scaffold_3.fasta
mv {params.project}/Ref/Fagopyrum-4-56655744*.fasta {params.project}/Ref/scaffold_4.fasta
mv {params.project}/Ref/Fagopyrum-5-53883329*.fasta {params.project}/Ref/scaffold_5.fasta
mv {params.project}/Ref/Fagopyrum-6-52287906*.fasta {params.project}/Ref/scaffold_6.fasta
mv {params.project}/Ref/Fagopyrum-7-51545819*.fasta {params.project}/Ref/scaffold_7.fasta
mv {params.project}/Ref/agopyrum-8-49982843*.fasta {params.project}/Ref/scaffold_8.fasta
mv {params.project}/Ref/scaffolds.fasta {params.project}/Ref/scaffold_9.fasta

######################################################################
#Chr_Merge:
#change rule to:
chmod +755  ../gff3sort/gff3sort.pl

 for i in {{2,3,4,5,6,7,8,9}} 
            do
            tail -n +4 {params.project}/Post_Maker_Files/MAKER_ORF_Filtered_{params.project}.scaffold_$i.AED_{wildcards.AED_filter}.gff3 >> {params.project}/FINAL_ANNOTATION/FINAL_{params.project}.AED_{wildcards.AED_filter}.gff3
            tail -n +4 {params.project}/Post_Maker_Files/predicted_Pseudogenes_MAKER_ORF_Filtered_{params.project}.scaffold_$i.AED_{wildcards.AED_filter}.gff3 >> {params.project}/FINAL_ANNOTATION/FINAL_{params.project}.AED_{wildcards.AED_filter}.gff3
            done		
		../gff3sort/gff3sort.pl --chr_order original {params.project}/FINAL_ANNOTATION/FINAL_{params.project}.AED_{wildcards.AED_filter}.gff3 > {params.project}/FINAL_ANNOTATION/FINAL_{params.project}.AED_{wildcards.AED_filter}.sorted.gff3


snakemake --snakefile Snakefile_Fagopyrum_tataricum --profile cc-slurm 



##############################################################################
#To Run:
tmux new-session -t Annotation_Fagopyrum
tmux attatch-session -t Annotation_Fagopyrum

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

snakemake --snakefile Snakefile_Fagopyrum_Rheum --profile cc-slurm 

snakemake --snakefile Snakefile_Fagopyrum_Rheum --profile cc-slurm 

snakemake --snakefile Snakefile_Fagopyrum_Rheum --profile cc-slurm 