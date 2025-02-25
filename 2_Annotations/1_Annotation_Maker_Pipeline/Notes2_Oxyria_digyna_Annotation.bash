############################################################
#Changes specific to Oxyria digyna
	#1. Renaming scaffolds
	#2. Make Assembly Chromosome splitter modifications
	#3. Change header for proper inputs (num chromosomes, file names)
	#4. Adjust join at end (from default CAP_Snakemake file)
	#5. Make Maker_Files_Polygon.zip which uses Rheum tangaticum instead of arabidopsis as protein file for homology
	#6. Modify inside of snakefile
	#6. Run the pipeline
############################################################
# Renaming scaffolds
grep ">" Oxy_dig_1.scaffolds_FINAL.final.review.fasta
#Very long list
>PGA_scaffold_1__1_contigs__length_77175000
>PGA_scaffold_2__2_contigs__length_80787722
>PGA_scaffold_3__1_contigs__length_76036264
>PGA_scaffold_4__7_contigs__length_70136240
>PGA_scaffold_5__1_contigs__length_73842794
>PGA_scaffold_6__1_contigs__length_73289246
>PGA_scaffold_7__3_contigs__length_88146246
>PGA_scaffold_8__1_contigs__length_3898474


module load StdEnv/2020 bioawk/1.0
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}'  Oxy_dig_1.scaffolds_FINAL.final.review.fasta  |\
awk -F '\t' '{printf("%d\t%s\n",length($2),$0);}' |\
sort -k1,1nr | cut -f 2- | tr "\t" "\n" > Oxyria.scaffold_sorted.fasta
bioawk -c fastx '{ print ">Oxy-" ++i  "-" length($seq) "\n" $seq}'  < Oxyria.scaffold_sorted.fasta > Oxyria.fasta
awk 'BEGIN{FS="[> ]"} /^>/{val=$2;next}  {print val,length($0)}'  Oxyria.fasta | head -n 23

grep ">" Oxyria.fasta

>Oxy-1-88146246
>Oxy-2-80787722
>Oxy-3-77175000
>Oxy-4-76036264
>Oxy-5-73842794
>Oxy-6-73289246
>Oxy-7-70136240
...
#######################################################
cp Assembly_Chr_splitter.R Assembly_Chr_splitter_Oxyria.R
nano Assembly_Chr_splitter_Oxyria.R

for (i in 1:7)
  {
  fasta_name = paste0(names(Assembly)[i],".fasta")
  writeXStringSet(Assembly[i],file=fasta_name)
}


print("Scaffold 8 contains all minor scaffolds")
fasta_name = "scaffolds.fasta"
writeXStringSet(Assembly[8:length(Assembly)],file=fasta_name)

########################################################
nano Snakefile_Oxyria
#Adjust first parameters:
CHRS = '1 2 3 4 5 6 7 8'.split()
PROJECT = "Oxyria"
REFERENCE = "Oxyria.fasta"
NANOPORE_FASTQ = "Nanopore.fastq"
AED_FILTER = '0.6 1'.split()


echo Assembly_Chr_splitter_Oxyria.R --args -f {input} 
echo Assembly split into Chromosomes
module load  StdEnv/2020 r/4.1.0
R --vanilla < Assembly_Chr_splitter_Oxyria.R --args -f {input} &&
mv Oxy-*-*.fasta {params.project}/Ref/
mv Oxyria/Ref/Oxy-1-88146246*.fasta Oxyria/Ref/scaffold_1.fasta
mv Oxyria/Ref/Oxy-2-80787722*.fasta Oxyria/Ref/scaffold_2.fasta
mv Oxyria/Ref/Oxy-3-77175000*.fasta Oxyria/Ref/scaffold_3.fasta
mv Oxyria/Ref/Oxy-4-76036264*.fasta Oxyria/Ref/scaffold_4.fasta
mv Oxyria/Ref/Oxy-5-73842794*.fasta Oxyria/Ref/scaffold_5.fasta
mv Oxyria/Ref/Oxy-6-73289246*.fasta Oxyria/Ref/scaffold_6.fasta
mv Oxyria/Ref/Oxy-7-70136240*.fasta Oxyria/Ref/scaffold_7.fasta
mv Oxyria/Ref/scaffolds.fasta Oxyria/Ref/scaffold_8.fasta

#########################################################################
#Adjust join at end:
"""
		cat {params.project}/Post_Maker_Files/predicted_Pseudogenes_MAKER_ORF_Filtered_{params.project}.scaffold_1.AED_{wildcards.AED_filter}.gff3 > {params.project}/FINAL_ANNOTATION/FINAL_{params.project}.AED_{wildcards.AED_filter}.gff3

        for i in {{2,3,4,5,6,7, 8}} 
            do
            tail -n +4 {params.project}/Post_Maker_Files/MAKER_ORF_Filtered_{params.project}.scaffold_$i.AED_{wildcards.AED_filter}.gff3 >> {params.project}/FINAL_ANNOTATION/FINAL_{params.project}.AED_{wildcards.AED_filter}.gff3
            tail -n +4 {params.project}/Post_Maker_Files/predicted_Pseudogenes_MAKER_ORF_Filtered_{params.project}.scaffold_$i.AED_{wildcards.AED_filter}.gff3 >> {params.project}/FINAL_ANNOTATION/FINAL_{params.project}.AED_{wildcards.AED_filter}.gff3
        	done

        ../gff3sort/gff3sort.pl --chr_order original {params.project}/FINAL_ANNOTATION/FINAL_{params.project}.AED_{wildcards.AED_filter}.gff3 > {params.project}/FINAL_ANNOTATION/FINAL_{params.project}.AED_{wildcards.AED_filter}.sorted.gff3
        """

#####################################################################
unzip Maker_Files.zip
cp -r Maker_Files Maker_Files_Polygon
cd Maker_Files_Polygon
wget https://figshare.com/ndownloader/files/39014222
mv 39014222 Rheum_tanguticum_Proteins.fa
cd ..
zip -r Maker_Files_Polygon.zip Maker_Files_Polygon

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

#########################################################################
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
snakemake --snakefile Snakefile_Oxyria_Rheum --profile cc-slurm 



