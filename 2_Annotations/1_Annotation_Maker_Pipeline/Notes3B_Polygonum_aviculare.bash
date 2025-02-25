########################################################################
#    Polygonum aviculare  Annotation specific modification notes:
#    Last updated: February 2024
#   
#    Specific changes to Snakemake pipeline as provided in Github to run on Polygonum aviculare 
#######################################################################

#Snakemake pipeline modifications: 
#make sure config.yaml matches parameters in Notes1
    #snakefile still needs modifications, and need to change parameters + chromosome splitting file but overall working until AED
    #Copy reference Fasta to Maker folder
    #Modify Chrsplitter script, to account for number of chromosomes
    #Modify Chrsplitter for proper copying
    #Running with maker until then 
    #adjust Rules for proper protein file
	
	
################################################
# data download
# from https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_934048045.1/
cd ~/projects/def-rieseber/Dryas_shared_data/Oxyria/Genomes_for_annotation
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/934/048/045/GCA_934048045.1_dcPolAvic1.1/GCA_934048045.1_dcPolAvic1.1_genomic.fna.gz
gunzip GCA_934048045.1_dcPolAvic1.1_genomic.fna.gz
mv GCA_934048045.1_dcPolAvic1.1_genomic.fna     Polygonum_aviculare_genome.fasta
 
#Set up:
cd ~/scratch/Annotation/
mkdir Polygonum_aviculare; cd Polygonum_aviculare

################################
# renaming scaffolds
# https://gist.github.com/darencard/e1933e544c9c96234e86d8cbccc709e0
# sort fasta by length
#https://www.biostars.org/p/153999/#154055

module load StdEnv/2020 bioawk/1.0
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}'  Polygonum_aviculare_genome.fasta  |\
awk -F '\t' '{printf("%d\t%s\n",length($2),$0);}' |\
sort -k1,1nr | cut -f 2- | tr "\t" "\n" > Polygonum_aviculare_sorted.fasta
bioawk -c fastx '{ print ">Polavi-" ++i "\n" $seq}'  < Polygonum_aviculare_sorted.fasta > Polygonum_aviculare.fasta
awk 'BEGIN{FS="[> ]"} /^>/{val=$2;next}  {print val,length($0)}'  Polygonum_aviculare.fasta | head -n 23

cp Polygonum_aviculare.fasta ../CAP_Snakemake/

##################################################
# Snakefile changes

# make a copy of the Snakefile_Oxyria
cp Snakefile_Oxyria Snakefile_Polygonum_aviculare

#--------------------------
#Adjust Assembly splitter file 
cp Assembly_Chr_splitter.R Assembly_Chr_splitter_Polavi.R
nano Assembly_Chr_splitter_Polyavi.R

for (i in 1:10)
  {
  fasta_name = paste0(names(Assembly)[i],".fasta")
  writeXStringSet(Assembly[i],file=fasta_name)
}

print("Scaffolds  contains all minor scaffolds")
fasta_name = "scaffolds.fasta"
writeXStringSet(Assembly[11:length(Assembly)],file=fasta_name)


#--------------------------
#Adjust first parameters:
nano Snakefile_Polygonum_aviculare

CHRS = '1 2 3 4 5 6 7 8 9 10 11'.split()
PROJECT = "Polavi" # use name that you changed scaffolds to above
REFERENCE = "Polygonum_aviculare.fasta"
NANOPORE_FASTQ = "Nanopore.fastq"
AED_FILTER = '0.6 1'.split()

...
# Chr_splitting: Split the Main Assembly in Chromosomes for easy handling
...
                R --vanilla < Assembly_Chr_splitter_{params.project}.R --args -f {input} &&
                mv {params.project}*.fasta {params.project}/Ref/
                mv scaffolds.fasta {params.project}/Ref/
                cd {params.project}/Ref/
                rename {params.project}- scaffold_ *
                mv scaffolds.fasta scaffold_11.fasta
                cd ../..

#--------------------------------
# Adjust join at the end

nano Snakefile_Polygonum_aviculare

                """
                cat {params.project}/Post_Maker_Files/predicted_Pseudogenes_MAKER_ORF_Filtered_{params.project}.scaffold_1.AED_{wildcards.AED_filter}.gff3 > {params.project}/FINAL_ANNOTATION/FINAL_{params.project}.AED_{wildcards.AED_filter}.gff3

                for i in {{2,3,4,5,6,7,8,9,10,11}} 
                        do
                        tail -n +4 {params.project}/Post_Maker_Files/MAKER_ORF_Filtered_{params.project}.scaffold_$i.AED_{wildcards.AED_filter}.gff3 >> {params.project}/FINAL_ANNOTATION/FINAL_{params.project}.AED_{wildcards.AED_filter}.gff3
                        tail -n +4 {params.project}/Post_Maker_Files/predicted_Pseudogenes_MAKER_ORF_Filtered_{params.project}.scaffold_$i.AED_{wildcards.AED_filter}.gff3 >> {params.project}/FINAL_ANNOTATION/FINAL_{params.project}.AED_{wildcards.AED_filter}.gff3
                        done

                ../gff3sort/gff3sort.pl --chr_order original {params.project}/FINAL_ANNOTATION/FINAL_{params.project}.AED_{wildcards.AED_filter}.gff3 > {params.project}/FINAL_ANNOTATION/FINAL_{params.project}.AED_{wildcards.AED_filter}.sorted.gff3
                """
				
				
##############################################################################
#Add correct protein file 
#Running with Rheum genome as protein file:
#Change init rule:
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
###################################################
#run

tmux new-session -t Annotation_Polygonum
tmux attach-session -t Annotation_Polygonum

module load  StdEnv/2020 r/4.1.0
export R_LIBS_USER=/home/celphin/R/x86_64-pc-linux-gnu-library/4.1.0/

R 
library(GenomicRanges)
library(Biostrings)
library(optparse)
library(GenomicFeatures)
library(Biostrings)
library(ORFik)
library(BSgenome)
library(rtracklayer)
q()


module load StdEnv/2020 intel/2020.1.217 metagenome-atlas/2.5.0 
snakemake --snakefile Snakefile_Polygonum_aviculare --profile cc-slurm 
