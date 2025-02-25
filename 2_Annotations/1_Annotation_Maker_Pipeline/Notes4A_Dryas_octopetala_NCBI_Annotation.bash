########################################################################
#Annotation of Dryas octopetala: Eight-petal Mountain-Avens, GCA_963921425.1
#Last updated: June 2024
#######################################################################

#------------------------
# get genome assemblies
cd /home/celphin/scratch/Annotation/CAP_Snakemake

cp ~/scratch/Dryas/Dryas_genomes/Dry-int_ragtag_output/ragtag.scaffold.fasta Genome_Dry-int-chr.fasta
cp ~/scratch/Dryas/Dryas_genomes/Dry-drumm_ragtag_output/ragtag.scaffold.fasta Genome_Dry-drumm-chr.fasta
cp ~/scratch/Dryas/Dryas_genomes/Dry-alask_ragtag_output/ragtag.scaffold.fasta Genome_Dry-alask-chr.fasta
cp ~/scratch/Dryas/Dryas_genomes/Dry-octo-H1_ragtag_output/ragtag.scaffold.fasta Genome_Dry-octo-H1-chr.fasta
cp ~/scratch/Dryas/Dryas_genomes/Dry-octo-H2_ragtag_output/ragtag.scaffold.fasta Genome_Dry-octo-H2-chr.fasta
cp ~/scratch/Dryas/Dryas_genomes/Dry-octo-H0.fasta Genome_Dry-octo-H0-chr.fasta

#----------------------------------
# make protein training file from other species in Rosaceae 
mkdir protein_training_files; cd protein_training_files

# download protein files from other genomes
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_019419815.1/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/019/419/815/GCF_019419815.1_Pyrus_bretschneideri_v1/GCF_019419815.1_Pyrus_bretschneideri_v1_protein.faa.gz
gunzip GCF_019419815.1_Pyrus_bretschneideri_v1_protein.faa.gz

# https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_916048215.2/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/916/048/215/GCF_916048215.2_drMalSylv7.2/GCF_916048215.2_drMalSylv7.2_protein.faa.gz
gunzip GCF_916048215.2_drMalSylv7.2_protein.faa.gz

# https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000346465.2/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/346/465/GCF_000346465.2_Prunus_persica_NCBIv2/GCF_000346465.2_Prunus_persica_NCBIv2_protein.faa.gz
gunzip GCF_000346465.2_Prunus_persica_NCBIv2_protein.faa.gz

# https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000184155.1/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/184/155/GCF_000184155.1_FraVesHawaii_1.0/GCF_000184155.1_FraVesHawaii_1.0_protein.faa.gz
gunzip GCF_000184155.1_FraVesHawaii_1.0_protein.faa.gz

# https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_958449725.1/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/958/449/725/GCF_958449725.1_drRosRugo1.1/GCF_958449725.1_drRosRugo1.1_protein.faa.gz
gunzip GCF_958449725.1_drRosRugo1.1_protein.faa.gz

# https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_933775445.1/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/933/775/445/GCF_933775445.1_drPotAnse1.1/GCF_933775445.1_drPotAnse1.1_protein.faa.gz
gunzip GCF_933775445.1_drPotAnse1.1_protein.faa.gz

cat *protein.faa > csa.trans.Protein.Rosaceae.fasta

#--------------------------
# zip protein file into Make_files
unzip Maker_Files.zip
mv  Maker_Files.zip  Maker_Files_orig.zip
cd Maker_Files

# add a copy of csa.trans.Protein.Rosaceae.fasta to the folder
cp ./protein_training_files/csa.trans.Protein.Rosaceae.fasta ./Maker_Files/csa.trans.Protein.Rosaceae.fa

# rezip
zip -r Maker_Files.zip Maker_Files

#-------------------------------
# make new Assembly Chr splitter using 9 chromosomes
cp Assembly_Chr_splitter_Polavi.R Assembly_Chr_splitter_Dryoct.R

for (i in 1:9)
  {
  fasta_name = paste0(names(Assembly)[i],".fasta")
  writeXStringSet(Assembly[i],file=fasta_name)
}
...10:...

#----------------------------
#Format chromosome names to 

cd /home/celphin/scratch/Annotation/CAP_Snakemake/genome_assemblies
module load StdEnv/2020 bioawk/1.0

awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}'  Genome_Dry-octo-H0-chr.fasta |\
awk -F '\t' '{printf("%d\t%s\n",length($2),$0);}' |\
sort -k1,1nr | cut -f 2- | tr "\t" "\n" |\
bioawk -c fastx '{ print ">DoctH0-" ++i "\n" $seq}' > Dry-octo-H0_test.fasta

grep ">" Genome_Dry-octo-H0-chr.fasta | head -n 15
grep ">" Dry-octo-H0_test.fasta | head -n 15

#-------------------------
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}'  Genome_Dry-int-chr.fasta |\
awk -F '\t' '{printf("%d\t%s\n",length($2),$0);}' |\
sort -k1,1nr | cut -f 2- | tr "\t" "\n" |\
bioawk -c fastx '{ print ">Dint-" ++i  "\n" $seq}' > Dry-int.fasta

grep ">" Genome_Dry-int-chr.fasta | head -n 15
grep ">" Dry-int.fasta | head -n 15

#-------------------------
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}'  Genome_Dry-drumm-chr.fasta |\
awk -F '\t' '{printf("%d\t%s\n",length($2),$0);}' |\
sort -k1,1nr | cut -f 2- | tr "\t" "\n" |\
bioawk -c fastx '{ print ">Ddru-" ++i  "\n" $seq}' > Dry-dru.fasta

grep ">" Genome_Dry-drumm-chr.fasta | head -n 15
grep ">" Dry-dru.fasta | head -n 15

#-------------------------
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Genome_Dry-alask-chr.fasta |\
awk -F '\t' '{printf("%d\t%s\n",length($2),$0);}' |\
sort -k1,1nr | cut -f 2- | tr "\t" "\n" |\
bioawk -c fastx '{ print ">Dala-" ++i  "\n" $seq}' > Dry-ala.fasta

grep ">" Genome_Dry-alask-chr.fasta | head -n 15
grep ">" Dry-ala.fasta | head -n 15

# move genomes to main folder
cp Dry-*.fasta ..

#----------------------------
# Make a copy of the template Snakefile
cp Snakefile_Polygonum_aviculare Snakefile_Dryas_octo

nano Snakefile_Dryas_octo

#In snake file replace:
# Oxyria.fasta -> Oxyria_ncbi.fasta
# project Oxyria: -> Oxyria_NCBI
# change link to protein training file in Snakefile
# change name of Assembly spliter R script

# Sample specifics
CHRS = '1 2 3 4 5 6 7 8 9 10'.split()
PROJECT = "Dryoct" # should match chr splitter file name
REFERENCE = "Dry-octo-H0.fasta"
...
expand("{Project}/Maker_Files/csa.trans.Protein.Rosaceae.fa",Project=PROJECT),
...
Protein_File=expand("{Project}/Maker_Files/csa.trans.Protein.Rosaceae.fa",Project=PROJ
...
echo Assembly_Chr_splitter_Dryoct.R --args -f {input}
...
mv scaffolds.fasta scaffold_10.fasta

# check all cases replaced
grep csa.trans.Protein. Snakefile_Dryas_octo
grep Assembly_Chr_splitter_ Snakefile_Dryas_octo

#--------------------------------
# Note that all names MUST match

#Scaffold names set above (e.g. DoctH0-)
#PROJECT ID
#Chr splitter R script name

#---------------------------------
# to run
tmux new-session -s Annotate_Dryocto
tmux attach-session -t Annotate_Dryocto

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
#q()

cd /home/celphin/scratch/Annotation/CAP_Snakemake
module load StdEnv/2020 intel/2020.1.217 metagenome-atlas/2.5.0 
snakemake --snakefile Snakefile_Dryas_octo --profile cc-slurm --latency-wait 30


#------------------------
(one of the commands exited with non-zero exit code; note that snakemake uses bash strict mode!)
#try
#bash Name_checker_pre.sh scaffold_7.fasta DoctH0
set +u
source ~/miniconda2/bin/activate EDTA
set -u
echo starting EDTA process on: scaffold_7.fasta
EDTA.pl --overwrite 0 --genome scaffold_7.fasta --sensitive 0 --anno 1 --evaluate 0 --threads 48 --force 1
conda deactivate
set +u
# bash Name_checker_post.sh scaffold_7.fasta DoctH0 scaffold_7.fasta.mod.EDTA.TEanno.gff3
# needed to remove this line of name check to get working

#---------------------------
# Add fix_nucelotides option

rule MAKER3:
        input:
                EST_File=rules.Init.output.Protein_File,
                reference=rules.Chr_splitting.output,
                Protein_File=rules.Init.output.Protein_File,
                Repeats_File=rules.EDTA_individual.output.repeats_file,
        output:
                "{Project}/Maker_Files/Chr{Chrs}/scaffold_{Chrs}.maker.output/scaffold_{Chrs}_master_da
tastore_index.log",
        params:
                project=PROJECT,
        threads: 48
        resources:
                mem_mb=187000, #4GB
                time=4000 #minutes
        shell:
                """
                BASEDIR=$PWD
                module load StdEnv/2020 gcc/9.3.0 openmpi/4.0.3 postgresql/13.2 maker/3.01.03 augustus/
3.4.0
                cd {params.project}/Maker_Files
                echo Creating Maker Files Chr: {wildcards.Chrs}
                mkdir -p Chr{wildcards.Chrs}
                cd Chr{wildcards.Chrs}
                maker -CTL
                cd ..

                sed -i 's|RepeatMasker=/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Compil
er/gcc9/repeatmasker/4.1.1/RepeatMasker|RepeatMasker=/home/celphin/scratch/Annotation/RepeatMasker/Repe
atMasker|g'  Chr{wildcards.Chrs}/maker_exe.ctl

                cp maker_opts.template Chr{wildcards.Chrs}/maker_opts.ctl
                sed -i 's/model_org=all/model_org=/g'  Chr{wildcards.Chrs}/maker_opts.ctl
                echo "########################" >> Chr{wildcards.Chrs}/maker_opts.ctl
                echo "#-----Custom Parameters (these are always required)" >> Chr{wildcards.Chrs}/maker
_opts.ctl
                echo genome=$BASEDIR/{params.project}/Ref/scaffold_{wildcards.Chrs}.fasta >> Chr{wildca
rds.Chrs}/maker_opts.ctl
                echo "est=$BASEDIR/{input.EST_File}" >> Chr{wildcards.Chrs}/maker_opts.ctl
                echo "protein=$BASEDIR/{input.Protein_File}" >> Chr{wildcards.Chrs}/maker_opts.ctl
                echo "repeat_protein=$BASEDIR/{input.Repeats_File}" >> Chr{wildcards.Chrs}/maker_opts.c
tl
                echo "augustus_species=arabidopsis" >> Chr{wildcards.Chrs}/maker_opts.ctl

                cd Chr{wildcards.Chrs}

                mpiexec --use-hwthread-cpus maker -fix_nucleotides

                module unload postgresql/13.2
                module unload openmpi/4.0.3
                module unload maker/3.01.03
                module unload augustus/3.4.0
                """

################################
maker -fix_nucleotides


# STATUS: Parsing control files...
# STATUS: Processing and indexing input FASTA files...
# ERROR: The nucleotide sequence file '/lustre04/scratch/celphin/Annotation/CAP_Snakemake/DoctH0/Maker_Files/csa.trans.Protein.Rosaceae.fa'
# contains either protein sequence or unsupported characters. Note the
# following nucleotides may be valid but are unsupported [RYKMSWBDHV]. This
# message is to get you to look at your input files, and verify that there
# is not a mistake. Both an explanation of the cause and a solution are
# indicated in this message. Do not post it to the mailing list. Please
# manually check/fix the file before continuing, or set -fix_nucleotides on
# the command line to automatically replace invalid nucleotides with 'N'.
# Invalid Character: 'P'

####################################
# After finished running - launch BUSCO

cp /home/celphin/scratch/Annotation/CAP_Snakemake/DoctH0/FINAL_ANNOTATION/* /home/celphin/scratch/Annotation/BUSCO
cp /home/celphin/scratch/Annotation/CAP_Snakemake/DoctH0/Ref/DoctH0_Main.fasta /home/celphin/scratch/Annotation/BUSCO
cd /home/celphin/scratch/Annotation/BUSCO

module load StdEnv/2020
module load gcc/9.3.0 python augustus hmmer blast+ metaeuk prodigal r bbmap

sbatch gff3_to_cds.sh FINAL_DoctH0.AED_0.6.gff3 DoctH0_Main.fasta Dryas-AED0.6-genes.fasta
sbatch gff3_to_cds.sh FINAL_Kelp.AED_1.sorted.gff3 DoctH0_Main.fasta Dryas-AED1-genes.fasta

#-----------------
tmux new-session -s BUSCO
tmux attach-session -t BUSCO

salloc -c10 --time 2:55:00 --mem 120000m --account def-rieseber

cd /home/celphin/scratch/Annotation/BUSCO
source ~/busco_env/bin/activate
module load StdEnv/2020
module load gcc/9.3.0 python augustus hmmer blast+ metaeuk prodigal r bbmap

# Dryas genome
busco --offline --in DoctH0_Main.fasta \
--out  BUSCO_Dryas_Assembly_Augustus  --lineage_dataset eudicots_odb10 --mode genome --cpu 10 \
--download_path ~/BUSCO_downloads/

INFO:

        --------------------------------------------------
        |Results from dataset eudicots_odb10              |
        --------------------------------------------------
        |C:98.0%[S:95.2%,D:2.8%],F:0.5%,M:1.5%,n:2326     |
        |2278   Complete BUSCOs (C)                       |
        |2214   Complete and single-copy BUSCOs (S)       |
        |64     Complete and duplicated BUSCOs (D)        |
        |12     Fragmented BUSCOs (F)                     |
        |36     Missing BUSCOs (M)                        |
        |2326   Total BUSCO groups searched               |
        --------------------------------------------------

# Dryas annotations 
busco --offline --in Dryas-AED0.6-genes.fasta \
--out  BUSCO_Dryas_Annotation0.6_Augustus  --lineage_dataset eudicots_odb10 --mode genome --cpu 10 \
--download_path ~/BUSCO_downloads/ 

INFO:

        --------------------------------------------------
        |Results from dataset eudicots_odb10              |
        --------------------------------------------------
        |C:91.8%[S:88.8%,D:3.0%],F:3.8%,M:4.4%,n:2326     |
        |2135   Complete BUSCOs (C)                       |
        |2065   Complete and single-copy BUSCOs (S)       |
        |70     Complete and duplicated BUSCOs (D)        |
        |88     Fragmented BUSCOs (F)                     |
        |103    Missing BUSCOs (M)                        |
        |2326   Total BUSCO groups searched               |
        --------------------------------------------------


busco --offline --in Dryas-AED1-genes.fasta \
--out  BUSCO_Dryas_Annotation1_Augustus  --lineage_dataset eudicots_odb10 --mode genome --cpu 10 \
--download_path ~/BUSCO_downloads/ 

 