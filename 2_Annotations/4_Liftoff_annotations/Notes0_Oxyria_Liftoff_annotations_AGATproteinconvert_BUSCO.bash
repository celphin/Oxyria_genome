#############################
# Liftoff to move annotations over
# Oxyria
# June 2024
# https://github.com/agshumate/Liftoff
##################################
# try liftoff to move Oxyria annotation over to other genomes

#####################################
# Transfer to Beluga from Cedar

# Ellesmere Hap1 and Hap2
# Protein mode = 86.0% AED=0.6: 82.9%; 
cp *.fasta ../../Annotation/liftoff/
cd /home/celphin/scratch/Annotation/liftoff

rename Oxy_draft_assembly Oxy_Elles *
rename .FINAL. . *

#---------------------------
# Svalbard_Hap1
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_029168935.1/
# Protein mode = 85.4% AED=0.6: 82.3%:

# Svalabrd Hap2
# https://springernature.figshare.com/collections/Whole-genome_sequencing_of_13_Arctic_plants_and_draft_genomes_of_Oxyria_digyna_and_Cochlearia_groenlandica/6965802/1

rename od. Oxy_Sval_ *
#------------------------
#DToL
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_964035995.1/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/964/035/995/GCA_964035995.1_dcOxyDigy1.hap1.1/GCA_964035995.1_dcOxyDigy1.hap1.1_genomic.fna.gz
gunzip GCA_964035995.1_dcOxyDigy1.hap1.1_genomic.fna.gz
mv GCA_964035995.1_dcOxyDigy1.hap1.1_genomic.fna DToL_h1.fasta

# https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_964036005.1/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/964/036/005/GCA_964036005.1_dcOxyDigy1.hap2.1/GCA_964036005.1_dcOxyDigy1.hap2.1_genomic.fna.gz
gunzip GCA_964036005.1_dcOxyDigy1.hap2.1_genomic.fna.gz
mv GCA_964036005.1_dcOxyDigy1.hap2.1_genomic.fna DToL_h2.fasta


########################################
# Run Liftoff to transfer annotations

# Liftoff
# https://github.com/agshumate/Liftoff
# cd ~ 
# ~/miniconda2/bin/conda create --name liftoff  
# source ~/miniconda2/bin/activate liftoff
# module load StdEnv/2020
# module load python/3.10.2
# module load scipy-stack/2022a
# module load gcc/9.3.0
# module load cuda/11.7
# module load hdf5/1.12.1
# module load parasail/2.5
# module load minimap2/2.24
# git clone https://github.com/agshumate/Liftoff liftoff 
# cd /home/celphin/scratch/Annotation/liftoff
# pip install .

# conda deactivate

#-------------------
# to avoid conda 
# module load gcc python/3.10 minimap2 parasail
# virtualenv --no-download --clear ~/ENV && source ~/ENV/bin/activate
# pip install git+https://github.com/agshumate/Liftoff.git


#------------------
# To run

# usage: liftoff [-h] (-g GFF | -db DB) [-o FILE] [-u FILE] [-exclude_partial]
               # [-dir DIR] [-mm2_options =STR] [-a A] [-s S] [-d D] [-flank F]
               # [-V] [-p P] [-m PATH] [-f TYPES] [-infer_genes]
               # [-infer_transcripts] [-chroms TXT] [-unplaced TXT] [-copies]
               # [-sc SC] [-overlap O] [-mismatch M] [-gap_open GO]
               # [-gap_extend GE]
               # target reference

 # -overlap O          maximum fraction [0.0-1.0] of overlap allowed by 2
                      # features; by default O=0.1
  # -mismatch M         mismatch penalty in exons when finding best mapping; by
                      # default M=2
  # -gap_open GO        gap open penalty in exons when finding best mapping; by
                      # default GO=2
  # -gap_extend GE      gap extend penalty in exons when finding best mapping;
                      # by default GE=1
#-----------------------

tmux new-session -s Liftoff
tmux attach-session -t Liftoff

salloc -c40 --time 4:55:00 --mem 120000m --account def-rieseber

nano run_liftoff.sh
#!/bin/bash
#SBATCH --account=def-rieseber
#SBATCH --job-name=liftoff
#SBATCH --time=0-10:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=95000m

# conda
source ~/miniconda2/bin/activate liftoff
module load StdEnv/2020
module load python/3.10.2
module load scipy-stack/2022a
module load gcc/9.3.0
module load cuda/11.7
module load hdf5/1.12.1
module load parasail/2.5
module load minimap2/2.24

cd /home/celphin/scratch/Oxyria/GeneSpace/Oxyria_genomes/liftoff

#------------------------------
SPP_Hap=$1

# or no conda
#module load gcc python/3.10 minimap2 parasail
#source ~/ENV/bin/activate
#python -c "import liftoff; print(liftoff.__version__)"

liftoff \
-g ./data/FINAL_Oxyria.AED_0.6.sorted.gff3 -p 40 -o ./output/${SPP_Hap}_liftoffpolish.fasta \
-dir ${SPP_Hap}_intermed -u ${SPP_Hap}_unmapped_features.txt \
-infer_genes -copies -a 0.95 -s 0.95 -d 5.0 -flank 0.8 -polish \
./data/${SPP_Hap}.fasta ./data/Oxyria_Main.fasta

sbatch run_liftoff.sh DToL_h1 # done
sbatch run_liftoff.sh Oxy_Elles_Hap1
sbatch run_liftoff.sh Oxy_Elles_Hap2 # done
sbatch run_liftoff.sh Oxy_Sval_h1 # done
sbatch run_liftoff.sh Oxy_Sval_h2 # done
sbatch run_liftoff.sh DToL_h2 # done

#-----------------------
SPP_Hap=DToL_h1
liftoff \
-g FINAL_Oxyria.AED_0.6.sorted.gff3 -p 40 -o ${SPP_Hap}_liftoffpolish1.gff3 \
-dir ${SPP_Hap}_intermed -u ${SPP_Hap}_unmapped_features.txt \
-infer_genes -copies -a 0.95 -s 0.95 -d 5.0 -flank 0.8 -polish \
${SPP_Hap}.fasta Oxyria_Main.fasta
# done - each about 4 hours

SPP_Hap=Oxy_Elles_Hap1
liftoff \
-g FINAL_Oxyria.AED_0.6.sorted.gff3 -p 40 -o ${SPP_Hap}_liftoffpolish1.gff3 \
-dir ${SPP_Hap}_intermed -u ${SPP_Hap}_unmapped_features.txt \
-infer_genes -copies -a 0.95 -s 0.95 -d 5.0 -flank 0.8 -polish \
${SPP_Hap}.fasta Oxyria_Main.fasta

SPP_Hap=Oxy_Elles_Hap2 
liftoff \
-g FINAL_Oxyria.AED_0.6.sorted.gff3 -p 40 -o ${SPP_Hap}_liftoffpolish1.gff3 \
-dir ${SPP_Hap}_intermed -u ${SPP_Hap}_unmapped_features.txt \
-infer_genes -copies -a 0.95 -s 0.95 -d 5.0 -flank 0.8 -polish \
${SPP_Hap}.fasta Oxyria_Main.fasta

SPP_Hap=Oxy_Sval_h1 
liftoff \
-g FINAL_Oxyria.AED_0.6.sorted.gff3 -p 40 -o ${SPP_Hap}_liftoffpolish1.gff3 \
-dir ${SPP_Hap}_intermed -u ${SPP_Hap}_unmapped_features.txt \
-infer_genes -copies -a 0.95 -s 0.95 -d 5.0 -flank 0.8 -polish \
${SPP_Hap}.fasta Oxyria_Main.fasta

SPP_Hap=Oxy_Sval_h2 
liftoff \
-g FINAL_Oxyria.AED_0.6.sorted.gff3 -p 40 -o ${SPP_Hap}_liftoffpolish1.gff3 \
-dir ${SPP_Hap}_intermed -u ${SPP_Hap}_unmapped_features.txt \
-infer_genes -copies -a 0.95 -s 0.95 -d 5.0 -flank 0.8 -polish \
${SPP_Hap}.fasta Oxyria_Main.fasta

SPP_Hap=DToL_h2
liftoff \
-g FINAL_Oxyria.AED_0.6.sorted.gff3 -p 40 -o ${SPP_Hap}_liftoffpolish1.gff3 \
-dir ${SPP_Hap}_intermed -u ${SPP_Hap}_unmapped_features.txt \
-infer_genes -copies -a 0.95 -s 0.95 -d 5.0 -flank 0.8 -polish \
${SPP_Hap}.fasta Oxyria_Main.fasta

##########################################################
# try Liftoff tools
# https://github.com/agshumate/LiftoffTools

tmux new-session -s liftoff
tmux attach-session -t liftoff

cd /home/celphin/scratch/Oxyria/GeneSpace/Oxyria_genomes/liftoff

# cp /home/celphin/scratch/Annotation/liftoff/data/Oxyria/*.gff3 .
# cp /home/celphin/scratch/Annotation/liftoff/data/Oxyria/*.fasta .
# cp /home/celphin/scratch/Annotation/liftoff/output/Oxy*.fasta_polished .
# cp /home/celphin/scratch/Annotation/liftoff/output/DToL*.fasta_polished .

rename fasta_polished gff3 *

#-----------------
cd /home/celphin/scratch/Oxyria/GeneSpace/Oxyria_genomes/liftoff

source ~/miniconda2/bin/activate liftoff
module load StdEnv/2020
module load python/3.10.2
module load scipy-stack/2022a
module load gcc/9.3.0
module load cuda/11.7
module load hdf5/1.12.1
module load parasail/2.5
module load minimap2/2.24
export PATH=$(pwd)/mmseqs/bin/:$PATH

# Install
# pip install liftofftools
# done

# Also need
# https://github.com/soedinglab/MMseqs2
# static build with AVX2 (fastest)
#wget https://mmseqs.com/latest/mmseqs-linux-avx2.tar.gz; tar xvfz mmseqs-linux-avx2.tar.gz; export PATH=$(pwd)/mmseqs/bin/:$PATH


# Run
liftofftools all -r Oxyria_Main.fasta -t DToL_h1.fasta \
-rg FINAL_Oxyria.AED_0.6.sorted.gff3 \
-tg DToL_h1_liftoffpolish1.gff3 \
-dir DToL_h1_liftofftools
# all failed
# gffutils database build failed with UNIQUE constraint failed: features.id

liftofftools all -r Oxyria_Main.fasta -t Oxy_Sval_h1.fasta \
-rg FINAL_Oxyria.AED_0.6.sorted.gff3 \
-tg Oxy_Sval_h1_liftoffpolish.gff3 \
-dir Oxy_Sval_h1_liftofftools

liftofftools all -r Oxyria_Main.fasta -t Oxy_Sval_h2.fasta \
-rg FINAL_Oxyria.AED_0.6.sorted.gff3 \
-tg Oxy_Sval_h2_liftoffpolish.gff3 \
-dir Oxy_Sval_h2_liftofftools

liftofftools all -r Oxyria_Main.fasta -t DToL_h2.fasta \
-rg FINAL_Oxyria.AED_0.6.sorted.gff3 \
-tg DToL_h2_liftoffpolish.gff3  \
-dir DToL_h2_liftofftools

#---------------------------
# Try just variants detection

tmux new-session -s liftoff
tmux attach-session -t liftoff

# parasail_memalign: posix_memalign failed: Cannot allocate memory
salloc -c1 --time 2:55:00 --mem 120000m --account def-rieseber

cd /home/celphin/scratch/Oxyria/GeneSpace/Oxyria_genomes/liftoff/
source ~/miniconda2/bin/activate liftoff
module load StdEnv/2020
module load python/3.10.2
module load scipy-stack/2022a
module load gcc/9.3.0
module load cuda/11.7
module load hdf5/1.12.1
module load parasail/2.5
module load minimap2/2.24
export PATH=$(pwd)/mmseqs/bin/:$PATH

liftofftools variants -r Oxyria_Main.fasta -t DToL_h1.fasta \
-rg FINAL_Oxyria.AED_0.6.sorted.gff3 \
-tg DToL_h1_liftoffpolish1.gff3_polished \
-dir DToL_h1_liftofftools_VARIANTS
# gffutils database build failed with UNIQUE constraint failed: features.id


liftofftools variants -r Oxyria_Main.fasta -t Oxy_Sval_h1.fasta \
-rg FINAL_Oxyria.AED_0.6.sorted.gff3 \
-tg Oxy_Sval_h1_liftoffpolish.gff3 \
-dir Oxy_Sval_h1_liftofftools_VARIANTS
#gffutils database build failed with UNIQUE constraint failed: features.id 

liftofftools variants -r Oxyria_Main.fasta -t Oxy_Sval_h2.fasta \
-rg FINAL_Oxyria.AED_0.6.sorted.gff3 \
-tg Oxy_Sval_h2_liftoffpolish.gff3 \
-dir Oxy_Sval_h2_liftofftools_VARIANTS

liftofftools variants -r Oxyria_Main.fasta -t DToL_h2.fasta \
-rg FINAL_Oxyria.AED_0.6.sorted.gff3 \
-tg DToL_h2_liftoffpolish.gff3 \
-dir DToL_h2_liftofftools_VARIANTS


#####################################
cd /home/celphin/scratch/Oxyria/GeneSpace/Oxyria_genomes/liftoff/Oxy_Sval_h1_liftofftools_VARIANTS/
Spp_hap=Oxy_Sval_h1

cd /home/celphin/scratch/Oxyria/GeneSpace/Oxyria_genomes/liftoff/Oxy_Sval_h2_liftofftools_VARIANTS/
Spp_hap=Oxy_Sval_h2

cd /home/celphin/scratch/Oxyria/GeneSpace/Oxyria_genomes/liftoff/DToL_h1_liftofftools_VARIANTS/
Spp_hap=DToL_h1

cd /home/celphin/scratch/Oxyria/GeneSpace/Oxyria_genomes/liftoff/DToL_h2_liftofftools_VARIANTS/
Spp_hap=DToL_h2


grep "frameshift"  variant_effects >${Spp_hap}_frameshift.txt
wc -l ${Spp_hap}_frameshift.txt
#
grep "inframe_deletion"  variant_effects >${Spp_hap}_inframe_deletion.txt
wc -l ${Spp_hap}_inframe_deletion.txt
# 

grep "start_lost"  variant_effects >${Spp_hap}_start_lost.txt
wc -l ${Spp_hap}_start_lost.txt
#

grep "3'_truncated" variant_effects >${Spp_hap}_3_truncated.txt
wc -l ${Spp_hap}_3_truncated.txt
#

grep "identical" variant_effects >${Spp_hap}_identical.txt
wc -l ${Spp_hap}_identical.txt
#

grep "nonsynonymous"  variant_effects >${Spp_hap}_nonsynonymous.txt
wc -l ${Spp_hap}_nonsynonymous.txt
#

grep "NA"  variant_effects >${Spp_hap}_NA.txt
wc -l ${Spp_hap}_NA.txt
#

grep -v "nonsynonymous"  variant_effects | grep "synonymous"  >${Spp_hap}_synonymous.txt
wc -l ${Spp_hap}_synonymous.txt
#

grep "stop_gained" variant_effects >${Spp_hap}_stop_gained.txt
wc -l ${Spp_hap}_stop_gained.txt
# 

grep "5'_truncated" variant_effects >${Spp_hap}_5_truncated.txt
wc -l ${Spp_hap}_5_truncated.txt


grep "inframe_insertion"  variant_effects >${Spp_hap}_inframe_insertion.txt
wc -l ${Spp_hap}_inframe_insertion.txt
# 

# https://en.wikipedia.org/wiki/Ka/Ks_ratio
# dN/dS=6450/7492 = 0.86 # dominant purifying or stabilizing selection


##############################
# try https://github.com/gpertea/gffread
# to get to protein conversion

# Install
 # cd /home/celphin/scratch/Annotation
  # git clone https://github.com/gpertea/gffread
  # cd gffread
  # make release

cd /home/celphin/scratch/Oxyria/GeneSpace/Oxyria_genomes/liftoff/

/home/celphin/scratch/Annotation/gffread/gffread -h

# to run proteins
/home/celphin/scratch/Annotation/gffread/gffread -g DToL_h1.fasta DToL_h1.gff3 -y DToL_h1_pep.fasta
/home/celphin/scratch/Annotation/gffread/gffread -g DToL_h2.fasta DToL_h2.gff3 -y DToL_h2_pep.fasta
/home/celphin/scratch/Annotation/gffread/gffread -g Oxy_Sval_h1.fasta Oxy_Sval_h1.gff3 -y Oxy_Sval_h1_pep.fasta
/home/celphin/scratch/Annotation/gffread/gffread -g Oxy_Sval_h2.fasta Oxy_Sval_h2.gff3 -y Oxy_Sval_h2_pep.fasta
/home/celphin/scratch/Annotation/gffread/gffread -g Oxy_Elles_Hap1.fasta Oxy_Elles_Hap1.gff3 -y Oxy_Elles_Hap1_pep.fasta
/home/celphin/scratch/Annotation/gffread/gffread -g Oxy_Elles_Hap2.fasta Oxy_Elles_Hap2.gff3 -y Oxy_Elles_Hap2_pep.fasta

#-------------------------
# to run cds
/home/celphin/scratch/Annotation/gffread/gffread -g DToL_h1.fasta DToL_h1.gff3 -x DToL_h1_cds.fa
/home/celphin/scratch/Annotation/gffread/gffread -g DToL_h2.fasta DToL_h2.gff3 -x DToL_h2_cds.fa
/home/celphin/scratch/Annotation/gffread/gffread -g Oxy_Sval_h1.fasta Oxy_Sval_h1.gff3 -x Oxy_Sval_h1_cds.fa
/home/celphin/scratch/Annotation/gffread/gffread -g Oxy_Sval_h2.fasta Oxy_Sval_h2.gff3 -x Oxy_Sval_h2_cds.fa
/home/celphin/scratch/Annotation/gffread/gffread -g Oxy_Elles_Hap1.fasta Oxy_Elles_Hap1.gff3 -x Oxy_Elles_Hap1_cds.fa
/home/celphin/scratch/Annotation/gffread/gffread -g Oxy_Elles_Hap2.fasta Oxy_Elles_Hap2.gff3 -x Oxy_Elles_Hap2_cds.fa


##########################
# check BUSCO scores

tmux new-session -s BUSCO
tmux attach-session -t BUSCO

cd /home/celphin/scratch/Oxyria/GeneSpace/Oxyria_genomes/liftoff

salloc -c10 --time 2:55:00 --mem 120000m --account def-rieseber

module load StdEnv/2020
module load gcc/9.3.0 python augustus hmmer blast+ metaeuk prodigal r bbmap
source ~/busco_env/bin/activate

# Oxyria genomes
busco --offline --in Oxy_Elles_Hap1.fasta \
--out  BUSCO_Oxy_Elles_Hap1_Assembly  --lineage_dataset eudicots_odb10 --mode genome --cpu 10 \
--download_path ~/BUSCO_downloads/
# C:92.7%[S:86.0%,D:6.7%],F:2.0%,M:5.3%,n:2326

busco --offline --in Oxy_Elles_Hap1_pep.fasta \
--out  BUSCO_Oxy_Elles_Hap1_Ann  --lineage_dataset eudicots_odb10 --mode protein --cpu 10 \
--download_path ~/BUSCO_downloads/
# C:83.1%[S:77.1%,D:6.0%],F:2.4%,M:14.5%,n:2326

#---------------------------
busco --offline --in Oxy_Elles_Hap2.fasta \
--out  BUSCO_Oxy_Elles_Hap2_Assembly  --lineage_dataset eudicots_odb10 --mode genome --cpu 10 \
--download_path ~/BUSCO_downloads/
# C:93.4%[S:86.3%,D:7.1%],F:1.6%,M:5.0%,n:2326

busco --offline --in Oxy_Elles_Hap2_pep.fasta \
--out  BUSCO_Oxy_Elles_Hap2_Ann  --lineage_dataset eudicots_odb10 --mode protein --cpu 10 \
--download_path ~/BUSCO_downloads/
# C:83.7%[S:77.5%,D:6.2%],F:2.3%,M:14.0%,n:2326

#---------------------------
busco --offline --in Oxy_Sval_h1.fasta \
--out  BUSCO_Oxy_Sval_h1_Assembly  --lineage_dataset eudicots_odb10 --mode genome --cpu 10 \
--download_path ~/BUSCO_downloads/
# C:93.6%[S:86.4%,D:7.2%],F:1.5%,M:4.9%,n:2326

busco --offline --in Oxy_Sval_h1_pep.fasta \
--out  BUSCO_Oxy_Sval_h1_Ann  --lineage_dataset eudicots_odb10 --mode protein --cpu 10 \
--download_path ~/BUSCO_downloads/
# C:4.3%[S:4.1%,D:0.2%],F:1.8%,M:93.9%,n:2326

#---------------------------
busco --offline --in Oxy_Sval_h2.fasta \
--out  BUSCO_Oxy_Sval_h2_Assembly  --lineage_dataset eudicots_odb10 --mode genome --cpu 10 \
--download_path ~/BUSCO_downloads/
# C:93.5%[S:86.6%,D:6.9%],F:1.5%,M:5.0%,n:2326

busco --offline --in Oxy_Sval_h2_pep.fasta \
--out  BUSCO_Oxy_Sval_h2_Ann  --lineage_dataset eudicots_odb10 --mode protein --cpu 10 \
--download_path ~/BUSCO_downloads/
# C:0.9%[S:0.9%,D:0.0%],F:0.2%,M:98.9%,n:2326

#---------------------------
busco --offline --in DToL_h1.fasta \
--out  BUSCO_DToL_h1_Assembly  --lineage_dataset eudicots_odb10 --mode genome --cpu 10 \
--download_path ~/BUSCO_downloads/
# C:93.6%[S:86.5%,D:7.1%],F:1.7%,M:4.7%,n:2326

busco --offline --in DToL_h1_pep.fasta \
--out  BUSCO_DToL_h1_Ann  --lineage_dataset eudicots_odb10 --mode protein --cpu 10 \
--download_path ~/BUSCO_downloads/
# C:81.7%[S:76.4%,D:5.3%],F:2.4%,M:15.9%,n:2326

#---------------------------
busco --offline --in DToL_h2.fasta \
--out  BUSCO_DToL_h2_Assembly  --lineage_dataset eudicots_odb10 --mode genome --cpu 10 \
--download_path ~/BUSCO_downloads/
# C:93.3%[S:86.4%,D:6.9%],F:1.6%,M:5.1%,n:2326

busco --offline --in DToL_h2_pep.fasta \
--out  BUSCO_DToL_h2_Ann  --lineage_dataset eudicots_odb10 --mode protein --cpu 10 \
--download_path ~/BUSCO_downloads/
# C:81.3%[S:76.1%,D:5.2%],F:2.4%,M:16.3%,n:2326

#---------------------------
busco --offline --in Oxyria_Main.fasta \
--out  BUSCO_Oxyria_Main_Assembly  --lineage_dataset eudicots_odb10 --mode genome --cpu 10 \
--download_path ~/BUSCO_downloads/
# C:93.6%[S:85.9%,D:7.7%],F:1.5%,M:4.9%,n:2326

busco --offline --in Oxyria_Main_pep.fasta \
--out  BUSCO_Oxyria_Main_Ann  --lineage_dataset eudicots_odb10 --mode protein --cpu 10 \
--download_path ~/BUSCO_downloads/
# C:3.2%[S:3.1%,D:0.1%],F:0.7%,M:96.1%,n:2326















































############################
# OLD - AGAT error
######################################
# Run Genespace/OrthoFinder to compare genomes

mkdir /home/celphin/scratch/Oxyria/GeneSpace/Oxyria_genomes
cd /home/celphin/scratch/Oxyria/GeneSpace/Oxyria_genomes

#----------------------------
# Oxyrias

# copy over genome and gff3
cp /home/celphin/scratch/Annotation/liftoff/output/*.fasta_polished .
cp /home/celphin/scratch/Annotation/liftoff/data/*.fasta .
cp /home/celphin/scratch/Annotation/liftoff/data/FINAL_Oxyria.AED_0.6.sorted.gff3 .
rename _liftoffpolish.fasta_polished .gff3 *
mv FINAL_Oxyria.AED_0.6.sorted.gff3 Oxyria_Main.gff3

#-------------------
# need protein files

cd /home/celphin/scratch/Oxyria/GeneSpace/Oxyria_genomes

#---------------------
module load StdEnv/2020 python/3.10.2
source ~/gff3_env/bin/activate

gff3_to_fasta -g  DToL_h1.sorted.gff3 -f DToL_h1.fasta -st pep -d simple -o DToL_h1
# ERROR    Line 394117: FORMAT: Duplicate ID: "Oxyria_Chr600004222-RA:cds" in non-adjacent lines: 394109,394110,394113,394114
# -> OZ038365.1   Liftoff CDS     8886022 8886105 .       +       2       ID=Oxyria_Chr600004222-RA:cds;Alias=character(0);Parent=Oxyria_Chr600004222-RA;extra_copy_number=0
# ERROR    Line 394122: FORMAT: Duplicate ID: "Oxyria_Chr600004222-RA:cds" in non-adjacent lines: 394109,394110,394113,394114,394117,394118
# -> OZ038365.1   Liftoff CDS     8886334 8886698 .       +       2       ID=Oxyria_Chr600004222-RA:cds;Alias=character(0);Parent=Oxyria_Chr600004222-RA;extra_copy_number=0
# ERROR    Line 394124: FORMAT: Duplicate ID: "Oxyria_Chr600004221" in non-adjacent lines: 394119
# -> OZ038365.1   Liftoff gene    8887414 8888824 .       -       .       ID=Oxyria_Chr600004221;Name=Oxyria_Chr600004221;Alias=augustus_masked-Oxy-6-73289246-processed-gene-642.15;Parent=character(0);coverage=1.0;sequence_ID=0.999;valid_ORFs=0;extra_copy_number=0;copy_num_ID=Oxyria_Chr600004221_0

deactivate

#------------------------
# reformat with AGAT??
# https://agat.readthedocs.io/en/latest/tools/agat_sp_fix_features_locations_duplicated.html

# Load modules on graham or cedar (or use instructions for conda at https://github.com/NBISweden/AGAT if using your own computer)
module load StdEnv/2023
module load apptainer/1.2.4

# Install `AGAT` via `Apptainer` into a temp directory
apptainer pull docker://quay.io/biocontainers/agat:1.4.0--pl5321hdfd78af_0
apptainer run agat_1.4.0--pl5321hdfd78af_0.sif

# Kick the tires, was it successfully installed?
agat_convert_sp_gxf2gxf.pl --help | less # It should show tool specific information. Use 'q' to return the prompt
agat_<tab>  # tab-completion should result in a list of all the possible scripts; hit space to advance to next page, or 'q' to return the prompt
agat_convert_sp_gxf2gxf.pl --help

agat_sp_extract_sequences.pl -g DToL_h1.gff3 -f DToL_h1.fasta -t cds -p -o DToL_h1_pep.fasta
agat_sp_extract_sequences.pl -g DToL_h2.gff3 -f DToL_h2.fasta -t cds -p -o DToL_h2_pep.fasta

agat_sp_extract_sequences.pl -g Oxy_Elles_Hap1.gff3 -f Oxy_Elles_Hap1.fasta -t cds -p -o Oxy_Elles_Hap1_pep.fasta
agat_sp_extract_sequences.pl -g Oxy_Elles_Hap2.gff3 -f Oxy_Elles_Hap2.fasta -t cds -p -o Oxy_Elles_Hap2_pep.fasta

agat_sp_extract_sequences.pl -g Oxy_Sval_h1.gff3 -f Oxy_Sval_h1.fasta -t cds -p -o Oxy_Sval_h1_pep.fasta
agat_sp_extract_sequences.pl -g Oxy_Sval_h2.gff3 -f Oxy_Sval_h2.fasta -t cds -p -o Oxy_Sval_h2_pep.fasta

# WARNING: Problem ! The size of the extracted sequence 374 is different than the specified span: 373.
# That often occurs when the fasta file does not correspond to the annotation file. Or the index file comes from another fasta file which had the same name and haven't been removed.
# As last possibility your gff contains location errors (Already encountered for a Maker annotation)
# Supplement information: seq_id=h1tg000020l ; seq_id_correct=h1tg000020l ; start=5914495 ; end=5914867 ; h1tg000020l sequence length: 6049140 )

agat_sp_extract_sequences.pl -g Oxyria_Main.gff3 -f Oxyria_Main.fasta -t cds -p -o Oxyria_Main_pep.fasta

exit


####################################
# Run BUSCO

tmux new-session -s BUSCO
tmux attach-session -t BUSCO

cd /home/celphin/scratch/Oxyria/GeneSpace/Oxyria_genomes/genomes

salloc -c10 --time 2:55:00 --mem 120000m --account def-rieseber

module load StdEnv/2020
module load gcc/9.3.0 python augustus hmmer blast+ metaeuk prodigal r bbmap
source ~/busco_env/bin/activate

# Oxyria genomes

busco --offline --in Oxy_Elles_Hap1.fasta \
--out  BUSCO_Oxy_Elles_Hap1_Assembly  --lineage_dataset eudicots_odb10 --mode genome --cpu 10 \
--download_path ~/BUSCO_downloads/
# C:92.7%[S:86.0%,D:6.7%],F:2.0%,M:5.3%,n:2326

busco --offline --in Oxy_Elles_Hap1_pep.fasta \
--out  BUSCO_Oxy_Elles_Hap1_Ann  --lineage_dataset eudicots_odb10 --mode protein --cpu 10 \
--download_path ~/BUSCO_downloads/
# C:83.1%[S:77.1%,D:6.0%],F:2.4%,M:14.5%,n:2326

#---------------------------
busco --offline --in Oxy_Elles_Hap2.fasta \
--out  BUSCO_Oxy_Elles_Hap2_Assembly  --lineage_dataset eudicots_odb10 --mode genome --cpu 10 \
--download_path ~/BUSCO_downloads/
# C:93.4%[S:86.3%,D:7.1%],F:1.6%,M:5.0%,n:2326

busco --offline --in Oxy_Elles_Hap2_pep.fasta \
--out  BUSCO_Oxy_Elles_Hap2_Ann  --lineage_dataset eudicots_odb10 --mode protein --cpu 10 \
--download_path ~/BUSCO_downloads/
# C:83.7%[S:77.5%,D:6.2%],F:2.3%,M:14.0%,n:2326

#---------------------------
busco --offline --in Oxy_Sval_h1.fasta \
--out  BUSCO_Oxy_Sval_h1_Assembly  --lineage_dataset eudicots_odb10 --mode genome --cpu 10 \
--download_path ~/BUSCO_downloads/
# C:93.6%[S:86.4%,D:7.2%],F:1.5%,M:4.9%,n:2326

busco --offline --in Oxy_Sval_h1_pep.fasta \
--out  BUSCO_Oxy_Sval_h1_Ann  --lineage_dataset eudicots_odb10 --mode protein --cpu 10 \
--download_path ~/BUSCO_downloads/
# C:4.3%[S:4.1%,D:0.2%],F:1.8%,M:93.9%,n:2326

#---------------------------
busco --offline --in Oxy_Sval_h2.fasta \
--out  BUSCO_Oxy_Sval_h2_Assembly  --lineage_dataset eudicots_odb10 --mode genome --cpu 10 \
--download_path ~/BUSCO_downloads/
# C:93.5%[S:86.6%,D:6.9%],F:1.5%,M:5.0%,n:2326

busco --offline --in Oxy_Sval_h2_pep.fasta \
--out  BUSCO_Oxy_Sval_h2_Ann  --lineage_dataset eudicots_odb10 --mode protein --cpu 10 \
--download_path ~/BUSCO_downloads/
# C:0.9%[S:0.9%,D:0.0%],F:0.2%,M:98.9%,n:2326

#---------------------------
busco --offline --in DToL_h1.fasta \
--out  BUSCO_DToL_h1_Assembly  --lineage_dataset eudicots_odb10 --mode genome --cpu 10 \
--download_path ~/BUSCO_downloads/
# C:93.6%[S:86.5%,D:7.1%],F:1.7%,M:4.7%,n:2326

busco --offline --in DToL_h1_pep.fasta \
--out  BUSCO_DToL_h1_Ann  --lineage_dataset eudicots_odb10 --mode protein --cpu 10 \
--download_path ~/BUSCO_downloads/
# C:81.7%[S:76.4%,D:5.3%],F:2.4%,M:15.9%,n:2326

#---------------------------
busco --offline --in DToL_h2.fasta \
--out  BUSCO_DToL_h2_Assembly  --lineage_dataset eudicots_odb10 --mode genome --cpu 10 \
--download_path ~/BUSCO_downloads/
# C:93.3%[S:86.4%,D:6.9%],F:1.6%,M:5.1%,n:2326

busco --offline --in DToL_h2_pep.fasta \
--out  BUSCO_DToL_h2_Ann  --lineage_dataset eudicots_odb10 --mode protein --cpu 10 \
--download_path ~/BUSCO_downloads/
# C:81.3%[S:76.1%,D:5.2%],F:2.4%,M:16.3%,n:2326

#---------------------------
busco --offline --in Oxyria_Main.fasta \
--out  BUSCO_Oxyria_Main_Assembly  --lineage_dataset eudicots_odb10 --mode genome --cpu 10 \
--download_path ~/BUSCO_downloads/
# C:93.6%[S:85.9%,D:7.7%],F:1.5%,M:4.9%,n:2326

busco --offline --in Oxyria_Main_pep.fasta \
--out  BUSCO_Oxyria_Main_Ann  --lineage_dataset eudicots_odb10 --mode protein --cpu 10 \
--download_path ~/BUSCO_downloads/
# C:3.2%[S:3.1%,D:0.1%],F:0.7%,M:96.1%,n:2326

#---------------------------
# Try cds for liftoff annotation?
cd /home/celphin/scratch/Oxyria/GeneSpace/Oxyria_genomes/genomes

sbatch gff3_to_cds.sh Oxyria_Main.gff3 Oxyria_Main.fasta Oxyria_Main_cds.fasta

busco --offline --in Oxyria_Main_cds.fasta \
--out  BUSCO_Oxyria_Main_Ann_cds  --lineage_dataset eudicots_odb10 --mode genome --cpu 10 \
--download_path ~/BUSCO_downloads/
# C:83.3%[S:76.6%,D:6.7%],F:3.1%,M:13.6%,n:2326

sbatch gff3_to_cds.sh Oxy_Sval_h1.gff3 Oxy_Sval_h1.fasta Oxy_Sval_h1_cds.fasta

busco --offline --in Oxy_Sval_h1_cds.fasta \
--out  BUSCO_Oxy_Sval_h1_Ann_cds  --lineage_dataset eudicots_odb10 --mode genome --cpu 10 \
--download_path ~/BUSCO_downloads/
# C:87.1%[S:14.7%,D:72.4%],F:2.5%,M:10.4%,n:2326

sbatch gff3_to_cds.sh Oxy_Sval_h2.gff3 Oxy_Sval_h2.fasta Oxy_Sval_h2_cds.fasta

busco --offline --in Oxy_Sval_h2_cds.fasta \
--out  BUSCO_Oxy_Sval_h2_Ann_cds  --lineage_dataset eudicots_odb10 --mode genome --cpu 10 \
--download_path ~/BUSCO_downloads/
# C:87.2%[S:14.9%,D:72.3%],F:2.5%,M:10.3%,n:2326

#####################################
# make folder structure

mkdir DToL_h1
cp DToL_h1* DToL_h1/

mkdir DToL_h2
cp DToL_h2* DToL_h2/

mkdir Oxy_Elles_Hap1
cp Oxy_Elles_Hap1* Oxy_Elles_Hap1/

mkdir Oxy_Elles_Hap2
cp Oxy_Elles_Hap2* Oxy_Elles_Hap2/

mkdir Oxy_Sval_h1
cp Oxy_Sval_h1* Oxy_Sval_h1/

mkdir Oxy_Sval_h2
cp Oxy_Sval_h2* Oxy_Sval_h2/

mkdir Oxyria_Main
cp Oxyria_Main* Oxyria_Main/

rm ./*/*.fasta.index
rm ./*/*.agat.log
rm ./*/*1.fasta
rm ./*/*2.fasta

####################################
# run parse annotations
# https://rdrr.io/github/jtlovell/GENESPACE/man/parse_annotations.html

tmux new-session -s GeneSpace1
tmux attach-session -t GeneSpace1

module load StdEnv/2020 r/4.2.2 glpk/5.0

R

library(GENESPACE)

parsing_files_other <- function(SPP_Hap, GeneID, Wd){
  parsedPaths_other <- parse_annotations(
    rawGenomeRepo = "/home/celphin/scratch/Oxyria/GeneSpace/Oxyria_genomes", 
    genomeDirs = SPP_Hap,
    genomeIDs = SPP_Hap,
    gffString = "gff",
    faString = "fa",
    headerEntryIndex = 1, 
    gffIdColumn = GeneID,
    genespaceWd = "/home/celphin/scratch/Oxyria/GeneSpace/Oxyria_genomes")
}

parse_ncbi_other <- function(SPP_Hap, gffID, gffStrip){
  parse_ncbi(
    rawGenomeRepo="/home/celphin/scratch/Oxyria/GeneSpace/Oxyria_genomes",
    genomeDirs = SPP_Hap,
    genomeIDs = SPP_Hap,
    gffString = "gff",
    faString = "fa",
    genespaceWd="/home/celphin/scratch/Oxyria/GeneSpace/Oxyria_genomes",
    troubleShoot = FALSE
  )
}

parsing_files_other("DToL_h1", "ID")
# DToL_h1: n unique sequences = 31522, n matched to gff = 31522

parsing_files_other("DToL_h2", "ID")
# DToL_h2: n unique sequences = 35476, n matched to gff = 35476

parsing_files_other("Oxy_Elles_Hap2", "ID")
# Oxy_Elles_Hap2: n unique sequences = 35196, n matched to gff = 35196

parsing_files_other("Oxy_Sval_h1", "ID")
# Oxy_Sval_h1: n unique sequences = 33401, n matched to gff = 33401

parsing_files_other("Oxyria_Main", "ID")
# Oxyria_Main: n unique sequences = 37871, n matched to gff = 37871

parsing_files_other("Oxy_Elles_Hap1", "ID")
# Oxy_Elles_Hap1: n unique sequences = 38198, n matched to gff = 38198

parsing_files_other("Oxy_Sval_h2", "ID")
# Oxy_Sval_h2: n unique sequences = 5363, n matched to gff = 5363
# remove

###################################
# to run

tmux new-session -s GeneSpace1
tmux attach-session -t GeneSpace1

module load StdEnv/2020 python/3.11.5 scipy-stack/2021a
module load StdEnv/2020 java/13.0.2
module load StdEnv/2020 r/4.2.2 glpk/5.0

export PATH=$PATH:/home/celphin/scratch/Oxyria/GeneSpace/OrthoFinder_source/orthofinder.py
alias orthofinder='python /home/celphin/scratch/Oxyria/GeneSpace/OrthoFinder_source/orthofinder.py'

R

library(GENESPACE)
gpar <- init_genespace(
  wd = "/home/celphin/scratch/Oxyria/GeneSpace/Oxyria_genomes", 
  path2mcscanx = "/home/celphin/scratch/Oxyria/GeneSpace/MCScanX/")
out <- run_genespace(gsParam = gpar) 

Checking annotation files (.bed and peptide .fa):
        DToL_h1       : 31522 / 31522 geneIDs exactly match (PASS)
        DToL_h2       : 35476 / 35476 geneIDs exactly match (PASS)
        Oxy_Elles_Hap1: 38198 / 38198 geneIDs exactly match (PASS)
        Oxy_Elles_Hap2: 35196 / 35196 geneIDs exactly match (PASS)
        Oxy_Sval_h1   : 33401 / 33401 geneIDs exactly match (PASS)
        Oxyria_Main   : 37871 / 37871 geneIDs exactly match (PASS)


# Could not find a valid path to the orthofinder program from R. To run         
# orthofinder, ensure that the orthofinder program is in the $PATH, 
# then call the following from the shell: 
# orthofinder -f /home/celphin/scratch/Oxyria/GeneSpace/Oxyria_genomes/tmp -t 16 -a 1 -X -o /home/celphin/scratch/Oxyria/GeneSpace/Oxyria_genomes/orthofinder

# Error in run_orthofinder(gsParam = gsParam, verbose = TRUE) :
  # Once OrthoFinder has been run, re-call run_genespace
