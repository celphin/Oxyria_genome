#############################
# Liftoff to move annotations over
# Dryas
# August 2024
# https://github.com/agshumate/Liftoff
##################################
# try liftoff to move Dryas octo annotation over to other genomes

cd ~/scratch/Annotation/liftoff/data
cp ~/scratch/Dryas/Dryas_genomes/Dry-int_ragtag_output/ragtag.scaffold.fasta  Dry-int-chr.fasta
cp ~/scratch/Dryas/Dryas_genomes/Dry-drumm_ragtag_output/ragtag.scaffold.fasta  Dry-drumm-chr.fasta
cp ~/scratch/Dryas/Dryas_genomes/Dry-alask_ragtag_output/ragtag.scaffold.fasta  Dry-alask-chr.fasta
cp ~/scratch/Dryas/Dryas_genomes/Dry-octo-H1_ragtag_output/ragtag.scaffold.fasta  Dry-octo-H1-chr.fasta
cp ~/scratch/Dryas/Dryas_genomes/Dry-octo-H2_ragtag_output/ragtag.scaffold.fasta  Dry-octo-H2-chr.fasta

cd /home/celphin/scratch/Annotation/CAP_Snakemake/DoctH0/FINAL_ANNOTATION
cp FINAL_DoctH0.AED_0.6.sorted.gff3 ~/scratch/Annotation/liftoff/data/

cd /home/celphin/scratch/Annotation/CAP_Snakemake/DoctH0/Ref
cp DoctH0_Main.fasta ~/scratch/Annotation/liftoff/data/

# change names to avoid overlap
cp FINAL_DoctH0.AED_0.6.sorted.gff3 Dry-int_FINAL_DoctH0.AED_0.6.sorted.gff3
cp FINAL_DoctH0.AED_0.6.sorted.gff3 Dry-drumm_FINAL_DoctH0.AED_0.6.sorted.gff3
cp FINAL_DoctH0.AED_0.6.sorted.gff3 Dry-alask_FINAL_DoctH0.AED_0.6.sorted.gff3
cp FINAL_DoctH0.AED_0.6.sorted.gff3 Dry-octo-H1_FINAL_DoctH0.AED_0.6.sorted.gff3
cp FINAL_DoctH0.AED_0.6.sorted.gff3 Dry-octo-H2_FINAL_DoctH0.AED_0.6.sorted.gff3

cp DoctH0_Main.fasta Dry-int_DoctH0_Main.fasta
cp DoctH0_Main.fasta Dry-drumm_DoctH0_Main.fasta
cp DoctH0_Main.fasta Dry-alask_DoctH0_Main.fasta
cp DoctH0_Main.fasta Dry-octo-H1_DoctH0_Main.fasta
cp DoctH0_Main.fasta Dry-octo-H2_DoctH0_Main.fasta

cd ~/scratch/Annotation/liftoff/

nano run_Dryas_liftoff.sh
#!/bin/bash
#SBATCH --account=def-rieseber
#SBATCH --job-name=liftoff
#SBATCH --time=0-15:00:00
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

cd /home/celphin/scratch/Annotation/liftoff

SPP_Hap=$1

# or no conda
#module load gcc python/3.10 minimap2 parasail
#source ~/ENV/bin/activate
#python -c "import liftoff; print(liftoff.__version__)"

liftoff \
-g ./data/${SPP_Hap}_FINAL_DoctH0.AED_0.6.sorted.gff3 -p 40 -o ./output/${SPP_Hap}_liftoffpolishslurm.fasta \
-dir ${SPP_Hap}_intermed_slurm -u ${SPP_Hap}_unmapped_features_slurm.txt \
-infer_genes -copies -a 0.95 -s 0.95 -d 5.0 -flank 0.8 -polish \
./data/${SPP_Hap}-chr.fasta ./data/${SPP_Hap}_DoctH0_Main.fasta


#--------------------------
# run on all the Dryas genomes

sbatch run_Dryas_liftoff.sh Dry-int 
sbatch run_Dryas_liftoff.sh Dry-octo-H1 
sbatch run_Dryas_liftoff.sh Dry-octo-H2 
sbatch run_Dryas_liftoff.sh Dry-drumm 
sbatch run_Dryas_liftoff.sh Dry-alask 

# not working properly?? - need chromosome level assemblies


######################
# Try interactive runs

tmux new-session -s liftoff
tmux attach-session -t liftoff

salloc -c40 --time 12:55:00 --mem 120000m --account def-rieseber

source ~/miniconda2/bin/activate liftoff
module load StdEnv/2020
module load python/3.10.2
module load scipy-stack/2022a
module load gcc/9.3.0
module load cuda/11.7
module load hdf5/1.12.1
module load parasail/2.5
module load minimap2/2.24

cd /home/celphin/scratch/Annotation/liftoff

SPP_Hap=Dry-int 

liftoff \
-g ./data/${SPP_Hap}_FINAL_DoctH0.AED_0.6.sorted.gff3 -p 40 -o ./output/${SPP_Hap}_liftoffpolish.fasta \
-dir ${SPP_Hap}_intermed -u ${SPP_Hap}_unmapped_features.txt \
-infer_genes -copies -a 0.95 -s 0.95 -d 5.0 -flank 0.8 -polish \
./data/${SPP_Hap}.fasta ./data/${SPP_Hap}_DoctH0_Main.fasta

SPP_Hap=Dry-drumm 

liftoff \
-g ./data/${SPP_Hap}_FINAL_DoctH0.AED_0.6.sorted.gff3 -p 40 -o ./output/${SPP_Hap}_liftoffpolish.fasta \
-dir ${SPP_Hap}_intermed -u ${SPP_Hap}_unmapped_features.txt \
-infer_genes -copies -a 0.95 -s 0.95 -d 5.0 -flank 0.8 -polish \
./data/${SPP_Hap}.fasta ./data/${SPP_Hap}_DoctH0_Main.fasta


SPP_Hap=Dry-alask 

liftoff \
-g ./data/${SPP_Hap}_FINAL_DoctH0.AED_0.6.sorted.gff3 -p 40 -o ./output/${SPP_Hap}_liftoffpolish.fasta \
-dir ${SPP_Hap}_intermed -u ${SPP_Hap}_unmapped_features.txt \
-infer_genes -copies -a 0.95 -s 0.95 -d 5.0 -flank 0.8 -polish \
./data/${SPP_Hap}.fasta ./data/${SPP_Hap}_DoctH0_Main.fasta


SPP_Hap=Dry-octo-H1 

liftoff \
-g ./data/${SPP_Hap}_FINAL_DoctH0.AED_0.6.sorted.gff3 -p 40 -o ./output/${SPP_Hap}_liftoffpolish.fasta \
-dir ${SPP_Hap}_intermed -u ${SPP_Hap}_unmapped_features.txt \
-infer_genes -copies -a 0.95 -s 0.95 -d 5.0 -flank 0.8 -polish \
./data/${SPP_Hap}.fasta ./data/${SPP_Hap}_DoctH0_Main.fasta


SPP_Hap=Dry-octo-H2 

liftoff \
-g ./data/${SPP_Hap}_FINAL_DoctH0.AED_0.6.sorted.gff3 -p 40 -o ./output/${SPP_Hap}_liftoffpolish.fasta \
-dir ${SPP_Hap}_intermed -u ${SPP_Hap}_unmapped_features.txt \
-infer_genes -copies -a 0.95 -s 0.95 -d 5.0 -flank 0.8 -polish \
./data/${SPP_Hap}.fasta ./data/${SPP_Hap}_DoctH0_Main.fasta

########################################
# copy over Dryas liftoff annotations and genomes

cd /home/celphin/scratch/Dryas/Dryas_genomes/Liftoff
cp /home/celphin/scratch/Annotation/liftoff/output/Dry-*fast* .
rename fasta_polished gff3 *

cp /home/celphin/scratch/Annotation/liftoff/data/Dry-*fasta .
cp /home/celphin/scratch/Annotation/liftoff/data/Dry-*gff3 .

####################################################
# try Liftoff tools
# https://github.com/agshumate/LiftoffTools

tmux new-session -s liftoff
tmux attach-session -t liftoff

cd /home/celphin/scratch/Dryas/Dryas_genomes/Liftoff
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
liftofftools all -r Dry-alask_DoctH0_Main.fasta -t Dry-alask-chr.fasta \
-rg Dry-alask_FINAL_DoctH0.AED_0.6.sorted.gff3 \
-tg Dry-alask_liftoffpolishslurm.gff3 \
-dir Dry-alask_liftofftools

liftofftools all -r Dry-int_DoctH0_Main.fasta -t Dry-int-chr.fasta \
-rg Dry-int_FINAL_DoctH0.AED_0.6.sorted.gff3 \
-tg Dry-int_liftoffpolishslurm.gff3 \
-dir Dry-int_liftofftools

liftofftools all -r Dry-drumm_DoctH0_Main.fasta -t Dry-drumm-chr.fasta \
-rg Dry-drumm_FINAL_DoctH0.AED_0.6.sorted.gff3 \
-tg Dry-drumm_liftoffpolishslurm.gff3 \
-dir Dry-drumm_liftofftools
# failed 

# Analyzing synteny
# Extracting transcript sequences
# Analyzing protein-coding clusters
# Traceback (most recent call last):
  # File "/home/celphin/.local/bin/liftofftools", line 8, in <module>
    # sys.exit(main())
  # File "/home/celphin/.local/lib/python3.10/site-packages/liftofftools/liftofftools.py", l
# ine 34, in main
    # analyze_clusters.main(ref_proteins, target_proteins, ref_trans, target_trans, ref_db,
# target_db, args)
  # File "/home/celphin/.local/lib/python3.10/site-packages/liftofftools/cluster_analysis/an
# alyze_clusters.py", line 14, in main
    # run_coding_workflow(ref_proteins, ref_db, target_db, target_proteins,ref_trans, target
# _trans, args)
  # File "/home/celphin/.local/lib/python3.10/site-packages/liftofftools/cluster_analysis/an
# alyze_clusters.py", line 21, in run_coding_workflow
    # unmapped_coding = analyze_proteins(ref_proteins, ref_db, target_db, target_proteins, a
# rgs)
  # File "/home/celphin/.local/lib/python3.10/site-packages/liftofftools/cluster_analysis/an
# alyze_clusters.py", line 37, in analyze_proteins
    # ref_clusters, target_clusters = cluster(ref_proteins, target_proteins, ref_protein_cod
# ing_genes,
  # File "/home/celphin/.local/lib/python3.10/site-packages/liftofftools/cluster_analysis/an
# alyze_clusters.py", line 48, in cluster
    # target_tsv = build_clusters(target_seqs,  target_genes, target_output, args)
  # File "/home/celphin/.local/lib/python3.10/site-packages/liftofftools/cluster_analysis/an
# alyze_clusters.py", line 107, in build_clusters
    # longest_seqs = sequence_dict.get_longest_isoform_dict(gene_list)
  # File "/home/celphin/.local/lib/python3.10/site-packages/liftofftools/sequences.py", line
 # 53, in get_longest_isoform_dict
    # longest_isoform_dict[gene] = self[longest_isoform]
# KeyError: ''

# OY992843.1_RagTag       Liftoff mRNA    113092  120929  .       -       .       ID=DoctH0_Chr1000000265-RA;Name=DoctH0_Chr1000000265-RA;Alias=maker-DoctH0-11-exonerate_protein2genome-gene-0.333-mRNA-1;Parent=DoctH0_Chr1000000265;matches_ref_protein=False;valid_ORF=False;missing_start_codon=True;extra_copy_number=0
# OY992843.1_RagTag       Liftoff mRNA    113092  120929  .       -       .       ID=DoctH0_Chr1000000265-RA;Name=DoctH0_Chr1000000265-RA;Alias=maker-DoctH0-11-exonerate_protein2genome-gene-0.333-mRNA-1;Parent=DoctH0_Chr1000000265;extra_copy_number=0
# OY992843.1_RagTag       Liftoff mRNA    113092  120929  .       -       .       ID=DoctH0_Chr1000000265-RA;Name=DoctH0_Chr1000000265-RA;Alias=maker-DoctH0-11-exonerate_protein2genome-gene-0.333-mRNA-1;Parent=DoctH0_Chr1000000265;matches_ref_protein=False;valid_ORF=False;missing_start_codon=True;extra_copy_number=0
# OY992843.1_RagTag       Liftoff exon    113092  120929  .       -       .       ID=DoctH0_Chr1000000265-RA:5;Alias=character(0);Parent=DoctH0_Chr1000000265-RA;extra_copy_number=0
# OY992843.1_RagTag       Liftoff exon    113092  120929  .       -       .       ID=DoctH0_Chr1000000265-RA:5;Alias=character(0);Parent=DoctH0_Chr1000000265-RA;extra_copy_number=0
# OY992843.1_RagTag       Liftoff CDS     113092  120929  .       -       0       ID=DoctH0_Chr1000000265-RA:cds;Alias=character(0);Parent=DoctH0_Chr1000000265-RA;extra_copy_number=0
# OY992843.1_RagTag       Liftoff CDS     113092  120929  .       -       0       ID=DoctH0_Chr1000000265-RA:cds;Alias=character(0);Parent=DoctH0_Chr1000000265-RA;extra_copy_number=0
#---------------------------
# Try just variants detection

tmux new-session -s liftoff
tmux attach-session -t liftoff

# parasail_memalign: posix_memalign failed: Cannot allocate memory
salloc -c1 --time 2:55:00 --mem 120000m --account def-rieseber

cd /home/celphin/scratch/Dryas/Dryas_genomes/Liftoff
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

liftofftools variants -r Dry-alask_DoctH0_Main.fasta -t Dry-alask-chr.fasta \
-rg Dry-alask_FINAL_DoctH0.AED_0.6.sorted.gff3 \
-tg Dry-alask_liftoffpolishslurm.gff3 \
-dir Dry-alask_liftofftools_VARIANTS

liftofftools variants -r Dry-int_DoctH0_Main.fasta -t Dry-int-chr.fasta \
-rg Dry-int_FINAL_DoctH0.AED_0.6.sorted.gff3 \
-tg Dry-int_liftoffpolishslurm.gff3 \
-dir Dry-int_liftofftools_VARIANTS

liftofftools variants -r Dry-drumm_DoctH0_Main.fasta -t Dry-drumm-chr.fasta \
-rg Dry-drumm_FINAL_DoctH0.AED_0.6.sorted.gff3 \
-tg Dry-drumm_liftoffpolishslurm.gff3 \
-dir Dry-drumm_liftofftools_VARIANTS
# fails
# liftoff tools gffutils database build failed with UNIQUE constraint failed: features.id

liftofftools variants -r Dry-octo-H1_DoctH0_Main.fasta -t Dry-octo-H1-chr.fasta \
-rg Dry-octo-H1_FINAL_DoctH0.AED_0.6.sorted.gff3 \
-tg Dry-octo-H1_liftoffpolishslurm.gff3 \
-dir Dry-octo-H1_liftofftools_VARIANTS
# fails
# liftoff tools gffutils database build failed with UNIQUE constraint failed: features.id

liftofftools variants -r Dry-octo-H2_DoctH0_Main.fasta -t Dry-octo-H2-chr.fasta \
-rg Dry-octo-H2_FINAL_DoctH0.AED_0.6.sorted.gff3 \
-tg Dry-octo-H2_liftoffpolishslurm.gff3 \
-dir Dry-octo-H2_liftofftools_VARIANTS

#####################################
cd /home/celphin/scratch/Dryas/Dryas_genomes/Liftoff/Dry-alask_liftofftools_VARIANTS/
Spp_hap=Dry-alask

cd /home/celphin/scratch/Dryas/Dryas_genomes/Liftoff/Dry-int_liftofftools_VARIANTS/
Spp_hap=Dry-int

cd /home/celphin/scratch/Dryas/Dryas_genomes/Liftoff/Dry-drumm_liftofftools_VARIANTS/
Spp_hap=Dry-drumm

cd /home/celphin/scratch/Dryas/Dryas_genomes/Liftoff/Dry-octo-H1_liftofftools_VARIANTS/
Spp_hap=Dry-octo-H1

cd /home/celphin/scratch/Dryas/Dryas_genomes/Liftoff/Dry-octo-H2_liftofftools_VARIANTS/
Spp_hap=Dry-octo-H2

#-------------------
grep "frameshift"  variant_effects >${Spp_hap}_frameshift.txt
wc -l ${Spp_hap}_frameshift.txt
#4222 
#6734 Dry-int_frameshift.txt
# 6709 Dry-octo-H2_frameshift.txt

grep "inframe_deletion"  variant_effects >${Spp_hap}_inframe_deletion.txt
wc -l ${Spp_hap}_inframe_deletion.txt
#1465 
#2319 Dry-int_inframe_deletion.txt
# 2318 Dry-octo-H2_inframe_deletion.txt

grep "start_lost"  variant_effects >${Spp_hap}_start_lost.txt
wc -l ${Spp_hap}_start_lost.txt
#4540
#6054 Dry-int_start_lost.txt

grep "3'_truncated" variant_effects >${Spp_hap}_3_truncated.txt
wc -l ${Spp_hap}_3_truncated.txt
#1008
#1339 Dry-int_3_truncated.txt
#[1290 Dry-octo-H2_3_truncated.txt

grep "identical" variant_effects >${Spp_hap}_identical.txt
wc -l ${Spp_hap}_identical.txt
# 743
# 1267 Dry-int_identical.txt
# 1297 Dry-octo-H2_identical.txt

grep "nonsynonymous"  variant_effects >${Spp_hap}_nonsynonymous.txt
wc -l ${Spp_hap}_nonsynonymous.txt
# 6450
# 11663 Dry-int_nonsynonymous.txt
# 11770 Dry-octo-H2_nonsynonymous.txt

grep "NA"  variant_effects >${Spp_hap}_NA.txt
wc -l ${Spp_hap}_NA.txt
# 13759
# 25904 Dry-int_NA.txt
# 26063 Dry-octo-H2_NA.txt

grep -v "nonsynonymous"  variant_effects | grep "synonymous"  >${Spp_hap}_synonymous.txt
wc -l ${Spp_hap}_synonymous.txt
# 1042 Dry-alask_synonymous.txt
# 1969 Dry-int_synonymous.txt
# 1920 Dry-octo-H2_synonymous.txt

#--------
grep "stop_gained" variant_effects >${Spp_hap}_stop_gained.txt
wc -l ${Spp_hap}_stop_gained.txt
# 1107 Dry-octo-H2_stop_gained.txt

grep "5'_truncated" variant_effects >${Spp_hap}_5_truncated.txt
wc -l ${Spp_hap}_5_truncated.txt
# 98 Dry-octo-H2_5_truncated.txt
# 103 Dry-int_5_truncated.txt
# 60 Dry-alask_5_truncated.txt

grep "inframe_insertion"  variant_effects >${Spp_hap}_inframe_insertion.txt
wc -l ${Spp_hap}_inframe_insertion.txt
# 1486 Dry-octo-H2_inframe_insertion.txt

wc -l variant_effects
#65593 variant_effects

# D. alask
# about 16 % of genes have no nonsyn changes
#------------------------
# D. int 

#---------------------------

# synonymous - A point mutation in the target transcript that does not change the amino acid sequence.
# nonsynonymous - A point mutation in the target transcript that changes the amino acid sequence.
# inframe deletion - A deletion in the target transcript sequence that is of some length divisible by 3.
# inframe insertion - An insertion in the target transcript sequence that is of some length divisible by 3.
# start codon loss - A point mutation in the start codon of the target transcript sequence.
# 5' truncation - A deletion of the 5' end of the target transcript. 
# 3' truncation - A deletion of the 3' end of the target transcript. 
# frameshift - A deletion or insertion in the target transcript that is of some length not divisible by 3. 
# stop codon gained - A point mutation in the target transcript that introduces a premature stop codon. 


##############################
# try https://github.com/gpertea/gffread
# to get to protein conversion

# Install
 cd /home/celphin/scratch/Annotation
  git clone https://github.com/gpertea/gffread
  cd gffread
  make release

cd /home/celphin/scratch/Dryas/Dryas_genomes/Liftoff

/home/celphin/scratch/Annotation/gffread/gffread -h

# to run proteins
/home/celphin/scratch/Annotation/gffread/gffread -g Dry-alask-chr.fasta Dry-alask_liftoffpolishslurm.gff3 -y Dry-alask_proteins.fa
/home/celphin/scratch/Annotation/gffread/gffread -g Dry-int-chr.fasta Dry-int_liftoffpolishslurm.gff3 -y Dry-int_proteins.fa
/home/celphin/scratch/Annotation/gffread/gffread -g Dry-drumm-chr.fasta Dry-drumm_liftoffpolishslurm.gff3 -y Dry-drumm_proteins.fa
/home/celphin/scratch/Annotation/gffread/gffread -g Dry-octo-H1-chr.fasta Dry-octo-H1_liftoffpolishslurm.gff3 -y Dry-octo-H1_proteins.fa
/home/celphin/scratch/Annotation/gffread/gffread -g Dry-octo-H2-chr.fasta Dry-octo-H2_liftoffpolishslurm.gff3 -y Dry-octo-H2_proteins.fa
/home/celphin/scratch/Annotation/gffread/gffread -g DoctH0_Main.fasta FINAL_DoctH0.AED_0.6.sorted.gff3 -y Dry-octo-H0_proteins.fa
rename octo-H2 ajan *
#-------------------------
# to run proteins
/home/celphin/scratch/Annotation/gffread/gffread -g Dry-alask-chr.fasta Dry-alask_liftoffpolishslurm.gff3 -x Dry-alask_cds.fa
/home/celphin/scratch/Annotation/gffread/gffread -g Dry-int-chr.fasta Dry-int_liftoffpolishslurm.gff3 -x Dry-int_cds.fa
/home/celphin/scratch/Annotation/gffread/gffread -g Dry-drumm-chr.fasta Dry-drumm_liftoffpolishslurm.gff3 -x Dry-drumm_cds.fa
/home/celphin/scratch/Annotation/gffread/gffread -g Dry-octo-H1-chr.fasta Dry-octo-H1_liftoffpolishslurm.gff3 -x Dry-octo-H1_cds.fa
/home/celphin/scratch/Annotation/gffread/gffread -g Dry-ajan-chr.fasta Dry-ajan_liftoffpolishslurm.gff3 -x Dry-ajan_cds.fa
/home/celphin/scratch/Annotation/gffread/gffread -g DoctH0_Main.fasta FINAL_DoctH0.AED_0.6.sorted.gff3 -x Dry-octo-H0_cds.fa

##########################
# check BUSCO scores

tmux new-session -s BUSCO
tmux attach-session -t BUSCO

cd /home/celphin/scratch/Dryas/Dryas_genomes/Liftoff

salloc -c10 --time 2:55:00 --mem 120000m --account def-rieseber

module load StdEnv/2020
module load gcc/9.3.0 python augustus hmmer blast+ metaeuk prodigal r bbmap
source ~/busco_env/bin/activate

# Dryas genomes
busco --offline --in Dry-alask_proteins.fa \
--out  BUSCO_Dryasalask_liftoff_pep  --lineage_dataset eudicots_odb10 --mode protein --cpu 10 \
--download_path ~/BUSCO_downloads/
# C:76.9%[S:75.5%,D:1.4%],F:8.2%,M:14.9%,n:2326

busco --offline --in Dry-int_proteins.fa \
--out  BUSCO_Dryasint_liftoff_pep  --lineage_dataset eudicots_odb10 --mode protein --cpu 10 \
--download_path ~/BUSCO_downloads/
# C:78.0%[S:76.4%,D:1.6%],F:7.7%,M:14.3%,n:2326 

busco --offline --in Dry-drumm_proteins.fa \
--out  BUSCO_Dryasdrumm_liftoff_pep  --lineage_dataset eudicots_odb10 --mode protein --cpu 10 \
--download_path ~/BUSCO_downloads/
# C:79.2%[S:77.0%,D:2.2%],F:5.5%,M:15.3%,n:2326

busco --offline --in Dry-octo-H0_proteins.fa \
--out  BUSCO_Dryasoct-main_pep  --lineage_dataset eudicots_odb10 --mode protein --cpu 10 \
--download_path ~/BUSCO_downloads/
# C:90.5%[S:87.4%,D:3.1%],F:4.2%,M:5.3%,n:2326
























###################################
# OLD - need to check conversion
# make protein files for Dryas from gff

module load StdEnv/2020 python/3.10.2
source ~/gff3_env/bin/activate

cd Dryas_drummondii
cp /home/celphin/scratch/Annotation/liftoff/data/Dry-drumm*fasta .
gff3_to_fasta -g  Dryas_drummondii.gff -f Dry-drumm-chr.fasta -st pep -d simple -o Dryas_drummondii_pep

# does not work

#----------------------------
# Try AGAT

cd Dryas_integrifolia
cp /home/celphin/scratch/Annotation/liftoff/data/Dry-int*fasta .

cd Dryas_alaskensis
cp /home/celphin/scratch/Annotation/liftoff/data/Dry-alask*fasta .

# reformat with AGAT instead??
# https://agat.readthedocs.io/en/latest/tools/agat_sp_fix_features_locations_duplicated.html

# Load modules on graham or cedar (or use instructions for conda at https://github.com/NBISweden/AGAT if using your own computer)
module load StdEnv/2023
module load apptainer/1.2.4

# Install `AGAT` via `Apptainer` into a temp directory
apptainer pull docker://quay.io/biocontainers/agat:1.4.0--pl5321hdfd78af_0
apptainer run agat_1.4.0--pl5321hdfd78af_0.sif

# Kick the tires, was it successfully installed?
# agat_convert_sp_gxf2gxf.pl --help | less # It should show tool specific information. Use 'q' to return the prompt
# agat_<tab>  # tab-completion should result in a list of all the possible scripts; hit space to advance to next page, or 'q' to return the prompt
# agat_convert_sp_gxf2gxf.pl --help

agat_sp_manage_IDs.pl --gff Dry-drumm_liftoffpolishslurm.gff3 -p all -o Dry-drumm_liftoff_cleaned.gff3
agat_sp_extract_sequences.pl -g Dry-drumm_liftoff_cleaned.gff3 -f Dry-drumm-chr.fasta -t cds -o Dry-drumm-cds.fasta

agat_sp_extract_sequences.pl -g Dryas_integrifolia.gff -f Dry-int-chr.fasta -t cds  -o Dry-int-cds.fasta

agat_sp_extract_sequences.pl -g Dryas_alaskensis.gff -f Dry-alask-chr.fasta -t cds  -o Dry-alask-cds.fasta

exit
