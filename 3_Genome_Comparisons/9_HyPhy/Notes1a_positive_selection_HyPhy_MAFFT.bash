#########################
# Running HyPhy programs
# Jan 2026
############################

# follow instructions here
# https://github.com/siribi/RELAX-workflow

# other tutorial: https://bioinformaticsworkbook.org/dataAnalysis/ComparativeGenomics/Finding_Positive_Selection_With_Codeml.html#gsc.tab=0

###############################
# Steps
# Run RELAX and absrel
# Explore results

######################
# Narval2
tmux new-session -s total
tmux attach-session -t total

# Load modules
module load  StdEnv/2020  gcc/9.3.0  openmpi/4.0.3 hyphy/2.5.49

##########################
# label Arctic branches as foreground in each gene tree
# ((A,B{Foreground}),C,(D,E{Foreground}));

# check tree formats
cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees
more OG0011565_tree.txt

# (((((((TAV2_LOCUS4956:0.073246,
# g11872.t1:0.137256):0.004606,
# 106302034:0.077317):0.039058,(
# LOC17883196:0.049298,
# LOC9310397:0.078487):0.06136):0.021773,
# maker-lg2-snap-gene-43.209-mRNA-1:0.22099):0.06837,
# g18107.t1:0.296852):0.343242,(((
# LOC133722319:0.036233,(
# 101294573:0.037006,
# LOC126799020:0.08705):0.007216):0.232907,
# DoctH0_Chr600000802:0.12274):0.007824,((
# LOC125470047:0.016387,
# LOC126593297:0.006243):0.144526,
# LOC18773133:0.102506):0.061272):0.228821):0.100199,(((((
# RtaG0011728.1:3.23998e-06,
# RtaG0015792.1:0.003127):0.000655,
# RnoG0007103.1:0.002416):0.026212,
# Polavi_Chr600001874:0.061237):0.033667,
# Oxyria_NCBI_Chr700002490:0.203675):0.010222,(
# FT01Gene17735.t1:0.017668,
# FEHAP213745.t1:0.0169):0.101281):0.100199);

#-----------------------------------
# Check gene  names for Arctic species
cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/genomes/
grep ">" Dryas_octopetala.f*
# Dryas_octopetala.fna:>DoctH0_Chr1000000308-RA
# Dryas_octopetala.fna:>DoctH0_Chr1000000309-RA
# Dryas_octopetala.fna:>DoctH0_Chr1000000312-RA
# Dryas_octopetala.fna:>DoctH0_Chr1000000315-RA
# Dryas_octopetala.fna:>DoctH0_Chr1000000316-RA

grep ">" Draba_nivalis.f*
# Draba_nivalis.fna:>snap_masked-lg7-processed-gene-95.56-mRNA-1 transcript offset:0 AED:0.44 eAED:0.44 QI:0|-1|0|1|-1|1|1|0|300
# Draba_nivalis.fna:>maker-lg7-snap-gene-95.218-mRNA-1 transcript offset:0 AED:0.44 eAED:0.45 QI:0|0|0|1|0.5|0.66|3|0|334
# Draba_nivalis.fna:>snap_masked-lg7-processed-gene-95.97-mRNA-1 transcript offset:0 AED:0.44 eAED:0.45 QI:0|-1|0|1|-1|1|1|0|307
# Draba_nivalis.fna:>snap_masked-lg7-processed-gene-95.98-mRNA-1 transcript offset:72 AED:0.46 eAED:0.46 QI:72|0.5|0.33|1|1|1|3|0|111
# Draba_nivalis.fna:>maker-lg7-snap-gene-95.230-mRNA-1 transcript offset:75 AED:0.47 eAED:0.47 QI:75|1|1|1|1|1|4|248|101
# Draba_nivalis.fna:>snap_masked-lg7-processed-gene-95.55-mRNA-1 transcript offset:0 AED:0.47 eAED:0.47 QI:0|-1|0|1|-1|1|1|0|141
# Draba_nivalis.fna:>snap_masked-lg7-processed-gene-95.69-mRNA-1 transcript offset:0 AED:0.51 eAED:0.53 QI:0|0|0|0.75|1|1|4|0|645
# Draba_nivalis.fna:>snap_masked-lg7-processed-gene-95.91-mRNA-1 transcript offset:0 AED:0.58 eAED:0.58 QI:0|-1|0|1|-1|1|1|0|192

grep ">" Oxyria_digyna_H1.f*
# Oxyria_digyna_H1.fna:>Oxyria_NCBI_Chr300003907-RA
# Oxyria_digyna_H1.fna:>Oxyria_NCBI_Chr300003908-RA
# Oxyria_digyna_H1.fna:>Oxyria_NCBI_Chr300003909-RA
# Oxyria_digyna_H1.fna:>Oxyria_NCBI_Chr300003910-RA
# Oxyria_digyna_H1.fna:>Oxyria_NCBI_Chr300003911-RA

grep ">" Cochlearia_groenlandica.f*
# Cochlearia_groenlandica.fna:> g2528.t1
# Cochlearia_groenlandica.fna:> g2529.t1
# Cochlearia_groenlandica.fna:> g2530.t1
# Cochlearia_groenlandica.fa:>g6166.t1
# Cochlearia_groenlandica.fa:>g6167.t1
# Cochlearia_groenlandica.fa:>g6168.t1
# Cochlearia_groenlandica.fa:>g6169.t1
# Cochlearia_groenlandica.fa:>g6170.t1



#------------------------------------
# Add {Foreground} label to Arctic branches
cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees

# test
sed -E 's/(DoctH0_Chr[^:(),]+)(:)/\1{Foreground}\2/g' OG0011565_tree.txt

sed -E '
   s/(DoctH0_Chr[^:(),]+)(:)/\1{Foreground}\2/g;
   s/(Oxyria_NCBI_Chr[^:(),]+)(:)/\1{Foreground}\2/g;
   s/((maker|snap_masked)-[^:(),]+)(:)/\1{Foreground}\3/g;
   s/(g[0-9]+\.t[0-9]+)(:)/\1{Foreground}\2/g
   ' OG0011565_tree.txt

# (((((((TAV2_LOCUS4956:0.073246,
# g11872.t1{Foreground}:0.137256):0.004606,
# 106302034:0.077317):0.039058,(
# LOC17883196:0.049298,
# LOC9310397:0.078487):0.06136):0.021773,
# maker-lg2-snap-gene-43.209-mRNA-1{Foreground}:0.22099):0.06837,
# g18107.t1{Foreground}:0.296852):0.343242,(((
# LOC133722319:0.036233,(
# 101294573:0.037006,
# LOC126799020:0.08705):0.007216):0.232907,
# DoctH0_Chr600000802{Foreground}:0.12274):0.007824,((
# LOC125470047:0.016387,
# LOC126593297:0.006243):0.144526,
# LOC18773133:0.102506):0.061272):0.228821):0.100199,(((((
# RtaG0011728.1:3.23998e-06,
# RtaG0015792.1:0.003127):0.000655,
# RnoG0007103.1:0.002416):0.026212,
# Polavi_Chr600001874:0.061237):0.033667,
# Oxyria_NCBI_Chr700002490{Foreground}:0.203675):0.010222,(
# FT01Gene17735.t1:0.017668,
# FEHAP213745.t1:0.0169):0.101281):0.100199);

#------------------------
# run for all
cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees

for f in *tree.txt
do
    # Check if tree contains at least one Arctic gene
    if grep -Eq 'DoctH0_Chr|Oxyria_NCBI_Chr|maker-|snap_masked-|g[0-9]+\.t[0-9]+' "$f"
    then
        sed -E '
        s/(DoctH0_Chr[^:(),]+)(:)/\1{Foreground}\2/g;
        s/(Oxyria_NCBI_Chr[^:(),]+)(:)/\1{Foreground}\2/g;
        s/((maker|snap_masked)-[^:(),]+)(:)/\1{Foreground}\3/g;
        s/(g[0-9]+\.t[0-9]+)(:)/\1{Foreground}\2/g
        ' "$f" > "${f}_HyPhy.txt"
    fi
done

ls *_HyPhy.txt | wc -l
#19,152

#########################
# Move only those trees and data with Arctic species to new folder

cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees

mkdir -p Arctic_trees

for f in *_HyPhy.txt
do
    base=${f%_HyPhy.txt}
    cp ${base}* Arctic_trees/
done

# check
cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees
ls *_HyPhy.txt | wc -l
# 0 
ls *_tree.txt | wc -l
# 3993

grep "Doct" *_tree.txt
# none

cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees/Arctic_trees
ls *_HyPhy.txt | wc -l
# 19152
ls *_tree.txt | wc -l
# 19152

########################
# Run remove duplicates script test
#  download script

module load  StdEnv/2020  gcc/9.3.0  openmpi/4.0.3 hyphy/2.5.49

wget https://raw.githubusercontent.com/veg/hyphy-analyses/master/remove-duplicates/remove-duplicates.bf
chmod +x remove-duplicates.bf

cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees/Arctic_trees

# test
hyphy remove-duplicates.bf --msa OG0011565_tree.txt_pal2nal.fasta --tree OG0011565_tree.txt_HyPhy.txt --output OG0011565_tree.txt_unique.nxh

#####################
# Run absrel test
# https://hyphy.org/tutorials/CL-prompt-tutorial/
# Enter 1 for "Selection Analyses", and then 6 for "aBSREL"
# Enter 1 to select the Universal genetic code.
# Select a coding sequence alignment file: 

#---------------
# Automatic command line test

cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees/Arctic_trees

hyphy absrel  --alignment OG0011565_tree.txt_unique.nxh --branches FOREGROUND
# * DOCTH0_CHR600000802, p-value =  0.01770
# * DOCTH0_CHR600000802, p-value =  0.00557 # with foreground

hyphy relax  --alignment OG0011565_tree.txt_unique.nxh --branches FOREGROUND
# Check errors.log for execution error details.


##################################
# loop through orthogroups that include Arctic spp 
# run as slurm array

# get file list of non empty files
cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees/Arctic_trees

mkdir -p logs

tmux attach-session -t backup

for f in *_tree.txt; do
    msa="${f}_pal2nal.fasta"
    tree="${f}_HyPhy.txt"

    # Require both files
    [[ -f "$msa" && -f "$tree" ]] || continue

    # Require at least 2 sequences
    seq_count=$(grep -c "^>" "$msa")
    [[ $seq_count -ge 2 ]] || continue

    # Skip zero-length sequences (ignoring gaps)
    if awk '
        /^>/ {
            if (seq_len == 0 && NR > 1) exit 1
            seq_len = 0
            next
        }
        { gsub("-", ""); seq_len += length($0) }
        END { if (seq_len == 0) exit 1 }
    ' "$msa"
    then
        echo "$f"
    fi

done > filtered_tree_list.txt

# check count
wc -l filtered_tree_list.txt
# 16401 filtered_tree_list.txt

#--------------------------
# check commands
hyphy absrel --help

wc -l  filtered_tree_list.txt
#16401

# check max array size: 
 sacctmgr show assoc user=$USER format=User,Account,MaxSubmitJobs
      # User    Account MaxSubmit
# ---------- ---------- ---------
   # celphin def-cronk+      1000
   # celphin def-cronk+      1000
   # celphin def-henry+      1000
   # celphin def-henry+      1000
   # celphin def-nbl_c+      1000
   # celphin def-nbl_g+      1000
   # celphin def-riese+      1000
   # celphin def-riese+      1000

# use nano to import text
cat << EOF > absrel_array.sh
#!/bin/bash
#SBATCH --account=def-henryg
#SBATCH --time=0-12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=32G
#SBATCH --output=logs/absrel_%A_%a.out
#SBATCH --error=logs/absrel_%A_%a.err

module load StdEnv/2020 gcc/9.3.0 openmpi/4.0.3 hyphy/2.5.49

cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees/Arctic_trees

CHUNK=20

START=$(( (SLURM_ARRAY_TASK_ID - 1) * CHUNK + 1 ))
END=$(( START + CHUNK - 1 ))

TOTAL=$(wc -l < filtered_tree_list.txt)

if [ $END -gt $TOTAL ]; then
    END=$TOTAL
fi

echo "Processing lines $START to $END"

sed -n "${START},${END}p" filtered_tree_list.txt | while read f
do
    echo "Processing $f"

    msa="${f}_pal2nal.fasta"
    tree="${f}_HyPhy.txt"
    output="${f}_unique.nxh"

    if [[ -f "$output" ]]; then
        echo "Output exists for $f — skipping"
        continue
    fi

    hyphy remove-duplicates.bf \
        --msa "$msa" \
        --tree "$tree" \
        --output "$output"

    hyphy absrel \
        --alignment "$output" \
        --branches FOREGROUND

    echo "Finished $f"
done

EOF

chmod +x absrel_array.sh
dos2unix absrel_array.sh

sbatch --array=1-850%100 absrel_array.sh

ls *nxh
ls *json

###############################
# Checking results original

cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees/Arctic_trees

ls *ABSREL.json | wc -l
# 16328

##########################

################################
# Rerun with prank alignments

# get file list of non empty files
cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees/Arctic_trees

mkdir -p logs

for f in *_tree.txt; do
    msa="${f}_prank_pal2nal.fasta"
    tree="${f}_HyPhy.txt"

    # Require both files
    [[ -f "$msa" && -f "$tree" ]] || continue

    # Require at least 2 sequences
    seq_count=$(grep -c "^>" "$msa")
    [[ $seq_count -ge 2 ]] || continue

    # Skip zero-length sequences (ignoring gaps)
    if awk '
        /^>/ {
            if (seq_len == 0 && NR > 1) exit 1
            seq_len = 0
            next
        }
        { gsub("-", ""); seq_len += length($0) }
        END { if (seq_len == 0) exit 1 }
    ' "$msa"
    then
        echo "$f"
    fi

done > filtered_prank_tree_list.txt

# check count
wc -l filtered_prank_tree_list.txt
# 4053 

#-------------------------
# use nano to import text
cat << EOF > absrel_prank_array.sh
#!/bin/bash
#SBATCH --account=def-henryg
#SBATCH --time=0-6:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=32G
#SBATCH --output=logs/absrel_prank_%A_%a.out
#SBATCH --error=logs/absrel_prank_%A_%a.err

module load StdEnv/2020 gcc/9.3.0 openmpi/4.0.3 hyphy/2.5.49

cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees/Arctic_trees

CHUNK=10

START=$(( (SLURM_ARRAY_TASK_ID - 1) * CHUNK + 1 ))
END=$(( START + CHUNK - 1 ))

TOTAL=$(wc -l < filtered_prank_tree_list.txt)

if [ $END -gt $TOTAL ]; then
    END=$TOTAL
fi

echo "Processing lines $START to $END"

sed -n "${START},${END}p" filtered_prank_tree_list.txt | while read f
do
    echo "Processing $f"

    msa="${f}_prank_pal2nal.fasta"
    tree="${f}_HyPhy.txt"
    output="${f}_prank_unique.nxh"

    if [[ -f "$output" ]]; then
        echo "Output exists for $f — skipping"
        continue
    fi

    hyphy remove-duplicates.bf \
        --msa "$msa" \
        --tree "$tree" \
        --output "$output"

    hyphy absrel \
        --alignment "$output" \
        --branches FOREGROUND

    echo "Finished $f"
done

EOF

chmod +x absrel_prank_array.sh
dos2unix absrel_prank_array.sh

sbatch --array=1-406%100 absrel_prank_array.sh
# Submitted batch job 57303060

#--------------------
# Check
cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees/Arctic_trees

ls *_prank_unique.nxh | wc -l
ls *_prank_unique.nxh.ABSREL.json  | wc -l
# 4053
# 3961

#########################################