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
    mv ${base}* Arctic_trees/
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


#-----------------------------
# OR Run each orthogroup in a new bash script of its own

for f in *_tree.txt; do
cat << EOF > absrel_${f}.sh
#!/bin/bash
#SBATCH --account=def-rieseber
#SBATCH --time=0-01:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=18G

module load  StdEnv/2020  gcc/9.3.0  openmpi/4.0.3 hyphy/2.5.49

cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees/Arctic_trees

hyphy remove-duplicates.bf --msa "${f}_pal2nal.fasta" --tree "${f}_HyPhy.txt" --output "${f}_unique.nxh"
hyphy absrel --alignment "${f}_unique.nxh\"
#hyphy relax --alignment "${f}_unique.nxh"

EOF

chmod +x absrel_${f}.sh
sbatch absrel_${f}.sh

done 


###############################
# Checking results

cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees/Arctic_trees

ls *ABSREL.json | wc -l
# 16328

##########################
# Which orthogroups are demonstrating positive selection

# look for found **1** branches under selection among **39** tested

for f in *ABSREL.json; do
  jq -r --arg file "$f" \
    '$file + "\t" + (.["test results"]["positive test results"]|tostring) + "\t" + (.["test results"].tested|tostring)' "$f"
done >> Total_ABSREL_results.txt

wc -l Total_ABSREL_results.txt
#16254 Total_ABSREL_results.txt

# filter for only those with positive selection
awk '$2 != 0' Total_ABSREL_results.txt | sort -k2,2nr > ABSREL_nonzero_sorted.txt

wc -l ABSREL_nonzero_sorted.txt
# 4076 ABSREL_nonzero_sorted.txt

more ABSREL_nonzero_sorted.txt

# OG0015441_tree.txt_unique.nxh.ABSREL.json       6       10
# OG0002756_tree.txt_unique.nxh.ABSREL.json       5       12
# OG0004437_tree.txt_unique.nxh.ABSREL.json       4       5
# OG0005531_tree.txt_unique.nxh.ABSREL.json       4       6
# OG0005857_tree.txt_unique.nxh.ABSREL.json       4       20
# OG0012691_tree.txt_unique.nxh.ABSREL.json       4       18
# OG0012770_tree.txt_unique.nxh.ABSREL.json       4       4
# OG0013538_tree.txt_unique.nxh.ABSREL.json       4       12
# OG0014351_tree.txt_unique.nxh.ABSREL.json       4       8
# OG0015407_tree.txt_unique.nxh.ABSREL.json       4       10
# OG0015559_tree.txt_unique.nxh.ABSREL.json       4       10
# OG0016078_tree.txt_unique.nxh.ABSREL.json       4       7
# OG0017813_tree.txt_unique.nxh.ABSREL.json       4       7
# OG0018973_tree.txt_unique.nxh.ABSREL.json       4       5
# OG0020180_tree.txt_unique.nxh.ABSREL.json       4       5
# OG0020223_tree.txt_unique.nxh.ABSREL.json       4       5
# OG0020857_tree.txt_unique.nxh.ABSREL.json       4       5
# OG0020866_tree.txt_unique.nxh.ABSREL.json       4       5
# OG0022631_tree.txt_unique.nxh.ABSREL.json       4       4
# OG0022756_tree.txt_unique.nxh.ABSREL.json       4       4


more OG0001467_tree.txt_unique.nxh.ABSREL.json
     # "DOCTH0_CHR100000472":"test",
     # "DOCTH0_CHR600001113":"test",
     # "G10838_T1":"test",
     # "G21743_T1":"test",
     # "G29482_T1":"test",
     # "MAKER_LG2_SNAP_GENE_0_246_MRNA_1":"test",
     # "MAKER_LG6_SNAP_GENE_73_102_MRNA_1":"test",
     # "OXYRIA_NCBI_CHR300002584":"test",
     # "OXYRIA_NCBI_CHR300005439":"test",
     # "OXYRIA_NCBI_CHR500005514":"test",
     # "OXYRIA_NCBI_CHR500005515":"test",
     # "OXYRIA_NCBI_CHR600004621":"test",

    # "OXYRIA_NCBI_CHR500005515":{
       # "Baseline MG94xREV":0.009652149439954898,
       # "Baseline MG94xREV omega ratio":0.5730726983670508,
       # "Corrected P-value":1.984387007775146e-05,

    # "OXYRIA_NCBI_CHR500005514":{
       # "Baseline MG94xREV":0.01302162728218238,
       # "Baseline MG94xREV omega ratio":1.44760088866469,
       # "Corrected P-value":1.831867990631508e-12,

     # "DOCTH0_CHR100000472":{
       # "Baseline MG94xREV":0.05141586457768321,
       # "Baseline MG94xREV omega ratio":0.1372979321207887,
       # "Corrected P-value":0.04545651860136113,

# find genes that have a P-value less than 0.05

jq -r '
  .["branch attributes"]["0"]
  | to_entries[]
  | select(.value["Corrected P-value"] < 0.05)
  | "\(.key)\t\(.value["Corrected P-value"])"
' OG0015441_tree.txt_unique.nxh.ABSREL.json

# DOCTH0_CHR300006814     0.0004372556584134046
# DOCTH0_CHR400000398     0.02450561644103588
# DOCTH0_CHR400001723     0.04809666165521614
# DOCTH0_CHR500000225     0.004120137030874105
# DOCTH0_CHR800006070     0.02062170869551722
# DOCTH0_CHR900001109     0.0003764289048757696
# NODE10  null
# NODE14  null
# NODE16  null
# NODE4   null
# NODE9   null

# to run for all OG
while read file _; do
  jq -r '
    .["branch attributes"]["0"]
    | to_entries[]
    | select(.value["Corrected P-value"] < 0.05)
    | "'"$file"'\t\(.key)\t\(.value["Corrected P-value"])"
  ' "$file"
done <  ABSREL_nonzero_sorted.txt > all_significant_genes.tsv


###########################
# Look at GO terms for significant genes

# Join all_significant_genes.tsv with interproscan data
# Total_interproscan_output_edited3.tsv .

cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/Gene_ontology

cp /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees/Arctic_trees/all_significant_genes.tsv . 

# Check combined Interproscan file 
wc -l Total_interproscan_output_edited3.tsv
# 109117 Total_interproscan_output_edited3.tsv

awk -F'\t' '$NF != "null"' all_significant_genes.tsv > all_significant_genes_filtered.tsv

wc -l all_significant_genes_filtered.tsv
# 4909 all_significant_genes_filtered.tsv

# OG0004437_tree.txt_unique.nxh.ABSREL.json       DOCTH0_CHR600004047     0.00002416028697815875
# OG0004437_tree.txt_unique.nxh.ABSREL.json       MAKER_LG5_SNAP_GENE_54_112_MRNA_1       0
# OG0004437_tree.txt_unique.nxh.ABSREL.json       OXYRIA_NCBI_CHR300001420        0

# OG0005531_tree.txt_unique.nxh.ABSREL.json       G20420_T1       0.01611781191827699
# OG0005531_tree.txt_unique.nxh.ABSREL.json       OXYRIA_NCBI_CHR500003522        0.03522219318626585
# OG0005531_tree.txt_unique.nxh.ABSREL.json       SNAP_MASKED_LG5_PROCESSED_GENE_102_35_MRNA_1    0.0354711


# OG0004910_tree.txt_unique.nxh.ABSREL.json       DOCTH0_CHR900003875     0.03651363647255934
# OG0004910_tree.txt_unique.nxh.ABSREL.json       MAKER_LG2_SNAP_GENE_64_140_MRNA_1       0.005872506338468753
# OG0004910_tree.txt_unique.nxh.ABSREL.json       OXYRIA_NCBI_CHR700000069        0.04587116189169627

# OG0005895_tree.txt_unique.nxh.ABSREL.json       DOCTH0_CHR900005879     0.03160461924412528
# OG0005895_tree.txt_unique.nxh.ABSREL.json       MAKER_LG5_SNAP_GENE_93_221_MRNA_1       0.03741486493303303
# OG0005895_tree.txt_unique.nxh.ABSREL.json       OXYRIA_NCBI_CHR700003176        0.04260344430396246

#----------------------------
# Check which OG have more than one species


awk '
{
    og = $1
    gene = $2

    if (gene ~ /^OXYRIA/) species="Oxyria"
    else if (gene ~ /^DOCT/) species="Dryas"
    else if (gene ~ /^G[0-9]+_T/) species="Cochgroen"
    else species="Draba"

    key = og SUBSEP species
    seen[key] = 1
}

END {
    for (k in seen) {
        split(k, a, SUBSEP)
        og = a[1]
        species = a[2]

        if (!(og in species_count)) {
            species_count[og] = 0
            species_list[og] = species
        } else if (species_list[og] !~ species) {
            species_list[og] = species_list[og] "," species
        }

        species_count[og]++
    }

    for (og in species_count) {
        if (species_count[og] > 1) {
            print og "\t" species_count[og] "\t" species_list[og]
        }
    }
}
' all_significant_genes_filtered.tsv | sort > OGs_multi_species.tsv


sort -t$'\t' -k2,2nr OGs_multi_species.tsv > OGs_multi_species.sorted.tsv


OG0005531_tree.txt_unique.nxh.ABSREL.json       4       Cochgroen,Draba,Dryas,Oxyria
OG0004358_tree.txt_unique.nxh.ABSREL.json       3       Dryas,Cochgroen,Draba
OG0004437_tree.txt_unique.nxh.ABSREL.json       3       Dryas,Oxyria,Draba
OG0004910_tree.txt_unique.nxh.ABSREL.json       3       Dryas,Draba,Oxyria
OG0005895_tree.txt_unique.nxh.ABSREL.json       3       Dryas,Oxyria,Draba
OG0007028_tree.txt_unique.nxh.ABSREL.json       3       Dryas,Draba,Oxyria
OG0007601_tree.txt_unique.nxh.ABSREL.json       3       Dryas,Cochgroen,Draba
OG0008262_tree.txt_unique.nxh.ABSREL.json       3       Draba,Cochgroen,Dryas
OG0008518_tree.txt_unique.nxh.ABSREL.json       3       Dryas,Cochgroen,Oxyria
OG0008574_tree.txt_unique.nxh.ABSREL.json       3       Oxyria,Draba,Cochgroen
OG0009332_tree.txt_unique.nxh.ABSREL.json       3       Dryas,Cochgroen,Draba
OG0009474_tree.txt_unique.nxh.ABSREL.json       3       Draba,Oxyria,Dryas
OG0009836_tree.txt_unique.nxh.ABSREL.json       3       Oxyria,Cochgroen,Dryas
OG0009890_tree.txt_unique.nxh.ABSREL.json       3       Oxyria,Draba,Dryas
OG0009938_tree.txt_unique.nxh.ABSREL.json       3       Cochgroen,Dryas,Draba
OG0011085_tree.txt_unique.nxh.ABSREL.json       3       Cochgroen,Draba,Oxyria
OG0011114_tree.txt_unique.nxh.ABSREL.json       3       Oxyria,Dryas,Cochgroen
OG0011174_tree.txt_unique.nxh.ABSREL.json       3       Oxyria,Draba,Cochgroen
OG0012177_tree.txt_unique.nxh.ABSREL.json       3       Draba,Cochgroen,Dryas
OG0012361_tree.txt_unique.nxh.ABSREL.json       3       Draba,Oxyria,Dryas
OG0013353_tree.txt_unique.nxh.ABSREL.json       3       Dryas,Cochgroen,Draba
OG0013973_tree.txt_unique.nxh.ABSREL.json       3       Dryas,Cochgroen,Draba
OG0014460_tree.txt_unique.nxh.ABSREL.json       3       Draba,Dryas,Cochgroen


grep OG0005531 InterProscan_ABSREL_sig_genes.tsv

# OG0005531       DOCTH0_CHR200004771     0.00612228594200048     Dryas_octopetala_interproscan_output.tsv        IPR007034,IPR012948,IPR027417,IPR030387,IPR037875,IPR039761  C-terminal, N-terminal,AARP2CN,Bms1/Tsr1-type G domain,P-loop containing nucleoside triphosphate hydrolase,Ribosome biogenesis protein Bms1,Ribosome biogenesis protein Bms1/Tsr1,Ribosome biogenesis protein BMS1/TSR1    GO:0005525,GO:0005634,GO:0042254
# OG0005531       G20420_T1       0.016117811918277       NA      NA      NA      NA
# OG0005531       OXYRIA_NCBI_CHR500003522        0.0352221931862658      Oxyria_digyna_H1_interproscan_output.tsv        IPR007034,IPR012948,IPR027417,IPR030387,IPR037875,IPR039761  C-terminal, N-terminal,AARP2CN,Bms1/Tsr1-type G domain,P-loop containing nucleoside triphosphate hydrolase,Ribosome biogenesis protein Bms1,Ribosome biogenesis protein Bms1/Tsr1,Ribosome biogenesis protein BMS1/TSR1    GO:0005525,GO:0005634,GO:0042254
# OG0005531       SNAP_MASKED_LG5_PROCESSED_GENE_102_35_MRNA_1    0.0354711550003797      NA      NA      NA      NA

# Ribosome biogenesis protein

##########################
# Join files in R 

tmux attach-session -t total 

module load  StdEnv/2020 r/4.2.2
export R_LIBS_USER=/home/celphin/R/x86_64-pc-linux-gnu-library/4.1.0/

R
library(dplyr)
library(tidyr)

#-------------------------
# load GO ont data
# formatted Interproscan to have no duplicates of genes - one row per gene

path="/home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/Gene_ontology/"
Gene_ont_file <- "Total_interproscan_output_edited3.tsv"
gene_ont <- read.delim(paste0(path,"/", Gene_ont_file), header = TRUE, sep = "\t", na.strings = "-", colClasses = c("character", "character", "character", "character"))

nrow(gene_ont)
# [1] 109 116

colnames(gene_ont) <- c("spp", "gene", "INTPRO", "descrip", "GOterm")
length(unique(gene_ont$INTPRO))
# [1] 15 433

#-------------------------------
# Load significant genes list

path="/home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/Gene_ontology/"
sig_genes <- read.delim(paste0(path,"/all_significant_genes_filtered.tsv"), header = FALSE, sep = "\t")

nrow(sig_genes)
# 4909

colnames(sig_genes) <- c("orthogroup", "gene", "p-value")

head(sig_genes)

sig_genes$orthogroup <- sub("_tree\\.txt_unique\\.nxh\\.ABSREL\\.json$", "", sig_genes$orthogroup)

#-------------------------
# Fix capitalization change 

sig_genes <- sig_genes %>% mutate(gene = toupper(gene))
gene_ont <- gene_ont %>% mutate(gene = toupper(gene))

#----------------------------
# Join sig_genes with GO info

merged_data <- sig_genes %>% left_join(gene_ont, by = "gene")

head(merged_data, 10)
colnames(merged_data)

write.table(merged_data, "InterProscan_ABSREL_sig_genes.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# remove descrip column

merged_data_nodescrip <- merged_data[,-6]

head(merged_data_nodescrip, 100)

q()
n 

################################
cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/Gene_ontology/

more InterProscan_ABSREL_sig_genes.tsv # file made below

wc -l InterProscan_ABSREL_sig_genes.tsv

# Oxyria
grep OXYRIA InterProscan_ABSREL_sig_genes.tsv | wc -l 
927

# Dryas
grep DOCT InterProscan_ABSREL_sig_genes.tsv | wc -l 
1726

# Cochgroen
grep G*_T InterProscan_ABSREL_sig_genes.tsv | wc -l 
797

# all the rest - Draba



#----------------
# Look at specific genes

grep OG0004437 InterProscan_ABSREL_sig_genes.tsv
# Histidine kinase/HSP90-like ATPase superfamily,MICRORCHIDIA ATPase family,Morc  GO:0016887
# Epigenetic: heterochromatin formation and transcriptional silencing.

grep OG0005531 InterProscan_ABSREL_sig_genes.tsv
#C-terminal, N-terminal,AARP2CN,Bms1/Tsr1-type G domain,P-loop containing nucleoside triphosphate hydrolase,Ribosome biogenesis protein Bms1,Ribosome biogenesis protein Bms1/Tsr1,Ribosome biogenesis protein BMS1/TSR1  GO:0005525,GO:0005634,GO:0042254
# protein synthesis capacity - growth rates

grep OG0004910 InterProscan_ABSREL_sig_genes.tsv
# family 31,Domain of unknown function DUF4094,Glycosyl transferase        GO:0006486,GO:0016020,GO:0016758
# Add sugar moieties to proteins or lipids

grep OG0005895 InterProscan_ABSREL_sig_genes.tsv
# unknown

# get counts of each unique GO term
cut -f7 InterProscan_ABSREL_sig_genes.tsv | \
grep -v "^NA$" | \
tr ',' '\n' | \
grep -v "^$" | \
sort | \
uniq -c | \
sort -nr > GO_counts.txt

head GO_counts.txt
    238 GO:0005515 # Protein binding
    170 GO:0005524 # ATP binding
    151 GO:0016020 # Membrane
    110 GO:0003676 # Nucleic acid binding
    106 GO:0003677 # DNA binding
     89 GO:0006355 # Regulation of transcription
     83 GO:0006468 # Protein phosphorylation
     83 GO:0004672 # Protein kinase activity
     75 GO:0003723 # RNA binding
     51 GO:0008270 # Zinc ion binding
     48 GO:0055085
     48 GO:0003700
     47 GO:0003824
     45 GO:0016887
     44 GO:0005975
     38 GO:0020037
     35 GO:0005634
     33 GO:0016491
     30 GO:0004523
     26 GO:0005506
     25 GO:0046872
     25 GO:0005509
     24 GO:0016705
     24 GO:0006508
     24 GO:0004553
     21 GO:0009451 # Response to stress
     21 GO:0006952 # Defense response
     21 GO:0004497
     19 GO:0046983
     19 GO:0043565
     19 GO:0022857
     17 GO:0015074
     17 GO:0006886
     16 GO:0003735
     15 GO:0016757
     15 GO:0008017
     15 GO:0006486
     15 GO:0006412
     15 GO:0006364
     15 GO:0005840
     15 GO:0004842
     14 GO:0140359
     14 GO:0016567
     14 GO:0008168
     14 GO:0006281
     13 GO:0016192
     13 GO:0006979 # Response to oxidative stress
     13 GO:0006470
     13 GO:0004601
     13 GO:0004252

# stress, redox, membrane transport, and chromatin-related functions








# Need to compare to frequencies in whole genome - ErmineJ

######################################
# subset Interproscan info by spp - to get null model with all genes

tmux attach-session -t total 

module load  StdEnv/2020 r/4.2.2
export R_LIBS_USER=/home/celphin/R/x86_64-pc-linux-gnu-library/4.1.0/

R
library(dplyr)
library(tidyr)

#-------------------------
# load GO ont data
# formatted Interproscan to have no duplicates of genes - one row per gene

path="/home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/Gene_ontology/"
Gene_ont_file <- "Total_interproscan_output_edited3.tsv"
gene_ont <- read.delim(paste0(path,"/", Gene_ont_file), header = TRUE, sep = "\t", na.strings = "-", colClasses = c("character", "character", "character", "character"))

#-------------------------------------------
#subset Interproscan info by spp - to get null model with all genes

unique(gene_ont$spp)
# [1] "Arabis_alpina_interproscan_output.tsv"
# [2] "Cochlearia_groenlandica_interproscan_output.tsv"
# [3] "Draba_nivalis_interproscan_output.tsv"
# [4] "Dryas_octopetala_interproscan_output.tsv"
# [5] "Oxyria_digyna_H1_interproscan_output.tsv"
# [6] "Rheum_nobile_H0_interproscan_output.tsv"

# Oxyria
gene_ont_Oxydig <- gene_ont[which(gene_ont$spp=="Oxyria_digyna_H1_interproscan_output.tsv"),]
# format for ermineJ
# https://erminej.msl.ubc.ca/help/input-files/
# https://erminej.msl.ubc.ca/help/input-files/gene-annotations/
gene_ont_Oxydig1 <- gene_ont_Oxydig %>%
  mutate(gene2 = gene) %>%
  select(gene, gene2, descrip, GOterm)
# write output for each spp
write.table(gene_ont_Oxydig1, "Oxydig_GO_mappings.ermineJ.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gene_ont_Oxydig, "Oxydig_interproscan_edited.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

#------------------------------------------
# Drabaniv
gene_ont_Drabaniv <- gene_ont[which(gene_ont$spp=="Draba_nivalis_interproscan_output.tsv"),]
# format for ermineJ
# https://erminej.msl.ubc.ca/help/input-files/
# https://erminej.msl.ubc.ca/help/input-files/gene-annotations/
gene_ont_Drabaniv1 <- gene_ont_Drabaniv %>%
  mutate(gene2 = gene) %>%
  select(gene, gene2, descrip, GOterm)
# write output for each spp
write.table(gene_ont_Drabaniv1, "Drabaniv_GO_mappings.ermineJ.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gene_ont_Drabaniv, "Drabaniv_interproscan_edited.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

#------------------------------------------
# Dryasoct
gene_ont_Dryasoct <- gene_ont[which(gene_ont$spp=="Dryas_octopetala_interproscan_output.tsv"),]
# format for ermineJ
# https://erminej.msl.ubc.ca/help/input-files/
# https://erminej.msl.ubc.ca/help/input-files/gene-annotations/
gene_ont_Dryasoct1 <- gene_ont_Dryasoct %>%
  mutate(gene2 = gene) %>%
  select(gene, gene2, descrip, GOterm)
# write output for each spp
write.table(gene_ont_Dryasoct1, "Dryasoct_GO_mappings.ermineJ.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gene_ont_Dryasoct, "Dryasoct_interproscan_edited.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

#------------------------------------------
# Cochgro
gene_ont_Cochgro <- gene_ont[which(gene_ont$spp=="Cochlearia_groenlandica_interproscan_output.tsv"),]
# format for ermineJ
# https://erminej.msl.ubc.ca/help/input-files/
# https://erminej.msl.ubc.ca/help/input-files/gene-annotations/
gene_ont_Cochgro1 <- gene_ont_Cochgro %>%
  mutate(gene2 = gene) %>%
  select(gene, gene2, descrip, GOterm)
# write output for each spp
write.table(gene_ont_Cochgro1, "Cochgro_GO_mappings.ermineJ.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gene_ont_Cochgro, "Cochgro_interproscan_edited.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

q()
n

#------------------------------------
# Need to get list of GO terms that were in original ABSREL test



#--------------------------------
# Prepare gene score file across families (total) - contracted

for taxon in Oxydig Rheumnob Arabalp Dryasoct Drabaniv Cochgro ; \
do awk 'BEGIN{FS="\t"}{print $2,"0"}' ~/scratch/Oxyria/CAFE/enrichment_analysis/"$taxon"_interproscan_edited.tsv | sort -u > "$taxon"_total_totalcontracted_genesets ; done

for taxon in Oxydig Rheumnob Arabalp Dryasoct Drabaniv Cochgro ; \
do cat ~/scratch/Oxyria/CAFE/contracted/"$taxon"_total_totalcontracted_geneIDs.txt | \
sed 's/ .*$//g' | while read gene ; do sed -i "s/$gene 0/$gene 1/g" "$taxon"_total_totalcontracted_genesets ; done ; done

for taxon in Oxydig Rheumnob Arabalp Dryasoct Drabaniv Cochgro ; \
do sed -i 's/ /\t/g' "$taxon"_total_totalcontracted_genesets; done


##############################################
# ermineJ

# https://erminej.msl.ubc.ca/help/tutorials/running-an-analysis-ora/

# As of ErmineJ 3, when using the ‘ORA’ method you have the option to use a simple “hit list” of genes,
# rather than preparing a score file yourself (a “quick list”). Caution: If you use this feature, 
# the “non-hits” will be all the rest of the genes listed in your annotation file. That might not 
# be appropriate if the annotation file includes genes that were not assayed in your experiment. 
# This is most likely to be a problem if your annotation file is a list of all the genes in the genome

# Note I should switch total to just be the orthogroups shared by all 


#######################################################################################################################################
tmux new-session -s Enrichment
tmux attach-session -t Enrichment

salloc -c1 --time 3:00:00 --mem 120000m --account def-rieseber

cd ~/scratch/Oxyria/CAFE/enrichment_analysis

#Expanded : Oxyria
ERMINEJ_HOME=/home/celphin/ermineJ-3.2
export JAVA_HOME=/cvmfs/soft.computecanada.ca/easybuild/software/2020/Core/java/13.0.2/

module load java/13.0.2

#-------------------
# total families combined - rapidly expanded/contracted

for taxon in Oxydig Cochgro Dryasoct Drabaniv ; do $ERMINEJ_HOME/bin/ermineJ.sh \
-a "$taxon"_GO_mappings.ermineJ.txt \
-s "$taxon"_total_expanded_genesets \
-c /home/celphin/ermineJ.data/go.obo \
--genesOut -aspects BCM \
-o "$taxon"_totalfam_rapidly_expanded_genesets.ermine.results -y 5 -b ; done


###########################################