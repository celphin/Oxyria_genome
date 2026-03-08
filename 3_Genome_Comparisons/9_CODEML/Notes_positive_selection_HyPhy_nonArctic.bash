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
# Copy trees from Arctic folder back up to main to rerun with non-Arctic

tmux attach-session -t backup1
ls *_tree.txt | wc -l
#14602
ls *_tree.txt_pal2nal.fasta | wc -l
# 3993

cp OG*_tree.txt ..
cp OG*_tree.txt_pal2nal.fasta ..

cd ..
ls *_tree.txt | wc -l
ls *_tree.txt_pal2nal.fasta | wc -l

# 23145
# 23145


############################
# Check gene  names for closest non-Arctic relatives

cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/CDS
grep ">" Arabis_alpina.f*
# >AALP_AAs64319U000100
# >AALP_AAs49271U000100
# >AALP_AAs60671U000100
# >AALP_AAs48102U000100
# >AALP_AAs60109U000100
# >AALP_AAs52579U000100
# >AALP_AAs45409U000100
# >AALP_AAs64185U000100
# >AALP_AAs61717U000100
# >AALP_AAs47510U000100
# >AALP_AAs68684U000100
# >AALP_AAs73943U000100
# >AALP_AAs42227U000100

grep ">" Prunus_persica.f*
# >LOC18768756
# >LOC18767021
# >LOC18767694
# >LOC18766660
# >LOC18767366
# >LOC18768373
# >LOC18768715
# >LOC18766300
# >LOC18767169
# >LOC18767794
# >LOC18767676
# >LOC18767280
# >LOC18766388
# >LOC18768341
# >LOC18768949
# >LOC18767264
# >LOC18768104
# >LOC18767743
# >LOC18766280
# >LOC18766638
# >LOC18766399
# >LOC18766758
# >LOC18767356
# >LOC18768784
# >LOC109950649

grep ">" Rheum_nobile_H0.f*
# >RnoG0007793.1
# >RnoG0007792.1
# >RnoG0007791.1
# >RnoG0007790.1
# >RnoG0007789.1
# >RnoG0007788.1
# >RnoG0007787.1
# >RnoG0007786.1
# >RnoG0024551.1
# >RnoG0024550.1
# >RnoG0024549.1
# >RnoG0024548.1
# >RnoG0024547.1
# >RnoG0024546.1

grep ">" Arabidopsis_lyrata.f*
# >LOC110224666
# >LOC9326461
# >LOC9326463
# >LOC9326464
# >LOC9326465
# >LOC9326466
# >LOC9326467
# >LOC9329189
# >LOC9326468
# >LOC9326469
# >LOC9329190
# >LOC9326470
# >LOC9329192
# >LOC9326471
# >LOC9329193
# >LOC9329194
# >LOC9326472
# >LOC9326473
# >LOC9326474
# >LOC9326475
# >LOC9329195
# >LOC9326476

cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/CDS

grep ">" Arabis_alpina.f* | sed 's/>//' > foreground_genes.txt
grep ">" Prunus_persica.f* | sed 's/>//' >> foreground_genes.txt
grep ">" Rheum_nobile_H0.f* | sed 's/>//' >> foreground_genes.txt
grep ">" Arabidopsis_lyrata.f* | sed 's/>//' >> foreground_genes.txt


mv foreground_genes.txt /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees

cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees

sort -u foreground_genes.txt > foreground_genes_unique.txt

#------------------------------------
# Add {Foreground} label to Arctic branches

# Narval1
tmux new-session -s total
tmux attach-session -t total

cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees

for f in *tree.txt
do
    echo "Processing $f"

    awk '
    NR==FNR {fg[$1]=1; next}

    {
        line=$0
        out=""
        pos=1

        while (match(substr(line,pos), /[^:(),]+:/)) {
            start = pos + RSTART - 1
            end = start + RLENGTH - 2
            gene = substr(line, start, RLENGTH-1)

            out = out substr(line, pos, RSTART-1)

            if (gene in fg)
                out = out gene "{Foreground}:"
            else
                out = out gene ":"

            pos = start + RLENGTH
        }

        out = out substr(line,pos)
        print out
    }
    ' foreground_genes_unique.txt "$f" > "${f}_HyPhy_nonArctic.txt"

done

ls *_HyPhy_nonArctic.txt | wc -l
# 23145

################################
# Remove new trees without Foreground branches
grep -L "{Foreground}" *_HyPhy_nonArctic.txt | xargs rm

ls *_HyPhy_nonArctic.txt | wc -l
# 18779

##################################
# loop through orthogroups that include nonArctic spp 
# run as slurm array

# get file list of non empty files
cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees/

mkdir -p logs

for f in *_tree.txt; do
    msa="${f}_pal2nal.fasta"
    tree="${f}_HyPhy_nonArctic.txt"

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

done > filtered_tree_list_nonArctic.txt

# check count
wc -l filtered_tree_list_nonArctic.txt

# 16401 filtered_tree_list.txt (Arctic)
# 16176 filtered_tree_list_nonArctic.txt

#----------------------
# use nano to import text
cat << EOF > absrel_array_nonArctic.sh
#!/bin/bash
#SBATCH --account=def-henryg
#SBATCH --time=0-12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=32G
#SBATCH --output=logs/absrel_nonArc_%A_%a.out
#SBATCH --error=logs/absrel_nonArc_%A_%a.err

module load StdEnv/2020 gcc/9.3.0 openmpi/4.0.3 hyphy/2.5.49

cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees/nonArctic_trees

CHUNK=20

START=$(( (SLURM_ARRAY_TASK_ID - 1) * CHUNK + 1 ))
END=$(( START + CHUNK - 1 ))

TOTAL=$(wc -l < filtered_tree_list_nonArctic.txt)

if [ $END -gt $TOTAL ]; then
    END=$TOTAL
fi

echo "Processing lines $START to $END"

sed -n "${START},${END}p" filtered_tree_list_nonArctic.txt | while read f
do
    echo "Processing $f"

    msa="${f}_pal2nal.fasta"
    tree="${f}_HyPhy_nonArctic.txt"
    output="${f}_unique_nonArctic.nxh"

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

chmod +x absrel_array_nonArctic.sh
dos2unix absrel_array_nonArctic.sh

cp Arctic_trees/remove-duplicates.bf .

sbatch --array=1-809%100 absrel_array_nonArctic.sh
# Submitted batch job 57490810

###############################
# Checking results
cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees

ls *_unique_nonArctic.nxh |wc -l
# 16176
ls *_unique_nonArctic.nxh.ABSREL.json | wc -l 
# 16044

##########################
# Which orthogroups are demonstrating positive selection

# look for found **1** branches under selection among **39** tested

for f in *_unique_nonArctic.nxh.ABSREL.json; do
  jq -r --arg file "$f" \
    '$file + "\t" + (.["test results"]["positive test results"]|tostring) + "\t" + (.["test results"].tested|tostring)' "$f"
done >> Total_ABSREL_results_nonArctic.txt

wc -l Total_ABSREL_results_nonArctic.txt
#16254 Total_ABSREL_results.txt (Arctic)
# 15975 Total_ABSREL_results_nonArctic.txt

# filter for only those with positive selection
awk '$2 != 0' Total_ABSREL_results_nonArctic.txt | sort -k2,2nr > ABSREL_nonzero_sorted_nonArctic.txt

wc -l ABSREL_nonzero_sorted_nonArctic.txt
# 4076 ABSREL_nonzero_sorted.txt
# 2696 ABSREL_nonzero_sorted_nonArctic.txt

more ABSREL_nonzero_sorted_nonArctic.txt

mv ABSREL_nonzero_sorted_nonArctic.txt ABSREL_nonArctic_nonzero_sorted.txt

#########################
# Fix names for nonArctic below here
#---------------------
# Check for paralogs

echo -e "file\tgene\tcorrected_p\tfull_adaptive_model\tnonsyn_subs\tsyn_subs" > all_genes_nonArctic.tsv

while read file _; do
  jq -r '
    .["branch attributes"]["0"]
    | to_entries[]
    | [
        "'"$file"'",
        .key,
        (.value["Corrected P-value"] // "NA"),
        (.value["Full adaptive model"] // "NA"),
        (.value["Full adaptive model (non-synonymous subs/site)"] // "NA"),
        (.value["Full adaptive model (synonymous subs/site)"] // "NA")
      ] | @tsv
  ' "$file"
done < ABSREL_nonArctic_nonzero_sorted.txt >> all_genes_nonArctic.tsv

# Find only the orthogroups with a p-value greater than 0.05 and 
# only one gene from that species in that group

awk -F'\t' '
BEGIN {
    OFS="\t"
    # Define species
    species_list[1]="Rhuem"
    species_list[2]="LOC"
    species_list[3]="Arabis"
    species_list[4]="Other"
}
NR==1 {
    # Save header
    header = $0
    next
}
{
    # Assign species based on gene ID
    if ($2 ~ /^RNOG/) species="Rhuem"
    else if ($2 ~ /^LOC/) species="LOC"
    else if ($2 ~ /^AALP_/) species="Arabis"
    else species="Other"

    # Count all genes per OG per species (for exclusion)
    all_count[$1,species]++

    # Only keep significant genes for output
    if ($3 < 0.05) {
        genes[$1] = genes[$1] $0 "\n"
        species_per_OG[$1][species]=1
    }
}
END {
    single_copy_file="OG_single_copy_significant_nonArctic.tsv"
    gene_names_file="OG_single_copy_significant_gene_names_nonArctic.txt"
    uniqueOG_file="OG_single_copy_significant_uniqueOG_nonArctic.txt"
    species_counts_file="OG_single_copy_species_counts_nonArctic.txt"

    # Header for single-copy file
    print header > single_copy_file

    # Header for species count summary
    print "OG", "Num_species_with_sig_gene", "Species" > species_counts_file

    for (og in genes) {
        exclude=0
        # Exclude OGs if any species has >1 gene (total genes, not only significant)
        for (i=1;i<=4;i++) {
            s=species_list[i]
            if (all_count[og,s]>1) {exclude=1; break}
        }

        if (!exclude) {
            # Save OG lines
            printf "%s", genes[og] >> single_copy_file

            # Save gene names
            split(genes[og], lines, "\n")
            for (i in lines) if (lines[i] != "") {
                split(lines[i], fields, FS)
                print fields[2] >> gene_names_file
            }

            # Unique OG list
            print og >> uniqueOG_file

            # Species counts summary (fixed order)
            n=0; species_list_str=""
            for (j=1;j<=4;j++) {
                s=species_list[j]
                if (species_per_OG[og][s]) {
                    n++
                    species_list_str = species_list_str ? species_list_str "," s : s
                }
            }
            print og, n, species_list_str >> species_counts_file
        }
    }

    # Sort species counts descending by number of species
    system("head -n 1 " species_counts_file " > tmp_header && tail -n +2 " species_counts_file " | sort -k2,2nr >> tmp_header && mv tmp_header \
	OG_single_copy_species_counts_sorted_nonArctic.tsv")
}
' all_genes_nonArctic.tsv

more OG_single_copy_species_counts_sorted_nonArctic.tsv
# OG      Num_species_with_sig_gene       Species
# OG0018292_tree.txt_unique_nonArctic.nxh.ABSREL.json     1       Rhuem
# OG0018370_tree.txt_unique_nonArctic.nxh.ABSREL.json     1       LOC
# OG0019695_tree.txt_unique_nonArctic.nxh.ABSREL.json     1       LOC
# OG0022997_tree.txt_unique_nonArctic.nxh.ABSREL.json     1       Rhuem


#---------------------------------
awk -F'\t' '
BEGIN {
    OFS="\t"
    # Define species
    species_list[1]="Rhuem"
    species_list[2]="LOC"
    species_list[3]="Arabis"
    species_list[4]="Other"
}
NR==1 {
    # Save header
    header = $0
    next
}
{
    # Only keep significant genes
    if ($3 >= 0.05) next

    # Assign species based on gene ID
    if ($2 ~ /^RNOG/) species="Rhuem"
    else if ($2 ~ /^LOC/) species="LOC"
    else if ($2 ~ /^AALP_/) species="Arabis"
    else species="Other"

    # Count genes per OG per species
    count[$1,species]++
    # Store line per OG
    genes[$1] = genes[$1] $0 "\n"
    # Track species per OG
    species_per_OG[$1][species]=1
}
END {
    # Output single-copy significant genes
    single_copy_file="OG_paralogs_significant_nonArctic.tsv"
    gene_names_file="OG_paralogs_significant_gene_names_nonArctic.txt"
    uniqueOG_file="OG_paralogs_significant_uniqueOG_nonArctic.txt"
    species_counts_file="OG_paralogs_species_counts_nonArctic.txt"

    # Header for single-copy file
    print header > single_copy_file

    # Header for species count summary
    print "OG", "Num_species_with_sig_gene", "Species" > species_counts_file

    for (og in genes) {
        exclude=0
        # Check if any species has >1 gene in this OG
        for (i=1;i<=4;i++) {
            s=species_list[i]
            if (count[og,s]>1) {exclude=1; break}
        }
        if (!exclude) {
            # Save OG lines
            printf "%s", genes[og] >> single_copy_file

            # Save gene names
            split(genes[og], lines, "\n")
            for (i in lines) if (lines[i] != "") {
                split(lines[i], fields, FS)
                print fields[2] >> gene_names_file
            }

            # Unique OG list
            print og >> uniqueOG_file

            # Species counts summary
            n=0; species_list_str=""
            for (s in species_per_OG[og]) {
                n++
                species_list_str = species_list_str ? species_list_str "," s : s
            }
            print og, n, species_list_str >> species_counts_file
        }
    }

    # Sort species counts descending by number of species
    system("head -n 1 " species_counts_file " > tmp_header && tail -n +2 " species_counts_file " | sort -k2,2nr >> tmp_header && \
	mv tmp_header OG_paralogs_species_counts_sorted_nonArctic.tsv")
}
' all_genes_nonArctic.tsv

more OG_paralogs_species_counts_sorted_nonArctic.tsv

# OG      Num_species_with_sig_gene       Species
# OG0007401_tree.txt_unique_nonArctic.nxh.ABSREL.json     3       Rhuem,LOC,Arabis
# OG0013935_tree.txt_unique_nonArctic.nxh.ABSREL.json     3       Rhuem,LOC,Arabis
# OG0000510_tree.txt_unique_nonArctic.nxh.ABSREL.json     2       Rhuem,LOC
# OG0000579_tree.txt_unique_nonArctic.nxh.ABSREL.json     2       Rhuem,Arabis
# OG0000847_tree.txt_unique_nonArctic.nxh.ABSREL.json     2       LOC,Arabis

################################
# Get counts of each value

# total XX OGs tested
awk '{count[$2]++} END {for (val in count) print val, count[val]}' \
 OG_paralogs_species_counts_sorted_nonArctic.tsv
# 1 2358
# 2 114
# 3 2

awk '{count[$2]++} END {for (val in count) print val, count[val]}' \
 OG_single_copy_species_counts_sorted_nonArctic.tsv
# 1 4

# Because not splitting LOC genes so all look like paralogs

###################
# Per species counts

grep -c 'LOC' OG_single_copy_species_counts_sorted_nonArctic.tsv
grep -c 'Rheum' OG_single_copy_species_counts_sorted_nonArctic.tsv
grep -c 'Arabis' OG_single_copy_species_counts_sorted_nonArctic.tsv
grep -c 'Other' OG_single_copy_species_counts_sorted_nonArctic.tsv
# 2
# 0
# 0
# 0

grep -c 'LOC' OG_paralogs_species_counts_sorted_nonArctic.tsv
grep -c 'Rhuem' OG_paralogs_species_counts_sorted_nonArctic.tsv
grep -c 'Arabis' OG_paralogs_species_counts_sorted_nonArctic.tsv
grep -c 'Other' OG_paralogs_species_counts_sorted_nonArctic.tsv
# 1024
# 1088
# 480
# 0

# ~ 500 genes per species

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