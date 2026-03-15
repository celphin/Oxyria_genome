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



################################
# Rerun HyPhy with prank alignments

# get file list of non empty files
cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees/Arctic_trees

mkdir -p logs

for f in *_tree.txt; do
    msa="${f}_guidance_pal2nal.fasta"
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

done > filtered_guidance_tree_list.txt

# check count
wc -l filtered_guidance_tree_list.txt
# 4053 

##############################
# Test ABSREL, RELAX and BUSTED

module load StdEnv/2020 gcc/9.3.0 openmpi/4.0.3 hyphy/2.5.49

cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees/Arctic_trees

msa="OG0022767_guidance_pal2nal.fasta"
tree="OG0022767_tree.txt_HyPhy.txt"
output="OG0022767_guidance_unique.nxh"

hyphy remove-duplicates.bf \
	--msa "$msa" \
	--tree "$tree" \
	--output "$output"

hyphy absrel \
	--alignment "$output" \
	--branches FOREGROUND

# Must be more than one branch for these tests?
hyphy relax \
	--alignment "$output"
	--branches FOREGROUND

hyphy busted \
	--alignment "$output"
	--branches FOREGROUND

#-------------------------
# use nano to import text
cat << 'EOF' > absrel_guidance_array.sh
#!/bin/bash
#SBATCH --account=def-henryg
#SBATCH --time=0-6:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=32G
#SBATCH --output=logs/absrel_guidance_%A_%a.out
#SBATCH --error=logs/absrel_guidance_%A_%a.err

module load StdEnv/2020 gcc/9.3.0 openmpi/4.0.3 hyphy/2.5.49

cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees/Arctic_trees

CHUNK=10

START=$(( (SLURM_ARRAY_TASK_ID - 1) * CHUNK + 1 ))
END=$(( START + CHUNK - 1 ))

TOTAL=$(wc -l < filtered_guidance_tree_list.txt)

if [ $END -gt $TOTAL ]; then
    END=$TOTAL
fi

echo "Processing lines $START to $END"

sed -n "${START},${END}p" filtered_guidance_tree_list.txt | while read f
do
    echo "Processing $f"

    msa="${f}_guidance_pal2nal.fasta"
    tree="${f}_HyPhy.txt"
    output="${f}_guidance_unique.nxh"

    if [[ -f "$output" ]]; then
        echo "Output exists for $f — skipping"
        continue
    fi

    # Remove duplicate sequences
    hyphy remove-duplicates.bf \
        --msa "$msa" \
        --tree "$tree" \
        --output "$output"

    # Run ABSREL
    hyphy absrel \
        --alignment "$output" \
        --branches FOREGROUND

    # Count number of foreground branches
    FG_COUNT=$(grep -o "Foreground" "$tree" | wc -l)

    # Count number of sequences
    SEQ_COUNT=$(grep -c "^>" "$output")

    # Count Arctic species from tree
    SPECIES_COUNT=$(sed 's/[(),]/\n/g' "$tree" | grep -v ':' | grep -v '^$' | awk '
    {
    gene=$1
    if (gene ~ /^OXYRIA/) species="Oxyria"
    else if (gene ~ /^DOCT/) species="Dryas"
    else if (gene ~ /^G[0-9]+_T/) species="Cochgroen"
    else species="Draba"
    print species
    }' | sort -u | wc -l)

    # Run RELAX only if criteria met
    # >=2 foreground branches
    # >=2 Arctic species in the tree
    # >=5 sequences in the alignment
    if [[ "$FG_COUNT" -ge 2 && "$SPECIES_COUNT" -ge 2 && "$SEQ_COUNT" -ge 5 ]]; then
        echo "Running RELAX (FG=$FG_COUNT, species=$SPECIES_COUNT, seqs=$SEQ_COUNT)"
        hyphy relax \
            --alignment "$output" \
            --branches FOREGROUND
    else
        echo "Skipping RELAX (FG=$FG_COUNT, species=$SPECIES_COUNT, seqs=$SEQ_COUNT)"
    fi

    # Run BUSTED
    hyphy busted \
        --alignment "$output" \
        --branches FOREGROUND

    echo "Finished $f"
done

EOF

chmod +x absrel_guidance_array.sh
dos2unix absrel_guidance_array.sh

sbatch --array=1-406%100 absrel_guidance_array.sh


#--------------------
# Check
cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees/Arctic_trees

ls *_guidance_unique.nxh | wc -l
ls *_guidance_unique.nxh.ABSREL.json  | wc -l


#------------------------
for f in *_guidance_unique.nxh.ABSREL.json; do
  jq -r --arg file "$f" \
    '$file + "\t" + (.["test results"]["positive test results"]|tostring) + "\t" + (.["test results"].tested|tostring)' "$f"
done >> Total_ABSREL_guidance_results.txt


# filter for only those with positive selection
awk '$2 != 0' Total_ABSREL_guidance_results.txt | sort -k2,2nr > ABSREL_guidance_nonzero_sorted.txt

# get gene names
while read file _; do
  jq -r '
    .["branch attributes"]["0"]
    | to_entries[]
    | select(.value["Corrected P-value"] < 0.05)
    | "'"$file"'\t\(.key)\t\(.value["Corrected P-value"])"
  ' "$file"
done <  ABSREL_guidance_nonzero_sorted.txt > all_significant_genes_guidance.tsv

#------------------
# Link to gene functions
cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/Gene_ontology

cp /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees/Arctic_trees/all_significant_genes_guidance.tsv . 

awk -F'\t' '$NF != "null"' all_significant_genes_guidance.tsv > all_significant_genes_guidance_filtered.tsv


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
' all_significant_genes_guidance_filtered.tsv | sort > OGs_multi_species_guidance.tsv

sort -t$'\t' -k2,2nr OGs_multi_species_guidance.tsv > OGs_multi_species_guidance.sorted.tsv

# OG0005531_tree.txt_guidance_unique.nxh.ABSREL.json 3       Cochgroen,Oxyria,Draba
# OG0011114_tree.txt_guidance_unique.nxh.ABSREL.json 3       Oxyria,Dryas,Cochgroen
# OG0011456_tree.txt_guidance_unique.nxh.ABSREL.json 3       Dryas,Oxyria,Draba
# OG0012809_tree.txt_guidance_unique.nxh.ABSREL.json 3       Dryas,Draba,Cochgroen

#---------------------
# Check for paralogs

echo -e "file\tgene\tcorrected_p\tfull_adaptive_model\tnonsyn_subs\tsyn_subs" > all_genes_guidance.tsv

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
done < ABSREL_guidance_nonzero_sorted.txt >> all_genes_guidance.tsv

# Find only the orthogroups with a p-value greater than 0.05 and 
# only one gene from that species in that group

awk -F'\t' '
BEGIN {
    OFS="\t"
    # Define species
    species_list[1]="Dryas"
    species_list[2]="Oxyria"
    species_list[3]="Draba"
    species_list[4]="Coch"
}
NR==1 {
    # Save header
    header = $0
    next
}
{
    # Assign species based on gene ID
    if ($2 ~ /^DOCTH0_CHR/) species="Dryas"
    else if ($2 ~ /^OXYRIA_NCBI_CHR/) species="Oxyria"
    else if ($2 ~ /^MAKER_/ || $2 ~ /^SNAP_MASKED_/) species="Draba"
    else if ($2 ~ /^G[0-9]+_T[0-9]+/) species="Coch"
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
    single_copy_file="OG_single_copy_significant_guidance.tsv"
    gene_names_file="OG_single_copy_significant_gene_names_guidance.txt"
    uniqueOG_file="OG_single_copy_significant_uniqueOG_guidance.txt"
    species_counts_file="OG_single_copy_species_counts_guidance.txt"

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
	OG_single_copy_species_counts_sorted_guidance.tsv")
}
' all_genes_guidance.tsv

more OG_single_copy_species_counts_sorted_guidance.tsv

# OG      Num_species_with_sig_gene       Species
# OG0011114_tree.txt_guidance_unique.nxh.ABSREL.json 3       Dryas,Oxyria,Coch
# OG0011456_tree.txt_guidance_unique.nxh.ABSREL.json 3       Dryas,Oxyria,Draba
# OG0005895_tree.txt_guidance_unique.nxh.ABSREL.json 2       Dryas,Oxyria
# OG0007572_tree.txt_guidance_unique.nxh.ABSREL.json 2       Dryas,Coch
# OG0007884_tree.txt_guidance_unique.nxh.ABSREL.json 2       Dryas,Oxyria
# OG0008152_tree.txt_guidance_unique.nxh.ABSREL.json 2       Draba,Coch
# OG0008378_tree.txt_guidance_unique.nxh.ABSREL.json 2       Dryas,Draba
# OG0008499_tree.txt_guidance_unique.nxh.ABSREL.json 2       Dryas,Draba

#---------------------------------
awk -F'\t' '
BEGIN {
    OFS="\t"
    # Define species
    species_list[1]="Dryas"
    species_list[2]="Oxyria"
    species_list[3]="Draba"
    species_list[4]="Coch"
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
    if ($2 ~ /^DOCTH0_CHR/) species="Dryas"
    else if ($2 ~ /^OXYRIA_NCBI_CHR/) species="Oxyria"
    else if ($2 ~ /^MAKER_/ || $2 ~ /^SNAP_MASKED_/) species="Draba"
    else if ($2 ~ /^G[0-9]+_T[0-9]+/) species="Coch"
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
    single_copy_file="OG_paralogs_significant_guidance.tsv"
    gene_names_file="OG_paralogs_significant_gene_names_guidance.txt"
    uniqueOG_file="OG_paralogs_significant_uniqueOG_guidance.txt"
    species_counts_file="OG_paralogs_species_counts_guidance.txt"

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
	mv tmp_header OG_paralogs_species_counts_sorted_guidance.tsv")
}
' all_genes_guidance.tsv

more OG_paralogs_species_counts_sorted_guidance.tsv


# OG      Num_species_with_sig_gene       Species
# OG0005531_tree.txt_guidance_unique.nxh.ABSREL.json 3       Coch,Oxyria,Draba
# OG0011114_tree.txt_guidance_unique.nxh.ABSREL.json 3       Coch,Dryas,Oxyria
# OG0011456_tree.txt_guidance_unique.nxh.ABSREL.json 3       Dryas,Oxyria,Draba
# OG0012809_tree.txt_guidance_unique.nxh.ABSREL.json 3       Coch,Dryas,Draba
# OG0000886_tree.txt_guidance_unique.nxh.ABSREL.json 2       Coch,Dryas
# OG0000896_tree.txt_guidance_unique.nxh.ABSREL.json 2       Coch,Draba
# OG0001197_tree.txt_guidance_unique.nxh.ABSREL.json 2       Oxyria,Draba
# OG0001403_tree.txt_guidance_unique.nxh.ABSREL.json 2       Coch,Oxyria
# OG0001639_tree.txt_guidance_unique.nxh.ABSREL.json 2       Dryas,Oxyria
# OG0002655_tree.txt_guidance_unique.nxh.ABSREL.json 2       Coch,Dryas
# OG0002705_tree.txt_guidance_unique.nxh.ABSREL.json 2       Dryas,Oxyria
# OG0002842_tree.txt_guidance_unique.nxh.ABSREL.json 2       Dryas,Draba
# OG0003138_tree.txt_guidance_unique.nxh.ABSREL.json 2       Dryas,Draba
# OG0003752_tree.txt_guidance_unique.nxh.ABSREL.json 2       Dryas,Draba
# OG0003915_tree.txt_guidance_unique.nxh.ABSREL.json 2       Dryas,Draba
# OG0003943_tree.txt_guidance_unique.nxh.ABSREL.json 2       Dryas,Oxyria
# OG0003944_tree.txt_guidance_unique.nxh.ABSREL.json 2       Oxyria,Draba
# OG0003965_tree.txt_guidance_unique.nxh.ABSREL.json 2       Coch,Dryas


grep OG0011456 OG_paralogs_species_counts_sorted.tsv
#OG0011456_tree.txt_unique.nxh.ABSREL.json       1       Oxyria

grep OG0011456 OG_single_copy_species_counts_sorted.tsv
#OG0011456_tree.txt_unique.nxh.ABSREL.json       1       Oxyria

grep OG0012809  OG_paralogs_species_counts_sorted.tsv
# OG0012809_tree.txt_unique.nxh.ABSREL.json       1       Dryas

#-------------------------
# Get counts of each value

# total 16401 OGs tested
awk '{count[$2]++} END {for (val in count) print val, count[val]}' \
 OG_paralogs_species_counts_sorted.tsv
# 1 3382
# 2 339
# 3 21
# 4 1

awk '{count[$2]++} END {for (val in count) print val, count[val]}' \
 OG_single_copy_species_counts_sorted.tsv
# 1 1501
# 2 149
# 3 13


# total 4053 OG tested
awk '{count[$2]++} END {for (val in count) print val, count[val]}' \
 OG_paralogs_species_counts_sorted_guidance.tsv
# 1 2007
# 2 144
# 3 4

awk '{count[$2]++} END {for (val in count) print val, count[val]}' \
 OG_single_copy_species_counts_sorted_guidance.tsv
# 1 872
# 2 61
# 3 2

###################
# Per species counts

grep -c 'Coch' OG_single_copy_species_counts_sorted_guidance.tsv
grep -c 'Dryas' OG_single_copy_species_counts_sorted_guidance.tsv
grep -c 'Oxyria' OG_single_copy_species_counts_sorted_guidance.tsv
grep -c 'Draba' OG_single_copy_species_counts_sorted_guidance.tsv
# 128
# 437
# 172
# 263

grep -c 'Coch' OG_paralogs_species_counts_sorted_guidance.tsv
grep -c 'Dryas' OG_paralogs_species_counts_sorted_guidance.tsv
grep -c 'Oxyria' OG_paralogs_species_counts_sorted_guidance.tsv
grep -c 'Draba' OG_paralogs_species_counts_sorted_guidance.tsv
# 320
# 842
# 426
# 719

grep -c 'Coch' OG_single_copy_species_counts_sorted.tsv
grep -c 'Dryas' OG_single_copy_species_counts_sorted.tsv
grep -c 'Oxyria' OG_single_copy_species_counts_sorted.tsv
grep -c 'Draba' OG_single_copy_species_counts_sorted.tsv
# 280
# 748
# 314
# 496

grep -c 'Coch' OG_paralogs_species_counts_sorted.tsv
grep -c 'Dryas' OG_paralogs_species_counts_sorted.tsv
grep -c 'Oxyria' OG_paralogs_species_counts_sorted.tsv
grep -c 'Draba' OG_paralogs_species_counts_sorted.tsv
# 672
# 1449
# 751
# 1255


#---------------------
# random chance? - yes

# 3382×1+339×2+21×3+1×4=4127
# 16401×4=65604
# p=4127/65604≈0.063

# species	observed	expected
# 0	12658	12629
# 1	3382	3395
# 2	339	345
# 3	21	17
# 4	1	0.27

# p≈0.37

# p>0.05

#-------------------
# 2007×1+144×2+4×3=2007+288+12=2307
# 4053×4=16212
# p=2307/16212≈0.142

# species	observed	expected
# 1	2007	3185
# 2	144	790
# 3	4	89
# 4	0	4

# Expect more overlap by random chance
# p<10−16

# PRANK is removing alignment-driven false positives


##########################
