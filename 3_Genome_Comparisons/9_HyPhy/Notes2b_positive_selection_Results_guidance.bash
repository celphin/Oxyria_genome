#########################
# Exploring HyPhy results with guidance
# Jan 2026
############################

# other tutorial: https://bioinformaticsworkbook.org/dataAnalysis/ComparativeGenomics/Finding_Positive_Selection_With_Codeml.html#gsc.tab=0

###############################

# Check
cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees/Arctic_trees

ls *_guidance_unique.nxh | wc -l
ls *_guidance_unique.nxh.ABSREL.json  | wc -l
ls *_guidance_unique.nxh.BUSTED.json  | wc -l
ls *_guidance_unique.nxh.RELAX.json  | wc -l

# 10493
# 9869
# 10481
# 1


more logs/absrel_guidance_57878989_1.out



##########################
# Explore results

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
