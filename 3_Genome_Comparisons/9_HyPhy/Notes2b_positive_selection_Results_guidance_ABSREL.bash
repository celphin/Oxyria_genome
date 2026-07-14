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
# 8613


more logs/absrel_guidance_57878989_1.out

##########################
# Explore results

# Narval1
tmux attach-session -t total

for f in *_guidance_unique.nxh.ABSREL.json; do
  jq -r --arg file "$f" \
    '$file + "\t" + (.["test results"]["positive test results"]|tostring) + "\t" + (.["test results"].tested|tostring)' "$f"
done >> Total_ABSREL_guidance_results.txt

wc -l Total_ABSREL_guidance_results.txt
# 9832

# filter for only those with positive selection
awk '$2 != 0' Total_ABSREL_guidance_results.txt | sort -k2,2nr > ABSREL_guidance_nonzero_sorted.txt

wc -l ABSREL_guidance_nonzero_sorted.txt
# 6876 ABSREL_guidance_nonzero_sorted.txt

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

more OGs_multi_species_guidance.sorted.tsv

wc -l OGs_multi_species_guidance.sorted.tsv

awk '{count[$2]++} END {print "2:",count[2]; print "3:",count[3]; print "4:",count[4]}' OGs_multi_species_guidance.sorted.tsv

# Counts not accounting for duplicate genes in one species 
# 2: 1983
# 3: 1170
# 4: 449

#---------------------
# Check for paralogs

cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees/Arctic_trees

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
# OG0003941_guidance_unique.nxh.ABSREL.json       4       Dryas,Oxyria,Draba,Coch
# OG0004962_guidance_unique.nxh.ABSREL.json       4       Dryas,Oxyria,Draba,Coch
# OG0006045_guidance_unique.nxh.ABSREL.json       4       Dryas,Oxyria,Draba,Coch
# OG0006379_guidance_unique.nxh.ABSREL.json       4       Dryas,Oxyria,Draba,Coch
# OG0006531_guidance_unique.nxh.ABSREL.json       4       Dryas,Oxyria,Draba,Coch
# OG0007066_guidance_unique.nxh.ABSREL.json       4       Dryas,Oxyria,Draba,Coch
# OG0007141_guidance_unique.nxh.ABSREL.json       4       Dryas,Oxyria,Draba,Coch
# OG0007200_guidance_unique.nxh.ABSREL.json       4       Dryas,Oxyria,Draba,Coch
# OG0007238_guidance_unique.nxh.ABSREL.json       4       Dryas,Oxyria,Draba,Coch
# OG0007280_guidance_unique.nxh.ABSREL.json       4       Dryas,Oxyria,Draba,Coch
# OG0007461_guidance_unique.nxh.ABSREL.json       4       Dryas,Oxyria,Draba,Coch
# OG0007468_guidance_unique.nxh.ABSREL.json       4       Dryas,Oxyria,Draba,Coch
# OG0007469_guidance_unique.nxh.ABSREL.json       4       Dryas,Oxyria,Draba,Coch
# OG0007598_guidance_unique.nxh.ABSREL.json       4       Dryas,Oxyria,Draba,Coch
# OG0007891_guidance_unique.nxh.ABSREL.json       4       Dryas,Oxyria,Draba,Coch
# OG0007989_guidance_unique.nxh.ABSREL.json       4       Dryas,Oxyria,Draba,Coch


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
# OG0000661_guidance_unique.nxh.ABSREL.json       4       Coch,Dryas,Oxyria,Draba
# OG0000852_guidance_unique.nxh.ABSREL.json       4       Coch,Dryas,Oxyria,Draba
# OG0000947_guidance_unique.nxh.ABSREL.json       4       Coch,Dryas,Oxyria,Draba
# OG0001359_guidance_unique.nxh.ABSREL.json       4       Coch,Dryas,Oxyria,Draba
# OG0001642_guidance_unique.nxh.ABSREL.json       4       Coch,Dryas,Oxyria,Draba
# OG0001643_guidance_unique.nxh.ABSREL.json       4       Coch,Dryas,Oxyria,Draba
# OG0001648_guidance_unique.nxh.ABSREL.json       4       Coch,Dryas,Oxyria,Draba
# OG0001725_guidance_unique.nxh.ABSREL.json       4       Coch,Dryas,Oxyria,Draba
# OG0001787_guidance_unique.nxh.ABSREL.json       4       Coch,Dryas,Oxyria,Draba
# OG0001800_guidance_unique.nxh.ABSREL.json       4       Coch,Dryas,Oxyria,Draba
# OG0001831_guidance_unique.nxh.ABSREL.json       4       Coch,Dryas,Oxyria,Draba
# OG0002045_guidance_unique.nxh.ABSREL.json       4       Coch,Dryas,Oxyria,Draba
# OG0002344_guidance_unique.nxh.ABSREL.json       4       Coch,Dryas,Oxyria,Draba
# OG0002354_guidance_unique.nxh.ABSREL.json       4       Coch,Dryas,Oxyria,Draba
# OG0002370_guidance_unique.nxh.ABSREL.json       4       Coch,Dryas,Oxyria,Draba

#-------------------------
# Get counts of each value

awk '{count[$2]++} END {for (val in count) print val, count[val]}' \
 OG_paralogs_species_counts_sorted_guidance.tsv

# 1 3004
# 2 1573
# 3 768
# 4 214

awk '{count[$2]++} END {for (val in count) print val, count[val]}' \
 OG_single_copy_species_counts_sorted_guidance.tsv

# 1 1864
# 2 850
# 3 399
# 4 122


wc -l OG_single_copy_species_counts_sorted_guidance.tsv
wc -l OG_paralogs_species_counts_sorted_guidance.tsv
#3236 OG_single_copy_species_counts_sorted_guidance.tsv
#5560 OG_paralogs_species_counts_sorted_guidance.tsv

#---------------------------------

#MAFFT
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

# Only MAFFT sig genes tested
awk '{count[$2]++} END {for (val in count) print val, count[val]}' \
 OG_paralogs_species_counts_sorted_prank.tsv
# 1 2007
# 2 144
# 3 4

awk '{count[$2]++} END {for (val in count) print val, count[val]}' \
 OG_single_copy_species_counts_sorted_prank.tsv
# 1 872
# 2 61
# 3 2

########################
# Get list of OGs in all 4 spp
cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees/Arctic_trees

# Extract full rows with all 4 species
awk '{split($3,species,","); if(length(species)==4) print $0}' \
OG_paralogs_species_counts_sorted_guidance.tsv \
> OG_4species_rows.tsv

# Extract full rows with all 4 species
awk '{split($3,species,","); if(length(species)==4) print $0}' \
OG_single_copy_species_counts_sorted_guidance.tsv \
> single_copy_4species_rows.tsv


###################
# Per species counts

grep -c 'Coch' OG_single_copy_species_counts_sorted_guidance.tsv
grep -c 'Dryas' OG_single_copy_species_counts_sorted_guidance.tsv
grep -c 'Oxyria' OG_single_copy_species_counts_sorted_guidance.tsv
grep -c 'Draba' OG_single_copy_species_counts_sorted_guidance.tsv

# 1112
# 1583
# 1140
# 1414

#~1500 genes per spp

grep -c 'Coch' OG_paralogs_species_counts_sorted_guidance.tsv
grep -c 'Dryas' OG_paralogs_species_counts_sorted_guidance.tsv
grep -c 'Oxyria' OG_paralogs_species_counts_sorted_guidance.tsv
grep -c 'Draba' OG_paralogs_species_counts_sorted_guidance.tsv

# 2093
# 2619
# 2134
# 2464

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


###################################

############################
# Run again with spp lists
cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/CDS
grep ">" Draba_nivalis.fa | sed 's/>//' | awk '{print $0 "\tDraba"}' > gene_Arctic_species_map.txt
grep ">" Cochlearia_groenlandica.fa | sed 's/>//' | awk '{print $0 "\tCochlearia"}' >> gene_Arctic_species_map.txt
grep ">" Dryas_octopetala.fa | sed 's/>//' | awk '{print $0 "\tDryas"}' >> gene_Arctic_species_map.txt
grep ">" Oxyria_digyna_H1.fa | sed 's/>//' | awk '{print $0 "\tOxyria"}' >> gene_Arctic_species_map.txt

cp gene_Arctic_species_map.txt /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees/Arctic_trees

cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees/Arctic_trees

# Fix for odd gene names
awk -F'\t' '
BEGIN {
    OFS="\t"
}

# Normalize gene IDs (fix case + suffix differences)
function normalize(g) {
    g = toupper(g)
    sub(/_[0-9]+$/, ".1", g)
    return g
}

# Load gene→species map (FIRST file)
FNR==NR {
    key = normalize($1)
    gene2species[key] = $2
    next
}

# Save header
NR==1 {
    header = $0
    next
}

{
    # Only keep significant genes
    if ($3 >= 0.05) next

    gene_raw = $2
    gene = normalize(gene_raw)

    # Skip + log unmapped genes (only once)
    if (!(gene in gene2species)) {
        if (!(gene in seen_unmapped)) {
            print gene_raw >> "unmapped_genes.txt"
            seen_unmapped[gene] = 1
        }
        next
    }

    species = gene2species[gene]

    # Count genes per OG per species
    count[$1,species]++

    # Store full line per OG
    genes[$1] = genes[$1] $0 "\n"

    # Track species per OG
    species_per_OG[$1][species] = 1
}

END {
    single_copy_file="OG_paralogs_significant_guidance_nonArctic.tsv"
    gene_names_file="OG_paralogs_significant_gene_names_guidance_nonArctic.txt"
    uniqueOG_file="OG_paralogs_significant_uniqueOG_guidance_nonArctic.txt"
    species_counts_file="OG_paralogs_species_counts_guidance_nonArctic.txt"

    # Headers
    print header > single_copy_file
    print "OG", "Num_species_with_sig_gene", "Species" > species_counts_file

    for (og in genes) {
        exclude = 0

        # Exclude OGs where any species has >1 gene
        for (s in species_per_OG[og]) {
            if (count[og,s] > 1) {
                exclude = 1
                break
            }
        }

        if (!exclude) {
            # Write OG lines
            printf "%s", genes[og] >> single_copy_file

            # Extract gene names
            split(genes[og], lines, "\n")
            for (i in lines) {
                if (lines[i] != "") {
                    split(lines[i], fields, FS)
                    print fields[2] >> gene_names_file
                }
            }

            # Unique OG list
            print og >> uniqueOG_file

            # Species summary
            n = 0
            species_list_str = ""

            for (s in species_per_OG[og]) {
                n++
                species_list_str = species_list_str ? species_list_str "," s : s
            }

            print og, n, species_list_str >> species_counts_file
        }
    }

    # Sort species counts descending
    system("head -n 1 " species_counts_file " > tmp_header && " \
           "tail -n +2 " species_counts_file " | sort -k2,2nr >> tmp_header && " \
           "mv tmp_header OG_paralogs_species_counts_sorted_guidance2.tsv")
}
' gene_Arctic_species_map.txt all_genes_guidancec.tsv


more OG_paralogs_species_counts_sorted_guidance2.tsv


##########################
# Single copy
awk -F'\t' '
BEGIN {
    OFS="\t"
}

# Normalize gene IDs (fix case + suffix differences)
function normalize(g) {
    g = toupper(g)
    sub(/_[0-9]+$/, ".1", g)
    return g
}

# Load gene→species map (FIRST file)
FNR==NR {
    key = normalize($1)
    gene2species[key] = $2
    next
}

# Save header
NR==1 {
    header = $0
    next
}

{
    gene_raw = $2
    gene = normalize(gene_raw)

    # Skip + log unmapped genes (only once)
    if (!(gene in gene2species)) {
        if (!(gene in seen_unmapped)) {
            print gene_raw >> "unmapped_genes_singlecopy.txt"
            seen_unmapped[gene] = 1
        }
        next
    }

    species = gene2species[gene]

    # Count ALL genes per OG per species (for exclusion)
    all_count[$1,species]++

    # Only keep significant genes for output
    if ($3 < 0.05) {
        genes[$1] = genes[$1] $0 "\n"
        species_per_OG[$1][species] = 1
    }
}

END {
    single_copy_file="OG_single_copy_significant_guidance_nonArctic.tsv"
    gene_names_file="OG_single_copy_significant_gene_names_guidance_nonArctic.txt"
    uniqueOG_file="OG_single_copy_significant_uniqueOG_guidance_nonArctic.txt"
    species_counts_file="OG_single_copy_species_counts_guidance_nonArctic.txt"

    # Header for single-copy file
    print header > single_copy_file

    # Header for species count summary
    print "OG", "Num_species_with_sig_gene", "Species" > species_counts_file

    for (og in genes) {
        exclude=0

        # Exclude OGs if any species has >1 gene (total genes)
        for (s in all_count) {
            split(s, keyparts, SUBSEP)
            if (keyparts[1] == og && all_count[s] > 1) {
                exclude=1
                break
            }
        }

        if (!exclude) {
            # Save OG lines
            printf "%s", genes[og] >> single_copy_file

            # Save gene names
            split(genes[og], lines, "\n")
            for (i in lines) {
                if (lines[i] != "") {
                    split(lines[i], fields, FS)
                    print fields[2] >> gene_names_file
                }
            }

            # Unique OG list
            print og >> uniqueOG_file

            # Species counts summary
            n=0
            species_list_str=""

            for (s in species_per_OG[og]) {
                n++
                species_list_str = species_list_str ? species_list_str "," s : s
            }

            print og, n, species_list_str >> species_counts_file
        }
    }

    # Sort species counts descending
    system("head -n 1 " species_counts_file " > tmp_header && " \
           "tail -n +2 " species_counts_file " | sort -k2,2nr >> tmp_header && " \
           "mv tmp_header OG_single_copy_species_counts_sorted_guidance2.tsv")
}
' gene_Arctic_species_map.txt all_genes_guidance.tsv

more OG_single_copy_species_counts_sorted_guidance2.tsv


awk '{count[$2]++} END {for (val in count) print val, count[val]}' \
 OG_paralogs_species_counts_sorted_guidance2.tsv

awk '{count[$2]++} END {for (val in count) print val, count[val]}' \
 OG_single_copy_species_counts_sorted_guidance2.tsv
