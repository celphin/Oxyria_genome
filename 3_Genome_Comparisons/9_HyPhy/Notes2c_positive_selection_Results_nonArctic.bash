#########################
# Running HyPhy programs
# Jan 2026
############################

# follow instructions here
# https://github.com/siribi/RELAX-workflow

# other tutorial: https://bioinformaticsworkbook.org/dataAnalysis/ComparativeGenomics/Finding_Positive_Selection_With_Codeml.html#gsc.tab=0

##########################
cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees/Arctic_trees

ls *_guidance_unique_nonArctic.nxh | wc -l
ls *_guidance_unique_nonArctic.nxh.ABSREL.json  | wc -l
ls *_guidance_unique_nonArctic.nxh.BUSTED.json  | wc -l
ls *_guidance_unique_nonArctic.nxh.RELAX.json  | wc -l

# 8926

##########################
# Which orthogroups are demonstrating positive selection

# Narval1
tmux new-session -s total
tmux attach-session -t total

# look for found **1** branches under selection among **39** tested

for f in *_guidance_unique_nonArctic.nxh.ABSREL.json; do
  jq -r --arg file "$f" \
    '$file + "\t" + (.["test results"]["positive test results"]|tostring) + "\t" + (.["test results"].tested|tostring)' "$f"
done >> Total_ABSREL_results_guidance_nonArctic.txt

wc -l Total_ABSREL_results_guidance_nonArctic.txt

# 8849 Total_ABSREL_results_guidance_nonArctic.txt

#-------------------------------
# WRONG column 2 is a count of positive tests
# For Arctic ABSREL results
awk '$2 < 0.05 {print $1}' Total_ABSREL_guidance_results.txt | sort | uniq > Arctic_sig_OGs.txt
awk '$2 < 0.05 {print $1}' Total_ABSREL_results_guidance_nonArctic.txt | sort | uniq > nonArctic_sig_OGs.txt

more Arctic_sig_OGs.txt
# OG0000911_guidance_unique.nxh.ABSREL.json
# OG0000944_guidance_unique.nxh.ABSREL.json
# OG0000949_guidance_unique.nxh.ABSREL.json
# OG0000953_guidance_unique.nxh.ABSREL.json

more nonArctic_sig_OGs.txt
# OG0000042_guidance_unique_nonArctic.nxh.ABSREL.json
# OG0000395_guidance_unique_nonArctic.nxh.ABSREL.json
# OG0000398_guidance_unique_nonArctic.nxh.ABSREL.json
# OG0000405_guidance_unique_nonArctic.nxh.ABSREL.json
# OG0000471_guidance_unique_nonArctic.nxh.ABSREL.json

# reformat
cut -d'_' -f1 Arctic_sig_OGs.txt | sort | uniq > Arctic_sig_IDs.txt
cut -d'_' -f1 nonArctic_sig_OGs.txt | sort | uniq > nonArctic_sig_IDs.txt

wc -l Arctic_sig_IDs.txt
wc -l nonArctic_sig_IDs.txt
# 2956 Arctic_sig_IDs.txt
# 3120 nonArctic_sig_IDs.txt

# Arctic only OG
grep -Fxv -f nonArctic_sig_IDs.txt Arctic_sig_IDs.txt > Arctic_only_OGs.txt

# Shared OG
grep -Fx -f nonArctic_sig_IDs.txt Arctic_sig_IDs.txt > shared_OGs.txt

wc -l Arctic_only_OGs.txt
wc -l shared_OGs.txt

# 1650 Arctic_only_OGs.txt
# 1306 shared_OGs.txt

##########################
# Overlap with those in all 4

cut -f1 OG_4species_rows.tsv | sed 's/_.*//' | sort | uniq > OG_4species_IDs.txt

grep -Fx -f nonArctic_sig_IDs.txt OG_4species_IDs.txt > nonArctic_4species_shared.txt
grep -Fxv -f nonArctic_sig_IDs.txt OG_4species_IDs.txt > 4species_not_in_nonArctic.txt

wc -l 4species_not_in_nonArctic.txt
wc -l nonArctic_4species_shared.txt

# 194 4species_not_in_nonArctic.txt
# 20 nonArctic_4species_shared.txt

wc -l nonArctic_sig_IDs.txt
wc -l OG_4species_IDs.txt

3120 nonArctic_sig_IDs.txt
214 OG_4species_IDs.txt

# total ~8000 OGs

# Expected (3120/8000)*214 = 84 but we only see 20
# 4-species OGs are underrepresented among nonArctic-selected genes

# in R to get P-value

matrix <- matrix(c(20, 3100, 194, 4686), nrow=2)
fisher.test(matrix, alternative="less")


##############################
# Check gene and spp counts

# filter for only those with positive selection
awk '$2 != 0' Total_ABSREL_results_guidance_nonArctic.txt | sort -k2,2nr > ABSREL_nonzero_sorted_guidance_nonArctic.txt

wc -l ABSREL_nonzero_sorted_guidance_nonArctic.txt
# 4076 ABSREL_nonzero_sorted.txt
# 6876 ABSREL_guidance_nonzero_sorted.txt # 1000 more in the Arctic again
# 5729 ABSREL_nonzero_sorted_guidance_nonArctic.txt

more ABSREL_nonzero_sorted_guidance_nonArctic.txt
# OG0000088_guidance_unique_nonArctic.nxh.ABSREL.json     12      43
# OG0001285_guidance_unique_nonArctic.nxh.ABSREL.json     10      12
# OG0001636_guidance_unique_nonArctic.nxh.ABSREL.json     10      11
# OG0000231_guidance_unique_nonArctic.nxh.ABSREL.json     9       30
# OG0000596_guidance_unique_nonArctic.nxh.ABSREL.json     9       15
# OG0000652_guidance_unique_nonArctic.nxh.ABSREL.json     9       17
# OG0000705_guidance_unique_nonArctic.nxh.ABSREL.json     9       11
# OG0000708_guidance_unique_nonArctic.nxh.ABSREL.json     9       18
# OG0000859_guidance_unique_nonArctic.nxh.ABSREL.json     8       14
# OG0000985_guidance_unique_nonArctic.nxh.ABSREL.json     8       16
# OG0001359_guidance_unique_nonArctic.nxh.ABSREL.json     8       11
# OG0001831_guidance_unique_nonArctic.nxh.ABSREL.json     8       10
# OG0001983_guidance_unique_nonArctic.nxh.ABSREL.json     8       8
# OG0000769_guidance_unique_nonArctic.nxh.ABSREL.json     7       13
# OG0001600_guidance_unique_nonArctic.nxh.ABSREL.json     7       11
# OG0001788_guidance_unique_nonArctic.nxh.ABSREL.json     7       9

more ABSREL_guidance_nonzero_sorted.txt
# OG0000151_guidance_unique.nxh.ABSREL.json       12      35
# OG0000415_guidance_unique.nxh.ABSREL.json       11      23
# OG0000708_guidance_unique.nxh.ABSREL.json       11      14
# OG0000721_guidance_unique.nxh.ABSREL.json       11      15
# OG0000231_guidance_unique.nxh.ABSREL.json       10      27
# OG0000652_guidance_unique.nxh.ABSREL.json       10      16
# OG0001278_guidance_unique.nxh.ABSREL.json       10      13
# OG0001312_guidance_unique.nxh.ABSREL.json       10      12
# OG0001628_guidance_unique.nxh.ABSREL.json       10      11
# OG0002844_guidance_unique.nxh.ABSREL.json       10      10
# OG0002874_guidance_unique.nxh.ABSREL.json       10      10
# OG0001127_guidance_unique.nxh.ABSREL.json       9       14
# OG0001330_guidance_unique.nxh.ABSREL.json       9       13
# OG0001710_guidance_unique.nxh.ABSREL.json       9       12
# OG0002108_guidance_unique.nxh.ABSREL.json       9       9
# OG0002225_guidance_unique.nxh.ABSREL.json       9       11
# OG0002709_guidance_unique.nxh.ABSREL.json       9       9
# OG0000238_guidance_unique.nxh.ABSREL.json       8       26
# OG0000369_guidance_unique.nxh.ABSREL.json       8       23
# OG0000596_guidance_unique.nxh.ABSREL.json       8       15
# OG0000912_guidance_unique.nxh.ABSREL.json       8       13

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


mv ABSREL_nonzero_sorted_guidance_nonArctic.txt ABSREL_guidance_nonArctic_nonzero_sorted.txt

#########################
# Fix names for guidance_nonArctic below here
#---------------------
# Check for paralogs

echo -e "file\tgene\tcorrected_p\tfull_adaptive_model\tnonsyn_subs\tsyn_subs" > all_genes_guidance_nonArctic.tsv

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
done < ABSREL_guidance_nonArctic_nonzero_sorted.txt >> all_genes_guidance_nonArctic.tsv

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
	OG_single_copy_species_counts_sorted_guidance_nonArctic.tsv")
}
' all_genes_guidance_nonArctic.tsv

more OG_single_copy_species_counts_sorted_guidance_nonArctic.tsv
OG      Num_species_with_sig_gene       Species
OG0003385_guidance_unique_nonArctic.nxh.ABSREL.json     1       Rhuem
OG0004280_guidance_unique_nonArctic.nxh.ABSREL.json     1       Arabis
OG0005239_guidance_unique_nonArctic.nxh.ABSREL.json     1       Rhuem
OG0010327_guidance_unique_nonArctic.nxh.ABSREL.json     1       LOC
OG0014878_guidance_unique_nonArctic.nxh.ABSREL.json     1       Rhuem
OG0015142_guidance_unique_nonArctic.nxh.ABSREL.json     1       Arabis
OG0015699_guidance_unique_nonArctic.nxh.ABSREL.json     1       LOC
OG0016698_guidance_unique_nonArctic.nxh.ABSREL.json     1       Rhuem
OG0016885_guidance_unique_nonArctic.nxh.ABSREL.json     1       Rhuem
OG0017941_guidance_unique_nonArctic.nxh.ABSREL.json     1       Rhuem


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
    single_copy_file="OG_paralogs_significant_guidance_nonArctic.tsv"
    gene_names_file="OG_paralogs_significant_gene_names_guidance_nonArctic.txt"
    uniqueOG_file="OG_paralogs_significant_uniqueOG_guidance_nonArctic.txt"
    species_counts_file="OG_paralogs_species_counts_guidance_nonArctic.txt"

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
	mv tmp_header OG_paralogs_species_counts_sorted_guidance_nonArctic.tsv")
}
' all_genes_guidance_nonArctic.tsv

more OG_paralogs_species_counts_sorted_guidance_nonArctic.tsv

# OG      Num_species_with_sig_gene       Species
# OG0007401_tree.txt_unique_guidance_nonArctic.nxh.ABSREL.json     3       Rhuem,LOC,Arabis
# OG0013935_tree.txt_unique_guidance_nonArctic.nxh.ABSREL.json     3       Rhuem,LOC,Arabis
# OG0000510_tree.txt_unique_guidance_nonArctic.nxh.ABSREL.json     2       Rhuem,LOC
# OG0000579_tree.txt_unique_guidance_nonArctic.nxh.ABSREL.json     2       Rhuem,Arabis
# OG0000847_tree.txt_unique_guidance_nonArctic.nxh.ABSREL.json     2       LOC,Arabis

################################
# Get counts of each value

# total XX OGs tested
awk '{count[$2]++} END {for (val in count) print val, count[val]}' \
 OG_paralogs_species_counts_sorted_guidance_nonArctic2.tsv

# non Arctic
# 1 2926
# 2 1418
# 3 546
# 4 131

# Arctic
# 1 3004
# 2 1573
# 3 768
# 4 214

###################
# Per species counts

# Because not splitting LOC genes so all look like paralogs
grep -c 'LOC' OG_single_copy_species_counts_sorted_guidance_nonArctic.tsv
grep -c 'Rheum' OG_single_copy_species_counts_sorted_guidance_nonArctic.tsv
grep -c 'Arabis' OG_single_copy_species_counts_sorted_guidance_nonArctic.tsv
grep -c 'Other' OG_single_copy_species_counts_sorted_guidance_nonArctic.tsv
# 2
# 0
# 0
# 0

grep -c 'LOC' OG_paralogs_species_counts_sorted_guidance_nonArctic.tsv
grep -c 'Rhuem' OG_paralogs_species_counts_sorted_guidance_nonArctic.tsv
grep -c 'Arabis' OG_paralogs_species_counts_sorted_guidance_nonArctic.tsv

# 3002 - two spp in category
# 1244
# 1443

# ~ 1500 genes per species


############################
# Run again with spp lists
cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/CDS
grep ">" Arabis_alpina.f* | sed 's/>//' | awk '{print $0 "\tArabis"}' > gene_species_map.txt
grep ">" Prunus_persica.f* | sed 's/>//' | awk '{print $0 "\tPrunus"}' >> gene_species_map.txt
grep ">" Rheum_nobile_H0.f* | sed 's/>//' | awk '{print $0 "\tRhuem"}' >> gene_species_map.txt
grep ">" Arabidopsis_lyrata.f* | sed 's/>//' | awk '{print $0 "\tAlyrata"}' >> gene_species_map.txt

cp gene_species_map.txt /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees/Arctic_trees

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
           "mv tmp_header OG_paralogs_species_counts_sorted_guidance_nonArctic.tsv")
}
' gene_species_map.txt all_genes_guidance_nonArctic.tsv


more OG_paralogs_species_counts_sorted_guidance_nonArctic.tsv


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
           "mv tmp_header OG_single_copy_species_counts_sorted_guidance_nonArctic.tsv")
}
' gene_species_map.txt all_genes_guidance_nonArctic.tsv

more OG_single_copy_species_counts_sorted_guidance_nonArctic.tsv


#-----------------
# How much overlap?

awk '{count[$2]++} END {for (val in count) print val, count[val]}' \
 OG_paralogs_species_counts_sorted_guidance_nonArctic.tsv
# 1 2926
# 2 1417
# 3 543
# 4 133

awk '{count[$2]++} END {for (val in count) print val, count[val]}' \
 OG_single_copy_species_counts_sorted_guidance_nonArctic.tsv
# 1 1949
# 2 913
# 3 318
# 4 87

#------------------
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


##########################
# Check overlap again

cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees/Arctic_trees

# Extract full rows with all 4 species
awk '{split($3,species,","); if(length(species)==4) print $0}' \
OG_paralogs_species_counts_sorted_guidance_nonArctic.tsv \
> OG_4species_rows_nonArctic.tsv
# OG0000640_guidance_unique_nonArctic.nxh.ABSREL.json     4       Alyrata,Rhuem,Prunus,Arabis
# OG0002006_guidance_unique_nonArctic.nxh.ABSREL.json     4       Alyrata,Rhuem,Prunus,Arabis
# OG0002108_guidance_unique_nonArctic.nxh.ABSREL.json     4       Alyrata,Rhuem,Prunus,Arabis
# OG0002162_guidance_unique_nonArctic.nxh.ABSREL.json     4       Alyrata,Rhuem,Prunus,Arabis
# OG0002882_guidance_unique_nonArctic.nxh.ABSREL.json     4       Alyrata,Rhuem,Prunus,Arabis
# OG0003081_guidance_unique_nonArctic.nxh.ABSREL.json     4       Alyrata,Rhuem,Prunus,Arabis
# OG0003271_guidance_unique_nonArctic.nxh.ABSREL.json     4       Alyrata,Rhuem,Prunus,Arabis
# OG0003390_guidance_unique_nonArctic.nxh.ABSREL.json     4       Alyrata,Rhuem,Prunus,Arabis
# OG0003413_guidance_unique_nonArctic.nxh.ABSREL.json     4       Alyrata,Rhuem,Prunus,Arabis
# OG0003423_guidance_unique_nonArctic.nxh.ABSREL.json     4       Alyrata,Rhuem,Prunus,Arabis
# OG0003727_guidance_unique_nonArctic.nxh.ABSREL.json     4       Alyrata,Rhuem,Prunus,Arabis
# OG0003929_guidance_unique_nonArctic.nxh.ABSREL.json     4       Alyrata,Rhuem,Prunus,Arabis

# Extract full rows with all 3 species
awk '{split($3,species,","); if(length(species)==4) print $0}' \
OG_single_copy_species_counts_sorted_guidance_nonArctic.tsv \
> single_copy_4species_rows_nonArctic.tsv


#####################
# Check overlap with Arctic

cut -f1 OG_paralogs_species_counts_sorted_guidance_nonArctic.tsv | sort > nonArctic_list.txt
cut -f1 OG_paralogs_species_counts_sorted_guidance.tsv | sort > arctic_list.txt 

sed 's/_.*//' arctic_list.txt | sort > arctic.og
sed 's/_.*//' nonArctic_list.txt | sort > nonArctic.og

comm -12 nonArctic.og arctic.og > overlap_OGs.txt

wc -l arctic.og
wc -l nonArctic.og
wc -l overlap_OGs.txt

# 5560 arctic.og
# 5020 nonArctic.og
# 3308 overlap_OGs.txt

#------------------
# OGs in all Arctic and nonArctic four spp
cut -f1 OG_4species_rows_nonArctic.tsv | sort > nonArctic_list4.txt
cut -f1 OG_4species_rows.tsv | sort > arctic_list4.txt 

sed 's/_.*//' arctic_list4.txt | sort > arctic4.og
sed 's/_.*//' nonArctic_list4.txt | sort > nonArctic4.og

comm -12 nonArctic4.og arctic4.og > overlap4_OGs.txt

wc -l arctic4.og
wc -l nonArctic4.og
wc -l overlap4_OGs.txt

# 214 arctic4.og
# 133 nonArctic4.og
# 16 overlap4_OGs.txt


comm -12 nonArctic.og arctic4.og > overlap_OGs4.txt

wc -l overlap_OGs4.txt
# 164 overlap_OGs4.txt


##########################
# Check GO terms

# copy over file to GO analyses
cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/Gene_ontology
cp /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees/Arctic_trees/single_copy_4species_rows.tsv . 

# Step 1: Extract OG IDs
cut -f1 OG_4species_rows.tsv > OG_4species_IDs.txt

# Step 2: Get genes with p-value < 0.05 for those OGs
awk -F'\t' 'NR==FNR {og[$1]; next} $1 in og && $3 < 0.05 {print $1, $2, $3}' OFS='\t' OG_4species_IDs.txt baseline_guidance_all_genes.tsv > genes_filtered.txt

# Step 3: Join with InterProScan data and assign species
awk -F'\t' 'NR==FNR {g[$2]=$1"\t"$3; next} 
     $2 in g {
         split(g[$2], ogp, "\t")
         og = ogp[1]
         pval = ogp[2]
         gene = $2
         
         # assign species
         species=""
         if (gene ~ /^DOCTH0_CHR/) species="Dryas"
         else if (gene ~ /^OXYRIA_NCBI_CHR/) species="Oxyria"
         else if (gene ~ /^MAKER_/ || gene ~ /^SNAP_MASKED_/) species="Draba"
         else if (gene ~ /^G[0-9]+_T[0-9]+/) species="Coch"
         else next

         # print OG, p-value, gene, species, description, GO
         print og, pval, gene, species, $10, $11
     }' OFS='\t' genes_filtered.txt InterProscan_guidance_all_genes.tsv \
     > genes_GO_description_4species_filtered.tsv

more genes_GO_description_4species_filtered.tsv

# extract GO terms
cut -f6 genes_GO_description_4species_filtered.tsv > GO_terms_4spp.txt
# split
cat GO_terms_4spp.txt | tr ',' '\n' > GO_terms_4spp_split.txt
# count occurances
sort GO_terms_4spp_split.txt | uniq -c | sort -nr > GO_terms_4spp_counts.txt
