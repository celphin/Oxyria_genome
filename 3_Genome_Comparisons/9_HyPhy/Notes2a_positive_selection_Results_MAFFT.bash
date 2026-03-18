#########################
# Exploring HyPhy results for alignments without guidance
# Jan 2026
############################

# follow instructions here
# https://github.com/siribi/RELAX-workflow

# other tutorial: https://bioinformaticsworkbook.org/dataAnalysis/ComparativeGenomics/Finding_Positive_Selection_With_Codeml.html#gsc.tab=0

###############################
# Steps
# Explore results from ABSERL, RELAX and BUSTED

######################
# Narval2
tmux new-session -s total
tmux attach-session -t total

##########################
# Which orthogroups are demonstrating positive selection

cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees/Arctic_trees

# look for found **1** branches under selection among **39** tested

for f in *ABSREL.json; do
  jq -r --arg file "$f" \
    '$file + "\t" + (.["test results"]["positive test results"]|tostring) + "\t" + (.["test results"].tested|tostring)' "$f"
done >> Total_ABSREL_results.txt

wc -l Total_ABSREL_results.txt
#16254 Total_ABSREL_results.txt

# baseline genes tested

echo -e "file\tgene\tcorrected_p\tfull_adaptive_model\tnonsyn_subs\tsyn_subs" > baseline_all_genes.tsv

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
done < Total_ABSREL_results.txt >> baseline_all_genes.tsv

#------------------
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


##################################
# Look at all the data for each gene

echo -e "file\tgene\tcorrected_p\tfull_adaptive_model\tnonsyn_subs\tsyn_subs" > all_genes.tsv

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
done < ABSREL_nonzero_sorted.txt >> all_genes.tsv

#------------------------------------
cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees/Arctic_trees

more all_genes.tsv
# file    gene    corrected_p     full_adaptive_model     nonsyn_subs     syn_subs
# OG0015441_tree.txt_unique.nxh.ABSREL.json       DOCTH0_CHR300006814     0.0004372556584134046   0.3516174924510933      0.3516111809241831      0.000006311526910256882
# OG0015441_tree.txt_unique.nxh.ABSREL.json       DOCTH0_CHR400000398     0.02450561644103588     81.70688001451325       81.61750677619105       0.08937323832213598
# OG0015441_tree.txt_unique.nxh.ABSREL.json       DOCTH0_CHR400001723     0.04809666165521614     0.7219684042673827      0.7218549227459655      0.0001134815214158162
# OG0015441_tree.txt_unique.nxh.ABSREL.json       DOCTH0_CHR500000225     0.004120137030874105    0.8490730084796092      0.8277188893013476      0.02135411917826048
# OG0015441_tree.txt_unique.nxh.ABSREL.json       DOCTH0_CHR600006291     0.8251896945627508      0.009922310372463279    0.009922310372463279    1E-10
# OG0015441_tree.txt_unique.nxh.ABSREL.json       DOCTH0_CHR700004842     0.7288055648788709      0.04929395054839377     0.03978031196134621     0.00951363858704757
# OG0015441_tree.txt_unique.nxh.ABSREL.json       DOCTH0_CHR800005359     1       0.05346671284950748     0.0333851610557265      0.02008155179378096
# OG0015441_tree.txt_unique.nxh.ABSREL.json       DOCTH0_CHR800006070     0.02062170869551722     1.588682372205865       1.588431437175698       0.0002509350301664349
# OG0015441_tree.txt_unique.nxh.ABSREL.json       DOCTH0_CHR900000114     0.2253137037661492      0.2083161911023486      0.2083161911023486      1E-10
# OG0015441_tree.txt_unique.nxh.ABSREL.json       DOCTH0_CHR900001109     0.0003764289048757696   0.7361764732201124      0.7288597604512717      0.007316712768843253
# OG0015441_tree.txt_unique.nxh.ABSREL.json       NODE10  NA      0.01747597313518143     0.008870791092959381    0.008605182042222015
# OG0015441_tree.txt_unique.nxh.ABSREL.json       NODE14  NA      0.4242486888707904      0.424219634219517       0.00002905465127344884
# OG0015441_tree.txt_unique.nxh.ABSREL.json       NODE16  NA      0.1181149077389998      1E-10   0.1181149077389998
# OG0015441_tree.txt_unique.nxh.ABSREL.json       NODE4   NA      0.03224344260775407     0.009815819793640891    0.02242762281411316
# OG0015441_tree.txt_unique.nxh.ABSREL.json       NODE9   NA      0.08640187693671385     0.04310605799912858     0.04329581893758524
# OG0002756_tree.txt_unique.nxh.ABSREL.json       101292950       NA      0.04224188206173039     0.01762990408411761     0.02461197797761279
# OG0002756_tree.txt_unique.nxh.ABSREL.json       106305749       NA      0.05526377480547504     0.04129273247652632     0.0139710423289487
# OG0002756_tree.txt_unique.nxh.ABSREL.json       106334671       NA      0.05685218901253847     0.02946906506949114     0.02738312394304725
# OG0002756_tree.txt_unique.nxh.ABSREL.json       106341938       NA      0.09253514647253609     0.04748065343387353     0.04505449303866266
# OG0002756_tree.txt_unique.nxh.ABSREL.json       AT3G24870       NA      0.005669645352322288    0.003809069473775934    0.001860575878546351
# OG0002756_tree.txt_unique.nxh.ABSREL.json       AT3G24880       NA      0       1E-10   1E-10
# OG0002756_tree.txt_unique.nxh.ABSREL.json       DOCTH0_CHR800000942     1       0.02908991671106302     0.005860313054690892    0.02322960365637216
# OG0002756_tree.txt_unique.nxh.ABSREL.json       FEHAP220557_T1  NA      0.9064784805021109      0.9012718063757053      0.005206674126413438
# OG0002756_tree.txt_unique.nxh.ABSREL.json       FEHAP227832_T1  NA      0.03087265623218266     0.01346867434473869     0.01740398188744395
# OG0002756_tree.txt_unique.nxh.ABSREL.json       FEHAP231234_T1  NA      0.02354448752524851     0.01701435386090207     0.00653013366434643
# OG0002756_tree.txt_unique.nxh.ABSREL.json       FT01GENE13907_T1        NA      0.02154504456836176     0.01294554579143482     0.0085994987769269
# OG0002756_tree.txt_unique.nxh.ABSREL.json       FT01GENE23674_T1        NA      0.02863323899943923     0.00837517491645824     0.02025806408298101
# OG0002756_tree.txt_unique.nxh.ABSREL.json       FT01GENE32782_T1        NA      0.01538833081088953     0.007794354732824121    0.007593976078065408
# OG0002756_tree.txt_unique.nxh.ABSREL.json       G13110_T1       1       0.01835770530710042     0.009186379267072341    0.009171326040028087
# OG0002756_tree.txt_unique.nxh.ABSREL.json       G13114_T1       1       0.01449186701940357     0.01217486056482487     0.002317006454578714
# OG0002756_tree.txt_unique.nxh.ABSREL.json       LOC103927702    NA      0.008105510869024106    0.00600629644860563     0.002099214420418468
# OG0002756_tree.txt_unique.nxh.ABSREL.json       LOC103960911    NA      0.0120865200058613      0.009656158848092868    0.002430361157768429
# OG0002756_tree.txt_unique.nxh.ABSREL.json       LOC110226977    NA      0.0680983111685811      0.03691776620321979     0.0311805449653612
# OG0002756_tree.txt_unique.nxh.ABSREL.json       LOC126585513    NA      0.007461265348889675    0.00550374520590679     0.001957520142982876
# OG0002756_tree.txt_unique.nxh.ABSREL.json       LOC126620741    NA      0.009141231247546525    0.005782647630398593    0.003358583617147915
# OG0002756_tree.txt_unique.nxh.ABSREL.json       LOC126805521    NA      0.04329538530776133     0.01765926655930957     0.02563611874845179
# OG0002756_tree.txt_unique.nxh.ABSREL.json       LOC133737018    NA      0.02288372680383075     0.007900611097733148    0.01498311570609757
# OG0002756_tree.txt_unique.nxh.ABSREL.json       LOC17893131     NA      0.02131636694347428     0.01396059275693361     0.007355774186540692
# OG0002756_tree.txt_unique.nxh.ABSREL.json       LOC17893888     NA      0.02091798373128021     0.01116967683421007     0.00974830689707013
# OG0002756_tree.txt_unique.nxh.ABSREL.json       LOC18768734     NA      0.046874236352324       0.0202918041743501      0.026582432177974
# OG0002756_tree.txt_unique.nxh.ABSREL.json       LOC9319635      NA      0.01540298146271036     0.007664874806765295    0.007738106655945072
# OG0002756_tree.txt_unique.nxh.ABSREL.json       MAKER_LG4_SNAP_GENE_88_212_MRNA_1       1       0.05577679831128731     0.01421545659384497     0.04156134171744238
# OG0002756_tree.txt_unique.nxh.ABSREL.json       MAKER_LG6_SNAP_GENE_34_144_MRNA_1       1       0.05730609043011364     0.02761174304367743     0.02969434738643626
# OG0002756_tree.txt_unique.nxh.ABSREL.json       NODE10  NA      0       1E-10   1E-10
# OG0002756_tree.txt_unique.nxh.ABSREL.json       NODE12  NA      0.001828331719167163    0.001828331719167163    1E-10
# OG0002756_tree.txt_unique.nxh.ABSREL.json       NODE15  NA      0.002279967118069129    0.002279967118069129    1E-10
# OG0002756_tree.txt_unique.nxh.ABSREL.json       NODE18  NA      0.2792557097959788      0.2576989867345455      0.02155672306143047
# OG0002756_tree.txt_unique.nxh.ABSREL.json       NODE2   NA      0.150355362923261       0.1147208593710522      0.03563450355220867
# OG0002756_tree.txt_unique.nxh.ABSREL.json       NODE21  NA      0.06746963677716496     0.04380567175271882     0.02366396502444627
# OG0002756_tree.txt_unique.nxh.ABSREL.json       NODE22  NA      0.0548309739785904      0.0408688980758322      0.01396207590275822
# OG0002756_tree.txt_unique.nxh.ABSREL.json       NODE23  NA      0.05867933677356833     0.02473046773865651     0.03394886903491185
# OG0002756_tree.txt_unique.nxh.ABSREL.json       NODE28  NA      0.009074282400812419    0.00353019130788312     0.005544091092929301
# OG0002756_tree.txt_unique.nxh.ABSREL.json       NODE29  NA      0.02566795546498141     0.01783363595732118     0.007834319507660262
# OG0002756_tree.txt_unique.nxh.ABSREL.json       NODE3   NA      0.01975559073323658     0.009635425942293035    0.01012016479094355
# OG0002756_tree.txt_unique.nxh.ABSREL.json       NODE33  NA      0.04125102071653155     0.02030002957309187     0.0209509911434398
# OG0002756_tree.txt_unique.nxh.ABSREL.json       NODE34  NA      0.01892897222403469     0.005146874327040652    0.01378209789699406
# OG0002756_tree.txt_unique.nxh.ABSREL.json       NODE35  NA      0.01194543544482175     0.009132554999078531    0.002812880445743206
# OG0002756_tree.txt_unique.nxh.ABSREL.json       NODE37  NA      0.01850525105393068     0.01008877362997801     0.008416477423952658

#-------------------------
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
    single_copy_file="OG_single_copy_significant.tsv"
    gene_names_file="OG_single_copy_significant_gene_names.txt"
    uniqueOG_file="OG_single_copy_significant_uniqueOG.txt"
    species_counts_file="OG_single_copy_species_counts.txt"

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
    system("head -n 1 " species_counts_file " > tmp_header && tail -n +2 " species_counts_file " | sort -k2,2nr >> tmp_header && mv tmp_header OG_single_copy_species_counts_sorted.tsv")
}
' all_genes.tsv

more OG_single_copy_species_counts_sorted.tsv

# OG      Num_species_with_sig_gene       Species
# OG0005895_tree.txt_unique.nxh.ABSREL.json       3       Dryas,Oxyria,Draba
# OG0008262_tree.txt_unique.nxh.ABSREL.json       3       Dryas,Draba,Coch
# OG0008518_tree.txt_unique.nxh.ABSREL.json       3       Dryas,Oxyria,Coch
# OG0009474_tree.txt_unique.nxh.ABSREL.json       3       Dryas,Oxyria,Draba
# OG0009938_tree.txt_unique.nxh.ABSREL.json       3       Dryas,Draba,Coch
# OG0011085_tree.txt_unique.nxh.ABSREL.json       3       Oxyria,Draba,Coch
# OG0011114_tree.txt_unique.nxh.ABSREL.json       3       Dryas,Oxyria,Coch
# OG0011174_tree.txt_unique.nxh.ABSREL.json       3       Oxyria,Draba,Coch
# OG0012177_tree.txt_unique.nxh.ABSREL.json       3       Dryas,Draba,Coch
# OG0012361_tree.txt_unique.nxh.ABSREL.json       3       Dryas,Oxyria,Draba
# OG0013353_tree.txt_unique.nxh.ABSREL.json       3       Dryas,Draba,Coch
# OG0013973_tree.txt_unique.nxh.ABSREL.json       3       Dryas,Draba,Coch
# OG0014460_tree.txt_unique.nxh.ABSREL.json       3       Dryas,Draba,Coch

wc -l OG_single_copy_species_counts_sorted.tsv
# 1664 OG_single_copy_species_counts_sorted.tsv



###########################
# Also get list for those with paralogs

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
    single_copy_file="OG_paralogs_significant.tsv"
    gene_names_file="OG_paralogs_significant_gene_names.txt"
    uniqueOG_file="OG_paralogs_significant_uniqueOG.txt"
    species_counts_file="OG_paralogs_species_counts.txt"

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
	mv tmp_header OG_paralogs_species_counts_sorted.tsv")
}
' all_genes.tsv

more OG_paralogs_species_counts_sorted.tsv

wc -l OG_paralogs_species_counts_sorted.tsv
# 3744 OG_paralogs_species_counts_sorted.tsv

# OG      Num_species_with_sig_gene       Species
# OG0005531_tree.txt_unique.nxh.ABSREL.json       4       Coch,Dryas,Oxyria,Draba
# OG0004358_tree.txt_unique.nxh.ABSREL.json       3       Coch,Dryas,Draba
# OG0004910_tree.txt_unique.nxh.ABSREL.json       3       Dryas,Oxyria,Draba
# OG0005895_tree.txt_unique.nxh.ABSREL.json       3       Dryas,Oxyria,Draba
# OG0007028_tree.txt_unique.nxh.ABSREL.json       3       Dryas,Oxyria,Draba
# OG0007601_tree.txt_unique.nxh.ABSREL.json       3       Coch,Dryas,Draba
# OG0008262_tree.txt_unique.nxh.ABSREL.json       3       Coch,Dryas,Draba
# OG0008518_tree.txt_unique.nxh.ABSREL.json       3       Coch,Dryas,Oxyria
# OG0008574_tree.txt_unique.nxh.ABSREL.json       3       Coch,Oxyria,Draba
# OG0009332_tree.txt_unique.nxh.ABSREL.json       3       Coch,Dryas,Draba
# OG0009474_tree.txt_unique.nxh.ABSREL.json       3       Dryas,Oxyria,Draba
# OG0009836_tree.txt_unique.nxh.ABSREL.json       3       Coch,Dryas,Oxyria
# OG0009890_tree.txt_unique.nxh.ABSREL.json       3       Dryas,Oxyria,Draba
# OG0009938_tree.txt_unique.nxh.ABSREL.json       3       Coch,Dryas,Draba
# OG0011085_tree.txt_unique.nxh.ABSREL.json       3       Coch,Oxyria,Draba
# OG0011114_tree.txt_unique.nxh.ABSREL.json       3       Coch,Dryas,Oxyria
# OG0011174_tree.txt_unique.nxh.ABSREL.json       3       Coch,Oxyria,Draba
# OG0012177_tree.txt_unique.nxh.ABSREL.json       3       Coch,Dryas,Draba
# OG0012361_tree.txt_unique.nxh.ABSREL.json       3       Dryas,Oxyria,Draba
# OG0013353_tree.txt_unique.nxh.ABSREL.json       3       Coch,Dryas,Draba
# OG0013973_tree.txt_unique.nxh.ABSREL.json       3       Coch,Dryas,Draba
# OG0014460_tree.txt_unique.nxh.ABSREL.json       3       Coch,Dryas,Draba

more OG_single_copy_species_counts_sorted.tsv

# OG      Num_species_with_sig_gene       Species
# OG0005895_tree.txt_unique.nxh.ABSREL.json       3       Dryas,Oxyria,Draba
# OG0008262_tree.txt_unique.nxh.ABSREL.json       3       Dryas,Draba,Coch
# OG0008518_tree.txt_unique.nxh.ABSREL.json       3       Dryas,Oxyria,Coch
# OG0009474_tree.txt_unique.nxh.ABSREL.json       3       Dryas,Oxyria,Draba
# OG0009938_tree.txt_unique.nxh.ABSREL.json       3       Dryas,Draba,Coch
# OG0011085_tree.txt_unique.nxh.ABSREL.json       3       Oxyria,Draba,Coch
# OG0011114_tree.txt_unique.nxh.ABSREL.json       3       Dryas,Oxyria,Coch
# OG0011174_tree.txt_unique.nxh.ABSREL.json       3       Oxyria,Draba,Coch
# OG0012177_tree.txt_unique.nxh.ABSREL.json       3       Dryas,Draba,Coch
# OG0012361_tree.txt_unique.nxh.ABSREL.json       3       Dryas,Oxyria,Draba
# OG0013353_tree.txt_unique.nxh.ABSREL.json       3       Dryas,Draba,Coch
# OG0013973_tree.txt_unique.nxh.ABSREL.json       3       Dryas,Draba,Coch
# OG0014460_tree.txt_unique.nxh.ABSREL.json       3       Dryas,Draba,Coch

# Functions
# | OG ID     | Count | Key Domains / Annotations                                                | Functional Theme                |
# | --------- | ----- | ------------------------------------------------------------------------ | ------------------------------- |
# | OG0005531 | 4     | Bms1/Tsr1-type G domain, P-loop NTPase, Ribosome biogenesis protein Bms1 | Ribosome / RNA biogenesis       |
# | OG0008518 | 3     | Large ribosomal subunit protein uL5, conserved sites                     | Ribosome / RNA biogenesis       |
# | OG0009836 | 3     | RsmD, RNA methyltransferase, SAM-dependent methyltransferase             | RNA modification / Ribosome     |
# | OG0009938 | 3     | U3 small nucleolar RNA-associated protein 13, WD40 repeats               | Ribosome / RNA biogenesis       |
# | OG0004437 | 3     | MICRORCHIDIA ATPase, S5 domain, HSP90-like ATPase                        | ATPase / Cytoskeleton           |
# | OG0007028 | 3     | FH2 domain, Formin-like                                                  | Cytoskeleton / Actin nucleation |
# | OG0009332 | 3     | Major facilitator superfamily (MFS transporter)                          | Membrane transport              |
# | OG0011174 | 3     | 7TM GPCR, GCR1-cAMP receptor                                             | Membrane signaling              |
# | OG0012361 | 3     | Nonaspanin (TM9SF)                                                       | Membrane transport              |
# | OG0004910 | 3     | Glycosyltransferase, DUF4094                                             | Protein / Glycosylation         |
# | OG0009474 | 3     | Beta-propeller, RING/FYVE/PHD, Clathrin, Zinc finger                     | Protein interaction / Scaffold  |
# | OG0011085 | 3     | Tetratricopeptide repeat (TPR), CLU domain                               | Protein interaction / Scaffold  |
# | OG0011114 | 3     | Pentatricopeptide repeat (PPR), TPR-like                                 | Protein interaction / Scaffold  |
# | OG0009890 | 3     | BRCA2 helical domain, OB-fold, Nucleic acid-binding                      | Protein interaction / Scaffold  |
# | OG0012177 | 3     | LEA_2 subgroup, Late embryogenesis abundant protein                      | Stress / Development            |
# | OG0008574 | 3     | Tim23-like                                                               | Mitochondrial import            |
# | OG0004358 | 3     | Inositol monophosphatase (HAL2)                                          | Metabolism / Signaling          |


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


# OG0005531_tree.txt_unique.nxh.ABSREL.json       4       Cochgroen,Draba,Dryas,Oxyria
# OG0004358_tree.txt_unique.nxh.ABSREL.json       3       Dryas,Cochgroen,Draba
# OG0004437_tree.txt_unique.nxh.ABSREL.json       3       Dryas,Oxyria,Draba
# OG0004910_tree.txt_unique.nxh.ABSREL.json       3       Dryas,Draba,Oxyria
# OG0005895_tree.txt_unique.nxh.ABSREL.json       3       Dryas,Oxyria,Draba
# OG0007028_tree.txt_unique.nxh.ABSREL.json       3       Dryas,Draba,Oxyria
# OG0007601_tree.txt_unique.nxh.ABSREL.json       3       Dryas,Cochgroen,Draba
# OG0008262_tree.txt_unique.nxh.ABSREL.json       3       Draba,Cochgroen,Dryas
# OG0008518_tree.txt_unique.nxh.ABSREL.json       3       Dryas,Cochgroen,Oxyria
# OG0008574_tree.txt_unique.nxh.ABSREL.json       3       Oxyria,Draba,Cochgroen
# OG0009332_tree.txt_unique.nxh.ABSREL.json       3       Dryas,Cochgroen,Draba
# OG0009474_tree.txt_unique.nxh.ABSREL.json       3       Draba,Oxyria,Dryas
# OG0009836_tree.txt_unique.nxh.ABSREL.json       3       Oxyria,Cochgroen,Dryas
# OG0009890_tree.txt_unique.nxh.ABSREL.json       3       Oxyria,Draba,Dryas
# OG0009938_tree.txt_unique.nxh.ABSREL.json       3       Cochgroen,Dryas,Draba
# OG0011085_tree.txt_unique.nxh.ABSREL.json       3       Cochgroen,Draba,Oxyria
# OG0011114_tree.txt_unique.nxh.ABSREL.json       3       Oxyria,Dryas,Cochgroen
# OG0011174_tree.txt_unique.nxh.ABSREL.json       3       Oxyria,Draba,Cochgroen
# OG0012177_tree.txt_unique.nxh.ABSREL.json       3       Draba,Cochgroen,Dryas
# OG0012361_tree.txt_unique.nxh.ABSREL.json       3       Draba,Oxyria,Dryas
# OG0013353_tree.txt_unique.nxh.ABSREL.json       3       Dryas,Cochgroen,Draba
# OG0013973_tree.txt_unique.nxh.ABSREL.json       3       Dryas,Cochgroen,Draba
# OG0014460_tree.txt_unique.nxh.ABSREL.json       3       Draba,Dryas,Cochgroen


#######################################
# check functions for all genes with multiple overlaps
grep OG0005531 InterProscan_ABSREL_sig_genes.tsv

# OG0005531       DOCTH0_CHR200004771     0.00612228594200048     Dryas_octopetala_interproscan_output.tsv        IPR007034,IPR012948,IPR027417,IPR030387,IPR037875,IPR039761  C-terminal, N-terminal,AARP2CN,Bms1/Tsr1-type G domain,P-loop containing nucleoside triphosphate hydrolase,Ribosome biogenesis protein Bms1,Ribosome biogenesis protein Bms1/Tsr1,Ribosome biogenesis protein BMS1/TSR1    GO:0005525,GO:0005634,GO:0042254
# OG0005531       G20420_T1       0.016117811918277       NA      NA      NA      NA
# OG0005531       OXYRIA_NCBI_CHR500003522        0.0352221931862658      Oxyria_digyna_H1_interproscan_output.tsv        IPR007034,IPR012948,IPR027417,IPR030387,IPR037875,IPR039761  C-terminal, N-terminal,AARP2CN,Bms1/Tsr1-type G domain,P-loop containing nucleoside triphosphate hydrolase,Ribosome biogenesis protein Bms1,Ribosome biogenesis protein Bms1/Tsr1,Ribosome biogenesis protein BMS1/TSR1    GO:0005525,GO:0005634,GO:0042254
# OG0005531       SNAP_MASKED_LG5_PROCESSED_GENE_102_35_MRNA_1    0.0354711550003797      NA      NA      NA      NA

# Ribosome biogenesis protein


cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/Gene_ontology

awk -F'\t' '{
    split($1,a,"_tree");
    print a[1] "\t" $2 "\t" $3
}' OGs_multi_species.sorted.tsv > OG_counts_species.tsv


awk -F'\t' '
NR==FNR {
    count[$1]=$2;
    spp[$1]=$3;
    next
}
($1 in count) {
    print $1 "\t" count[$1] "\t" spp[$1] "\t" $0
}
' OG_counts_species.tsv InterProscan_ABSREL_sig_genes.tsv >  OG_full_annotation.tsv

cut -f1,2,9 OG_full_annotation.sorted.tsv > OG_description_count.tsv

sort OG_description_count.tsv | uniq > OG_description_count.unique.tsv

grep -v -P '\tNA$' OG_description_count.tsv | sort | uniq > OG_description_count.noNA.tsv

# Most of these genes fall into four broad categories:
# Ribosome and RNA biogenesis / modification – OG0005531, OG0008518, OG0009836, OG0009938
# ATPase/GTPase and cytoskeleton remodeling – OG0004437, OG0007028
# Membrane transport / signaling – OG0009332, OG0011174, OG0012361
# Protein interaction / structural repeat scaffolds – OG0009474, OG0011085, OG0011114, OG0009890

# Secondary themes include stress response, 
# glycosylation, and mitochondrial import.

# OG0012177 – LEA (Late embryogenesis abundant) protein
# OG0004910 – Glycosyltransferase (DUF4094)
# OG0009836 – SAM-dependent RNA methyltransferase


# | OG ID     | Count | Key Domains / Annotations                                                | Functional Theme                |
# | --------- | ----- | ------------------------------------------------------------------------ | ------------------------------- |
# | OG0005531 | 4     | Bms1/Tsr1-type G domain, P-loop NTPase, Ribosome biogenesis protein Bms1 | Ribosome / RNA biogenesis       |
# | OG0008518 | 3     | Large ribosomal subunit protein uL5, conserved sites                     | Ribosome / RNA biogenesis       |
# | OG0009836 | 3     | RsmD, RNA methyltransferase, SAM-dependent methyltransferase             | RNA modification / Ribosome     |
# | OG0009938 | 3     | U3 small nucleolar RNA-associated protein 13, WD40 repeats               | Ribosome / RNA biogenesis       |
# | OG0004437 | 3     | MICRORCHIDIA ATPase, S5 domain, HSP90-like ATPase                        | ATPase / Cytoskeleton           |
# | OG0007028 | 3     | FH2 domain, Formin-like                                                  | Cytoskeleton / Actin nucleation |
# | OG0009332 | 3     | Major facilitator superfamily (MFS transporter)                          | Membrane transport              |
# | OG0011174 | 3     | 7TM GPCR, GCR1-cAMP receptor                                             | Membrane signaling              |
# | OG0012361 | 3     | Nonaspanin (TM9SF)                                                       | Membrane transport              |
# | OG0004910 | 3     | Glycosyltransferase, DUF4094                                             | Protein / Glycosylation         |
# | OG0009474 | 3     | Beta-propeller, RING/FYVE/PHD, Clathrin, Zinc finger                     | Protein interaction / Scaffold  |
# | OG0011085 | 3     | Tetratricopeptide repeat (TPR), CLU domain                               | Protein interaction / Scaffold  |
# | OG0011114 | 3     | Pentatricopeptide repeat (PPR), TPR-like                                 | Protein interaction / Scaffold  |
# | OG0009890 | 3     | BRCA2 helical domain, OB-fold, Nucleic acid-binding                      | Protein interaction / Scaffold  |
# | OG0012177 | 3     | LEA_2 subgroup, Late embryogenesis abundant protein                      | Stress / Development            |
# | OG0008574 | 3     | Tim23-like                                                               | Mitochondrial import            |
# | OG0004358 | 3     | Inositol monophosphatase (HAL2)                                          | Metabolism / Signaling          |


# | OG ID     | Functional Category             | Evolution in Most Plants | Likely Arctic-Specific Positive Selection | Notes                                                           |
# | --------- | ------------------------------- | ------------------------ | ----------------------------------------- | --------------------------------------------------------------- |
# | OG0005531 | Ribosome / RNA biogenesis       | Highly conserved         | Likely                                    | Cold-stable ribosome assembly advantageous in Arctic plants     |
# | OG0008518 | Ribosome / RNA biogenesis       | Highly conserved         | Likely                                    | uL5 ribosomal protein may evolve subtle stabilizing mutations   |
# | OG0009836 | RNA modification / Ribosome     | Moderately conserved     | Likely                                    | RNA methyltransferases stabilize rRNA under cold stress         |
# | OG0009938 | Ribosome / RNA biogenesis       | Highly conserved         | Likely                                    | WD40 scaffold for snoRNA processing, supports ribosome assembly |
# | OG0004437 | ATPase / Cytoskeleton           | Moderately conserved     | Possible                                  | Morc ATPase may have adaptive variants in Arctic plants         |
# | OG0007028 | Cytoskeleton / Actin nucleation | Moderately conserved     | Likely                                    | Formins critical for actin dynamics at low temperatures         |
# | OG0009332 | Membrane transport              | Variable                 | Likely                                    | MFS transporter function affected by membrane fluidity          |
# | OG0011174 | Membrane signaling (GPCR)       | Variable                 | Possible                                  | GPCR adaptations may tune stress signaling pathways             |
# | OG0012361 | Membrane transport              | Moderately conserved     | Possible                                  | TM9SF function may be fine-tuned for cold environments          |
# | OG0004910 | Protein / Glycosylation         | Conserved                | Possible                                  | Glycosylation may stabilize proteins in cold stress             |
# | OG0009474 | Protein interaction / Scaffold  | Moderately conserved     | Possible                                  | Beta-propeller/RING may aid multi-protein complex stability     |
# | OG0011085 | Protein interaction / Scaffold  | Conserved                | Possible                                  | TPR scaffolds stabilize protein interactions under stress       |
# | OG0011114 | Protein interaction / Scaffold  | Conserved                | Possible                                  | PPR proteins involved in RNA editing, may adapt for cold        |
# | OG0009890 | Protein interaction / Scaffold  | Conserved                | Possible                                  | BRCA2-like scaffold may support DNA repair under stress         |
# | OG0012177 | Stress / Development            | Variable                 | Likely                                    | LEA proteins are classic cold/drought adaptation targets        |
# | OG0008574 | Mitochondrial import            | Highly conserved         | Possible                                  | Tim23-like ensures energy production in cold conditions         |
# | OG0004358 | Metabolism / Signaling          | Conserved                | Possible                                  | Inositol monophosphatase may help osmotic balance in cold       |

# Ribosome / RNA biogenesis genes 
# (OG0005531, OG0008518, OG0009836, OG0009938)
# Cold slows translation and RNA processing.
 # Mutations that stabilize ribosomes or RNA under 
 # cold stress increase fitness.
# Studies in Arctic and alpine plants (like Arabidopsis lyrata, 
# Brassica species) have shown ribosomal proteins and RNA 
# modification enzymes under positive selection in cold-adapted 
# populations.

# Cytoskeleton / ATPase genes (OG0004437, OG0007028)
# Maintaining actin dynamics and protein complexes at low 
# temperatures is critical. Adaptive changes in actin nucleators
 # and ATPases have been documented in extremophile plants.

# Membrane transport / signaling genes (OG0009332, OG0011174, OG0012361)
# Membrane fluidity and transporter function are strongly 
# influenced by temperature. Selection often favors alleles 
# that maintain transporter efficiency in freezing conditions.

# Stress / LEA proteins (OG0012177)
# LEA proteins are classic examples of cold and desiccation-adaptive 
# genes. Positive selection in LEA proteins is well-documented 
# in Arctic, alpine, and drought-tolerant plants.
# https://pmc.ncbi.nlm.nih.gov/articles/PMC11772582/
# https://pmc.ncbi.nlm.nih.gov/articles/PMC9031148/

########################
# Prank results from below

# OG0005531_tree.txt_prank_unique.nxh.ABSREL.json 3       Cochgroen,Oxyria,Draba
# OG0011114_tree.txt_prank_unique.nxh.ABSREL.json 3       Oxyria,Dryas,Cochgroen
# OG0011456_tree.txt_prank_unique.nxh.ABSREL.json 3       Dryas,Oxyria,Draba
# OG0012809_tree.txt_prank_unique.nxh.ABSREL.json 3       Dryas,Draba,Cochgroen

# How many are paralogs?

# Functions
# | OG ID     | Count | Key Domains / Annotations                                                | Functional Theme                |
# | --------- | ----- | ------------------------------------------------------------------------ | ------------------------------- |
# | OG0005531 | 4     | Bms1/Tsr1-type G domain, P-loop NTPase, Ribosome biogenesis protein Bms1 | Ribosome / RNA biogenesis       |
# | OG0011114 | 3     | Pentatricopeptide repeat (PPR), TPR-like                                 | Protein interaction / Scaffold  |


################################
# Rerun with prank alignments for just significant genes above

#------------------------
for f in *_prank_unique.nxh.ABSREL.json; do
  jq -r --arg file "$f" \
    '$file + "\t" + (.["test results"]["positive test results"]|tostring) + "\t" + (.["test results"].tested|tostring)' "$f"
done >> Total_ABSREL_prank_results.txt


# filter for only those with positive selection
awk '$2 != 0' Total_ABSREL_prank_results.txt | sort -k2,2nr > ABSREL_prank_nonzero_sorted.txt

# get gene names
while read file _; do
  jq -r '
    .["branch attributes"]["0"]
    | to_entries[]
    | select(.value["Corrected P-value"] < 0.05)
    | "'"$file"'\t\(.key)\t\(.value["Corrected P-value"])"
  ' "$file"
done <  ABSREL_prank_nonzero_sorted.txt > all_significant_genes_prank.tsv

#------------------
# Link to gene functions
cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/Gene_ontology

cp /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees/Arctic_trees/all_significant_genes_prank.tsv . 

awk -F'\t' '$NF != "null"' all_significant_genes_prank.tsv > all_significant_genes_prank_filtered.tsv


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
' all_significant_genes_prank_filtered.tsv | sort > OGs_multi_species_prank.tsv

sort -t$'\t' -k2,2nr OGs_multi_species_prank.tsv > OGs_multi_species_prank.sorted.tsv

# OG0005531_tree.txt_prank_unique.nxh.ABSREL.json 3       Cochgroen,Oxyria,Draba
# OG0011114_tree.txt_prank_unique.nxh.ABSREL.json 3       Oxyria,Dryas,Cochgroen
# OG0011456_tree.txt_prank_unique.nxh.ABSREL.json 3       Dryas,Oxyria,Draba
# OG0012809_tree.txt_prank_unique.nxh.ABSREL.json 3       Dryas,Draba,Cochgroen

#---------------------
# Check for paralogs

echo -e "file\tgene\tcorrected_p\tfull_adaptive_model\tnonsyn_subs\tsyn_subs" > all_genes_prank.tsv

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
done < ABSREL_prank_nonzero_sorted.txt >> all_genes_prank.tsv

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
    single_copy_file="OG_single_copy_significant_prank.tsv"
    gene_names_file="OG_single_copy_significant_gene_names_prank.txt"
    uniqueOG_file="OG_single_copy_significant_uniqueOG_prank.txt"
    species_counts_file="OG_single_copy_species_counts_prank.txt"

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
	OG_single_copy_species_counts_sorted_prank.tsv")
}
' all_genes_prank.tsv

more OG_single_copy_species_counts_sorted_prank.tsv

# OG      Num_species_with_sig_gene       Species
# OG0011114_tree.txt_prank_unique.nxh.ABSREL.json 3       Dryas,Oxyria,Coch
# OG0011456_tree.txt_prank_unique.nxh.ABSREL.json 3       Dryas,Oxyria,Draba
# OG0005895_tree.txt_prank_unique.nxh.ABSREL.json 2       Dryas,Oxyria
# OG0007572_tree.txt_prank_unique.nxh.ABSREL.json 2       Dryas,Coch
# OG0007884_tree.txt_prank_unique.nxh.ABSREL.json 2       Dryas,Oxyria
# OG0008152_tree.txt_prank_unique.nxh.ABSREL.json 2       Draba,Coch
# OG0008378_tree.txt_prank_unique.nxh.ABSREL.json 2       Dryas,Draba
# OG0008499_tree.txt_prank_unique.nxh.ABSREL.json 2       Dryas,Draba

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
    single_copy_file="OG_paralogs_significant_prank.tsv"
    gene_names_file="OG_paralogs_significant_gene_names_prank.txt"
    uniqueOG_file="OG_paralogs_significant_uniqueOG_prank.txt"
    species_counts_file="OG_paralogs_species_counts_prank.txt"

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
	mv tmp_header OG_paralogs_species_counts_sorted_prank.tsv")
}
' all_genes_prank.tsv

more OG_paralogs_species_counts_sorted_prank.tsv


# OG      Num_species_with_sig_gene       Species
# OG0005531_tree.txt_prank_unique.nxh.ABSREL.json 3       Coch,Oxyria,Draba
# OG0011114_tree.txt_prank_unique.nxh.ABSREL.json 3       Coch,Dryas,Oxyria
# OG0011456_tree.txt_prank_unique.nxh.ABSREL.json 3       Dryas,Oxyria,Draba
# OG0012809_tree.txt_prank_unique.nxh.ABSREL.json 3       Coch,Dryas,Draba
# OG0000886_tree.txt_prank_unique.nxh.ABSREL.json 2       Coch,Dryas
# OG0000896_tree.txt_prank_unique.nxh.ABSREL.json 2       Coch,Draba
# OG0001197_tree.txt_prank_unique.nxh.ABSREL.json 2       Oxyria,Draba
# OG0001403_tree.txt_prank_unique.nxh.ABSREL.json 2       Coch,Oxyria
# OG0001639_tree.txt_prank_unique.nxh.ABSREL.json 2       Dryas,Oxyria
# OG0002655_tree.txt_prank_unique.nxh.ABSREL.json 2       Coch,Dryas
# OG0002705_tree.txt_prank_unique.nxh.ABSREL.json 2       Dryas,Oxyria
# OG0002842_tree.txt_prank_unique.nxh.ABSREL.json 2       Dryas,Draba
# OG0003138_tree.txt_prank_unique.nxh.ABSREL.json 2       Dryas,Draba
# OG0003752_tree.txt_prank_unique.nxh.ABSREL.json 2       Dryas,Draba
# OG0003915_tree.txt_prank_unique.nxh.ABSREL.json 2       Dryas,Draba
# OG0003943_tree.txt_prank_unique.nxh.ABSREL.json 2       Dryas,Oxyria
# OG0003944_tree.txt_prank_unique.nxh.ABSREL.json 2       Oxyria,Draba
# OG0003965_tree.txt_prank_unique.nxh.ABSREL.json 2       Coch,Dryas


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
 OG_paralogs_species_counts_sorted_prank.tsv
# 1 2007
# 2 144
# 3 4

awk '{count[$2]++} END {for (val in count) print val, count[val]}' \
 OG_single_copy_species_counts_sorted_prank.tsv
# 1 872
# 2 61
# 3 2

###################
# Per species counts

grep -c 'Coch' OG_single_copy_species_counts_sorted_prank.tsv
grep -c 'Dryas' OG_single_copy_species_counts_sorted_prank.tsv
grep -c 'Oxyria' OG_single_copy_species_counts_sorted_prank.tsv
grep -c 'Draba' OG_single_copy_species_counts_sorted_prank.tsv
# 128
# 437
# 172
# 263

grep -c 'Coch' OG_paralogs_species_counts_sorted_prank.tsv
grep -c 'Dryas' OG_paralogs_species_counts_sorted_prank.tsv
grep -c 'Oxyria' OG_paralogs_species_counts_sorted_prank.tsv
grep -c 'Draba' OG_paralogs_species_counts_sorted_prank.tsv
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
