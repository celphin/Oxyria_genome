#########################
# Exploring HyPhy results with guidance
# Jan 2026
############################

# other tutorial: https://bioinformaticsworkbook.org/dataAnalysis/ComparativeGenomics/Finding_Positive_Selection_With_Codeml.html#gsc.tab=0

###############################

# Check
cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees/Arctic_trees

ls *_guidance_unique.nxh | wc -l
ls *_guidance_unique.nxh.RELAX.json  | wc -l

# 10493
# 8613

# Many empty
# Find non-empty results
find . -type f -name "*_guidance_unique.nxh.RELAX.json" -size +0 > RELAX_nonempty.txt

# OG0015506_guidance_unique.nxh.RELAX.json
# OG0015532_guidance_unique.nxh.RELAX.json
# OG0001631_guidance_unique.nxh.RELAX.json
# OG0001605_guidance_unique.nxh.RELAX.json
# OG0002382_guidance_unique.nxh.RELAX.json
# OG0002724_guidance_unique.nxh.RELAX.json
# OG0000858_guidance_unique.nxh.RELAX.json
# OG0000874_guidance_unique.nxh.RELAX.json
# OG0002324_guidance_unique.nxh.RELAX.json
# OG0000203_guidance_unique.nxh.RELAX.json
# OG0013543_guidance_unique.nxh.RELAX.json
# OG0013563_guidance_unique.nxh.RELAX.json
# OG0013244_guidance_unique.nxh.RELAX.json

more OG0015506_guidance_unique.nxh.RELAX.json


##########################
# Explore results

# Narval1
tmux attach-session -t total
cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees/Arctic_trees

jq -r '
  [
    input_filename,
    ."test results".LRT,
    ."test results"."p-value",
    ."test results"."relaxation or intensification parameter",
    (."tested"."0" | to_entries | map(select(.value=="Test") | .key) | join(",")),
    (."tested"."0" | to_entries | map(select(.value=="Reference") | .key) | join(","))
  ] | @tsv
' *_guidance_unique.nxh.RELAX.json > Total_RELAX_guidance_results.txt


#----------------------------
# filter for only those with positive selection

# intensified selection
awk '$3 < 0.05 && $4 > 1.2' Total_RELAX_guidance_results.txt > guidance_intensified.txt
wc -l guidance_intensified.txt

# relaxed selection
awk '$3 < 0.05 && $4 < 0.8' Total_RELAX_guidance_results.txt > guidance_relaxed.txt
wc -l guidance_relaxed.txt

# 62 guidance_intensified.txt
# 52 guidance_relaxed.txt

#--------------------------
# intensified selection
awk '$3 < 0.05 && $4 > 1.01' Total_RELAX_guidance_results.txt > guidance_intensified_full.txt
wc -l guidance_intensified_full.txt

# relaxed selection
awk '$3 < 0.05 && $4 < 0.99' Total_RELAX_guidance_results.txt > guidance_relaxed_full.txt
wc -l guidance_relaxed_full.txt

# 69 guidance_intensified_full.txt
# 54 guidance_relaxed_full.txt

#----------------------------
# Check which Arctic species are involved

awk -F'\t' '
BEGIN{OFS="\t"}
{
  split($5, genes, ",")
  delete seen

  for (i in genes) {
    g = genes[i]

    if (g ~ /^DOCTH0_CHR/) species="Dryas"
    else if (g ~ /^OXYRIA_NCBI_CHR/) species="Oxyria"
    else if (g ~ /^MAKER_/ || g ~ /^SNAP_MASKED_/) species="Draba"
    else if (g ~ /^G[0-9]+_T[0-9]+/) species="Coch"
    else continue

    seen[species] = 1
  }

  count = 0
  species_list=""
  for (s in seen) {
    count++
    species_list = species_list (species_list ? "," : "") s
  }

  print $1, count, species_list
}
' guidance_intensified.txt > intensified_species_counts.txt


sort -t $'\t' -k2,2nr intensified_species_counts.txt > intensified_species_counts_sorted.txt

more intensified_species_counts_sorted.txt

# OG0012366_guidance_unique.nxh.RELAX.json        4       Coch,Dryas,Oxyria,Draba
# OG0001914_guidance_unique.nxh.RELAX.json        3       Coch,Dryas,Oxyria
# OG0005808_guidance_unique.nxh.RELAX.json        3       Coch,Oxyria,Draba
# OG0000748_guidance_unique.nxh.RELAX.json        2       Coch,Oxyria
# OG0002336_guidance_unique.nxh.RELAX.json        2       Coch,Oxyria
# OG0002753_guidance_unique.nxh.RELAX.json        2       Dryas,Oxyria
# OG0003167_guidance_unique.nxh.RELAX.json        2       Coch,Oxyria
# OG0004606_guidance_unique.nxh.RELAX.json        2       Coch,Draba
# OG0005375_guidance_unique.nxh.RELAX.json        2       Dryas,Oxyria
# OG0008256_guidance_unique.nxh.RELAX.json        2       Coch,Oxyria
# OG0008931_guidance_unique.nxh.RELAX.json        2       Dryas,Draba
# OG0010248_guidance_unique.nxh.RELAX.json        2       Coch,Dryas
# OG0010271_guidance_unique.nxh.RELAX.json        2       Coch,Draba
# OG0011082_guidance_unique.nxh.RELAX.json        2       Coch,Draba
# OG0011160_guidance_unique.nxh.RELAX.json        2       Oxyria,Draba
# OG0012671_guidance_unique.nxh.RELAX.json        2       Dryas,Oxyria
# OG0013338_guidance_unique.nxh.RELAX.json        2       Dryas,Draba
# OG0014344_guidance_unique.nxh.RELAX.json        2       Dryas,Oxyria
# OG0015804_guidance_unique.nxh.RELAX.json        2       Coch,Draba
# OG0016396_guidance_unique.nxh.RELAX.json        2       Coch,Draba
# OG0016658_guidance_unique.nxh.RELAX.json        2       Coch,Draba
# OG0016933_guidance_unique.nxh.RELAX.json        2       Dryas,Oxyria
# OG0020127_guidance_unique.nxh.RELAX.json        2       Coch,Dryas
# OG0021301_guidance_unique.nxh.RELAX.json        2       Coch,Draba

#-----------------------
awk -F'\t' '
BEGIN{OFS="\t"}
{
  split($5, genes, ",")
  delete seen

  for (i in genes) {
    g = genes[i]

    if (g ~ /^DOCTH0_CHR/) species="Dryas"
    else if (g ~ /^OXYRIA_NCBI_CHR/) species="Oxyria"
    else if (g ~ /^MAKER_/ || g ~ /^SNAP_MASKED_/) species="Draba"
    else if (g ~ /^G[0-9]+_T[0-9]+/) species="Coch"
    else continue

    seen[species] = 1
  }

  count = 0
  species_list=""
  for (s in seen) {
    count++
    species_list = species_list (species_list ? "," : "") s
  }

  print $1, count, species_list
}
' guidance_relaxed.txt > relaxed_species_counts.txt

sort -t $'\t' -k2,2nr relaxed_species_counts.txt > relaxed_species_counts_sorted.txt

more relaxed_species_counts_sorted.txt

# OG0000199_guidance_unique.nxh.RELAX.json        3       Coch,Dryas,Oxyria
# OG0002019_guidance_unique.nxh.RELAX.json        3       Dryas,Oxyria,Draba
# OG0000916_guidance_unique.nxh.RELAX.json        2       Coch,Oxyria
# OG0001631_guidance_unique.nxh.RELAX.json        2       Dryas,Oxyria
# OG0001665_guidance_unique.nxh.RELAX.json        2       Dryas,Oxyria
# OG0004232_guidance_unique.nxh.RELAX.json        2       Dryas,Oxyria
# OG0005678_guidance_unique.nxh.RELAX.json        2       Dryas,Oxyria
# OG0009496_guidance_unique.nxh.RELAX.json        2       Coch,Draba
# OG0009755_guidance_unique.nxh.RELAX.json        2       Dryas,Oxyria
# OG0012681_guidance_unique.nxh.RELAX.json        2       Dryas,Oxyria
# OG0014019_guidance_unique.nxh.RELAX.json        2       Dryas,Oxyria
# OG0014472_guidance_unique.nxh.RELAX.json        2       Coch,Draba
# OG0014630_guidance_unique.nxh.RELAX.json        2       Oxyria,Draba
# OG0014873_guidance_unique.nxh.RELAX.json        2       Coch,Dryas
# OG0016108_guidance_unique.nxh.RELAX.json        2       Dryas,Oxyria
# OG0017264_guidance_unique.nxh.RELAX.json        2       Coch,Draba


##########################
# Search for these OG and genes in ABSREL results
cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/Gene_ontology

# in 3 or 4 of the 4 arctic species
# intensified 
# OG0012366_guidance_unique.nxh.RELAX.json        4       Coch,Dryas,Oxyria,Draba
# OG0001914_guidance_unique.nxh.RELAX.json        3       Coch,Dryas,Oxyria
# OG0005808_guidance_unique.nxh.RELAX.json        3       Coch,Oxyria,Draba

# relaxed
# OG0000199_guidance_unique.nxh.RELAX.json        3       Coch,Dryas,Oxyria
# OG0002019_guidance_unique.nxh.RELAX.json        3       Dryas,Oxyria,Draba

#-----------------------------------
grep OG0012366 baseline_guidance_all_genes.tsv

# MAKER_LG5_SNAP_GENE_85_258_MRNA_1 0.04
# G4015_T1 0.06
# DOCTH0_CHR900003261 0
# OXYRIA_NCBI_CHR100004345  1.7e-6

grep CHR900003261 Dryasoct.ermine.results
# nothing

grep CHR100004345 Oxydig.ermine.results
# carbon-carbon lyase activity #53 of them

grep G4015  Cochgro.ermine.results
#  carbon-carbon lyase activity # 43 of them

grep "MAKER_LG5_SNAP_GENE_85_258_MRNA_1" Drabaniv.ermine.results
# carbon-carbon lyase activity #  GO:0016830 # 45 of them

# all under positive selection
# Carbon-carbon lyase activity plays critical roles in metabolic 
# pathways like glycolysis

cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees/Arctic_trees

more OG0012366_tree.txt_HyPhy_relax.txt

# ((((FT01Gene23781.t1:0.02456,FEHAP231105.t1:0.016882):0.036459,
# (Oxyria_NCBI_Chr100004345{Foreground}:0.057383,
# (RtaG0041881.1:0.018107,RnoG0013552.1:0.007086):0.02
# 5078):0.021019):0.01198,Polavi_Chr100002344:0.164844):0.044558,
# (((LOC9325994:0.025893,LOC17897262:0.027542):0.006746,
# (((maker-lg5-snap-gene-85.258-mRNA-1{Foreground}:0.102749,
# TAV2_LOCUS552:0.026482):0.02024,AALP_AA1G133300:0.023597):0.010836,
# g4015.t1{Foreground}:0.040074)
# :0.001826):0.114449,(((LOC126795455:0.35078,(LOC133743809:0.015015,
# LOC133743744:0.021387):0.007143):0.03799,
# DoctH0_Chr900003261{Foreground}:0.100593):0.009759,
# ((LOC103952293:0.014215,LOC126619830:0.015099):0.062689
# ,LOC18777121:0.097959):0.009769):0.061913):0.044558);



#--------------------------------
grep OG0001914 baseline_guidance_all_genes.tsv
OXYRIA_NCBI_CHR800003891
# many more oxyria genes
G23814_T1 
DOCTH0_CHR800004810
# None under positive selection

grep Chr800003891 Total_interproscan_output_edited3.tsv 
grep Chr800004810 Total_interproscan_output_edited3.tsv 
grep g23814.t1 Total_interproscan_output_edited3.tsv 

# Unknown function

#--------------------------------
grep OG0005808 baseline_guidance_all_genes.tsv
OXYRIA_NCBI_CHR100004659 # under positive selection
OXYRIA_NCBI_CHR600001741 # under positive selection
MAKER_LG1_SNAP_GENE_67_134_MRNA_1
G29267_T1

grep OXYRIA_NCBI_CHR600001741  Oxydig.ermine.results
grep OXYRIA_NCBI_CHR100004659 Oxydig.ermine.results

grep Chr600001741 Total_interproscan_output_edited3.tsv 
grep g29267.t1 Total_interproscan_output_edited3.tsv 

# Unknown function

#--------------------------------
grep OG0000199 baseline_guidance_all_genes.tsv
# OG0001048_guidance_unique.nxh.ABSREL.json       RnoG0000199_1   NA      0.3341087961810071      0.3341087961810071 1e-10

grep RnoG0000199.1 Total_interproscan_output_edited3.tsv 
# "FKBP-type peptidyl-prolyl cis-trans isomerase domain,
# Peptidyl-prolyl cis-trans isomerase domain superfamily,
# Tetratricopeptide repeat,Tetratricopeptide-like helical 
# domain superfamily"    "GO:0003755,GO:0005515"

# Protein folding chaperone

#--------------------------------
grep OG0002019 baseline_guidance_all_genes.tsv

#Unknown

##############################
# Arctic adaptation involves tightening selection on metabolic processes
#  while allowing greater flexibility in protein folding machinery

#----------------------------
cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees/Arctic_trees

wc -l relaxed_species_counts_sorted.txt
# 52
wc -l intensified_species_counts_sorted.txt
# 62

awk '{count[$2]++} END {for (val in count) print val, count[val]}' \
 relaxed_species_counts_sorted.txt
# 1 36
# 2 14
# 3 2

awk '{count[$2]++} END {for (val in count) print val, count[val]}' \
 intensified_species_counts_sorted.txt
# 1 38
# 2 21
# 3 2
# 4 1
