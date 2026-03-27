#########################
# Running Gene Ontology enrichment on guidance results
# Mar 2026
############################

# Narval1
tmux new-session -s total
tmux attach-session -t total

cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/Gene_ontology

cp /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees/Arctic_trees/all_genes_guidance.tsv . 

more all_genes_guidance.tsv
# file    gene    corrected_p     full_adaptive_model     nonsyn_subs     syn_subs
# OG0000042_tree.txt_unique.nxh.ABSREL.json       106297419       NA      0.0554937291668885      1E-10   0.0554937291668885
# OG0000042_tree.txt_unique.nxh.ABSREL.json       106302014       NA      0.07021517831434533     1E-10   0.07021517831434533
# OG0015441_tree.txt_unique.nxh.ABSREL.json       DOCTH0_CHR800006070     0.02062170869551722     1.588682372205865       1.588431437175698       0.0002509350301664349
# OG0015441_tree.txt_unique.nxh.ABSREL.json       DOCTH0_CHR900000114     0.2253137037661492      0.2083161911023486      0.2083161911023486      1E-10
# OG0015441_tree.txt_unique.nxh.ABSREL.json       DOCTH0_CHR900001109     0.0003764289048757696   0.7361764732201124      0.7288597604512717      0.007316712768843253
# OG0015441_tree.txt_unique.nxh.ABSREL.json       NODE10  NA      0.01747597313518143     0.008870791092959381    0.008605182042222015
# OG0015441_tree.txt_unique.nxh.ABSREL.json       NODE14  NA      0.4242486888707904      0.424219634219517       0.00002905465127344884
# OG0015441_tree.txt_unique.nxh.ABSREL.json       NODE16  NA      0.1181149077389998      1E-10   0.1181149077389998
# OG0015441_tree.txt_unique.nxh.ABSREL.json       NODE4   NA      0.03224344260775407     0.009815819793640891    0.02242762281411316
# OG0015441_tree.txt_unique.nxh.ABSREL.json       NODE9   NA      0.08640187693671385     0.04310605799912858     0.04329581893758524
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

mv all_genes_guidance.tsv baseline_guidance_all_genes.tsv

more Total_interproscan_output_edited3.tsv


##########################
# Join files in R 


module load StdEnv/2020 r/4.2.2

export R_LIBS_USER=$HOME/R/x86_64-pc-linux-gnu-library/4.2
mkdir -p $R_LIBS_USER

R

library(dplyr)
library(tidyr)
library(stringr)

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

unique(gene_ont$spp)

#-------------------------------
# Load all genes list

path="/home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/Gene_ontology/"
all_genes <- read.delim(paste0(path,"/baseline_guidance_all_genes.tsv"), header = TRUE, sep = "\t")

nrow(all_genes)
# 626 371

colnames(all_genes) 

head(all_genes)

all_genes$orthogroup <- sub("_tree\\.txt_unique\\.nxh\\.ABSREL\\.json$", "", all_genes$file)

#-------------------------
# Fix capitalization change 
# Also fix . to _ and - to _

all_genes <- all_genes %>%
  mutate(gene = gene %>%
           toupper() %>%
           str_replace_all("[.-]", "_"))

gene_ont <- gene_ont %>%
  mutate(gene = gene %>%
           toupper() %>%
           str_replace_all("[.-]", "_"))


#----------------------------
# Join sig_genes with GO info

merged_data <- all_genes %>% left_join(gene_ont, by = "gene")

head(merged_data, 10)
colnames(merged_data)

write.table(merged_data, "InterProscan_guidance_all_genes.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# remove descrip column
merged_data_nodescrip <- merged_data[,-6]
head(merged_data_nodescrip, 100)

# remove rows with no GO terms
merged_data_filtered <- merged_data %>%
  filter(!is.na(GOterm))

head(merged_data_filtered, 10)

gene_ont <- merged_data_filtered

#-------------------------------------------
# Need to compare to frequencies in whole genome - ErmineJ
# Need to get list of GO terms that were in original ABSREL test

#subset Interproscan info by spp - to get null model with all genes

unique(gene_ont$spp)
# [2] "Cochlearia_groenlandica_interproscan_output.tsv"
# [3] "Draba_nivalis_interproscan_output.tsv"
# [4] "Dryas_octopetala_interproscan_output.tsv"
# [5] "Oxyria_digyna_H1_interproscan_output.tsv"

###########################
# format for ermineJ
# https://erminej.msl.ubc.ca/help/input-files/
# https://erminej.msl.ubc.ca/help/input-files/gene-annotations/

# Make mapping file and genesets

# Oxyria
gene_ont_Oxydig <- gene_ont[which(gene_ont$spp=="Oxyria_digyna_H1_interproscan_output.tsv"),]
write.table(gene_ont_Oxydig, "Oxydig_interproscan_guidance.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

gene_ont_Oxydig1 <- gene_ont_Oxydig %>%
  mutate(gene2 = gene) %>%
  select(gene, gene2, descrip, GOterm)
write.table(gene_ont_Oxydig1, "Oxydig_GO_mappings_guidance.ermineJ.txt", sep = "\t", quote = FALSE, row.names = FALSE)

gene_ont_Oxydig2 <- gene_ont_Oxydig %>%
  select(gene, corrected_p)
write.table(gene_ont_Oxydig2, "Oxydig_genescores_guidance", sep = "\t", quote = FALSE, col.names=TRUE, row.names = FALSE)

#------------------------------------------
# Drabaniv
gene_ont_Drabaniv <- gene_ont[which(gene_ont$spp=="Draba_nivalis_interproscan_output.tsv"),]
write.table(gene_ont_Drabaniv, "Drabaniv_interproscan_guidance.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

gene_ont_Drabaniv1 <- gene_ont_Drabaniv %>%
  mutate(gene2 = gene) %>%
  select(gene, gene2, descrip, GOterm)
write.table(gene_ont_Drabaniv1, "Drabaniv_GO_mappings_guidance.ermineJ.txt", sep = "\t", quote = FALSE, row.names = FALSE)

gene_ont_Drabaniv2 <- gene_ont_Drabaniv %>%
  select(gene, corrected_p)
write.table(gene_ont_Drabaniv2, "Drabaniv_genescores_guidance", sep = "\t", quote = FALSE, col.names=TRUE, row.names = FALSE)

#------------------------------------------
# Dryasoct
gene_ont_Dryasoct <- gene_ont[which(gene_ont$spp=="Dryas_octopetala_interproscan_output.tsv"),]
write.table(gene_ont_Dryasoct, "Dryasoct_interproscan_guidance.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

gene_ont_Dryasoct1 <- gene_ont_Dryasoct %>%
  mutate(gene2 = gene) %>%
  select(gene, gene2, descrip, GOterm)
write.table(gene_ont_Dryasoct1, "Dryasoct_GO_mappings_guidance.ermineJ.txt", sep = "\t", quote = FALSE, row.names = FALSE)

gene_ont_Dryasoct2 <- gene_ont_Dryasoct %>%
  select(gene, corrected_p,)
write.table(gene_ont_Dryasoct2, "Dryasoct_genescores_guidance", sep = "\t", quote = FALSE, col.names=TRUE, row.names = FALSE)

#------------------------------------------
# Cochgro
gene_ont_Cochgro <- gene_ont[which(gene_ont$spp=="Cochlearia_groenlandica_interproscan_output.tsv"),]
write.table(gene_ont_Cochgro, "Cochgro_interproscan_guidance.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

gene_ont_Cochgro1 <- gene_ont_Cochgro %>%
  mutate(gene2 = gene) %>%
  select(gene, gene2, descrip, GOterm)
write.table(gene_ont_Cochgro1, "Cochgro_GO_mappings_guidance.ermineJ.txt", sep = "\t", quote = FALSE, row.names = FALSE)

gene_ont_Cochgro2 <- gene_ont_Cochgro %>%
  select(gene, corrected_p)
write.table(gene_ont_Cochgro2, "Cochgro_genescores_guidance", sep = "\t", quote = FALSE, col.names=TRUE, row.names = FALSE)

##################################
# Remove gene families that have paralogs in Arctic species

# Oxyria
single_gene_ont_Oxydig <- gene_ont_Oxydig %>%
  group_by(orthogroup) %>%
  filter(n_distinct(gene) == 1) %>%
  ungroup()

single_gene_ont_Oxydig1 <- single_gene_ont_Oxydig %>%
  mutate(gene2 = gene) %>%
  select(gene, gene2, descrip, GOterm)
write.table(single_gene_ont_Oxydig1, "single_Oxydig_GO_mappings_guidance.ermineJ.txt", sep = "\t", quote = FALSE, row.names = FALSE)

single_gene_ont_Oxydig2 <- single_gene_ont_Oxydig %>%
  select(gene, corrected_p)
write.table(single_gene_ont_Oxydig2, "single_Oxydig_genescores_guidance", sep = "\t", quote = FALSE, col.names=TRUE, row.names = FALSE)

#------------------------------------------
# Drabaniv
single_gene_ont_Drabaniv <- gene_ont_Drabaniv %>%
  group_by(orthogroup) %>%
  filter(n_distinct(gene) == 1) %>%
  ungroup()

single_gene_ont_Drabaniv1 <- single_gene_ont_Drabaniv %>%
  mutate(gene2 = gene) %>%
  select(gene, gene2, descrip, GOterm)
write.table(single_gene_ont_Drabaniv1, "single_Drabaniv_GO_mappings_guidance.ermineJ.txt", sep = "\t", quote = FALSE, row.names = FALSE)

single_gene_ont_Drabaniv2 <- single_gene_ont_Drabaniv %>%
  select(gene, corrected_p)
write.table(single_gene_ont_Drabaniv2, "single_Drabaniv_genescores_guidance", sep = "\t", quote = FALSE, col.names=TRUE, row.names = FALSE)

#------------------------------------------
# Dryasoct
single_gene_ont_Dryasoct <- gene_ont_Dryasoct %>%
  group_by(orthogroup) %>%
  filter(n_distinct(gene) == 1) %>%
  ungroup()

single_gene_ont_Dryasoct1 <- single_gene_ont_Dryasoct %>%
  mutate(gene2 = gene) %>%
  select(gene, gene2, descrip, GOterm)
write.table(single_gene_ont_Dryasoct1, "single_Dryasoct_GO_mappings_guidance.ermineJ.txt", sep = "\t", quote = FALSE, row.names = FALSE)

single_gene_ont_Dryasoct2 <- single_gene_ont_Dryasoct %>%
  select(gene, corrected_p,)
write.table(single_gene_ont_Dryasoct2, "single_Dryasoct_genescores_guidance", sep = "\t", quote = FALSE, col.names=TRUE, row.names = FALSE)

#------------------------------------------
# Cochgro
single_gene_ont_Cochgro <- gene_ont_Cochgro %>%
  group_by(orthogroup) %>%
  filter(n_distinct(gene) == 1) %>%
  ungroup()

single_gene_ont_Cochgro1 <- single_gene_ont_Cochgro %>%
  mutate(gene2 = gene) %>%
  select(gene, gene2, descrip, GOterm)
write.table(single_gene_ont_Cochgro1, "single_Cochgro_GO_mappings_guidance.ermineJ.txt", sep = "\t", quote = FALSE, row.names = FALSE)

single_gene_ont_Cochgro2 <- single_gene_ont_Cochgro %>%
  select(gene, corrected_p)
write.table(single_gene_ont_Cochgro2, "single_Cochgro_genescores_guidance", sep = "\t", quote = FALSE, col.names=TRUE, row.names = FALSE)


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
# narval1
tmux new-session -s Enrichment
tmux attach-session -t Enrichment

# https://erminej.msl.ubc.ca/help/tutorials/running-an-analysis-ora/

salloc -c1 --time 3:00:00 --mem 120000m --account def-rieseber

cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/Gene_ontology

ERMINEJ_HOME=/home/celphin/ermineJ-3.2
export JAVA_HOME=/cvmfs/soft.computecanada.ca/easybuild/software/2020/Core/java/13.0.2/

module load StdEnv/2020 java/13.0.2

#-------------------
# Scores are p-values

for taxon in Oxydig Cochgro Dryasoct Drabaniv ; do $ERMINEJ_HOME/bin/ermineJ.sh \
-a "$taxon"_GO_mappings_guidance.ermineJ.txt \
-s "$taxon"_genescores_guidance \
-c /home/celphin/ermineJ.data/go.obo \
--genesOut --logTrans -n ORA \
-o "$taxon"_ORA_guidance.ermine.results -y 5 --mtc FDR ; done


#------------------------
# Run again removing any genes that are paralogs in Arctic species
for taxon in Oxydig Cochgro Dryasoct Drabaniv ; do $ERMINEJ_HOME/bin/ermineJ.sh \
-a single_"$taxon"_GO_mappings_guidance.ermineJ.txt \
-s single_"$taxon"_genescores_guidance \
-c /home/celphin/ermineJ.data/go.obo \
--genesOut --logTrans -n ORA \
-o single_"$taxon"_ORA_guidance.ermine.results -y 5 --mtc FDR ; done


###########################################
# Explore data

# Extract GO ID 
more Oxydig_ORA_guidance.ermine.results

# Broad themes
# Oxyria glycosylation
# Dryas DNA repair
# Draba defense Response
# Coch Oxidioreductase 

#-----------------------------

awk -F'\t' '$1=="!" {print $2,$3,$6,$7}' OFS='\t' Oxydig_ORA_guidance.ermine.results \
| sort -t $'\t' -k4,4g | head -n 30

lipid binding   GO:0008289      38.000  7.983e-04
COPII vesicle coat      GO:0030127      7.000   2.693e-03
positive regulation of growth   GO:0045927      7.000   2.693e-03
aminoacyltransferase activity   GO:0016755      32.000  3.357e-03
ubiquitin-like protein transferase activity     GO:0019787      32.000  3.357e-03
lipid kinase activity   GO:0001727      13.000  3.927e-03
endoplasmic reticulum to Golgi vesicle-mediated transport       GO:0006888     11.000   4.241e-03
regulation of response to stress        GO:0080134      5.000   4.243e-03
transcription corepressor activity      GO:0003714      5.000   4.243e-03
cell cycle process      GO:0022402      17.000  4.798e-03
ubiquitin-protein transferase activity  GO:0004842      29.000  4.967e-03
chromatin organization  GO:0006325      12.000  7.8e-03
clathrin binding        GO:0030276      7.000   8.59e-03
GTPase binding  GO:0051020      10.000  8.893e-03
small GTPase binding    GO:0031267      10.000  8.893e-03
phosphatidylinositol kinase activity    GO:0052742      9.000   9.203e-03
cytoskeletal motor activity     GO:0003774      15.000  9.703e-03
microtubule-based movement      GO:0007018      15.000  9.703e-03
microtubule motor activity      GO:0003777      15.000  9.703e-03
regulation of response to stimulus      GO:0048583      15.000  9.703e-03
phospholipid binding    GO:0005543      23.000  0.01098403
Golgi vesicle transport GO:0048193      20.000  0.01160645
aspartic-type endopeptidase activity    GO:0004190      13.000  0.01219857
aspartic-type peptidase activity        GO:0070001      13.000  0.01219857
UDP-glycosyltransferase activity        GO:0008194      33.000  0.01222347
chromatin remodeling    GO:0006338      10.000  0.0166458
monoatomic anion transport      GO:0006820      10.000  0.0166458
phosphatidylinositol phospholipase C activity   GO:0004435      5.000   0.01834862
phospholipase C activity        GO:0004629      5.000   0.01834862
protein-DNA complex organization        GO:0071824      13.000  0.01956287


awk -F'\t' '$1=="!" {print $2,$3,$6,$7}' OFS='\t' Dryasoct_ORA_guidance.ermine.results \
| sort -t $'\t' -k4,4g | head -n 30
ubiquitin-protein transferase activity  GO:0004842      39.000  5.006e-05
aminoacyltransferase activity   GO:0016755      41.000  5.568e-05
ubiquitin-like protein transferase activity     GO:0019787      40.000  9.168e-05
transcription regulator activity        GO:0140110      46.000  1.698e-04
chromosome organization GO:0051276      21.000  1.983e-04
phosphatidylinositol binding    GO:0035091      16.000  3.404e-04
cellular response to stress     GO:0033554      41.000  6.299e-04
regulation of transcription by RNA polymerase II        GO:0006357      19.000 1.477e-03
DNA-binding transcription factor activity       GO:0003700      31.000  1.5e-03
ubiquitin protein ligase activity       GO:0061630      15.000  1.721e-03
cellular response to stimulus   GO:0051716      42.000  2.021e-03
ubiquitin-like protein ligase activity  GO:0061659      16.000  2.056e-03
nucleus GO:0005634      52.000  2.095e-03
nuclear protein-containing complex      GO:0140513      45.000  2.635e-03
DNA recombination       GO:0006310      13.000  2.678e-03
polysaccharide binding  GO:0030247      10.000  2.96e-03
phospholipid binding    GO:0005543      20.000  3.07e-03
sequence-specific DNA binding   GO:0043565      20.000  3.07e-03
DNA damage response     GO:0006974      35.000  3.409e-03
endocytosis     GO:0006897      8.000   3.976e-03
import into cell        GO:0098657      8.000   3.976e-03
enzyme activator activity       GO:0008047      16.000  4.24e-03
double-strand break repair      GO:0006302      6.000   4.252e-03
ubiquitin binding       GO:0043130      6.000   4.252e-03
ubiquitin-like protein binding  GO:0032182      6.000   4.252e-03
DNA repair      GO:0006281      34.000  5.066e-03
lipid binding   GO:0008289      33.000  7.448e-03
cell cycle process      GO:0022402      16.000  7.999e-03
molecular function activator activity   GO:0140677      17.000  8.39e-03
oxidoreductase activity, acting on the CH-CH group of donors    GO:0016627     25.000   8.889e-03

awk -F'\t' '$1=="!" {print $2,$3,$6,$7}' OFS='\t' Cochgro_ORA_guidance.ermine.results \
| sort -t $'\t' -k4,4g | head -n 30
endoplasmic reticulum to Golgi vesicle-mediated transport       GO:0006888      15.000  1.622e-03
guanosine tetraphosphate metabolic process      GO:0015969      5.000   3.039e-03
purine ribonucleoside bisphosphate metabolic process    GO:0034035      5.000   3.039e-03
N-acetyltransferase activity    GO:0008080      11.000  4.752e-03
regulation of phosphate metabolic process       GO:0019220      6.000   4.876e-03
regulation of phosphorus metabolic process      GO:0051174      6.000   4.876e-03
DNA-binding transcription factor activity       GO:0003700      41.000  6.453e-03
N-acyltransferase activity      GO:0016410      13.000  6.581e-03
sequence-specific DNA binding   GO:0043565      26.000  9.638e-03
carbon-nitrogen ligase activity, with glutamine as amido-N-donor        GO:0016884      7.000  0.01379173
regulation of developmental process     GO:0050793      7.000   0.01379173
branched-chain amino acid biosynthetic process  GO:0009082      6.000   0.01430692
sulfur compound binding GO:1901681      6.000   0.01430692
transcription coactivator activity      GO:0003713      6.000   0.01430692
rRNA metabolic process  GO:0016072      18.000  0.01497825
anatomical structure development        GO:0048856      9.000   0.02060988
branched-chain amino acid metabolic process     GO:0009081      9.000   0.02060988
translational initiation        GO:0006413      8.000   0.02400244
cysteine-type peptidase activity        GO:0008234      17.000  0.03351397
regulation of transcription by RNA polymerase II        GO:0006357      17.000  0.03351397
positive regulation of DNA-templated transcription      GO:0045893      9.000   0.03419093
positive regulation of RNA biosynthetic process GO:1902680      9.000   0.03419093
positive regulation of RNA metabolic process    GO:0051254      9.000   0.03419093
Golgi apparatus GO:0005794      5.000   0.03495789
negative regulation of cell communication       GO:0010648      5.000   0.03495789
negative regulation of signaling        GO:0023057      5.000   0.03495789
negative regulation of signal transduction      GO:0009968      5.000   0.03495789
regulation of multicellular organismal development      GO:2000026      5.000   0.03495789
regulation of multicellular organismal process  GO:0051239      5.000   0.03495789
regulation of post-embryonic development        GO:0048580      5.000   0.03495789


awk -F'\t' '$1=="!" {print $2,$3,$6,$7}' OFS='\t' Drabaniv_ORA_guidance.ermine.results \
| sort -t $'\t' -k4,4g | head -n 30
DNA-binding transcription factor activity       GO:0003700      50.000  9.63e-05
cytoskeletal motor activity     GO:0003774      19.000  2.406e-03
microtubule-based movement      GO:0007018      19.000  2.406e-03
microtubule motor activity      GO:0003777      19.000  2.406e-03
tubulin binding GO:0015631      31.000  2.819e-03
GTPase activator activity       GO:0005096      13.000  3.001e-03
lipid binding   GO:0008289      36.000  3.211e-03
microtubule-based process       GO:0007017      35.000  3.578e-03
GTPase regulator activity       GO:0030695      19.000  4.001e-03
nucleoside-triphosphatase regulator activity    GO:0060589      19.000  4.001e-03
MCM complex     GO:0042555      5.000   5.193e-03
microtubule binding     GO:0008017      27.000  6.118e-03
catalytic activity, acting on DNA       GO:0140097      37.000  6.906e-03
phosphatidylinositol binding    GO:0035091      16.000  8.32e-03
DNA duplex unwinding    GO:0032508      6.000   8.89e-03
DNA geometric change    GO:0032392      6.000   8.89e-03
negative regulation of response to stimulus     GO:0048585      6.000   8.89e-03
peptidyl-arginine methylation   GO:0018216      6.000   8.89e-03
peptidyl-arginine modification  GO:0018195      6.000   8.89e-03
protein methyltransferase activity      GO:0008276      14.000  9.761e-03
ATP-dependent activity, acting on DNA   GO:0008094      25.000  0.01066891
DNA replication initiation      GO:0006270      7.000   0.01101737
N-methyltransferase activity    GO:0008170      11.000  0.01176182
positive regulation of DNA-templated transcription      GO:0045893      11.000  0.01176182
positive regulation of RNA biosynthetic process GO:1902680      11.000  0.01176182
positive regulation of RNA metabolic process    GO:0051254      11.000  0.01176182
phosphatidylinositol phosphate binding  GO:1901981      9.000   0.01234408
macromolecule methylation       GO:0043414      16.000  0.01320162
cytoskeletal protein binding    GO:0008092      35.000  0.01445144
methylation     GO:0032259      19.000  0.01450707



##############################
# Witout paralogs

awk -F'\t' '$1=="!" {print $2,$3,$6,$7}' OFS='\t' single_Oxydig_ORA_guidance.ermine.results \
| sort -t $'\t' -k4,4g | head -n 30
lipid binding   GO:0008289      17.000  6.641e-04
organelle organization  GO:0006996      33.000  3.742e-03
chromatin organization  GO:0006325      9.000   4.657e-03
protein-DNA complex organization        GO:0071824      10.000  6.038e-03
vesicle tethering complex       GO:0099023      10.000  6.038e-03
chromatin remodeling    GO:0006338      7.000   7.016e-03
positive regulation of biological process       GO:0048518      11.000  7.141e-03
clathrin binding        GO:0030276      5.000   8.832e-03
COPII vesicle coat      GO:0030127      5.000   8.832e-03
regulation of response to stress        GO:0080134      5.000   8.832e-03
vesicle coat    GO:0030120      5.000   8.832e-03
membrane coat   GO:0030117      8.000   9.978e-03
protein-containing complex binding      GO:0044877      8.000   9.978e-03
phospholipid binding    GO:0005543      9.000   0.01221389
Golgi vesicle transport GO:0048193      12.000  0.01555654
vesicle-mediated transport      GO:0016192      29.000  0.01812314
aminoacyltransferase activity   GO:0016755      21.000  0.02006525
ubiquitin-like protein transferase activity     GO:0019787      21.000  0.02006525
endoplasmic reticulum to Golgi vesicle-mediated transport       GO:0006888      8.000   0.02416611
UDP-glycosyltransferase activity        GO:0008194      12.000  0.02750306
ubiquitin-protein transferase activity  GO:0004842      19.000  0.03214733
positive regulation of growth   GO:0045927      5.000   0.03587561
nucleus GO:0005634      30.000  0.03939835
GTPase activity GO:0003924      14.000  0.04136307
regulation of metabolic process GO:0019222      46.000  0.04157287
post-translational protein modification GO:0043687      20.000  0.04263685
protein modification by small protein conjugation or removal    GO:0070647      20.000  0.04263685
inorganic anion transmembrane transporter activity      GO:0015103      6.000   0.04298312
phosphatidylinositol binding    GO:0035091      6.000   0.04298312
regulation of growth    GO:0040008      6.000   0.04298312


awk -F'\t' '$1=="!" {print $2,$3,$6,$7}' OFS='\t' single_Dryasoct_ORA_guidance.ermine.results \
| sort -t $'\t' -k4,4g | head -n 30
ubiquitin-protein transferase activity  GO:0004842      32.000  2.137e-03
aminoacyltransferase activity   GO:0016755      34.000  2.703e-03
nuclear protein-containing complex      GO:0140513      38.000  3.903e-03
ubiquitin-like protein transferase activity     GO:0019787      33.000  3.975e-03
enzyme activator activity       GO:0008047      12.000  8.549e-03
double-strand break repair      GO:0006302      6.000   0.01001149
ubiquitin binding       GO:0043130      6.000   0.01001149
ubiquitin-like protein binding  GO:0032182      6.000   0.01001149
lipid binding   GO:0008289      21.000  0.01348488
phosphatidylinositol binding    GO:0035091      11.000  0.01501765
ubiquitin protein ligase activity       GO:0061630      12.000  0.01979063
chromosome organization GO:0051276      16.000  0.02135807
vesicle-mediated transport      GO:0016192      44.000  0.02146158
cell cycle process      GO:0022402      13.000  0.02442122
molecular function activator activity   GO:0140677      13.000  0.02442122
ubiquitin-like protein ligase activity  GO:0061659      13.000  0.02442122
catalytic complex       GO:1902494      53.000  0.02519215
phospholipid binding    GO:0005543      15.000  0.03295927
transcription regulator activity        GO:0140110      34.000  0.04105949
damaged DNA binding     GO:0003684      6.000   0.04222762
endocytosis     GO:0006897      6.000   0.04222762
import into cell        GO:0098657      6.000   0.04222762
positive regulation of macromolecule metabolic process  GO:0010604      9.000   0.04416595
positive regulation of metabolic process        GO:0009893      9.000   0.04416595
nitrogen compound transport     GO:0071705      45.000  0.04832816
post-translational protein modification GO:0043687      32.000  0.05268816
protein modification by small protein conjugation or removal    GO:0070647      32.000  0.05268816
DNA recombination       GO:0006310      10.000  0.05338942
catalytic activity, acting on a tRNA    GO:0140101      22.000  0.05445519
cellular response to stress     GO:0033554      34.000  0.05462646


awk -F'\t' '$1=="!" {print $2,$3,$6,$7}' OFS='\t' single_Cochgro_ORA_guidance.ermine.results \
| sort -t $'\t' -k4,4g | head -n 30
DNA-binding transcription factor activity       GO:0003700      20.000  2.055e-03
transcription regulator activity        GO:0140110      25.000  3.222e-03
cellular homeostasis    GO:0019725      6.000   8.661e-03
nuclear transport       GO:0051169      6.000   8.661e-03
RNA modification        GO:0009451      41.000  0.01337098
vesicle tethering complex       GO:0099023      10.000  0.02162536
endoplasmic reticulum protein-containing complex        GO:0140534      5.000   0.02165055
intracellular chemical homeostasis      GO:0055082      5.000   0.02165055
nucleocytoplasmic transport     GO:0006913      5.000   0.02165055
ATP-dependent chromatin remodeler activity      GO:0140658      8.000   0.04435361
homeostatic process     GO:0042592      7.000   0.04836528
endoplasmic reticulum to Golgi vesicle-mediated transport       GO:0006888      6.000   0.05192653
intracellular monoatomic ion homeostasis        GO:0006873      4.000   0.0528149
molecular carrier activity      GO:0140104      4.000   0.0528149
regulation of developmental process     GO:0050793      4.000   0.0528149
supramolecular complex  GO:0099080      4.000   0.0528149
protein serine/threonine phosphatase activity   GO:0004722      9.000   0.06433396
establishment of localization in cell   GO:0051649      25.000  0.06452869
intracellular transport GO:0046907      25.000  0.06452869
protein dimerization activity   GO:0046983      11.000  0.07365658
regulation of DNA-templated transcription       GO:0006355      32.000  0.0789449
regulation of RNA biosynthetic process  GO:2001141      32.000  0.0789449
regulation of RNA metabolic process     GO:0051252      32.000  0.0789449
chemical homeostasis    GO:0048878      6.000   0.09231386
phosphatidylinositol binding    GO:0035091      6.000   0.09231386
cysteine-type peptidase activity        GO:0008234      9.000   0.0957229
arginine metabolic process      GO:0006525      4.000   0.11514361
nucleoside triphosphate biosynthetic process    GO:0009142      4.000   0.11514361
purine nucleoside triphosphate biosynthetic process     GO:0009145      4.000   0.11514361
purine ribonucleoside triphosphate biosynthetic process GO:0009206      4.000   0.11514361


awk -F'\t' '$1=="!" {print $2,$3,$6,$7}' OFS='\t' single_Drabaniv_ORA_guidance.ermine.results \
| sort -t $'\t' -k4,4g | head -n 30
transcription regulator activity        GO:0140110      37.000  6.751e-04
DNA-binding transcription factor activity       GO:0003700      29.000  1.565e-03
lipid metabolic process GO:0006629      52.000  2.269e-03
MCM complex     GO:0042555      5.000   0.01086875
negative regulation of response to stimulus     GO:0048585      5.000   0.01086875
cytoskeletal motor activity     GO:0003774      14.000  0.0136802
microtubule-based movement      GO:0007018      14.000  0.0136802
microtubule motor activity      GO:0003777      14.000  0.0136802
mRNA metabolic process  GO:0016071      19.000  0.01523803
N-methyltransferase activity    GO:0008170      9.000   0.01659767
protein transport       GO:0015031      30.000  0.01726122
cellular lipid metabolic process        GO:0044255      34.000  0.02006616
DNA duplex unwinding    GO:0032508      6.000   0.02010238
DNA geometric change    GO:0032392      6.000   0.02010238
establishment of protein localization   GO:0045184      31.000  0.02266444
sequence-specific DNA binding   GO:0043565      15.000  0.02408547
DNA replication initiation      GO:0006270      7.000   0.02684042
nucleocytoplasmic transport     GO:0006913      7.000   0.02684042
intracellular protein transport GO:0006886      24.000  0.0296072
cellular macromolecule localization     GO:0070727      31.000  0.03024463
macromolecule localization      GO:0033036      31.000  0.03024463
protein localization    GO:0008104      31.000  0.03024463
cellular localization   GO:0051641      42.000  0.0303193
RNA catabolic process   GO:0006401      8.000   0.03153893
establishment of protein localization to organelle      GO:0072594      11.000  0.03797651
phospholipid biosynthetic process       GO:0008654      11.000  0.03797651
protein localization to organelle       GO:0033365      11.000  0.03797651
lipid binding   GO:0008289      13.000  0.0387272
lipid catabolic process GO:0016042      5.000   0.04324914
establishment of localization in cell   GO:0051649      34.000  0.04422629


#########################
# Present in both datasets

# trafficking + membranes + lipid interactions
# secretion
# protein sorting
# stress adaptation (especially in plants)
# endomembrane system activity
# vesicle trafficking and membrane dynamics

# Protein modification + ubiquitin system
# active protein quality control
# signaling rewiring
# stress or developmental regulation

# Chromatin + nuclear regulation
# environmental response
# developmental shifts
# stress adaptation



# Chromatin/transcription changes → alter gene expression
# Ubiquitin system → removes/replaces proteins
# Vesicle trafficking → redistributes proteins/lipids

# Broad groups
# Vesicle trafficking & membrane systems
# Protein modification (ubiquitin system)
# Nuclear regulation (chromatin + transcription)
# Transport & localization
# Lipid biology
# RNA & homeostasis (Draba and Coch)


# | Module                  | Oxydig | Dryasoct | Cochgro | Drabaniv |
# | ----------------------- | ------ | -------- | ------- | -------- |
# | Vesicle trafficking     | ✅✅     | ✅        | ✅       | ✅        |
# | Ubiquitin/protein mod   | ✅      | ✅✅       | ⚪       | ✅        |
# | Chromatin/transcription | ✅      | ✅        | ✅✅      | ✅✅       |
# | Transport/localization  | ✅      | ✅        | ✅       | ✅        |
# | Lipid biology           | ✅      | ✅        | ⚪       | ✅✅       |
# | RNA/homeostasis         | ⚪      | ⚪        | ✅       | ⚪        |


# The strong recurrence of vesicle-related terms across all datasets is unusual and meaningful

# The ubiquitin system enrichment suggests active regulation, not just passive changes

# The transcription + chromatin signals confirm upstream control



# Membrane lipids adjusted
# → maintain fluidity

# Vesicle trafficking optimized
# → proteins reach correct locations

# Misfolded proteins cleared (ubiquitin)
# → prevent damage accumulation

# Transcription rewired
# → respond to seasonal/stress signals

# DNA repair enhanced
# → maintain integrity under stress

# # All target cellular robustness systems

# Linked in functions - all high network centrality
# lipid composition → affects membranes
# membranes → affect vesicle trafficking
# trafficking → controls protein localization
# ubiquitin → regulates protein turnover
# transcription → regulates all of the above

# See transcription terms but these are more heavily annotated
# ubiquitin ligases are large gene families and likely to appear under selection

##################
# Metalbolism - fitting RELAX results
aminoacyltransferase activity
glycosyltransferase activity
metabolic process regulation


# ######################
# Next steps?

# test for convergent enrichment across species (same GO terms repeatedly selected)

# identify shared genes under selection within these GO categories


# REVIGO

# clustering GO terms (semantic similarity)

# collapsing into GO slim categories

# plotting enrichment across datasets (heatmap)














################################
# Older

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
