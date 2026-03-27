#########################
# Running Gene Ontology enrichment
# Mar 2026
############################

# Narval1
tmux new-session -s total
tmux attach-session -t total

cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/Gene_ontology

cp /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/orthofinder/Results_Aug18/Gene_Trees/Arctic_trees/baseline_all_genes.tsv . 

more baseline_all_genes.tsv
# file    gene    corrected_p     full_adaptive_model     nonsyn_subs     syn_subs
# OG0000042_tree.txt_unique.nxh.ABSREL.json       101295001       NA      0.06709438122731588     1E-10   0.06709438122731588
# OG0000042_tree.txt_unique.nxh.ABSREL.json       101295289       NA      0.04783298929132371     1E-10   0.04783298929132371
# OG0000042_tree.txt_unique.nxh.ABSREL.json       101295580       NA      0       1E-10   1E-10
# OG0000042_tree.txt_unique.nxh.ABSREL.json       101296064       NA      0.05069316386452072     1E-10   0.05069316386452072
# OG0000042_tree.txt_unique.nxh.ABSREL.json       101300007       NA      0.1077841536586944      0.05235768716482594     0.05542646649386831
# OG0000042_tree.txt_unique.nxh.ABSREL.json       101300449       NA      0       1E-10   1E-10
# OG0000042_tree.txt_unique.nxh.ABSREL.json       101300916       NA      0.05020476492605291     0.05020476492605291     1E-10
# OG0000042_tree.txt_unique.nxh.ABSREL.json       101302942       NA      0.06005795878196007     1E-10   0.06005795878196007
# OG0000042_tree.txt_unique.nxh.ABSREL.json       101303099       NA      0.05928206485260853     0.05928206485260853     1E-10
# OG0000042_tree.txt_unique.nxh.ABSREL.json       101304970       NA      0       1E-10   1E-10
# OG0000042_tree.txt_unique.nxh.ABSREL.json       106296244       NA      0.1225864624578918      0.05023313020707065     0.07235333225082106
# OG0000042_tree.txt_unique.nxh.ABSREL.json       106297419       NA      0.0554937291668885      1E-10   0.0554937291668885
# OG0000042_tree.txt_unique.nxh.ABSREL.json       106302014       NA      0.07021517831434533     1E-10   0.07021517831434533
# OG0000042_tree.txt_unique.nxh.ABSREL.json       106303396       NA      0       1E-10   1E-10
# OG0000042_tree.txt_unique.nxh.ABSREL.json       106305074       NA      0       1E-10   1E-10
# OG0000042_tree.txt_unique.nxh.ABSREL.json       106305316       NA      0.09755979421286293     0.09755979421286293     1E-10
# OG0000042_tree.txt_unique.nxh.ABSREL.json       106315413       NA      0.1872692664976272      1E-10   0.1872692664976272
# OG0015441_tree.txt_unique.nxh.ABSREL.json       DOCTH0_CHR300006814     0.0004372556584134046   0.3516174924510933      0.3516111809241831      0.000006311526   910256882
# OG0015441_tree.txt_unique.nxh.ABSREL.json       DOCTH0_CHR400000398     0.02450561644103588     81.70688001451325       81.61750677619105       0.089373238322   13598
# OG0015441_tree.txt_unique.nxh.ABSREL.json       DOCTH0_CHR400001723     0.04809666165521614     0.7219684042673827      0.7218549227459655      0.000113481521  415862
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
all_genes <- read.delim(paste0(path,"/baseline_all_genes.tsv"), header = TRUE, sep = "\t")

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

write.table(merged_data, "InterProscan_ABSREL_all_genes.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

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
write.table(gene_ont_Oxydig, "Oxydig_interproscan_edited.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

gene_ont_Oxydig1 <- gene_ont_Oxydig %>%
  mutate(gene2 = gene) %>%
  select(gene, gene2, descrip, GOterm)
write.table(gene_ont_Oxydig1, "Oxydig_GO_mappings.ermineJ.txt", sep = "\t", quote = FALSE, row.names = FALSE)

gene_ont_Oxydig2 <- gene_ont_Oxydig %>%
  select(gene, corrected_p)
write.table(gene_ont_Oxydig2, "Oxydig_genescores", sep = "\t", quote = FALSE, col.names=TRUE, row.names = FALSE)

#------------------------------------------
# Drabaniv
gene_ont_Drabaniv <- gene_ont[which(gene_ont$spp=="Draba_nivalis_interproscan_output.tsv"),]
write.table(gene_ont_Drabaniv, "Drabaniv_interproscan_edited.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

gene_ont_Drabaniv1 <- gene_ont_Drabaniv %>%
  mutate(gene2 = gene) %>%
  select(gene, gene2, descrip, GOterm)
write.table(gene_ont_Drabaniv1, "Drabaniv_GO_mappings.ermineJ.txt", sep = "\t", quote = FALSE, row.names = FALSE)

gene_ont_Drabaniv2 <- gene_ont_Drabaniv %>%
  select(gene, corrected_p)
write.table(gene_ont_Drabaniv2, "Drabaniv_genescores", sep = "\t", quote = FALSE, col.names=TRUE, row.names = FALSE)

#------------------------------------------
# Dryasoct
gene_ont_Dryasoct <- gene_ont[which(gene_ont$spp=="Dryas_octopetala_interproscan_output.tsv"),]
write.table(gene_ont_Dryasoct, "Dryasoct_interproscan_edited.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

gene_ont_Dryasoct1 <- gene_ont_Dryasoct %>%
  mutate(gene2 = gene) %>%
  select(gene, gene2, descrip, GOterm)
write.table(gene_ont_Dryasoct1, "Dryasoct_GO_mappings.ermineJ.txt", sep = "\t", quote = FALSE, row.names = FALSE)

gene_ont_Dryasoct2 <- gene_ont_Dryasoct %>%
  select(gene, corrected_p,)
write.table(gene_ont_Dryasoct2, "Dryasoct_genescores", sep = "\t", quote = FALSE, col.names=TRUE, row.names = FALSE)

#------------------------------------------
# Cochgro
gene_ont_Cochgro <- gene_ont[which(gene_ont$spp=="Cochlearia_groenlandica_interproscan_output.tsv"),]
write.table(gene_ont_Cochgro, "Cochgro_interproscan_edited.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

gene_ont_Cochgro1 <- gene_ont_Cochgro %>%
  mutate(gene2 = gene) %>%
  select(gene, gene2, descrip, GOterm)
write.table(gene_ont_Cochgro1, "Cochgro_GO_mappings.ermineJ.txt", sep = "\t", quote = FALSE, row.names = FALSE)

gene_ont_Cochgro2 <- gene_ont_Cochgro %>%
  select(gene, corrected_p)
write.table(gene_ont_Cochgro2, "Cochgro_genescores", sep = "\t", quote = FALSE, col.names=TRUE, row.names = FALSE)

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
write.table(single_gene_ont_Oxydig1, "single_Oxydig_GO_mappings.ermineJ.txt", sep = "\t", quote = FALSE, row.names = FALSE)

single_gene_ont_Oxydig2 <- single_gene_ont_Oxydig %>%
  select(gene, corrected_p)
write.table(single_gene_ont_Oxydig2, "single_Oxydig_genescores", sep = "\t", quote = FALSE, col.names=TRUE, row.names = FALSE)

#------------------------------------------
# Drabaniv
single_gene_ont_Drabaniv <- gene_ont_Drabaniv %>%
  group_by(orthogroup) %>%
  filter(n_distinct(gene) == 1) %>%
  ungroup()

single_gene_ont_Drabaniv1 <- single_gene_ont_Drabaniv %>%
  mutate(gene2 = gene) %>%
  select(gene, gene2, descrip, GOterm)
write.table(single_gene_ont_Drabaniv1, "single_Drabaniv_GO_mappings.ermineJ.txt", sep = "\t", quote = FALSE, row.names = FALSE)

single_gene_ont_Drabaniv2 <- single_gene_ont_Drabaniv %>%
  select(gene, corrected_p)
write.table(single_gene_ont_Drabaniv2, "single_Drabaniv_genescores", sep = "\t", quote = FALSE, col.names=TRUE, row.names = FALSE)

#------------------------------------------
# Dryasoct
single_gene_ont_Dryasoct <- gene_ont_Dryasoct %>%
  group_by(orthogroup) %>%
  filter(n_distinct(gene) == 1) %>%
  ungroup()

single_gene_ont_Dryasoct1 <- single_gene_ont_Dryasoct %>%
  mutate(gene2 = gene) %>%
  select(gene, gene2, descrip, GOterm)
write.table(single_gene_ont_Dryasoct1, "single_Dryasoct_GO_mappings.ermineJ.txt", sep = "\t", quote = FALSE, row.names = FALSE)

single_gene_ont_Dryasoct2 <- single_gene_ont_Dryasoct %>%
  select(gene, corrected_p,)
write.table(single_gene_ont_Dryasoct2, "single_Dryasoct_genescores", sep = "\t", quote = FALSE, col.names=TRUE, row.names = FALSE)

#------------------------------------------
# Cochgro
single_gene_ont_Cochgro <- gene_ont_Cochgro %>%
  group_by(orthogroup) %>%
  filter(n_distinct(gene) == 1) %>%
  ungroup()

single_gene_ont_Cochgro1 <- single_gene_ont_Cochgro %>%
  mutate(gene2 = gene) %>%
  select(gene, gene2, descrip, GOterm)
write.table(single_gene_ont_Cochgro1, "single_Cochgro_GO_mappings.ermineJ.txt", sep = "\t", quote = FALSE, row.names = FALSE)

single_gene_ont_Cochgro2 <- single_gene_ont_Cochgro %>%
  select(gene, corrected_p)
write.table(single_gene_ont_Cochgro2, "single_Cochgro_genescores", sep = "\t", quote = FALSE, col.names=TRUE, row.names = FALSE)


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

# https://erminej.msl.ubc.ca/help/tutorials/running-an-analysis-ora/

salloc -c1 --time 3:00:00 --mem 120000m --account def-rieseber

cd /home/celphin/scratch/Oxyria_Positive_Selection_Test/Total_genomes/Gene_ontology

ERMINEJ_HOME=/home/celphin/ermineJ-3.2
export JAVA_HOME=/cvmfs/soft.computecanada.ca/easybuild/software/2020/Core/java/13.0.2/

module load StdEnv/2020 java/13.0.2

#-------------------
# Scores are p-values

for taxon in Oxydig Cochgro Dryasoct Drabaniv ; do $ERMINEJ_HOME/bin/ermineJ.sh \
-a "$taxon"_GO_mappings.ermineJ.txt \
-s "$taxon"_genescores \
-c /home/celphin/ermineJ.data/go.obo \
--genesOut --logTrans -n ORA \
-o "$taxon"_ORA.ermine.results -y 5 --mtc FDR ; done


#------------------------
# Run again removing any genes that are paralogs in Arctic species
for taxon in Oxydig Cochgro Dryasoct Drabaniv ; do $ERMINEJ_HOME/bin/ermineJ.sh \
-a single_"$taxon"_GO_mappings.ermineJ.txt \
-s single_"$taxon"_genescores \
-c /home/celphin/ermineJ.data/go.obo \
--genesOut --logTrans -n ORA \
-o single_"$taxon"_ORA.ermine.results -y 5 --mtc FDR ; done






###########################################
# Explore data

# Extract GO ID 
more Oxydig.ermine.results

# Broad themes
# Oxyria glycosylation
# Dryas DNA repair
# Draba defense Response
# Coch Oxidioreductase 

#-----------------------------

awk -F'\t' '$1=="!" {print $2,$3,$6,$7}' OFS='\t' Oxydig.ermine.results \
| sort -t $'\t' -k4,4g | head -n 50

# H4/H2A histone acetyltransferase complex        GO:0043189      4.000   2.188e-05
# H4 histone acetyltransferase complex    GO:1902562      4.000   2.188e-05
# NuA4 histone acetyltransferase complex  GO:0035267      4.000   2.188e-05
# histone acetyltransferase complex       GO:0000123      5.000   4.663e-05
# acetyltransferase complex       GO:1902493      5.000   9.612e-05
# protein acetyltransferase complex       GO:0031248      5.000   9.612e-05
# protein import  GO:0017038      3.000   1.739e-04
# cytoskeleton organization       GO:0007010      7.000   4.041e-04
# 1,4-alpha-glucan branching enzyme activity      GO:0003844      2.000   4.201e-03
# palmitoyltransferase activity   GO:0016409      3.000   5.888e-03
# protein polyubiquitination      GO:0000209      2.000   6.214e-03
# polysaccharide biosynthetic process     GO:0000271      4.000   9.541e-03
# macromolecule glycosylation     GO:0043413      6.000   0.0109247
# protein glycosylation   GO:0006486      6.000   0.0109247
# 4-alpha-hydroxytetrahydrobiopterin dehydratase activity GO:0008124      2.000   0.01128283
# diol biosynthetic process       GO:0034312      2.000   0.01128283
# diol metabolic process  GO:0034311      2.000   0.01128283
# tetrahydrobiopterin biosynthetic process        GO:0006729      2.000   0.01128283
# tetrahydrobiopterin metabolic process   GO:0046146      2.000   0.01128283
# pteridine-containing compound metabolic process GO:0042558      3.000   0.01183288
# glycosylation   GO:0070085      6.000   0.01474199
# microtubule binding     GO:0008017      6.000   0.01856283
# polyol biosynthetic process     GO:0046173      2.000   0.02126416
# transcription factor TFIID complex      GO:0005669      2.000   0.02126416
# glucan biosynthetic process     GO:0009250      3.000   0.02437415
# transferase complex     GO:1990234      5.000   0.02749037
# energy reserve metabolic process        GO:0006112      2.000   0.02933932
# glycogen biosynthetic process   GO:0005978      2.000   0.02933932
# glycogen metabolic process      GO:0005977      2.000   0.02933932
# phosphate ion transmembrane transport   GO:0035435      2.000   0.02933932
# carbohydrate biosynthetic process       GO:0016051      4.000   0.02961807
# ATP-dependent activity, acting on RNA   GO:0008186      4.000   0.03313381
# RNA helicase activity   GO:0003724      4.000   0.03313381
# inorganic anion transmembrane transport GO:0098661      2.000   0.03376353
# phosphate ion transport GO:0006817      2.000   0.03376353
# pteridine-containing compound biosynthetic process      GO:0042559      2.000   0.03376353
# polysaccharide metabolic process        GO:0005976      5.000   0.03742848
# phosphate transmembrane transporter activity    GO:0005315      2.000   0.03842885
# lipoate biosynthetic process    GO:0009107      2.000   0.04332332
# lipoate metabolic process       GO:0009106      2.000   0.04332332
# lipoate synthase activity       GO:0016992      2.000   0.04332332
# alcohol biosynthetic process    GO:0046165      2.000   0.04843536
# rRNA processing GO:0006364      4.000   0.06169119

awk -F'\t' '$1=="!" {print $2,$3,$6,$7}' OFS='\t' Dryasoct.ermine.results \
| sort -t $'\t' -k4,4g | head -n 30
# ATP-dependent activity, acting on DNA   GO:0008094      9.000   1.358e-03
# transferase activity, transferring sulphur-containing groups    GO:0016782      4.000   1.766e-03
# DNA helicase activity   GO:0003678      4.000   2.333e-03
# DNA recombination       GO:0006310      5.000   2.616e-03
# helicase activity       GO:0004386      8.000   2.924e-03
# rRNA metabolic process  GO:0016072      7.000   5.694e-03
# rRNA processing GO:0006364      7.000   5.694e-03
# double-strand break repair      GO:0006302      3.000   6.002e-03
# DNA replication GO:0006260      6.000   8.173e-03
# ATP-dependent DNA damage sensor activity        GO:0140664      4.000   0.01004508
# DNA damage sensor activity      GO:0140612      4.000   0.01004508
# double-stranded DNA binding     GO:0003690      7.000   0.01396006
# acylglycerol O-acyltransferase activity GO:0016411      2.000   0.02109263
# indolalkylamine biosynthetic process    GO:0046219      2.000   0.02109263
# indole-containing compound biosynthetic process GO:0042435      2.000   0.02109263
# MCM complex     GO:0042555      2.000   0.02109263
# NAD+ binding    GO:0070403      2.000   0.02109263
# trehalose metabolic process     GO:0005991      2.000   0.02109263
# tryptophan biosynthetic process GO:0000162      2.000   0.02109263
# sequence-specific double-stranded DNA binding   GO:1990837      4.000   0.02334864
# lipid kinase activity   GO:0001727      3.000   0.02348876
# palmitoyltransferase activity   GO:0016409      3.000   0.02770338
# double-strand break repair via homologous recombination GO:0000724      2.000   0.02876431
# double-stranded RNA binding     GO:0003725      2.000   0.02876431
# primary amine oxidase activity  GO:0008131      2.000   0.02876431
# quinone binding GO:0048038      2.000   0.02876431
# NAD binding     GO:0051287      5.000   0.02900502
# positive regulation of macromolecule metabolic process  GO:0010604      5.000   0.03435476
# positive regulation of metabolic process        GO:0009893      5.000   0.03435476
# ubiquitin protein ligase activity       GO:0061630      5.000   0.03724389
# S-adenosylmethionine-dependent methyltransferase activity       GO:0008757      6.000   0.04585428
# L-ascorbic acid binding GO:0031418      2.000   0.04680006
# oxidoreductase activity, acting on the CH-NH2 group of donors, oxygen as acceptor       GO:0016641      2.000   0.04680006
# chromosome organization GO:0051276      4.000   0.0481783
# dicarboxylic acid metabolic process     GO:0043648      3.000   0.05429115
# DNA replication initiation      GO:0006270      2.000   0.05699948
# mismatched DNA binding  GO:0030983      2.000   0.05699948
# phosphatidylinositol kinase activity    GO:0052742      2.000   0.05699948
# amine metabolic process GO:0009308      4.000   0.06702383
# protein processing      GO:0016485      3.000   0.06737438
# water-soluble vitamin biosynthetic process      GO:0042364      3.000   0.06737438
# amylase activity        GO:0016160      2.000   0.06788529
# malate dehydrogenase activity   GO:0016615      2.000   0.06788529
# mismatch repair GO:0006298      2.000   0.06788529

awk -F'\t' '$1=="!" {print $2,$3,$6,$7}' OFS='\t' Cochgro.ermine.results \
| sort -t $'\t' -k4,4g | head -n 30
# sexual reproduction     GO:0019953      2.000   2.509e-03
# regulation of defense response to fungus        GO:1900150      2.000   4.602e-03
# 3-hydroxyisobutyryl-CoA hydrolase activity      GO:0003860      2.000   8.808e-03
# regulation of response to biotic stimulus       GO:0002831      2.000   8.808e-03
# regulation of response to external stimulus     GO:0032101      2.000   8.808e-03
# short-chain fatty acyl-CoA dehydrogenase activity       GO:0016937      2.000   8.808e-03
# acyl-CoA dehydrogenase activity GO:0003995      2.000   0.01047874
# oxidoreductase activity, acting on the CH-CH group of donors, with a flavin as acceptor GO:0052890      2.000   0.01047874
# regulation of defense response  GO:0031347      2.000   0.01047874
# pseudouridine synthase activity GO:0009982      2.000   0.01420075
# pseudouridine synthesis GO:0001522      2.000   0.01420075
# regulation of response to stress        GO:0080134      2.000   0.02306845
# water-soluble vitamin biosynthetic process      GO:0042364      2.000   0.02816006
# vitamin biosynthetic process    GO:0009110      2.000   0.03365589
# water-soluble vitamin metabolic process GO:0006767      2.000   0.03654772
# vitamin metabolic process       GO:0006766      2.000   0.04260441
# transaminase activity   GO:0008483      2.000   0.05232947
# transferase activity, transferring nitrogenous groups   GO:0016769      2.000   0.05232947
# vacuolar transport      GO:0007034      2.000   0.05573088
# ATP-dependent chromatin remodeler activity      GO:0140658      2.000   0.06275754
# intramolecular transferase activity     GO:0016866      2.000   0.06275754


awk -F'\t' '$1=="!" {print $2,$3,$6,$7}' OFS='\t' Drabaniv.ermine.results \
| sort -t $'\t' -k4,4g | head -n 30

# methylenetetrahydrofolate dehydrogenase [NAD(P)+] activity      GO:0004486      4.000   5.085e-05
# methylenetetrahydrofolate dehydrogenase (NADP+) activity        GO:0004488      4.000   5.085e-05
# oxidoreductase activity, acting on the CH-NH group of donors, NAD or NADP as acceptor   GO:0016646      4.000   6.301e-04
# protein-containing complex binding      GO:0044877      7.000   5.304e-03
# oxidoreductase activity, acting on the CH-NH group of donors    GO:0016645      4.000   6.136e-03
# membrane organization   GO:0061024      4.000   7.083e-03
# aminopeptidase activity GO:0004177      3.000   0.01335124
# intermembrane lipid transfer    GO:0120009      2.000   0.01708648
# lipid transfer activity GO:0120013      2.000   0.01708648
# oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor   GO:0016616      6.000   0.02083865
# COPII vesicle coat      GO:0030127      2.000   0.02233409
# protein polyubiquitination      GO:0000209      2.000   0.02233409
# regulation of response to biotic stimulus       GO:0002831      2.000   0.02815217
# regulation of response to external stimulus     GO:0032101      2.000   0.02815217
# malate dehydrogenase activity   GO:0016615      2.000   0.03450193
# regulation of defense response  GO:0031347      2.000   0.03450193
# serine family amino acid metabolic process      GO:0009069      3.000   0.03787732
# cysteine metabolic process      GO:0006534      2.000   0.0413464
# translational elongation        GO:0006414      3.000   0.04189395
# oxidoreductase activity, acting on CH-OH group of donors        GO:0016614      6.000   0.04447004
# microtubule-based movement      GO:0007018      4.000   0.04856391
# microtubule motor activity      GO:0003777      4.000   0.04856391
# small-subunit processome        GO:0032040      2.000   0.0486503
# isoprenoid metabolic process    GO:0006720      4.000   0.05187672
# ATP-dependent chromatin remodeler activity      GO:0140658      3.000   0.05515992
# lipid transporter activity      GO:0005319      2.000   0.05638005
# cytoskeletal motor activity     GO:0003774      4.000   0.06253469


# Themes
| Functional Category                          | Example GO Terms (from your results)                                                                                                                           | Species Showing Enrichment          | Biological Interpretation                                                                                                                   |
| -------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------- |
| **Chromatin regulation & DNA maintenance**   | histone acetyltransferase complex, chromatin remodeler activity, DNA helicase activity, DNA recombination, double-strand break repair, chromosome organization | Oxytropis, Dryas, Draba, Cochlearia | Suggests adaptation involving **genome stability and transcriptional regulation**, likely responding to UV damage and environmental stress. |
| **Redox metabolism & oxidative stress**      | oxidoreductase activity, malate dehydrogenase activity, acyl-CoA dehydrogenase activity, tetrahydrobiopterin metabolism, NAD/NADP binding                      | Oxytropis, Dryas, Cochlearia, Draba | Indicates selection on enzymes involved in **reactive oxygen species management and metabolic balance** under cold stress.                  |
| **Stress and defense signaling**             | regulation of defense response, response to external stimulus, response to biotic stimulus, regulation of response to stress                                   | Cochlearia, Draba                   | Suggests adaptation in **stress signaling networks**, possibly balancing growth and stress tolerance in extreme climates.                   |
| **Membrane and lipid metabolism**            | lipid transfer activity, palmitoyltransferase activity, acylglycerol acyltransferase activity, membrane organization                                           | Oxytropis, Dryas, Draba             | Membrane lipid modification helps maintain **membrane fluidity and stability at low temperatures**.                                         |
| **Cytoskeleton and intracellular transport** | cytoskeleton organization, microtubule binding, microtubule motor activity, COPII vesicle coat                                                                 | Oxytropis, Draba                    | Cytoskeletal systems may be under selection to maintain **cell structure and transport processes during cold stress**.                      |
| **Carbohydrate and energy metabolism**       | glycogen biosynthesis, glucan biosynthesis, polysaccharide metabolism, trehalose metabolism                                                                    | Oxytropis, Dryas                    | Carbohydrate pathways may contribute to **energy storage and cryoprotection** in freezing conditions.                                       |
| **RNA and ribosome processing**              | rRNA processing, rRNA metabolic process, RNA helicase activity                                                                                                 | Oxytropis, Dryas, Draba             | Suggests selection on **protein synthesis machinery**, potentially optimizing growth during short Arctic seasons.                           |

# Shared terms
| GO Term                                     | GO ID      | Species Showing Enrichment |
| ------------------------------------------- | ---------- | -------------------------- |
| rRNA processing                             | GO:0006364 | Oxytropis, Dryas           |
| palmitoyltransferase activity               | GO:0016409 | Oxytropis, Dryas           |
| protein polyubiquitination                  | GO:0000209 | Oxytropis, Draba           |
| regulation of defense response              | GO:0031347 | Cochlearia, Draba          |
| regulation of response to external stimulus | GO:0032101 | Cochlearia, Draba          |
| regulation of response to biotic stimulus   | GO:0002831 | Cochlearia, Draba          |
| water-soluble vitamin biosynthetic process  | GO:0042364 | Dryas, Cochlearia          |
| malate dehydrogenase activity               | GO:0016615 | Dryas, Draba               |
| ATP-dependent chromatin remodeler activity  | GO:0140658 | Cochlearia, Draba          |

stress response regulation
chromatin remodeling
redox metabolism
membrane modification


##############################
# Witout paralogs

awk -F'\t' '$1=="!" {print $2,$3,$6,$7}' OFS='\t' single_Oxydig.ermine.results \
| sort -t $'\t' -k4,4g | head -n 30

# almitoyltransferase activity   GO:0016409      2.000   7.676e-03
# protein transport       GO:0015031      5.000   0.01436896
# intracellular protein transport GO:0006886      4.000   0.01786147
# establishment of protein localization   GO:0045184      5.000   0.0196577
# cellular macromolecule localization     GO:0070727      5.000   0.02167118
# macromolecule localization      GO:0033036      5.000   0.02167118
# protein localization    GO:0008104      5.000   0.02167118
# intracellular protein-containing complex        GO:0140535      4.000   0.02732024
# nitrogen compound transport     GO:0071705      5.000   0.04109043
# DNA-templated transcription initiation  GO:0006352      2.000   0.04856101
# transcription coregulator activity      GO:0003712      2.000   0.05879127


awk -F'\t' '$1=="!" {print $2,$3,$6,$7}' OFS='\t' single_Dryasoct.ermine.results \
| sort -t $'\t' -k4,4g | head -n 30
# ATP-dependent activity, acting on DNA   GO:0008094      9.000   3.008e-03
# DNA recombination       GO:0006310      5.000   3.168e-03
# DNA helicase activity   GO:0003678      4.000   4.505e-03
# helicase activity       GO:0004386      7.000   5.537e-03
# DNA replication GO:0006260      6.000   6.392e-03
# NAD binding     GO:0051287      5.000   6.773e-03
# sequence-specific double-stranded DNA binding   GO:1990837      4.000   7.286e-03
# double-strand break repair      GO:0006302      3.000   9.95e-03
# transferase activity, transferring sulphur-containing groups    GO:0016782      3.000   9.95e-03
# sequence-specific DNA binding   GO:0043565      9.000   0.01029736
# DNA repair      GO:0006281      10.000  0.01151728
# catalytic activity, acting on DNA       GO:0140097      10.000  0.01329212
# rRNA metabolic process  GO:0016072      6.000   0.01357375
# rRNA processing GO:0006364      6.000   0.01357375
# ATP-dependent DNA damage sensor activity        GO:0140664      4.000   0.01577975
# DNA damage sensor activity      GO:0140612      4.000   0.01577975
# DNA damage response     GO:0006974      10.000  0.01632977
# NAD+ binding    GO:0070403      2.000   0.02048215
# small-subunit processome        GO:0032040      2.000   0.02048215
# double-stranded DNA binding     GO:0003690      6.000   0.02070334
# palmitoyltransferase activity   GO:0016409      3.000   0.02133431
# MCM complex     GO:0042555      2.000   0.02976737
# protein heterodimerization activity     GO:0046982      3.000   0.03164541
# oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor   GO:0016616      5.000   0.03256232
# positive regulation of macromolecule metabolic process  GO:0010604      4.000   0.03705679
# positive regulation of metabolic process        GO:0009893      4.000   0.03705679
# transcription cis-regulatory region binding     GO:0000976      3.000   0.03761326
# transcription regulatory region nucleic acid binding    GO:0001067      3.000   0.03761326
# double-strand break repair via homologous recombination GO:0000724      2.000   0.04038304
# double-stranded RNA binding     GO:0003725      2.000   0.04038304

awk -F'\t' '$1=="!" {print $2,$3,$6,$7}' OFS='\t' single_Cochgro.ermine.results \
| sort -t $'\t' -k4,4g | head -n 30
# phosphatidylinositol binding    GO:0035091      2.000   0.02180979
# ATP-dependent chromatin remodeler activity      GO:0140658      2.000   0.02417668
# phospholipid binding    GO:0005543      2.000   0.02920863
# anatomical structure development        GO:0048856      2.000   0.03745958
# lyase activity  GO:0016829      3.000   0.04476684
# 1-phosphatidylinositol binding  GO:0005545      1.000   0.0627273
# arginine N-methyltransferase activity   GO:0016273      1.000   0.0627273
# carbon-sulfur lyase activity    GO:0016846      1.000   0.0627273
# clathrin binding        GO:0030276      1.000   0.0627273
# clathrin coat assembly  GO:0048268      1.000   0.0627273
# clathrin-coated vesicle GO:0030136      1.000   0.0627273
# coated vesicle  GO:0030135      1.000   0.0627273
# cytoplasmic vesicle     GO:0031410      1.000   0.0627273
# intracellular vesicle   GO:0097708      1.000   0.0627273
# phosphatase regulator activity  GO:0019208      1.000   0.0627273
# primary alcohol metabolic process       GO:0034308      1.000   0.0627273
# protein-arginine N-methyltransferase activity   GO:0016274      1.000   0.0627273
# root development        GO:0048364      1.000   0.0627273
# system development      GO:0048731      1.000   0.0627273
# thiamine-containing compound biosynthetic process       GO:0042724      1.000   0.0627273
# thiamine-containing compound metabolic process  GO:0042723      1.000   0.0627273
# thiamine metabolic process      GO:0006772      1.000   0.0627273
# tissue development      GO:0009888      1.000   0.0627273
# vesicle GO:0031982      1.000   0.0627273
# developmental process   GO:0032502      2.000   0.06661252


awk -F'\t' '$1=="!" {print $2,$3,$6,$7}' OFS='\t' single_Drabaniv.ermine.results \
| sort -t $'\t' -k4,4g | head -n 30

# phospholipid binding    GO:0005543      4.000   2.806e-03
# membrane organization   GO:0061024      3.000   4.387e-03
# GTPase activity GO:0003924      5.000   5.654e-03
# guanyl nucleotide binding       GO:0019001      6.000   7.542e-03
# intermembrane lipid transfer    GO:0120009      2.000   7.951e-03
# lipid transfer activity GO:0120013      2.000   7.951e-03
# vesicle coat    GO:0030120      2.000   0.01169998
# lipid transporter activity      GO:0005319      2.000   0.01606919
# lipid transport GO:0006869      3.000   0.01917756
# GTP binding     GO:0005525      5.000   0.02429545
# guanyl ribonucleotide binding   GO:0032561      5.000   0.02429545
# membrane coat   GO:0030117      2.000   0.02651539
# dicarboxylic acid metabolic process     GO:0043648      2.000   0.03251983
# lipid binding   GO:0008289      4.000   0.0408587
# alcohol metabolic process       GO:0006066      2.000   0.04592196
# polygalacturonase activity      GO:0004650      2.000   0.04592196
# endoplasmic reticulum to Golgi vesicle-mediated transport       GO:0006888      2.000   0.05325639
# FAD binding     GO:0071949      2.000   0.05325639
# microtubule-based movement      GO:0007018      3.000   0.06004565
# microtubule motor activity      GO:0003777      3.000   0.06004565
# phosphatidylinositol binding    GO:0035091      2.000   0.06097326
# ligase activity, forming carbon-nitrogen bonds  GO:0016879      2.000   0.06904443
# cytoskeletal motor activity     GO:0003774      3.000   0.06992031



| Theme                              | Oxydig               | Dryasoct                 | Cochgro            | Drabaniv         |
| ---------------------------------- | -------------------- | ------------------------ | ------------------ | ---------------- |
| Protein/macromolecule localization | ✅                    |                          |                    |                  |
| Vesicle-mediated transport         |                      |                          | ✅                  | ✅                |
| Lipid/membrane binding & transport |                      |                          | ✅                  | ✅                |
| DNA replication & repair           |                      | ✅                        |                    |                  |
| Enzymatic activity (modifications) | palmitoyltransferase | DNA-related transferases | methyltransferases | GTPases, ligases |

# Shared GO terms
| GO Term                                    | GO ID      | Species Showing Enrichment |
| ------------------------------------------ | ---------- | -------------------------- |
| palmitoyltransferase activity              | GO:0016409 | Oxytropis, Dryas           |
| rRNA processing                            | GO:0006364 | Dryas                      |
| phospholipid binding                       | GO:0005543 | Cochlearia, Draba          |
| phosphatidylinositol binding               | GO:0035091 | Cochlearia, Draba          |
| ATP-dependent chromatin remodeler activity | GO:0140658 | Cochlearia                 |
| microtubule motor activity                 | GO:0003777 | Draba                      |
| microtubule-based movement                 | GO:0007018 | Draba                      |

# Shared functional categories

| Functional Category                        | Example GO Terms                                                                  | Species                             |
| ------------------------------------------ | --------------------------------------------------------------------------------- | ----------------------------------- |
| DNA repair & genome maintenance            | DNA recombination, DNA helicase activity, DNA repair, double-strand break repair  | Dryas                               |
| Chromatin regulation                       | ATP-dependent chromatin remodeler activity, transcription coregulator activity    | Cochlearia, Oxytropis               |
| Membrane & lipid interactions              | phospholipid binding, phosphatidylinositol binding, palmitoyltransferase activity | Oxytropis, Dryas, Cochlearia, Draba |
| Vesicle trafficking & protein localization | protein transport, vesicle coat, ER→Golgi transport                               | Oxytropis, Draba, Cochlearia        |
| Cytoskeleton & intracellular movement      | microtubule motor activity, microtubule-based movement                            | Draba                               |

#########################
# Present in both datasets


| Functional Category                          | All Genes (with paralogs)                                                          | Single-Copy Genes Only                                               | Interpretation                                                                                                     |
| -------------------------------------------- | ---------------------------------------------------------------------------------- | -------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------ |
| **Chromatin regulation & transcription**     | Histone acetyltransferase complexes, chromatin remodeling, transcription complexes | Chromatin remodeler activity, transcription coregulator activity     | Selection on **gene expression regulation** appears in both analyses → likely a robust signal.                     |
| **DNA repair & genome maintenance**          | DNA helicase activity, recombination, double-strand break repair, DNA replication  | DNA recombination, helicase activity, DNA repair (mainly in *Dryas*) | Suggests adaptation to **DNA damage from UV and environmental stress**.                                            |
| **Redox metabolism**                         | Oxidoreductase activity, dehydrogenases, vitamin metabolism                        | Some oxidoreductase activity but less prominent                      | Redox metabolism appears stronger when paralogs are included, suggesting some **gene family expansion effects**.   |
| **Membrane & lipid biology**                 | Lipid metabolism, palmitoyltransferase activity, membrane organization             | Palmitoyltransferase activity, phospholipid binding, lipid transport | **Very consistent signal across both datasets**, indicating strong selection on **membrane systems**.              |
| **Stress & defense regulation**              | Regulation of defense response, response to stress                                 | Less prominent                                                       | Stress signaling pathways may involve **expanded gene families**, explaining reduced signal after paralog removal. |
| **Cytoskeleton & intracellular movement**    | Cytoskeleton organization, microtubule binding                                     | Microtubule motor activity, microtubule-based movement               | Suggests adaptation affecting **cell structure and intracellular transport**.                                      |
| **Protein localization & vesicle transport** | Protein import, macromolecule localization                                         | Protein transport, vesicle coat, ER→Golgi transport                  | Indicates selection on **intracellular trafficking systems**.                                                      |
| **Carbohydrate & energy metabolism**         | Polysaccharide metabolism, glycogen biosynthesis, trehalose metabolism             | Mostly absent                                                        | Likely influenced by **duplicated metabolic genes**.                                                               |

# Cold, stress and UV impact 
# DNA repair, oxidative stress, and membrane/lipids

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
