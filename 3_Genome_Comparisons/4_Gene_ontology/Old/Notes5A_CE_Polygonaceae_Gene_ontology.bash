#############################
# Oxyria and other Polygonaceae Gene Ontology Enrichment
# Feb 2024
##############################
#Note: see other notes labeled as Notes5 for Go ontology analysis
#######################################

# Tutorial https://github.com/elsemikk/Willisornis_Genome_Assembly/blob/master/4.2_Gene_Family_Expansions.md
# GO enrichment and plotting

# What genes are in these orthogroups?

cd /lustre04/scratch/celphin/Oxyria/CAFE/expanded

# Oxydig_expandedfams.sig

#get a fasta of the expanded families from the OrthoFinder output
for taxon in Oxydig ; do cat ../analysis/"$taxon"_sig_changes | cut -f 1 | \
while read Orthogroup ; do \
cp /lustre04/scratch/celphin/Oxyria/CAFE/orthofinder/Results_Mar04/Orthogroup_Sequences/"$Orthogroup".fa . 
done ; done

#-------------------------------------
#Get a list of the Oxydig genes in these orthogroups
for taxon in Oxydig ; do cat ../analysis/"$taxon"_sig_changes | \
cut -f 1 | while read Orthogroup ; do grep "Oxyria" "$Orthogroup".fa >> "$taxon"_expanded_geneIDs.txt ; done ; done

for taxon in Oxydig ; do sed -i 's/>//g' "$taxon"_expanded_geneIDs.txt ; done
#for taxon in Oxydig ; do sed -i 's/-RA//g' "$taxon"_expanded_geneIDs.txt ; done

#-----------------------------------
#human readable
for taxon in Oxydig ; do cat ../analysis/"$taxon"_sig_changes | cut -f 1 | \
while read Orthogroup ; do echo "$Orthogroup" >> human_"$taxon"_expanded_geneIDs.txt 
grep "Oxyria" "$Orthogroup".fa >> human_"$taxon"_expanded_geneIDs.txt ; done ; done

for taxon in Oxydig ; do sed -i 's/>//g' human_"$taxon"_expanded_geneIDs.txt ; done
#for taxon in Oxydig ; do sed -i 's/-RA//g' human_"$taxon"_expanded_geneIDs.txt ; done

######################################
# Enrichment analysis
# https://erminej.msl.ubc.ca/help/tutorials/running-an-analysis-ora/

#-----------------------------
# copy over the annotation file
# Annotation file This is a tab-delimited file .txt file. The first line is a header. 
# The first column is gene names. The second column is a duplicate of the first column in our case. 
# The third column is for a description. I will make it another duplicate for now. 
# The fourth column is a list of GO terms (eg GO:12345,GO2345,GO5643).

mkdir enrichment_analysis
cd /lustre04/scratch/celphin/Oxyria/CAFE/enrichment_analysis

# copy over the annotation file for Oxyria
cp /lustre04/scratch/celphin/Oxyria/Polygonaceae_Genomes_Annotations/MS_CE_Oxyria_digyna/Oxyria_AED0.6_interproscan.tsv /lustre04/scratch/celphin/Oxyria/CAFE/enrichment_analysis

#extract GO terms
# first column is gene name 
# 14th column is GOterms

awk -v FS="\t" '{print $1 "\t" $1 "\t" $13 "\t" $14}'  Oxyria_AED0.6_interproscan.tsv  > GO_mappings.ermineJ.txt

cat > temp_header.txt
gene	gene	description	GO
#ctrl-d
 
cat temp_header.txt GO_mappings.ermineJ.txt > temp && mv temp GO_mappings.ermineJ.txt
rm temp_header.txt
sed -i 's/ /-/g' GO_mappings.ermineJ.txt

#-------------------------------
# copy over the Gene Score File
# The gene score file will tell ErmineJ which genes are in our chosen set and which are not. 
# There must be a 1 line header (the first line will be skipped). There will be two columns: 
# the first has our gene names, and the second has the gene score. 
# In this case it will be a binary score of "NA" (not on our list) or "1" (on our list).
# As a refinement, we could put the p-value instead of a 1 vs 0.

cd /lustre04/scratch/celphin/Oxyria/CAFE/enrichment_analysis

#this is a bad  way
for taxon in Oxydig ; do awk 'BEGIN{FS="\t"}{print $2,"0"}' ~/Willisornis/functional_ann/Gene_mappings.txt > "$taxon"_expanded_genesets ; done
for taxon in Oxydig ; do cat /home/0_GENOMES1/0_WEIRLAB_GENOMES_CHROMIUMX/Willisornis/Cafe/expanded/"$taxon"_expanded_geneIDs.txt | \
sed 's/ .*$//g' | while read gene ; do sed -i "s/$gene 0/$gene 1/g" "$taxon"_expanded_genesets ; done ; done

#make tab separated instead of space separated
for taxon in Oxydig ; do sed -i 's/ /\t/g' "$taxon"_expanded_genesets ; done

#repeat for contracted families
#awk 'BEGIN{FS="\t"}{print $2,"0"}' ~/Willisornis/functional_ann/Gene_mappings.txt > contracted_genesets
#cat /home/0_GENOMES1/0_WEIRLAB_GENOMES_CHROMIUMX/Willisornis/Cafe/contracted/contracted_geneIDs.txt | while read gene ; do sed -i "s/$gene 0/$gene 1/g" contracted_genesets ; done
#sed -i 's/ /\t/g' contracted_genesets

####################################
# install ermine

cd /lustre04/scratch/celphin/Oxyria/CAFE/enrichment_analysis
wget https://home.pavlab.msl.ubc.ca/ermineJ/distributions/ermineJ-3.2-src.tar.gz

tar -xvf ermineJ-3.2-src.tar.gz

#---------------------------------------
# Run ermine!
# Ignore the data profile files, we do not have this type of data.

cd /lustre04/scratch/celphin/Oxyria/CAFE/enrichment_analysis

ERMINEJ_HOME=/home/0_PROGRAMS/ermineJ-3.1.2/
for taxon in Oxydig ; do /lustre04/scratch/celphin/Oxyria/CAFE/enrichment_analysis/ermineJ-3.2/src/bin/ermineJ.sh \
-a GO_mappings.ermineJ.txt \
-s "$taxon"_expanded_genesets \
-c /home/0_PROGRAMS/ermineJ-3.1.2/go_daily-termdb.rdf-xml.gz \
-o "$taxon"_expanded_genesets.ermine.results -y 5 -b ; done

#-------------------------
#contracted families
#/home/0_PROGRAMS/ermineJ-3.1.2/bin/ermineJ.sh -a GO_mappings.ermineJ.txt -s contracted_genesets -c /home/0_PROGRAMS/ermineJ-3.1.2/go_daily-termdb.rdf-xml.gz -o contracted_genesets.ermine.results -y 5 -b

#    -y: min size of gene sets to consider
#    -b: scores go from small to big

#--------------------------
# Now a ranked list regardless of significance

cd /home/0_GENOMES1/0_WEIRLAB_GENOMES_CHROMIUMX/Willisornis/Cafe/analysis
cut -f 1,41 ../error_model_1/Base_change.tab > orthogroup_expansiveness
cut -f 1,2 /home/0_GENOMES1/0_Genetic_resources/orthologues/PasserinesPlus/OrthoFinder/Results_Jan20/Orthogroups/Orthogroups.tsv > orthogroup2gene.dict

#----------------------------
# Need to replace Orthogroup IDs with gene IDs:

cat > orthogroup2gene_scores.py

import os
file1="orthogroup_expansiveness" #This file contains the list of orthologs and their scores
file2='orthogroup2gene.dict' #This file contains the list of gene names that belong to each ortholog name
mainpath=os.getcwd()
fileI=[x.strip().split('\t') for x in open(os.path.join(mainpath,file1),'r').readlines()]
fileII=[x.strip().split('\t') for x in open(os.path.join(mainpath,file2),'r').readlines()]
count=0
#for i in range(len(fileI)):
#    if count<100:
#        if i!=0:#The first row in the table are the headers, which we don't need
#            print(fileI[i])
#    else:
#        break
#    count+=1
first=True
finalstringfileII=[]
for h in fileII:
    if not first:
        try:
            for name in h[1].split(', '):
                finalstringfileII.append([h[0],name]) #Put the names of each gene into a new row with its ortholog
        except:
            pass
    else:
        pass
    first=False

finalstring=''
for r in finalstringfileII:
    for e in fileI:
        if r[0] == e[0]:
            finalstring+=r[1]+'\t'+e[1]+'\n'

export=open(os.path.join(mainpath,"output.txt"),"w")
export.write(finalstring)
export.close()

EOF

#-----------------------------------
#run python script
python orthogroup2gene_scores.py
mv output.txt genescores_expansiveness.txt
sed -i 's/-RA//g' genescores_expansiveness.txt
sed -i 's/+//g' genescores_expansiveness.txt


#----------------------------------
cd /home/0_GENOMES1/0_WEIRLAB_GENOMES_CHROMIUMX/Willisornis/Cafe/enrichment_analysis
ERMINEJ_HOME=/home/0_PROGRAMS/ermineJ-3.1.2/

/home/0_PROGRAMS/ermineJ-3.1.2/bin/ermineJ.sh -a GO_mappings.ermineJ.txt \
-s ../analysis/genescores_expansiveness.txt -c /home/0_PROGRAMS/ermineJ-3.1.2/go_daily-termdb.rdf-xml.gz \
-o expansiveness_genesets.ermine.results -y 5 -b
#Gene set resampling
/home/0_PROGRAMS/ermineJ-3.1.2/bin/ermineJ.sh -a GO_mappings.ermineJ.txt \
-s ../analysis/genescores_expansiveness.txt -c /home/0_PROGRAMS/ermineJ-3.1.2/go_daily-termdb.rdf-xml.gz \
-o expansiveness_genesets.ermine.GSR.results -y 5 -b -n GSR


######################################
# Visualization
# I am going to use R to visualize the results on a phylogeny. 
# I would like to show the number of genes within orthogroups at each tip of the phylogeny.


















############################################
# OLDER - Jan 2024
################################################
# Gene ontology enrichment
# https://rpubs.com/jrgonzalezISGlobal/enrichment
# Note should also look here for RNAseq sorting in Dryas

#---------------------------
#  change Interproscan output to GO annotation file for R
# from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7952955/

cd ~/projects/def-rieseber/Dryas_shared_data/Oxyria/
mkdir GO_enrichment; cd GO_enrichment

wget https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7952955/bin/BioProtoc-11-03-3912-s001.zip
unzip BioProtoc-11-03-3912-s001.zip

#-------------------
# Interproscan data
cp ~/projects/def-rieseber/Dryas_shared_data/Oxyria/Polygonaceae_Genomes_Annotations/MS_CE_Fagopyrum_tataricum/Fagopyrum_AED0.6_interproscan.tsv .
# 274408, 79%
cp ~/projects/def-rieseber/Dryas_shared_data/Oxyria/Polygonaceae_Genomes_Annotations/MS_CE_Oxyria_digyna/Oxyria_AED0.6_interproscan.tsv .
# 351100, 83%
cp ~/projects/def-rieseber/Dryas_shared_data/Oxyria/Polygonaceae_Genomes_Annotations/Rheum_tangaticum/Rheum_tangaticum_interproscan_out.tsv .
# 268361, 93%
cp ~/projects/def-rieseber/Dryas_shared_data/Oxyria/Polygonaceae_Genomes_Annotations/Polygonum_aviculare/Polavi_AED0.6_interproscan.tsv .
# 265885, 83%

#----------------------
# Column number of GO terms
head -n 1 Oxyria_AED0.6_interproscan.tsv | awk '{print NF}'
15

#--------------------
# Substitute INPUT_FILE.tsv by the name of your input file.
# [Optional] Line 2: Change the name of the output file (OUTPUT_FILE_GOs.txt).
# Line 12: Based on your annotation file, set the column (N) where the GO ID numbers were reported (15).

nano GO_retriever.py

input_fl = "Fagopyrum_tat_AED0.6_interproscan.tsv"
output = open("Fagopyrum_tat_GO.txt", "w")

with open(input_fl) as fl:
for line in fl:
line = line.strip()
cols = line.split("\t")
ID = cols[0]

# search for go_term on the Nth column
try:
go_term = cols[13]
except:
continue

# loop through each go_term and write one per line
go_terms = go_term.split("|")
for go in go_terms:
go_num = go.lstrip("GO:")
output.write(ID + " = " + go_num + "\n")

# A pipeline for non-model organisms for de novo transcriptome assembly, annotation, and gene ontology analysis using open tools: case study with Scots pine

# Gustavo T. Duarte, Polina Yu. Volkova, Stanislav A. Geras�kin


#-------------------------------
module load StdEnv/2020 python/3.11.5

python GO_retriever_Fagopyrum_tat.py Fagopyrum_AED0.6_interproscan.tsv

python GO_retriever_Oxyria.py Oxyria_AED0.6_interproscan.tsv

python GO_retriever_Rheum_tan.py Rheum_tangaticum_interproscan_out.tsv

python GO_retriever_Polavi.py Polavi_AED0.6_interproscan.tsv

#-------------------------
# change "=" to tab and "-" to NA

sed -i 's/=/\t/g' Oxyria_GOs.txt
sed -i 's/-/NA/g' Oxyria_GOs.txt

#------------------------------
# R packages for GO enrichment

tmux new-session -s GO
tmux attach-session -t GO

module load StdEnv/2020
module load r/4.2.2

R

#------------------------
# Install libraries

# clusterProfiler
# https://rpubs.com/jrgonzalezISGlobal/enrichment

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("clusterProfiler")

#----------------------
# GOCompare
# https://github.com/ccsosa/GOCompare

#CRAN
install.packages("GOCompare")
#Alternative: GitHub
library(devtools)
remotes::install_github("ccsosa/GOCompare")

################################
# load data

Fagopyrum_tat_GO <- as.data.frame(utils::read.table("Fagopyrum_tat_GO.txt", check.names = FALSE))

Oxyria_GO <- as.data.frame(utils::read.table("Oxyria_GOs.txt", check.names = FALSE))

Polavi_GO <- as.data.frame(utils::read.table("Polavi_GOs.txt", check.names = FALSE))

Rheum_tan_GO <- as.data.frame(utils::read.table("Rheum_tan_GOs.txt", check.names = FALSE))


#-------------------------------
# From 
# https://onlinelibrary.wiley.com/doi/full/10.1111/1755-0998.13280
# As a first analysis of gene content, we annotated Pfam domains (El-Gebali et al., 2018) for the 
# predicted genes of each assembly using InterProScan (Jones et al., 2014). Pfam domains were quantified
# for each species, and domains with a Z-score above 1.96 or below �1.96 in D. nivalis were considered 
# significantly enriched or contracted, respectively.

########################
# Maybe try: https://scienceparkstudygroup.github.io/rna-seq-lesson/07-functional-enrichment/index.html#12-over-representation-analysis-ora






#####################################
#------------------------------
# https://bioconductor.org/packages/release/bioc/manuals/clusterProfiler/man/clusterProfiler.pdf
# https://scienceparkstudygroup.github.io/rna-seq-lesson/07-functional-enrichment/index.html#12-over-representation-analysis-ora 
# https://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-go.html

library(clusterProfiler)
# GSEA()
# enricher()

Oxyria_genes <- unique(as.factor(Oxyria_GO[,2]))[c(1:300)]

# https://rdrr.io/bioc/clusterProfiler/man/enricher.html
Oxydig_enrich <- enricher(gene = Oxyria_genes, pAdjustMethod ="BH", qvalueCutoff = 0.05, TERM2GENE=Oxyria_GO)


#--------------
# this is only needed if some genes can be enriched
Oxyria_genes_sort <- sort(Oxyria_genes, decreasing = TRUE)

# https://rdrr.io/bioc/clusterProfiler/man/GSEA.html
Oxydig_GSEA <- GSEA(gene = Oxyria_genes_sort, pAdjustMethod ="BH", qvalueCutoff = 0.2, TERM2GENE=Oxyria_GO)

preparing geneSet collections...
--> Expected input gene ID: 160,16491,16851,8418,16811,15995
Error in check_gene_id(geneList, geneSets) :
  --> No gene can be mapped....


#------------------------------------
write_delim(x = as.data.frame(Oxydig_enrich), 
            path = "go_results_Oxyria.tsv", 
            delim = "\t")

ora_analysis_bp_simplified <- clusterProfiler::simplify(ora_analysis_bp) 
ora_analysis_bp_simplified@result[1:5,1:8]
dotplot(ora_analysis_bp_simplified)

#-----------------------------
tab.go <- as.data.frame(ans.go)
tab.go<- subset(tab.go, Count>5)
tab.go[1:5, 1:6]

library(enrichplot)
p1 <- barplot(ans.dis, showCategory=10)
p1


#--------------------------------
# Try topGO
# https://bioconductor.org/packages/devel/bioc/vignettes/topGO/inst/doc/topGO.pdf


#---------------------
# Try GOCompare

library(GOCompare)

compareGOspecies()