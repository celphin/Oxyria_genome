#######################################################################################################################################
#WRONG do not follow (May 2024)
##########################################################################
#1 Install Phenogram:
wget https://ritchielab.org/files/RL_software/ruby_install.sh
wget https://ritchielab.org/files/RL_software/pheno_gram.rb
sh ruby_install.sh
ruby pheno_gram.rb

#######################################################################################################################################
#2 For Oxyria Rheum nobile intersection :  Get list of siginificant genes and their go terms,  as well as start end positions and their start/end points 
mkdir OxyRheum_Phenogram; cd OxyRheum_Phenogram
awk '$2 == "1" {print $1}' ~/scratch/Oxyria/CAFE/enrichment_analysis/Oxydig_Rheumnob_expanded_genesets | grep "Oxy" > OxyRheum_expanded_genes
awk '$2 == "1" {print $1}' ~/scratch/Oxyria/CAFE/enrichment_analysis/Oxydig_Rheumnob_contracted_genesets | grep "Oxy" > OxyRheum_contracted_genes

grep -f OxyRheum_expanded_genes  ~/scratch/Oxyria/CAFE/enrichment_analysis/Oxydig_GO_mappings.ermineJ.txt  | cut -f1,4 > OxyRheum_expanded_goterms
grep -f OxyRheum_contracted_genes ~/scratch/Oxyria/CAFE/enrichment_analysis/Oxydig_GO_mappings.ermineJ.txt  | cut -f1,4 > OxyRheum_contracted_goterms

grep -f OxyRheum_expanded_genes ../FINAL_Oxyria.AED_0.6.sorted.gff3 | cut -f3,4,9 | grep "gene"   > OxyRheum_expanded_genes_annotation
grep -f OxyRheum_contracted_genes ../FINAL_Oxyria.AED_0.6.sorted.gff3 |cut -f3,4,9 | grep "gene" > OxyRheum_contracted_genes_annotation

cd ..
#######################################################################################################################################
#3 For Oxyria:  Get list of siginificant genes and their go terms,  as well as start end positions and their start/end points 
mkdir Oxydig_Phenograms; cd Oxydig_Phenograms
awk '$2 == "1" {print $1}' ~/scratch/Oxyria/CAFE/enrichment_analysis/Oxydig_expanded_genesets | grep "Oxy" > Oxydig_expanded_genes
awk '$2 == "1" {print $1}' ~/scratch/Oxyria/CAFE/enrichment_analysis/Oxydig_contracted_genesets | grep "Oxy" > Oxydig_contracted_genes

grep -f Oxydig_expanded_genes  ~/scratch/Oxyria/CAFE/enrichment_analysis/Oxydig_GO_mappings.ermineJ.txt  | cut -f1,4 > Oxydig_expanded_goterms
grep -f Oxydig_contracted_genes ~/scratch/Oxyria/CAFE/enrichment_analysis/Oxydig_GO_mappings.ermineJ.txt  | cut -f1,4 > Oxydig_contracted_goterms

#
grep -f Oxydig_expanded_genes ../FINAL_Oxyria.AED_0.6.sorted.gff3 | cut -f3,4,9 | grep "gene"   > Oxydig_expanded_genes_annotation
grep -f Oxydig_contracted_genes ../FINAL_Oxyria.AED_0.6.sorted.gff3 |cut -f3,4,9 | grep "gene" > Oxydig_contracted_genes_annotation
cd .. 
######################################################################################################################################
#4 Make Genome file:
grep ">" Oxyria_Main.fasta | head -n 7 > Oxyria_Genome.txt ##Now add header ID Size, and reformat
######################################################################################################################################
module load  StdEnv/2020 r/4.1.0
export R_LIBS_USER=/home/msandler/R/x86_64-pc-linux-gnu-library/4.1.0/

#Follow along Oxyria_Phenogram_Prep.R in R
######################################################################################################################################
#Remove any Chromosome 8 
awk '$1 != 8' "Oxydig_Phenograms/Oxydig_expanded.txt" > Oxydig_Phenograms/Oxydig_expanded_input.txt
awk '$1 != 8' "Oxydig_Phenograms/Oxydig_contracted.txt" > "Oxydig_Phenograms/Oxydig_contracted_input.txt"
awk '$1 != 8' "OxyRheum_Phenogram/OxyRheum_expanded.txt" > "OxyRheum_Phenogram/OxyRheum_expanded_input.txt"
awk '$1 != 8' "OxyRheum_Phenogram/OxyRheum_contracted.txt" > "OxyRheum_Phenogram/OxyRheum_contracted_input.txt"
######################################################################################################################################
#Run phenogram
ruby pheno_gram.rb -i Oxydig_Phenograms/Oxydig_expanded_input.txt -g Oxyria_Genome.txt -t "Expanded Orthogroup Associated GO-terms for Oxyria digyna" -o Oxydig_Phenograms/oxyria_expanded -f jpg
ruby pheno_gram.rb -i Oxydig_Phenograms/Oxydig_contracted_input.txt -g Oxyria_Genome.txt -t "Contracted Orthogroup Associated GO-terms for Oxyria digyna" -o  Oxydig_Phenograms/oxyria_contracted -f jpg

ruby pheno_gram.rb -i OxyRheum_Phenogram/OxyRheum_expanded_input.txt -g Oxyria_Genome.txt -t "Expanded Orthogroup for Rheum nobile and Oxyria digyna Associated GO-terms for Oxyria digyna" -o OxyRheum_Phenogram/oxyriarheum_expanded -f jpg
ruby pheno_gram.rb -i OxyRheum_Phenogram/OxyRheum_contracted_input.txt -g Oxyria_Genome.txt -t "Contracted Orthogroup Associated for Rheum nobile and Oxyria digyna GO-terms for Oxyria digyna" -o  OxyRheum_Phenogram/oxyriarheum_contracted -f jpg

xdg-open oxyria_expanded.jpg
xdg-open oxyria_contracted.jpg
######################################################################################################################################
#Copy to locals
scp -v msandler@graham.computecanada.ca:/home/msandler/scratch/Oxyria/Phenograms/*/*.jpg .