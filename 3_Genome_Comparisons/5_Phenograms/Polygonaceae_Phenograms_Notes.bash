##########################################################
#January 2024
#Phenograms for 4 polygonaceae genomes
##########################################################

wget http://current.geneontology.org/ontology/subsets/goslim_plant.obo
wget https://current.geneontology.org/ontology/external2go/interpro2go

wget https://ritchielab.org/files/RL_software/ruby_install.sh
wget https://ritchielab.org/files/RL_software/pheno_gram.rb

sh ruby_install.sh
ruby pheno_gram.rb
##########################################################
#Create genome files: Chr #, chromsome length:
#grep ">" Oxyria_Main.fasta | head -n 8 > Oxyria_Genome.txt #Now add header ID Size, and reformat
#grep ">" Polygonum_Main.fasta | head -n 10 > Polygonum_Genome.txt #Now add header ID Size, and reformat
#grep ">" Rheum_Main.fasta | head -n 8 > Rheum_Genome.txt #Now add header ID Size, and reformat


grep ">" Fagopyrum_Main.fasta | head -n 8 > Fagopyrum_Genome.txt #Now add header ID Size, and reformat

##########################################################
#cut -f1,14 interproscan_oxyria_v5_63_95.tsv| grep "GO" | grep -v "Chr08">  Oxy_goterm_file.tsv 
cut -f1,14 Fagopyrum_AED0.6_interproscan.tsv | grep "GO" | grep -v "Chr9">  Fagopyrum_goterm_file.tsv
#cut -f1,14 interproscan_oxyria_v5_63_95.tsv| grep "GO" | grep -v "Chr08">  Oxy_goterm_file.tsv
#cut -f1,14 interproscan_oxyria_v5_63_95.tsv| grep "GO" | grep -v "Chr08">  Oxy_goterm_file.tsv

##########################################################
module load StdEnv/2020
module load r/4.2.1
export R_LIBS_USER=/home/msandler/R/x86_64-pc-linux-gnu-library/4.2.1/
Rscript Polygonaceae_Phenogram_Prep.R  FINAL_Fagopyrum.AED_0.6.sorted.gff3 Fagopyrum_goterm_file.tsv Fagopyrum



#############################################################


ruby pheno_gram.rb -i Oxyria_Defenses.txt -g Oxyria_Genome.txt -t "Oxyria Defenses Go Plot" -o oxyria_defenses -f jpg
ruby pheno_gram.rb -i Oxyria_Hormone.txt -g Oxyria_Genome.txt -t "Oxyria Responses to Hormones Go Plot" -o oxyria_hormone -f jpg
ruby pheno_gram.rb -i Oxyria_Oxidative_Stress.txt -g Oxyria_Genome.txt -t "Oxyria Responses to Oxidative Stress Go Plot" -o oxyria_oxidative_stress -f jpg
ruby pheno_gram.rb -i Oxyria_Light.txt -g Oxyria_Genome.txt -t "Oxyria Light Responses Go Plot" -o oxyria_light -f jpg
ruby pheno_gram.rb -i Oxyria_hum_heat.txt -g Oxyria_Genome.txt -t "Oxyria Water and Heat Related Responses Go Plot" -o oxyria_hum_heat -f jpg
ruby pheno_gram.rb -i Oxyria_Ion.txt -g Oxyria_Genome.txt -t "Oxyria Response to Ions Go Plot" -o oxyria_ion -f jpg
ruby pheno_gram.rb -i Oxyria_DNA_Damage.txt -g Oxyria_Genome.txt -t "Oxyria DNA Damage Response Go Plot" -o oxyria_dna_dam -f jpg
ruby pheno_gram.rb -i Oxyria_ER.txt -g Oxyria_Genome.txt -t "Oxyria Endplasmic Reticulum Stress Go Plot" -o oxyria_er -f jpg

xdg-open example.png
ruby pheno_gram.rb -i Fagopyrum_Defenses.txt -g Fagopyrum_Genome.txt -t "Fagopyrum Defenses Go Plot" -o Fagopyrum_defenses -f jpg
ruby pheno_gram.rb -i Fagopyrum_Hormone.txt -g Fagopyrum_Genome.txt -t "Fagopyrum Responses to Hormones Go Plot" -o Fagopyrum_hormone -f jpg
ruby pheno_gram.rb -i Fagopyrum_Oxidative_Stress.txt -g Fagopyrum_Genome.txt -t "Fagopyrum Responses to Oxidative Stress Go Plot" -o Fagopyrum_oxidative_stress -f jpg
ruby pheno_gram.rb -i Fagopyrum_Light.txt -g Fagopyrum_Genome.txt -t "Fagopyrum Light Responses Go Plot" -o Fagopyrum_light -f jpg
ruby pheno_gram.rb -i Fagopyrum_hum_heat.txt -g Fagopyrum_Genome.txt -t "Fagopyrum Water and Heat Related Responses Go Plot" -o Fagopyrum_hum_heat -f jpg
ruby pheno_gram.rb -i Fagopyrum_Ion.txt -g Fagopyrum_Genome.txt -t "Fagopyrum Response to Ions Go Plot" -o Fagopyrum_ion -f jpg
ruby pheno_gram.rb -i Fagopyrum_DNA_Damage.txt -g Fagopyrum_Genome.txt -t "Fagopyrum DNA Damage Response Go Plot" -o Fagopyrum_dna_dam -f jpg
ruby pheno_gram.rb -i Fagopyrum_ER.txt -g Fagopyrum_Genome.txt -t "Fagopyrum Endplasmic Reticulum Stress Go Plot" -o Fagopyrum_er -f jpg