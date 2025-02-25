#Last modified January 2024
#Mapping gff3s to fasta
    #Download gff3 toolkit
    #Rheum nobile
########################################################################
#Download gff3 toolkit
#https://gff3toolkit.readthedocs.io/en/latest/gff3_to_fasta.html
module load python/3.10.2
virtualenv ~/gff3_env

source ~/gff3_env/bin/activate
pip install gff3tool
cd gff3tool

########################################################################
#Rheum nobile run
cd ~/scratch
mkdir gff3_to_proteins
gff3_to_fasta -g R_nobile.gff3 -f R_nobile.fasta -st pep -d simple -o R_nobile_proteins
cp R_nobile_proteins_pep.fa ~/projects/def-rieseber/Dryas_shared_data/Oxyria/Polygonaceae_Genomes_Annotations/Rheum_nobile/R_nobile_proteins.fasta


##########################################################################
#Rumex hastatulus scratch (not working)

grep "gene" R_hastalus.gff3 | wc -l
42993

grep "CDS" R_hastalus.gff3 | wc -l
183608

grep "mRNA" R_hastalus.gff3 | wc -l
43786

grep -v "intron" R_hastalus.gff3 > R_hastalus_cleaned.gff3

gff3_QC -g R_hastalus.gff3 -f R_hastalus.fasta -o report.txt -s statistic.txt
gff3_to_fasta -g R_hastalus_cleaned.gff3 -f R_hastalus.fasta -st pep -d simple -o R_hastalus_cleaned
gff3_to_fasta -g R_hastalus.gff3 -f R_hastalus_linkage.fasta -st pep -d simple -o R_hastalus
cp R_hastalus_pep.fa ~/projects/def-rieseber/Dryas_shared_data/Oxyria/Polygonaceae_Genomes_Annotations/Rumex_hastalus/R_hastalus_proteins.fasta


gff3_to_fasta -g R_hastatulus_cleaned.gff3 -f R_hastalus_linkage.fasta -st pep -d simple -o R_hastatulus

#########################################################################