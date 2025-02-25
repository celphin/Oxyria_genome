##############################
# Macro-synteny of Dryas genomes
# Mummer
# June 2024
#####################################

###########################
# Run Mummer to compare/align genomes (?) 
cd ~/scratch/Dryas/Dryas_genomes/
mkdir Mummer; cd Mummer
mkdir output

#----------------------------
# remove short scaffolds
cd ~/scratch/Dryas/Dryas_genomes/Mummer
module load seqkit/2.3.1
seqkit seq -m 1000000 ~/scratch/Dryas/Dryas_genomes/Dry-int_ragtag_output/ragtag.scaffold.fasta > Dry-int-chr.fasta
seqkit seq -m 5000000 ~/scratch/Dryas/Dryas_genomes/Dry-drumm_ragtag_output/ragtag.scaffold.fasta > Dry-drumm-chr.fasta
seqkit seq -m 5000000 ~/scratch/Dryas/Dryas_genomes/Dry-alask_ragtag_output/ragtag.scaffold.fasta > Dry-alask-chr.fasta
seqkit seq -m 5000000 ~/scratch/Dryas/Dryas_genomes/Dry-octo-H1_ragtag_output/ragtag.scaffold.fasta > Dry-octo-H1-chr.fasta
seqkit seq -m 5000000 ~/scratch/Dryas/Dryas_genomes/Dry-octo-H2_ragtag_output/ragtag.scaffold.fasta > Dry-octo-H2-chr.fasta
seqkit seq -m 5000000 ~/scratch/Dryas/Dryas_genomes/Dry-octo-H0.fasta > Dry-octo-H0-chr.fasta


###########################
# rename scaffolds

module load StdEnv/2020 bioawk/1.0
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}'  Dry-int-chr.fasta  |\
awk -F '\t' '{printf("%d\t%s\n",length($2),$0);}' |\
sort -k1,1nr | cut -f 2- | tr "\t" "\n" |\
bioawk -c fastx '{ print ">Dint-" ++i "-" length($seq) "\n" $seq}'  > Dint_names.fasta 

awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}'  Dry-drumm-chr.fasta  |\
awk -F '\t' '{printf("%d\t%s\n",length($2),$0);}' |\
sort -k1,1nr | cut -f 2- | tr "\t" "\n" |\
bioawk -c fastx '{ print ">Ddru-" ++i "-" length($seq) "\n" $seq}'  > Ddru_names.fasta 

awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}'  Dry-alask-chr.fasta  |\
awk -F '\t' '{printf("%d\t%s\n",length($2),$0);}' |\
sort -k1,1nr | cut -f 2- | tr "\t" "\n" |\
bioawk -c fastx '{ print ">Dala-" ++i "-" length($seq) "\n" $seq}'  > Dala_names.fasta 

awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}'  Dry-octo-H0-chr.fasta  |\
awk -F '\t' '{printf("%d\t%s\n",length($2),$0);}' |\
sort -k1,1nr | cut -f 2- | tr "\t" "\n" |\
bioawk -c fastx '{ print ">Doct-" ++i "-" length($seq) "\n" $seq}'  > Doct_names.fasta 

awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}'  Dry-octo-H1-chr.fasta  |\
awk -F '\t' '{printf("%d\t%s\n",length($2),$0);}' |\
sort -k1,1nr | cut -f 2- | tr "\t" "\n" |\
bioawk -c fastx '{ print ">Doct-" ++i "-" length($seq) "\n" $seq}'  > Doct1_names.fasta 

grep ">" *_names.fasta

Dala_names.fasta:>Dala-1-27876026
Dala_names.fasta:>Dala-2-27382827
Dala_names.fasta:>Dala-3-24975665
Dala_names.fasta:>Dala-4-22505024
Dala_names.fasta:>Dala-5-21923467
Dala_names.fasta:>Dala-6-20333848
Dala_names.fasta:>Dala-7-18775577
Dala_names.fasta:>Dala-8-7904532
Dala_names.fasta:>Dala-9-7463978
Ddru_names.fasta:>Ddru-1-32137890
Ddru_names.fasta:>Ddru-2-26519491
Ddru_names.fasta:>Ddru-3-26263223
Ddru_names.fasta:>Ddru-4-25316817
Ddru_names.fasta:>Ddru-5-21834029
Ddru_names.fasta:>Ddru-6-20569620
Ddru_names.fasta:>Ddru-7-20106494
Ddru_names.fasta:>Ddru-8-19617937
Ddru_names.fasta:>Ddru-9-19018337
Dint_names.fasta:>Dint-1-28381354
Dint_names.fasta:>Dint-2-27139721
Dint_names.fasta:>Dint-3-25140901
Dint_names.fasta:>Dint-4-22705339
Dint_names.fasta:>Dint-5-21599655
Dint_names.fasta:>Dint-6-21410073
Dint_names.fasta:>Dint-7-18488805
Dint_names.fasta:>Dint-8-7421606
Dint_names.fasta:>Dint-9-6866929
Doct_names.fasta:>Doct-1-33681925
Doct_names.fasta:>Doct-2-33033793
Doct_names.fasta:>Doct-3-28467151
Doct_names.fasta:>Doct-4-26579288
Doct_names.fasta:>Doct-5-22578703
Doct_names.fasta:>Doct-6-21929824
Doct_names.fasta:>Doct-7-21730711
Doct_names.fasta:>Doct-8-21219003
Doct_names.fasta:>Doct-9-20035207



#####################
# Mummer
#-l Set the minimum length of a single exact match (20) - leave at 20, set 100 for Oxyria
#-c Sets the minimum length of a cluster of matches (65) - leave at 65, set to 2000 Oxyria
#-g Set the maximum gap between two adjacent matches in a cluster (90) - leave at 90

cd ~/scratch/Dryas/Dryas_genomes/Mummer

cat << EOF > mummer_Dint_Doct.sh
#!/bin/bash
#SBATCH --account=def-rieseber
#SBATCH --time=00:50:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=191000M

module load StdEnv/2020
module load python/3.10.2
module load minimap2/2.24
module load scipy-stack/2022a
module load gcc/9.3.0
module load r-bundle-bioconductor/3.12 # includes plotly
module load seqtk/1.3
module load StdEnv/2020
module load mummer/4.0.0beta2

cd ~/scratch/Dryas/Dryas_genomes/Mummer

nucmer -t 40 -c 2000 -l 100 -p ./output/Dint_Doct Dint_names.fasta  Doct_names.fasta
delta-filter -l 10000 -q -r ./output/Dint_Doct.delta > Dint_Doct_filter.delta
mummerplot Dint_Doct_filter.delta -R Dint_names.fasta -Q Doct_names.fasta  --png -p out_Dint_Doct

EOF

sbatch mummer_Dint_Doct.sh

#------------------
salloc -c40 --time 2:55:00 --mem 190000m --account def-rieseber

module load StdEnv/2020
module load python/3.10.2
module load minimap2/2.24
module load scipy-stack/2022a
module load gcc/9.3.0
module load r-bundle-bioconductor/3.12 # includes plotly
module load seqtk/1.3
module load StdEnv/2020
module load mummer/4.0.0beta2

cd ~/scratch/Dryas/Dryas_genomes/Mummer

nucmer -t 40 -c 2000 -l 100 -p ./output/Dint_Doct Dint_names.fasta  Doct_names.fasta
delta-filter -l 10000 -q -r ./output/Dint_Doct.delta > Dint_Doct_filter.delta
mummerplot Dint_Doct_filter.delta -R Dint_names.fasta -Q Doct_names.fasta  --png -p out_Dint_Doct

nucmer -t 40 -c 2000 -l 100 -p ./output/Dala_Doct Dala_names.fasta  Doct_names.fasta
delta-filter -l 10000 -q -r ./output/Dala_Doct.delta > Dala_Doct_filter.delta
mummerplot Dala_Doct_filter.delta -R Dala_names.fasta -Q Doct_names.fasta  --png -p out_Dala_Doct

nucmer -t 40 -c 2000 -l 100 -p ./output/Ddru_Doct Ddru_names.fasta  Doct_names.fasta
delta-filter -l 10000 -q -r ./output/Ddru_Doct.delta > Ddru_Doct_filter.delta
mummerplot Ddru_Doct_filter.delta -R Ddru_names.fasta -Q Doct_names.fasta  --png -p out_Ddru_Doct

nucmer -t 40 -c 2000 -l 100 -p ./output/Doct1_Doct Doct1_names.fasta  Doct_names.fasta
delta-filter -l 10000 -q -r ./output/Doct1_Doct.delta > Doct1_Doct_filter.delta
mummerplot Doct1_Doct_filter.delta -R Doct1_names.fasta -Q Doct_names.fasta  --png -p out_Doct1_Doct

