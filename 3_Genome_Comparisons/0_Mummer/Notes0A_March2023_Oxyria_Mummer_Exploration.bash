# Macro-synteny 

# Oxyria and Rumex comparisons (2n=10)
# Oxyria and Fagopyrum tataricum 2n = 2x = 16

# https://docs.google.com/document/d/1D7qtl-i80d0k_P7l1g9uWYixXlmdvVfdc1Z0LACRpp8/edit#

#-----------------------------------
# make synteny plot
tmux new-session -s Oxyria
tmux attach-session -t Oxyria

cd /home/celphin/scratch/
mkdir Oxyria/; cd Oxyria

cd /home/celphin/scratch/Oxyria
cp /home/celphin/scratch/repeats/input_chromosomes/Oxyria/Oxy_dig_1.scaffolds_FINAL.final.review.fasta .
cp /home/celphin/scratch/repeats/input_chromosomes/Oxyria/Oxyria_RT_2.fasta .
cp /home/celphin/scratch/repeats/input_chromosomes/Fagopyrum/Fagopyrum_2.fasta .
cp /home/celphin/scratch/repeats/input_chromosomes/Polygonum/Polygonum_2.fasta .

mv Oxy_dig_1.scaffolds_FINAL.final.review.fasta Oxyria_genome.fasta 

mkdir Mummer; cd Mummer
mkdir output

#----------
# keep only chromosomes
# remove sequences less than 1Mbp
cd /home/celphin/scratch/Oxyria
module load seqkit/2.3.1
seqkit seq -m 5000000 ./total_genomes/Fagopyrum_2.fasta > Fagopyrum.fasta
seqkit seq -m 5000000 ./total_genomes/Oxyria_genome.fasta > Oxyria.fasta
seqkit seq -m 5000000 ./total_genomes/Oxyria_RT_2.fasta > Oxyria_RT.fasta
seqkit seq -m 5000000 ./total_genomes/Polygonum_2.fasta > Polygonum.fasta

###########################
# rename scaffolds

module load StdEnv/2020 bioawk/1.0
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}'  Fagopyrum.fasta  |\
awk -F '\t' '{printf("%d\t%s\n",length($2),$0);}' |\
sort -k1,1nr | cut -f 2- | tr "\t" "\n" |\
bioawk -c fastx '{ print ">Fago-" ++i "-" length($seq) "\n" $seq}'  > Fagopyrum_names.fasta 

awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}'  Oxyria.fasta  |\
awk -F '\t' '{printf("%d\t%s\n",length($2),$0);}' |\
sort -k1,1nr | cut -f 2- | tr "\t" "\n" |\
bioawk -c fastx '{ print ">Oxy-" ++i "-" length($seq) "\n" $seq}'  > Oxyria_names.fasta 

awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}'  Oxyria_RT.fasta  |\
awk -F '\t' '{printf("%d\t%s\n",length($2),$0);}' |\
sort -k1,1nr | cut -f 2- | tr "\t" "\n" |\
bioawk -c fastx '{ print ">OxyRT-" ++i "-" length($seq) "\n" $seq}'  > Oxyria_RT_names.fasta 

awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}'  Polygonum.fasta  |\
awk -F '\t' '{printf("%d\t%s\n",length($2),$0);}' |\
sort -k1,1nr | cut -f 2- | tr "\t" "\n" |\
bioawk -c fastx '{ print ">Poly-" ++i "-" length($seq) "\n" $seq}'  > Polygonum_names.fasta 

grep ">" *_names.fasta


#####################
# Mummer
#-l Set the minimum length of a single exact match (20) - leave at 20, set 100 for Oxyria
#-c Sets the minimum length of a cluster of matches (65) - leave at 65, set to 2000 Oxyria
#-g Set the maximum gap between two adjacent matches in a cluster (90) - leave at 90

cd /home/celphin/scratch/Oxyria/Mummer

cat << EOF > mummer_Oxyria_Oxyria_RT.sh
#!/bin/bash
#SBATCH --account=def-rieseber
#SBATCH --time=00:50:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128000M

module load StdEnv/2020
module load python/3.10.2
module load minimap2/2.24
module load scipy-stack/2022a
module load gcc/9.3.0
module load r-bundle-bioconductor/3.12 # includes plotly
module load seqtk/1.3
module load StdEnv/2020
module load mummer/4.0.0beta2

cd /home/celphin/scratch/Oxyria/Mummer

nucmer -t 32 -c 2000 -l 100 -p ./output/Oxyria_Oxyria_RT Oxyria_RT_names.fasta Oxyria_names.fasta
delta-filter -l 10000 -q -r ./output/Oxyria_Oxyria_RT.delta > Oxyria_Oxyria_RT_filter.delta
mummerplot Oxyria_Oxyria_RT_filter.delta -R Oxyria_RT_names.fasta -Q Oxyria_names.fasta  --png -p out_Oxyria_RT 

EOF

sbatch mummer_Oxyria_Oxyria_RT.sh

#-------------------

cat << EOF > mummer_Oxyria_Fagopyrum.sh
#!/bin/bash
#SBATCH --account=def-rieseber
#SBATCH --time=00:50:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128000M

module load StdEnv/2020
module load python/3.10.2
module load minimap2/2.24
module load scipy-stack/2022a
module load gcc/9.3.0
module load r-bundle-bioconductor/3.12 # includes plotly
module load seqtk/1.3
module load StdEnv/2020
module load mummer/4.0.0beta2

cd /home/celphin/scratch/Oxyria/Mummer

nucmer -t 32 -p ./output/Oxyria_Fagopyrum Fagopyrum_names.fasta Oxyria_names.fasta
delta-filter -l 1000 -q -r ./output/Oxyria_Fagopyrum.delta > Oxyria_Fagopyrum_filter.delta
mummerplot Oxyria_Fagopyrum_filter.delta -R Fagopyrum_names.fasta -Q Oxyria_names.fasta --png -p out_Fagopyrum

EOF

sbatch mummer_Oxyria_Fagopyrum.sh

#-------------------

cat << EOF > mummer_Oxyria_Polygonum.sh
#!/bin/bash
#SBATCH --account=def-rieseber
#SBATCH --time=00:50:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128000M

module load StdEnv/2020
module load python/3.10.2
module load minimap2/2.24
module load scipy-stack/2022a
module load gcc/9.3.0
module load r-bundle-bioconductor/3.12 # includes plotly
module load seqtk/1.3
module load StdEnv/2020
module load mummer/4.0.0beta2

cd /home/celphin/scratch/Oxyria/Mummer

nucmer -t 32 -p ./output/Oxyria_Polygonum Polygonum_names.fasta Oxyria_names.fasta
delta-filter -l 1000 -q -r ./output/Oxyria_Polygonum.delta > Oxyria_Polygonum_filter.delta
mummerplot Oxyria_Polygonum_filter.delta -R Polygonum_names.fasta -Q Oxyria_names.fasta -p out_Polygonum --png 

EOF

sbatch mummer_Oxyria_Polygonum.sh

#------------------------------
cat << EOF > mummer_Fagopyrum_Polygonum.sh
#!/bin/bash
#SBATCH --account=def-rieseber
#SBATCH --time=00:50:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128000M

module load StdEnv/2020
module load python/3.10.2
module load minimap2/2.24
module load scipy-stack/2022a
module load gcc/9.3.0
module load r-bundle-bioconductor/3.12 # includes plotly
module load seqtk/1.3
module load StdEnv/2020
module load mummer/4.0.0beta2

cd /home/celphin/scratch/Oxyria/Mummer

nucmer -t 32 -p ./output/Fagopyrum_Polygonum Polygonum_names.fasta Fagopyrum_names.fasta
delta-filter -l 1000 -q -r ./output/Fagopyrum_Polygonum.delta > Fagopyrum_Polygonum_filter.delta
mummerplot Fagopyrum_Polygonum_filter.delta -R Polygonum_names.fasta -Q Fagopyrum_names.fasta -p out_Poly_Fago --png 

EOF

sbatch mummer_Fagopyrum_Polygonum.sh

##########################
# plotting

#-----------------
# edit the .gp files to make dots smaller
nano out.gp

set style line 1  lt 1 lw 0.3 pt 6 ps 0.1
set style line 2  lt 3 lw 0.3 pt 6 ps 0.1
set style line 3  lt 2 lw 0.3 pt 6 ps 0.1


############################
# Other plotting??
#--------------------
# https://schneebergerlab.github.io/syri/install.html

#-------------------------
# https://github.com/ctb/2020-stacked-dot-plots

#-------------------------
# plot?
# https://github.com/jmonlong/Hippocamplus/issues/2
# https://github.com/jmonlong/Hippocamplus/blob/master/content/post/2018-06-09-ClusterEqualSize.Rmd
# https://elki-project.github.io/tutorial/same-size_k_means


module load StdEnv/2020 r/4.1.2
R

ggplot(mumgp.filt.diag, aes(x=rs, xend=re, y=qs, yend=qe, colour=strand)) +
  geom_segment(show.legend=FALSE) + geom_point(alpha=0.09) + theme_bw() + 
  facet_grid(qid~rid, scales='free', space='free', switch='both') +
  guides(colour=guide_legend(override.aes=list(alpha=1))) + 
  theme(strip.text.y=element_text(angle=180, size=10),
        strip.text.x=element_text(size=10),
        strip.background=element_blank(),
        legend.position=c(1,-.03), legend.justification=c(1,1),
        legend.direction='horizontal',
        axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        panel.spacing=unit(0, 'cm')) +
  xlab('reference sequence') + ylab('assembly') + scale_colour_brewer(palette='Set1')
  
  
  #########################################
  # to plot a single chromosomes
  
# example

cat << EOF > mummer_Oxyria_Chr7.sh
cat << EOF > mummer_Oxyria_Polygonum.sh
#!/bin/bash
#SBATCH --account=def-rieseber
#SBATCH --time=03:00:00
#SBATCH --time=02:50:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=187G
#SBATCH --cpus-per-task=32
#SBATCH --mem=128000M

module load StdEnv/2020
module load python/3.10.2
@ -322,187 +147,33 @@ module load seqtk/1.3
module load StdEnv/2020
module load mummer/4.0.0beta2

cd /home/celphin/scratch/Oxyria/Mummer

nucmer -t 48 -p ./output/Oxyria_Rumex_Chr7 Rumex_genome_chr.fasta /home/celphin/scratch/repeats/Oxyria/chromosome_files/Oxyria_H0_Chr7.fasta
delta-filter -q -r ./output/Oxyria_Rumex_Chr7.delta > Oxyria_Rumex_Chr7_filter.delta
mummerplot Oxyria_Rumex_Chr7_filter.delta -R Rumex_genome_chr.fasta -Q /home/celphin/scratch/repeats/Oxyria/chromosome_files/Oxyria_H0_Chr7.fasta -p out_Chr7

EOF

sbatch mummer_Oxyria_Chr1.sh
sbatch mummer_Oxyria_Chr2.sh
sbatch mummer_Oxyria_Chr3.sh
sbatch mummer_Oxyria_Chr4.sh
sbatch mummer_Oxyria_Chr5.sh
sbatch mummer_Oxyria_Chr6.sh
sbatch mummer_Oxyria_Chr7.sh