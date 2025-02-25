############################
# Look for types of duplications from CAFE
# July  2024
# https://github.com/qiao-xin/DupGen_finder
##############################

# Install

cd /home/celphin/scratch/Oxyria  # or any directory of your choice
git clone https://github.com/qiao-xin/DupGen_finder.git
cd DupGen_finder
make
chmod 775 DupGen_finder.pl
chmod 775 DupGen_finder-unique.pl
chmod 775 set_PATH.sh
source set_PATH.sh

# test
cd /home/celphin/scratch/Oxyria/DupGen_finder
DupGen_finder.pl

# Usage: DupGen_finder.pl -i data_directory -t target_species -c outgroup_species -o output_directory
# #####################
# Optional:
# -a 1 or 0(are segmental duplicates ancestral loci or not? default: 1, yes)
# -d number_of_genes(maximum distance to call proximal, default: 10)
# #####################
# The following are optional MCScanX parameters:
# -k match_score(cutoff score of collinear blocks for MCScanX, default: 50)
# -g gap_penalty(gap penalty for MCScanX, default: -1)
# -s match_size(number of genes required to call a collinear block for MCScanX, default: 5)
# -e e_value(alignment significance for MCScanX, default: 1e-05)
# -m max_gaps(maximum gaps allowed for MCScanX, default: 25)
# -w overlap_window(maximum distance in terms of gene number, to collapse BLAST matches for MCScanX, default: 5)

###############################
# Preparing input files

    # For the target genome in which gene duplicaiton modes will be classified, please prepare two input files:
        # target_species.gff, a gene position file for the target species, following a tab-delimited format. For example, "Ath.gff".
        # target_species.blast, a blastp output file (-outfmt 6) for the target species (self-genome comparison). For example, "Ath.blast".

    # For the outgroup genome, please prepare two input files:
        # [target_species]_[outgroup_species].gff, a gene position file for the target_species and outgroup_species, following a tab-delimited format.
        # [target_species]_[outgroup_species].blast, a blastp output file (-outfmt 6) between the target and outgroup species (cross-genome comparison).

# blast files
cd /home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/orthofinder/Results_Aug17/WorkingDirectory

more SpeciesIDs.txt
# 0: Arabidopsis_lyrata.fa
# 1: Arabidopsis_thaliana.fa
*# 2: Arabis_alpina.fa
# 3: Argentina_anserina.fa
# 4: Brassica_oleracea.fa
# 5: Capsella_rubella.fa
*# 6: Draba_nivalis.fa
*# 7: Dryas_octopetala.fa
# 8: Fagopyrum_escelentum_H2.fa
# 9: Fagopyrum_tataricum_H1.fa
# 10: Fragaria_vesca.fa
# 11: Malus_sylvestris.fa
*# 12: Oxyria_digyna_H1.fa
# 13: Polygunum_aviculare_H0.fa
# 14: Prunus_persica.fa
# 15: Pyrus_bretschneideri.fa
# 16: Rheum_nobile_H0.fa
# 17: Rheum_tangaticum_H0.fa
# 18: Rosa_rugosa.fa
# 19: Thlaspi_arvense.fa

more SpeciesTree_rooted_ids.txt
# http://etetoolkit.org/treeview/
(((13:0.156996,(12:0.136637,(17:0.0440109,16:0.0354046)0.885208:0.0630934)0.591091:0.0341106)0.57481
4:0.0347272,(9:0.0360394,8:0.037885)0.974872:0.103084)0.966019:0.093549,((((5:0.0516886,(0:0.0229793
,1:0.077621)0.801828:0.02683)0.797259:0.0314861,(6:0.0862001,2:0.0831463)0.731582:0.0326396)0.220731
:0.0158126,(4:0.0837281,19:0.104684)0.346374:0.0201726)0.968304:0.215788,((18:0.0480735,(3:0.0780935
,10:0.0594828)0.50771:0.0187218)0.954026:0.0798625,(7:0.143586,(14:0.0789428,(11:0.0172918,15:0.0191
209)0.981725:0.0817642)0.717304:0.0315978)0.626785:0.029253)0.962593:0.114578)0.966019:0.093549);

# Brass
# 2/6 and 5
cp Blast2_2.txt.gz /home/celphin/scratch/Oxyria/DupGen_finder/data
cp Blast6_6.txt.gz /home/celphin/scratch/Oxyria/DupGen_finder/data
cp Blast2_5.txt.gz /home/celphin/scratch/Oxyria/DupGen_finder/data
cp Blast6_5.txt.gz /home/celphin/scratch/Oxyria/DupGen_finder/data

# Rose
# 7 and 18
cp Blast7_7.txt.gz /home/celphin/scratch/Oxyria/DupGen_finder/data
cp Blast7_18.txt.gz /home/celphin/scratch/Oxyria/DupGen_finder/data

# Poly
# 12 and 13
cp Blast12_12.txt.gz /home/celphin/scratch/Oxyria/DupGen_finder/data
cp Blast12_13.txt.gz /home/celphin/scratch/Oxyria/DupGen_finder/data

# gene IDs
cp SequenceIDs.txt /home/celphin/scratch/Oxyria/DupGen_finder/data

#--------------------
cd /home/celphin/scratch/Oxyria/DupGen_finder/data
gunzip *.gz
rename .txt .blast *

#-------------
# gff3 files
cd /home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/bed

cp Dryas_octopetala.bed /home/celphin/scratch/Oxyria/DupGen_finder/data
cp Draba_nivalis.bed /home/celphin/scratch/Oxyria/DupGen_finder/data
cp Arabis_alpina.bed /home/celphin/scratch/Oxyria/DupGen_finder/data
cp Oxyria_digyna_H1.bed /home/celphin/scratch/Oxyria/DupGen_finder/data

cp Rosa_rugosa.bed /home/celphin/scratch/Oxyria/DupGen_finder/data
cp Capsella_rubella.bed /home/celphin/scratch/Oxyria/DupGen_finder/data
cp Polygunum_aviculare_H0.bed /home/celphin/scratch/Oxyria/DupGen_finder/data

#-------------------
cd /home/celphin/scratch/Oxyria/DupGen_finder/data
more SequenceIDs.txt

mv Dryas_octopetala.bed 7.gff
mv Draba_nivalis.bed 6.gff
mv Arabis_alpina.bed 2.gff
mv Oxyria_digyna_H1.bed 12.gff
mv Rosa_rugosa.bed 18.gff
mv Capsella_rubella.bed 5.gff
mv Polygunum_aviculare_H0.bed 13.gff

#---------------------
# fix names of geneids
grep Polavi_Chr100000050 *

13.gff:Polavi-1 283812  285985  Polavi_Chr100000050
13.gff:Polavi-10        10467046        10467352        Polavi_Chr1000000500
13.gff:Polavi-10        10483672        10488354        Polavi_Chr1000000501
13.gff:Polavi-10        10497423        10497729        Polavi_Chr1000000503
13.gff:Polavi-10        10565760        10566259        Polavi_Chr1000000505
13.gff:Polavi-10        10568349        10568682        Polavi_Chr1000000506
13.gff:Polavi-10        10569363        10573843        Polavi_Chr1000000507
13.gff:Polavi-10        10582034        10583504        Polavi_Chr1000000508
13.gff:Polavi-10        10583762        10584690        Polavi_Chr1000000509
SequenceIDs.txt:13_46: Polavi_Chr100000050
SequenceIDs.txt:13_3895: Polavi_Chr1000000500
SequenceIDs.txt:13_3896: Polavi_Chr1000000501
SequenceIDs.txt:13_3897: Polavi_Chr1000000503
SequenceIDs.txt:13_3898: Polavi_Chr1000000505
SequenceIDs.txt:13_3899: Polavi_Chr1000000506
SequenceIDs.txt:13_3900: Polavi_Chr1000000507
SequenceIDs.txt:13_3901: Polavi_Chr1000000508
SequenceIDs.txt:13_3902: Polavi_Chr1000000509

grep 13_3895 *
BLAST12_13.blast:12_23650       13_3895 41.2    97      56      1       1       97      1       96 3.81e-18 72.4
BLAST12_13.blast:12_23651       13_3895 40.6    96      57      0       1       96      1       96 3.44e-21 80.1
SequenceIDs.txt:13_3895: Polavi_Chr1000000500

#---------------------------------
# left_join in R

# join # ids to gff files
# reorder columns
# Species_abbrev-Chr_ID	gene_ID	start_position	end_position

tmux new-session -s DupGen
tmux attach-session -t DupGen

module load StdEnv/2023
module load r/4.4.0
 
R 

# https://dplyr.tidyverse.org/reference/mutate-joins.html
install.packages("dplyr")

# read in sequnece ids
inpath="/home/celphin/scratch/Oxyria/DupGen_finder/data/"
SequenceIDs<-as.matrix(utils::read.table(paste0(inpath, "SequenceIDs.txt"), header = FALSE, sep = ":", check.names = FALSE))
head(SequenceIDs)
colnames(SequenceIDs) <-  c("genenum", "geneID")
SequenceIDs <- trimws(SequenceIDs)
head(SequenceIDs)

  genenum   geneID
1     0_0 LOC9330760
2     0_1 LOC9328081
3     0_2 LOC9325393
4     0_3 LOC9325394
5     0_4 LOC9328083
6     0_5 LOC9328084

#------------------
# read in gff files
 # 13.gff
 # 18.gff
 # 2.gff
 # 5.gff
 # 6.gff
 # 7.gff

gff13 <- as.matrix(utils::read.table(paste0(inpath, "12.gff"), header = FALSE, check.names = FALSE))
colnames(gff13) <- c("chrom", "start", "end", "geneID")
head(gff13)
gff13 <- trimws(gff13)
head(gff13)

     # chrom start   end              geneID
# 1 Polavi-1 14594 18403 Polavi_Chr100000001
# 2 Polavi-1 19174 20875 Polavi_Chr100000002
# 3 Polavi-1 21594 39613 Polavi_Chr100000003
# 4 Polavi-1 52961 58022 Polavi_Chr100000004
# 5 Polavi-1 60547 69571 Polavi_Chr100000005
# 6 Polavi-1 70144 72735 Polavi_Chr100000006


# types issue
gff13 <-  as.data.frame(gff13)
SequenceIDs <-  as.data.frame(SequenceIDs)

typeof(gff13)
typeof(SequenceIDs)

# left_join
gff13IDs <- dplyr::left_join(gff13, SequenceIDs)

     # chrom start   end              geneID genenum
# 1 Polavi-1 14594 18403 Polavi_Chr100000001    13_0
# 2 Polavi-1 19174 20875 Polavi_Chr100000002    13_1
# 3 Polavi-1 21594 39613 Polavi_Chr100000003    13_2
# 4 Polavi-1 52961 58022 Polavi_Chr100000004    13_3
# 5 Polavi-1 60547 69571 Polavi_Chr100000005    13_4
# 6 Polavi-1 70144 72735 Polavi_Chr100000006    13_5


# change order of columns
gff13IDs_ord <- gff13IDs[,c(1,5,2,3)]
head(gff13IDs_ord)
     # chrom genenum start   end
# 1 Polavi-1    13_0 14594 18403
# 2 Polavi-1    13_1 19174 20875
# 3 Polavi-1    13_2 21594 39613
# 4 Polavi-1    13_3 52961 58022
# 5 Polavi-1    13_4 60547 69571
# 6 Polavi-1    13_5 70144 72735

write.table(gff13IDs_ord, file="12.gff3", sep="\t", col.names = F, row.names = F, quote=FALSE)

###########################
mkdir orig_format
mv *.gff orig_format/
rename gff3 gff *

mv BLAST12_12.blast 12.blast
mv BLAST12_13.blast 12_13.blast

mv BLAST2_2.blast 2.blast
mv BLAST2_5.blast 2_5.blast

mv BLAST7_18.blast 7_18.blast
mv BLAST7_7.blast 7.blast

mv BLAST6_5.blast 6_5.blast
mv BLAST6_6.blast 6.blast

# join gff files

cat 12.gff 13.gff >> 12_13.gff
cat 2.gff 5.gff >> 2_5.gff
cat 7.gff 18.gff >> 7_18.gff
cat 6.gff 5.gff >> 6_5.gff

##############################
# Running

DupGen_finder.pl -i /home/celphin/scratch/Oxyria/DupGen_finder/data \
-t 12 -c 13 \
-o /home/celphin/scratch/Oxyria/DupGen_finder/output

# Reading BLAST file and pre-processing
# Generating BLAST list
# 371368 matches imported (148161 discarded)
# 11883 pairwise comparisons
# 1243 alignments generated
# Pairwise collinear blocks written to /home/celphin/scratch/Oxyria/DupGen_finder/output/12.collinearity [25.589 seconds elapsed]
# Reading BLAST file and pre-processing
# Generating BLAST list
# 459761 matches imported (0 discarded)
# 3134 pairwise comparisons
# 1094 alignments generated
# Pairwise collinear blocks written to /home/celphin/scratch/Oxyria/DupGen_finder/output/12_13.collinearity [19.656 seconds elapsed]


DupGen_finder.pl -i /home/celphin/scratch/Oxyria/DupGen_finder/data \
-t 2 -c 5 \
-o /home/celphin/scratch/Oxyria/DupGen_finder/output

DupGen_finder.pl -i /home/celphin/scratch/Oxyria/DupGen_finder/data \
-t 7 -c 18 \
-o /home/celphin/scratch/Oxyria/DupGen_finder/output

DupGen_finder.pl -i /home/celphin/scratch/Oxyria/DupGen_finder/data \
-t 6 -c 5 \
-o /home/celphin/scratch/Oxyria/DupGen_finder/output

# need to make more files to run the reverse

####################
cd /home/celphin/scratch/Oxyria/DupGen_finder/output

more 12.pairs.stats - Oxyria_digyna_H1.fa
Types   NO. of gene pairs
WGD-pairs       9779
TD-pairs        1595
PD-pairs        1321
TRD-pairs       16104
DSD-pairs       27726

more 6.pairs.stats - Draba_nivalis.fa
Types   NO. of gene pairs
WGD-pairs       4643
TD-pairs        2058
PD-pairs        1968
TRD-pairs       7105
DSD-pairs       24429

more 7.pairs.stats - Dryas_octopetala.fa
Types   NO. of gene pairs
WGD-pairs       18345
TD-pairs        2457
PD-pairs        2505
TRD-pairs       10313
DSD-pairs       31112

more 2.pairs.stats - Arabis_alpina.fa
Types   NO. of gene pairs
WGD-pairs       2310
TD-pairs        1338
PD-pairs        560
TRD-pairs       4728
DSD-pairs       14408

#########################
nano Arctic_plant_DupGen_results.csv

Spp,Types,Numb_gene_pairs
Oxyria_digyna,WGD-pairs,9779
Oxyria_digyna,TD-pairs,1595
Oxyria_digyna,PD-pairs,1321
Oxyria_digyna,TRD-pairs,16104
Oxyria_digyna,DSD-pairs,27726
Draba_nivalis,WGD-pairs,4643
Draba_nivalis,TD-pairs,2058
Draba_nivalis,PD-pairs,1968
Draba_nivalis,TRD-pairs,7105
Draba_nivalis,DSD-pairs,24429
Dryas_octopetala,WGD-pairs,18345
Dryas_octopetala,TD-pairs,2457
Dryas_octopetala,PD-pairs,2505
Dryas_octopetala,TRD-pairs,10313
Dryas_octopetala,DSD-pairs,31112
Arabis_alpina,WGD-pairs,2310
Arabis_alpina,TD-pairs,1338
Arabis_alpina,PD-pairs,560
Arabis_alpina,TRD-pairs,4728
Arabis_alpina,DSD-pairs,14408


##################################
cd /home/celphin/scratch/Oxyria/DupGen_finder/output

tmux new-session -s DupGen
tmux attach-session -t DupGen

module load StdEnv/2023
module load r/4.4.0
 
R 

library(ggplot2)

# https://datatricks.co.uk/multiple-bar-charts-in-r

data  <- read.csv("/home/celphin/scratch/Oxyria/DupGen_finder/output/Arctic_plant_DupGen_results.csv",sep=",")

theme <- theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               axis.line.x = element_line(color="#919191", size = 0.1),
               axis.line.y = element_line(color="#919191", size = 0.1),
               axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)
               )

grDevices::png(file = "Arctic_plant_DupGen_results.png", width = 1000, height = 800)
gg <- ggplot(data)+ 
    geom_bar(aes(x = Spp, y = Numb_gene_pairs, fill = Types), position = "stack", stat = "identity")+ 
    theme
print(gg)
grDevices::dev.off()

grDevices::pdf(file = "Arctic_plant_DupGen_results.pdf")
gg <- ggplot(data)+ 
    geom_bar(aes(x = Spp, y = Numb_gene_pairs, fill = Types), position = "stack", stat = "identity")+ 
    theme
print(gg)
grDevices::dev.off()

###############################
# Fraction of wgd pairs that are on different chromosomes

cd /home/celphin/scratch/Oxyria/DupGen_finder/output

#!/bin/bash

# Initialize counters
total_rows=0
different_count=0

# Read the file line by line, skipping the header
{
    read # skip header
    while read -r col1 id1 col2 id2 evalue; do
        # Extract the DoctH0-1 part from both columns
        id1=$(echo "$id1" | grep -o 'DoctH0-[0-9]*')
        id2=$(echo "$id2" | grep -o 'DoctH0-[0-9]*')

        # Only compare if both IDs are found
        if [[ -n "$id1" && -n "$id2" ]]; then
            # Increment the total row counter
            total_rows=$((total_rows + 1))

            # Print the current comparison for debugging
            #echo "Comparing: $id1 vs $id2"

            # Check if the extracted identifiers are different
            if [[ "$id1" != "$id2" ]]; then
                different_count=$((different_count + 1))
            fi
        fi
    done
} < 7.wgd.pairs

# Calculate the fraction
if [[ $total_rows -gt 0 ]]; then
    fraction=$(echo "scale=4; $different_count / $total_rows" | bc)
    echo "Fraction of rows with different chromosomes: $fraction"
else
    echo "No valid rows in the file."
fi

# .9089

#############################
# calculate Ka Ks and plot for duplciate genes
# https://github.com/qiao-xin/Scripts_for_GB

cd /home/celphin/scratch/Oxyria/DupGen_finder/

mkdir Ka_Ks; cd Ka_Ks

git clone https://github.com/qiao-xin/Scripts_for_GB.git

# NOTE need CDS 
cp /home/celphin/scratch/Dryas/Dryas_genomes/Liftoff/Dry-octo-H0_cds.fa ./data/

# Dependancies
/home/celphin/scratch/Oxyria/wgdi/pal2nal.v14/pal2nal.pl

# https://sourceforge.net/projects/kakscalculator2/
wget https://sourceforge.net/projects/kakscalculator2/files/KaKs_Calculator2.0.tar.gz/download
mv download KaKs_Calculator2.0.tar.gz
tar -xvzf KaKs_Calculator2.0.tar.gz
chmod -R  +755 *

export PATH=/home/celphin/scratch/Oxyria/DupGen_finder/Ka_Ks/KaKs_Calculator2.0/bin/Linux:$PATH
# test
KaKs_Calculator

#----------------
# run gene conversion

module load StdEnv/2020 gcc/9.3.0 bioperl/1.7.8 mafft/7.471 paml/4.9j clustal-omega/1.2.4

export PATH=/home/celphin/scratch/Oxyria/DupGen_finder/Ka_Ks/KaKs_Calculator2.0/bin/Linux:$PATH
export PATH=/home/celphin/scratch/Oxyria/wgdi/pal2nal.v14/:$PATH

perl ./Ka_Ks/Scripts_for_GB/calculate_Ka_Ks_pipeline/calculate_Ka_Ks_pipe.pl -d ./data/Dry-octo-H0_cds.fa -g ./output/7.wgd.pairs -o ./output/Dryas.td.kaks

# sh: KaKs_Calculator: command not found
# readline() on closed filehandle IN at ./Ka_Ks/Scripts_for_GB/calculate_Ka_Ks_pipeline/calculate_Ka_Ks_pipe.pl line 87.
# sed: can't read ./output/Dryas.td.kaks.KKC.format: No such file or directory

#---------------------
# try gene conversion

cd output
cp ../data/Dry-octo-H0_cds.fa 7.cds
perl ../Ka_Ks/Scripts_for_GB/detect_gene_conversion/GeConScan.pl 7 18 r 7.wgd.pairs

# MSG: Could not read file '18.cds': No such file or directory
# mv Rosa_rugosa.bed 18.gff
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_958449725.1/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/958/449/725/GCF_958449725.1_drRosRugo1.1/GCF_958449725.1_drRosRugo1.1_cds_from_genomic.fna.gz
gunzip GCF_958449725.1_drRosRugo1.1_cds_from_genomic.fna.gz
mv  GCF_958449725.1_drRosRugo1.1_cds_from_genomic.fna   18.cds

#--------------------------
# try again 
cd /home/celphin/scratch/Oxyria/DupGen_finder/output
cp 18* /home/celphin/scratch/Oxyria/DupGen_finder/Ka_Ks/Scripts_for_GB/detect_gene_conversion/
cp 7* /home/celphin/scratch/Oxyria/DupGen_finder/Ka_Ks/Scripts_for_GB/detect_gene_conversion/
cd /home/celphin/scratch/Oxyria/DupGen_finder/Ka_Ks/Scripts_for_GB/detect_gene_conversion/

module load StdEnv/2020 gcc/9.3.0 bioperl/1.7.8 mafft/7.471 paml/4.9j clustal-omega/1.2.4

export PATH=/home/celphin/scratch/Oxyria/DupGen_finder/Ka_Ks/KaKs_Calculator2.0/bin/Linux:$PATH
export PATH=/home/celphin/scratch/Oxyria/wgdi/pal2nal.v14/:$PATH
export PATH=/home/celphin/scratch/Oxyria/DupGen_finder/Ka_Ks/Scripts_for_GB/detect_gene_conversion:$PATH

perl GeConScan.pl 7 18 r 7.wgd.pairs

# ------------- EXCEPTION: Bio::Root::Exception -------------
# MSG: Could not read file 'cds/7.7.fasta': No such file or directory
# STACK: Error::throw
# STACK: Bio::Root::Root::throw /cvmfs/soft.computecanada.ca/easybuild/software/2020/Core/bioperl/1.7.8/lib/perl5/site_perl/5.30.2/Bio/Root/Root.pm:449
# STACK: Bio::Root::IO::_initialize_io /cvmfs/soft.computecanada.ca/easybuild/software/2020/Core/bioperl/1.7.8/lib/perl5/site_perl/5.30.2/Bio/Root/IO.pm:272
# STACK: Bio::SeqIO::_initialize /cvmfs/soft.computecanada.ca/easybuild/software/2020/Core/bioperl/1.7.8/lib/perl5/site_perl/5.30.2/Bio/SeqIO.pm:508
# STACK: Bio::SeqIO::fasta::_initialize /cvmfs/soft.computecanada.ca/easybuild/software/2020/Core/bioperl/1.7.8/lib/perl5/site_perl/5.30.2/Bio/SeqIO/fasta.pm:88
# STACK: Bio::SeqIO::new /cvmfs/soft.computecanada.ca/easybuild/software/2020/Core/bioperl/1.7.8/lib/perl5/site_perl/5.30.2/Bio/SeqIO.pm:384
# STACK: Bio::SeqIO::new /cvmfs/soft.computecanada.ca/easybuild/software/2020/Core/bioperl/1.7.8/lib/perl5/site_perl/5.30.2/Bio/SeqIO.pm:430
# STACK: main::main calculate.K.quartet.genes.CV.BP.all.pl:73
# STACK: calculate.K.quartet.genes.CV.BP.all.pl:44
