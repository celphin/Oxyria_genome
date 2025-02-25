#############################
# SYRI to compare structural variations
# Oxyria
# July 2024
# https://github.com/agshumate/Liftoff
##################################

# try SYRI
# https://github.com/schneebergerlab/syri
# https://schneebergerlab.github.io/syri/

# Install

# Prereqs
    # Python >=3.8 and the following packages: Cython-0.29.23, numpy-1.21.2, scipy-1.6.2, 
	# pandas-1.2.4, python-igraph-0.9.1, psutil-5.8.0, pysam-0.16.0.1, and matplotlib-3.3.4
    # C/C++ compiler: g++

module load StdEnv/2020 python/3.9 scipy-stack/2022a
virtualenv ~/syri
source ~/syri/bin/activate
pip install --upgrade pip --no-index
pip install git+https://github.com/schneebergerlab/syri.git
pip install git+https://github.com/schneebergerlab/plotsr.git
pip install pysam

pip list

syri -h
# works

######################
# to run 
# https://schneebergerlab.github.io/syri/pipeline.html

# Perform whole genome alignment
# Using minimap2 for generating alignment. Any other whole genome alignment tool can also be used.
#-----------------------
# remove all short scaffolds
# make copy of genomes
cd /home/celphin/scratch/Oxyria/GeneSpace/Oxyria_genomes/genomes/syri

module load seqkit/2.3.1
seqkit seq -m 5000000 Oxyria_Main.fasta > Oxyria_Main_short.fasta
seqkit seq -m 5000000 Oxy_Elles_Hap1.fasta > Oxy_Elles_Hap1_short.fasta
seqkit seq -m 5000000 Oxy_Elles_Hap2.fasta > Oxy_Elles_Hap2_short.fasta

seqkit seq -m 5000000 Oxy_Sval_h1.fasta > Oxy_Sval_h1_short.fasta
seqkit seq -m 5000000  Oxy_Sval_h2.fasta > Oxy_Sval_h2_short.fasta

seqkit seq -m 5000000 DToL_h1.fasta > DToL_h1_short.fasta
seqkit seq -m 5000000 DToL_h2.fasta > DToL_h2_short.fasta

##############################
module load bioawk/1.0
bioawk -c fastx '{ print $name, length($seq) }' < Oxy_Sval_h1_short.fasta
# h1tg000001l     79472951
# h1tg000002l     45063795
# h1tg000004l     16686886
# h1tg000005l     56138638
# h1tg000006l     36875860
# h1tg000007l     22272060
# h1tg000008l     17445183
# h1tg000009l     52776258
# h1tg000010l     37061391
# h1tg000011l     34330788
# h1tg000012l     9152510
# h1tg000013l     29747821
# h1tg000015l     19157742
# h1tg000017l     8009034
# h1tg000018l     10794573
# h1tg000019l     20250416
# h1tg000020l     6048833
# h1tg000022l     17901664
# h1tg000024l     14293109

bioawk -c fastx '{ print $name, length($seq) }' < Oxy_Elles_Hap2_short2.fasta
# HiC_scaffold_1  78501445
# HiC_scaffold_2  78458516
# HiC_scaffold_3  72720770
# HiC_scaffold_4  70815000
# HiC_scaffold_5  74731636
# HiC_scaffold_6  72112163
# HiC_scaffold_7  86137797

bioawk -c fastx '{ print $name, length($seq) }' < Oxy_Sval_h2_short.fasta
# h2tg000001l     19009535
# h2tg000002l     5917825
# h2tg000003l     45019965
# h2tg000005l     29386052
# h2tg000006l     5466857
# h2tg000007l     23490958
# h2tg000008l     55831360
# h2tg000009l     36798760
# h2tg000010l     5734725
# h2tg000011l     22121999
# h2tg000012l     32318812
# h2tg000013l     17392007
# h2tg000014l     36980249
# h2tg000015l     53952182
# h2tg000016l     28400861
# h2tg000017l     9140748
# h2tg000018l     10987281
# h2tg000020l     17151209
# h2tg000022l     8016925
# h2tg000024l     6091349
# h2tg000025l     15974674
# h2tg000028l     17399715
# h2tg000032l     14187580

bioawk -c fastx '{ print $name, length($seq) }' < DToL_h1_short.fasta
# OZ038361.1      88077321
# OZ038362.1      79981515
# OZ038363.1      78332497
# OZ038364.1      77668887
# OZ038365.1      76803065
# OZ038366.1      75968744
# OZ038367.1      74436958

bioawk -c fastx '{ print $name, length($seq) }' < DToL_h2_short.fasta
# CAXIWB010000078.1       86622132
# CAXIWB010000079.1       78290123
# CAXIWB010000080.1       76877250
# CAXIWB010000081.1       76726383
# CAXIWB010000082.1       75562310
# CAXIWB010000083.1       73468461
# CAXIWB010000084.1       73286991


##############################
# need to RagTag the Sval genomes
salloc -c40 --time 2:55:00 --mem 190000m --account def-rieseber

cd /home/celphin/scratch/Oxyria/GeneSpace/Oxyria_genomes/genomes/syri
source ~/RagTag/bin/activate

module load StdEnv/2020
module load python/3.10.2
module load minimap2/2.24
module load scipy-stack/2022a

pip list

# Try running
/home/celphin/scratch/RagTag/RagTag-master/ragtag.py \
scaffold Oxyria_Main.fasta Oxy_Sval_h1.fasta \
-t 39 -u -o ./RagTag_out_Sval1/

/home/celphin/scratch/RagTag/RagTag-master/ragtag.py \
scaffold Oxyria_Main.fasta Oxy_Sval_h2.fasta \
-t 39 -u -o ./RagTag_out_Sval2/

###############################
# Make short
module load StdEnv/2023 seqkit/2.5.1

RagTag_out_Sval1
cp  ragtag.scaffold.fasta ../Oxy_Sval_h1RT.fasta
cd ..
seqkit seq -m 5000000 Oxy_Sval_h1RT.fasta > Oxy_Sval_h1_short.fasta
bioawk -c fastx '{ print $name, length($seq) }' < Oxy_Sval_h1_short.fasta
# Oxy-1-88146246_RagTag   86582034
# Oxy-2-80787722_RagTag   79714091
# Oxy-3-77175000_RagTag   79472951
# Oxy-4-76036264_RagTag   78410798
# Oxy-5-73842794_RagTag   73303751
# Oxy-6-73289246_RagTag   72361354
# Oxy-7-70136240_RagTag   76064323

cd RagTag_out_Sval2
cp  ragtag.scaffold.fasta ../Oxy_Sval_h2RT.fasta
cd ..
seqkit seq -m 5000000  Oxy_Sval_h2RT.fasta > Oxy_Sval_h2_short.fasta
bioawk -c fastx '{ print $name, length($seq) }' < Oxy_Sval_h2_short.fasta
# Oxy-1-88146246_RagTag   88360298
# Oxy-2-80787722_RagTag   79729308
# Oxy-3-77175000_RagTag   76449106
# Oxy-4-76036264_RagTag   78099446
# Oxy-5-73842794_RagTag   71352700
# Oxy-6-73289246_RagTag   71762064
# Oxy-7-70136240_RagTag   74444656

#################################
# Rename contigs of all genomes

#awk '/^>/{print ">Oxyria_Main" ++i; next}{print}' < Oxyria_Main_short.fasta
#awk '/^>/{print ">Elles_Hap1" ++i; next}{print}' < Oxy_Elles_Hap1_short.fasta
#awk '/^>/{print ">Sval_h1" ++i; next}{print}' < Oxy_Sval_h1_short.fasta

awk '/^>/{print ">Elles_Hap1" ++i; next}{print}' < Oxy_Elles_Hap1_short3.fasta > Oxy_Elles_Hap1_short4.fasta
awk '/^>/{print ">Sval_h2" ++i; next}{print}' < Oxy_Sval_h2_short.fasta > Oxy_Sval_h2_short1.fasta
awk '/^>/{print ">DToL_h1" ++i; next}{print}' < DToL_h1_short.fasta > DToL_h1_short1.fasta
awk '/^>/{print ">DToL_h2" ++i; next}{print}' < DToL_h2_short.fasta > DToL_h2_short1.fasta

mv Oxy_Sval_h2_short1.fasta Oxy_Sval_h2_short.fasta
mv DToL_h1_short1.fasta DToL_h1_short.fasta
mv DToL_h2_short1.fasta DToL_h2_short.fasta

##################################
# Minimap2

tmux new-session -s syri3
tmux attach-session -t syri3

salloc -c40 --time 2:55:00 --mem 191000M --account def-rieseber

cd /home/celphin/scratch/Oxyria/GeneSpace/Oxyria_genomes/genomes/syri
module load StdEnv/2020 python/3.9 scipy-stack/2022a
module load minimap2/2.24

minimap2 -ax asm5 -t 40 --eqx Oxyria_Main_short.fasta Oxy_Elles_Hap1_short.fasta > Main_Elles1_out.sam

minimap2 -ax asm5 -t 40 --eqx Oxy_Elles_Hap1_short2.fasta Oxy_Elles_Hap2_short.fasta > Elles1_Elles2_out.sam

minimap2 -ax asm5 -t 40 --eqx Oxy_Elles_Hap2_short2.fasta Oxy_Sval_h1_short.fasta > Elles2_Sval1_out.sam

minimap2 -ax asm5 -t 40 --eqx Oxy_Sval_h1_short2.fasta Oxy_Sval_h2_short.fasta > Sval1_Sval2_out.sam

minimap2 -ax asm5 -t 40 --eqx Oxy_Sval_h2_short2.fasta DToL_h1_short.fasta > Sval2_DToL1_out.sam

minimap2 -ax asm5 -t 40 --eqx DToL_h1_short2.fasta DToL_h2_short.fasta > DToL1_DToL2_out.sam

################################
# Remap

salloc -c40 --time 2:55:00 --mem 191000M --account def-rieseber

cd /home/celphin/scratch/Oxyria/GeneSpace/Oxyria_genomes/genomes/syri
module load StdEnv/2020 python/3.9 scipy-stack/2022a
module load minimap2/2.24

minimap2 -ax asm5 -t 10 --eqx Oxyria_Main_short.fasta Oxy_Elles_Hap1_short2.fasta > Main_Elles1_out2.sam

minimap2 -ax asm5 -t 10 --eqx Oxy_Elles_Hap1_short5.fasta Oxy_Elles_Hap2_short2.fasta > Elles1_Elles2_out3.sam

minimap2 -ax asm5 -t 10 --eqx Oxy_Elles_Hap2_short2.fasta Oxy_Sval_h1_short2.fasta > Elles2_Sval1_out2.sam

minimap2 -ax asm5 -t 10 --eqx Oxy_Sval_h1_short2.fasta Oxy_Sval_h2_short2.fasta > Sval1_Sval2_out2.sam

minimap2 -ax asm5 -t 10 --eqx Oxy_Sval_h2_short2.fasta DToL_h1_short2.fasta > Sval2_DToL1_out2.sam

minimap2 -ax asm5 -t 10 --eqx DToL_h1_short2.fasta DToL_h2_short2.fasta > DToL1_DToL2_out2.sam

####################################
# Ryn SyRI
# Elles Hap 1
cd /home/celphin/scratch/Oxyria/GeneSpace/Oxyria_genomes/genomes/syri

mkdir Main_Elles1_out
mv Main_Elles1_out.sam Main_Elles1_out
cd Main_Elles1_out

#---------------
module load StdEnv/2020 python/3.9 scipy-stack/2022a
module load minimap2/2.24
source ~/syri/bin/activate

syri -c Main_Elles1_out.sam -r ../Oxyria_Main_short.fasta -q ../Oxy_Elles_Hap1_short.fasta -k -F S

# Reading Coords - WARNING - Chromosomes IDs do not match.
# Reading Coords - WARNING - Matching them automatically. For each reference genome, most similar query genome will be selected. Check mapids.txt for mapping used.
# Reading Coords - WARNING - Reference chromosome Oxy-2-80787722 has high fraction of inverted alignments with its homologous chromosome in the query genome (HiC_scaffold_5). Ensure that same chromosome-strands are being compared in the two genomes, as different strand can result in unexpected errors.
# Reading Coords - WARNING - Reference chromosome Oxy-4-76036264 has high fraction of inverted alignments with its homologous chromosome in the query genome (HiC_scaffold_1). Ensure that same chromosome-strands are being compared in the two genomes, as different strand can result in unexpected errors.
# Reading Coords - WARNING - Reference chromosome Oxy-6-73289246 has high fraction of inverted alignments with its homologous chromosome in the query genome (HiC_scaffold_2). Ensure that same chromosome-strands are being compared in the two genomes, as different strand can result in unexpected errors.
# Reading Coords - WARNING - Reference chromosome Oxy-7-70136240 has high fraction of inverted alignments with its homologous chromosome in the query genome (HiC_scaffold_4). Ensure that same chromosome-strands are being compared in the two genomes, as different strand can result in unexpected errors.
# syri.Oxy-2-80787722 - ERROR - No syntenic region found for chromosome: Oxy-2-80787722. This is potentially caused by the two assemblies having different strands for this chromosomes. Reverse complement the chromosome to ensure that the same strands are analysed. Exiting.
# syri.Oxy-4-76036264 - ERROR - No syntenic region found for chromosome: Oxy-4-76036264. This is potentially caused by the two assemblies having different strands for this chromosomes. Reverse complement the chromosome to ensure that the same strands are analysed. Exiting.
# syri.Oxy-6-73289246 - ERROR - No syntenic region found for chromosome: Oxy-6-73289246. This is potentially caused by the two assemblies having different strands for this chromosomes. Reverse complement the chromosome to ensure that the same strands are analysed. Exiting.
# syri.Oxy-7-70136240 - ERROR - No syntenic region found for chromosome: Oxy-7-70136240. This is potentially caused by the two assemblies having different strands for this chromosomes. Reverse complement the chromosome to ensure that the same strands are analysed. Exiting.

#----
# Align genomes to each other - remove backwards chromosomes
# Manually rearrange

# Elles1
# syri.Oxy-2-80787722 
# syri.Oxy-4-76036264 
# syri.Oxy-6-73289246 
# syri.Oxy-7-70136240

# Oxy-1-88146246  HiC_scaffold_7
# Oxy-2-80787722  HiC_scaffold_5
# Oxy-3-77175000  HiC_scaffold_6
# Oxy-4-76036264  HiC_scaffold_1
# Oxy-5-73842794  HiC_scaffold_3
# Oxy-6-73289246  HiC_scaffold_2
# Oxy-7-70136240  HiC_scaffold_4

# split the fasta by chromosome
mkdir chromosome_files; cd chromosome_files
awk '/^>/ { file=substr($1,2) ".fasta" } { print > file }' ../../Oxy_Elles_Hap1_short.fasta
find -type f -size -5000000c -delete

# reverse complement the sequence 
seqkit seq -t DNA --reverse --complement HiC_scaffold_5.fasta > RHiC_scaffold_5.fasta
seqkit seq -t DNA --reverse --complement HiC_scaffold_1.fasta > RHiC_scaffold_1.fasta
seqkit seq -t DNA --reverse --complement HiC_scaffold_2.fasta > RHiC_scaffold_2.fasta
seqkit seq -t DNA --reverse --complement HiC_scaffold_4.fasta > RHiC_scaffold_4.fasta

cat HiC_scaffold_7.fasta  \
RHiC_scaffold_5.fasta  \
HiC_scaffold_6.fasta  \
RHiC_scaffold_1.fasta  \
HiC_scaffold_3.fasta \
RHiC_scaffold_2.fasta \
RHiC_scaffold_4.fasta > ../../Oxy_Elles_Hap1_short2.fasta

cd ../..

#--------
# Rerun syri

mv Main_Elles1_out2.sam Main_Elles1_out
cd Main_Elles1_out
mkdir txt
mv *.txt txt

module load StdEnv/2020 python/3.9 scipy-stack/2022a
module load minimap2/2.24
source ~/syri/bin/activate

syri -c Main_Elles1_out2.sam -r ../Oxyria_Main_short.fasta -q ../Oxy_Elles_Hap1_short2.fasta -k -F S

# Reading Coords - WARNING - Chromosomes IDs do not match.
# Reading Coords - WARNING - Matching them automatically. For each reference genome, most similar query genome will be selected. Check mapids.txt for mapping used.

#----------
# Output
more syri.summary

#Variation_type Count   Length_ref      Length_qry
Syntenicregions        374     451340179       441587762
Inversions      265     27528133        26560023
Translocations  195     7177247 6553530
Duplications(reference)        39      383625  -
Duplications(query)    278     -       1726194
Notaligned(reference) 855     53073235        -
Notaligned(query)     1055    -       34770628
SNPs    310180  310180  310180
Insertions      22643   -       1482208
Deletions       22100   1522759 -
Copygains       39      -       200282
Copylosses      50      179949  -
Highlydiverged 3824    119572015       108246914
Tandemrepeats  11      4643    3655


####################################
# Ryn SyRI
# Elles H2
cd /home/celphin/scratch/Oxyria/GeneSpace/Oxyria_genomes/genomes/syri

mkdir Elles1_Elles2_out
mv Elles1_Elles2_out.sam Elles1_Elles2_out
cd Elles1_Elles2_out

#---------------
# Elles Hap 2
module load StdEnv/2020 python/3.9 scipy-stack/2022a
module load minimap2/2.24
source ~/syri/bin/activate

syri -c Elles1_Elles2_out.sam -r ../Oxy_Elles_Hap1_short2.fasta -q ../Oxy_Elles_Hap2_short.fasta -k -F S

# Reading Coords - WARNING - Chromosomes IDs do not match.
# Reading Coords - WARNING - Matching them automatically. For each reference genome, most similar query genome will be selected. Check mapids.txt for mapping used.
# Reading Coords - WARNING - Reference chromosome HiC_scaffold_3  (Elles_Hap24). 

deactivate
#----
# Align genomes to each other - remove backwards chromosomes
# Manually rearrange
# Elles_Hap2
HiC_scaffold_1
HiC_scaffold_3

# split the fasta by chromosome
mkdir chromosome_files; cd chromosome_files
awk '/^>/ { file=substr($1,2) ".fasta" } { print > file }' ../../Oxy_Elles_Hap1_short2.fasta
find -type f -size -5000000c -delete

# reverse complement the sequence 
seqkit seq -t DNA --reverse --complement HiC_scaffold_6.fasta > RHiC_scaffold_6.fasta
seqkit seq -t DNA --reverse --complement HiC_scaffold_3.fasta > RHiC_scaffold_3.fasta
seqkit seq -t DNA --reverse --complement HiC_scaffold_4.fasta > RHiC_scaffold_4.fasta

cat HiC_scaffold_1.fasta \
HiC_scaffold_2.fasta \
HiC_scaffold_3.fasta \
RHiC_scaffold_4.fasta \
HiC_scaffold_5.fasta \
HiC_scaffold_6.fasta \
HiC_scaffold_7.fasta \
> ../../Oxy_Elles_Hap1_short3.fasta

cd ../..

#--------
# Rerun syri

mv Elles1_Elles2_out2.sam Elles1_Elles2_out
cd Elles1_Elles2_out
mkdir txt
mv *.txt txt

module load StdEnv/2020 python/3.9 scipy-stack/2022a
module load minimap2/2.24
source ~/syri/bin/activate

syri -c Elles1_Elles2_out2.sam -r ../Oxy_Elles_Hap1_short4.fasta -q ../Oxy_Elles_Hap2_short2.fasta -k -F S

# Reading Coords - WARNING - Chromosomes IDs do not match.
# Reading Coords - WARNING - Matching them automatically. For each reference genome, most similar query genome will be selected. Check mapids.txt for mapping used.
# Reading Coords - WARNING - Reference chromosome Elles_Hap13 has high fraction of inverted alignments with its homologous chromosome in the query genome (HiC_scaffold_4). Ensure that same chromosome-strands are being compared in the two genomes, as different strand can result in unexpected errors.
# Reading Coords - WARNING - Reference chromosome Elles_Hap16 has high fraction of inverted alignments with its homologous chromosome in the query genome (HiC_scaffold_1). Ensure that same chromosome-strands are being compared in the two genomes, as different strand can result in unexpected errors.

#--------------------------------
# split the fasta by chromosome
mkdir chromosome_files2; cd chromosome_files2
awk '/^>/ { file=substr($1,2) ".fasta" } { print > file }' ../../Oxy_Elles_Hap1_short4.fasta
find -type f -size -5000000c -delete

# reverse complement the sequence 
seqkit seq -t DNA --reverse --complement Elles_Hap13.fasta > RElles_Hap13.fasta
seqkit seq -t DNA --reverse --complement Elles_Hap16.fasta > RElles_Hap16.fasta

cat Elles_Hap11.fasta \
Elles_Hap12.fasta \
RElles_Hap13.fasta \
Elles_Hap14.fasta \
Elles_Hap15.fasta \
RElles_Hap16.fasta \
Elles_Hap17.fasta \
> ../../Oxy_Elles_Hap1_short5.fasta

cd ../..

#-------------------------
mv Elles1_Elles2_out3.sam Elles1_Elles2_out
cd Elles1_Elles2_out
mkdir txt
mv *.txt txt

module load StdEnv/2020 python/3.9 scipy-stack/2022a
module load minimap2/2.24
source ~/syri/bin/activate

syri -c Elles1_Elles2_out3.sam -r ../Oxy_Elles_Hap1_short5.fasta -q ../Oxy_Elles_Hap2_short2.fasta -k -F S

# Reading Coords - WARNING - Chromosomes IDs do not match.
# Reading Coords - WARNING - Matching them automatically. For each reference genome, most similar query genome will be selected. Check mapids.txt for mapping used.

#----------
# Output

#Variation_type Count   Length_ref      Length_qry
Syntenicregions        117     62648856        64052875
Inversions      66      80203296        132227344
Translocations  1624    166450139       171232184
Duplications(reference)        1020    72130179        -
Duplications(query)    827     -       29820077
Notaligned(reference) 2594    171078153       -
Notaligned(query)     2479    -       157501083
SNPs    442598  442598  442598
Insertions      31532   -       1985645
Deletions       32156   2087590 -
Copygains       67      -       466450
Copylosses      74      280514  -
Highlydiverged 4411    155339687       213465333
Tandemrepeats  9       2052    1904


#Variation_type Count   Length_ref      Length_qry
Syntenicregions        119     62629888        63965689
Inversions      63      75695337        125180850
Translocations  1661    189229507       195388331
Duplications(reference)        895     74234541        -
Duplications(query)    772     -       30850300
Not aligned(reference) 2448    160204700       -
Not aligned(query)     2529    -       157152110
SNPs    475360  475360  475360
Insertions      34340   -       2269775
Deletions       34930   2311549 -
Copygains       66      -       443752
Copylosses      75      312992  -
Highlydiverged 4845    162515983       219406948
Tandemrepeats  11      2244    2068


#Variation_type Count   Length_ref      Length_qry
Syntenicregions        751     382159136       386886576
Inversions      328     35689088        36255983
Translocations  551     6283225 6256819
Duplications(reference)        110     1051100 -
Duplications(query)    496     -       2496246
Notaligned(reference) 1688    86966130        -
Notaligned(query)     2061    -       102144401
SNPs    682326  682326  682326
Insertions      47390   -       2855694
Deletions       47962   2987118 -
Copygains       93      -       348119
Copylosses      103     463725  -
Highlydiverged 8087    220102009       225610661
Tandemrepeats  21      4388    4880


########################
# Sval H1

cd /home/celphin/scratch/Oxyria/GeneSpace/Oxyria_genomes/genomes/syri

mkdir Elles2_Sval1_out
mv Elles2_Sval1_out.sam Elles2_Sval1_out
cd Elles2_Sval1_out

#---------------
# Elles Hap 2
module load StdEnv/2020 python/3.9 scipy-stack/2022a
module load minimap2/2.24
source ~/syri/bin/activate

syri -c Elles2_Sval1_out.sam -r ../Oxy_Elles_Hap2_short2.fasta -q ../Oxy_Sval_h1_short.fasta -k -F S

# 2024-07-26 13:17:13,224 - Reading Coords - WARNING - syri:134 - Reference chromosome HiC_scaffold_1
# 2024-07-26 13:17:13,249 - Reading Coords - WARNING - syri:134 - Reference chromosome HiC_scaffold_3
# 2024-07-26 13:17:13,261 - Reading Coords - WARNING - syri:134 - Reference chromosome HiC_scaffold_4

# HiC_scaffold_1  Oxy-3-77175000_RagTag
# HiC_scaffold_2  Oxy-2-80787722_RagTag
# HiC_scaffold_3  Oxy-7-70136240_RagTag
# HiC_scaffold_4  Oxy-5-73842794_RagTag
# HiC_scaffold_5  Oxy-4-76036264_RagTag
# HiC_scaffold_6  Oxy-6-73289246_RagTag
# HiC_scaffold_7  Oxy-1-88146246_RagTag

deactivate

#----
# Align genomes to each other - remove backwards chromosomes
# Manually rearrange

# Sval_h1
HiC_scaffold_1 
HiC_scaffold_3 
HiC_scaffold_4

# split the fasta by chromosome
mkdir chromosome_files; cd chromosome_files
awk '/^>/ { file=substr($1,2) ".fasta" } { print > file }' ../../Oxy_Sval_h1_short.fasta
find -type f -size -5000000c -delete

# reverse complement the sequence 
seqkit seq -t DNA --reverse --complement Oxy-3-77175000_RagTag.fasta > ROxy-3-77175000_RagTag.fasta
seqkit seq -t DNA --reverse --complement Oxy-7-70136240_RagTag.fasta > ROxy-7-70136240_RagTag.fasta
seqkit seq -t DNA --reverse --complement Oxy-5-73842794_RagTag.fasta > ROxy-5-73842794_RagTag.fasta


cat ROxy-3-77175000_RagTag.fasta \
Oxy-2-80787722_RagTag.fasta \
ROxy-7-70136240_RagTag.fasta \
ROxy-5-73842794_RagTag.fasta \
Oxy-4-76036264_RagTag.fasta \
Oxy-6-73289246_RagTag.fasta \
Oxy-1-88146246_RagTag.fasta \
> ../../Oxy_Sval_h1_short2.fasta

cd ../..

#--------
# Rerun syri

mv Elles2_Sval1_out2.sam Elles2_Sval1_out
cd Elles2_Sval1_out
mkdir txt
mv *.txt txt

module load StdEnv/2020 python/3.9 scipy-stack/2022a
module load minimap2/2.24
source ~/syri/bin/activate

syri -c Elles2_Sval1_out2.sam -r ../Oxy_Elles_Hap2_short2.fasta -q ../Oxy_Sval_h1_short2.fasta -k -F S

# Reading Coords - WARNING - Chromosomes IDs do not match.
# Reading Coords - WARNING - Matching them automatically. For each reference genome, most similar query genome will be selected. Check mapids.txt for mapping used.


#----------
# Output

#Variation_type Count   Length_ref      Length_qry
Syntenicregions        662     430842577       433438977
Inversions      202     39647360        41118918
Translocations  468     3623474 3491714
Duplications(reference)        204     1538832 -
Duplications(query)    297     -       1239060
Notaligned(reference) 1405    58500459        -
Notaligned(query)     1577    -       66611189
SNPs    1106783 1106783 1106783
Insertions      74247   -       2653306
Deletions       76496   2641069 -
Copygains       203     -       900201
Copylosses      171     737949  -
Highlydiverged 14751   350403040       354149634
Tandemrepeats  30      12356   9953


###########################
# Sval H2

cd /home/celphin/scratch/Oxyria/GeneSpace/Oxyria_genomes/genomes/syri

mkdir Sval1_Sval2_out
mv Sval1_Sval2_out.sam Sval1_Sval2_out
cd Sval1_Sval2_out

#---------------
module load StdEnv/2020 python/3.9 scipy-stack/2022a
module load minimap2/2.24
source ~/syri/bin/activate

syri -c Sval1_Sval2_out.sam -r ../Oxy_Sval_h1_short2.fasta -q ../Oxy_Sval_h2_short.fasta -k -F S

deactivate
#----
# Align genomes to each other - remove backwards chromosomes
# Manually rearrange

# Sval_h2
# syri.Oxy-3-77175000_RagTag - ERROR - No syntenic region found for chromosome: Oxy-3-77175000_RagTag. This is potentially caused by the two assemblies having different strands for this chromosomes. Reverse complement the chromosome to ensure that the same strands are analysed. Exiting.
# syri.Oxy-5-73842794_RagTag - ERROR - No syntenic region found for chromosome: Oxy-5-73842794_RagTag. This is potentially caused by the two assemblies having different strands for this chromosomes. Reverse complement the chromosome to ensure that the same strands are analysed. Exiting.
# syri.Oxy-7-70136240_RagTag - ERROR - No syntenic region found for chromosome: Oxy-7-70136240_RagTag. This is potentially caused by the two assemblies having different strands for this chromosomes. Reverse complement the chromosome to ensure that the same strands are analysed. Exiting.

# Oxy-1-88146246_RagTag   Sval_h21
# Oxy-2-80787722_RagTag   Sval_h22
*# Oxy-3-77175000_RagTag   Sval_h23
# Oxy-4-76036264_RagTag   Sval_h24
*# Oxy-5-73842794_RagTag   Sval_h25
# Oxy-6-73289246_RagTag   Sval_h26
*# Oxy-7-70136240_RagTag   Sval_h27

# split the fasta by chromosome
mkdir chromosome_files; cd chromosome_files
awk '/^>/ { file=substr($1,2) ".fasta" } { print > file }' ../../Oxy_Sval_h2_short.fasta
find -type f -size -5000000c -delete

# reverse complement the sequence 
seqkit seq -t DNA --reverse --complement Sval_h23.fasta > RSval_h23.fasta
seqkit seq -t DNA --reverse --complement Sval_h25.fasta > RSval_h25.fasta
seqkit seq -t DNA --reverse --complement Sval_h27.fasta > RSval_h27.fasta

cat Sval_h21.fasta \
Sval_h22.fasta \
RSval_h23.fasta \
Sval_h24.fasta \
RSval_h25.fasta \
Sval_h26.fasta \
RSval_h27.fasta \
> ../../Oxy_Sval_h2_short2.fasta

cd ../..

#--------
# Rerun syri

mv Sval1_Sval2_out2.sam Sval1_Sval2_out
cd Sval1_Sval2_out
mkdir txt
mv *.txt txt

module load StdEnv/2020 python/3.9 scipy-stack/2022a
module load minimap2/2.24
source ~/syri/bin/activate

syri -c Sval1_Sval2_out2.sam -r ../Oxy_Sval_h1_short2.fasta -q ../Oxy_Sval_h2_short2.fasta -k -F S

# Reading Coords - WARNING - Chromosomes IDs do not match.
# Reading Coords - WARNING - Matching them automatically. For each reference genome, most similar query genome will be selected. Check mapids.txt for mapping used.

#----------
# Output

#Variation_type Count   Length_ref      Length_qry
Syntenicregions        69      529217537       522199908
Inversions      39      8572798 7946526
Translocations  17      2041894 2062708
Duplications(reference)        5       11100   -
Duplications(query)    66      -       510997
Notaligned(reference) 113     7002547 -
Notaligned(query)     165     -       7536145
SNPs    275336  275336  275336
Insertions      21718   -       1953385
Deletions       21658   2771239 -
Copygains       25      -       159803
Copylosses      40      312377  -
Highlydiverged 3492    125664619       119011434
Tandemrepeats  4       785     1408



########################################
# DToL H1
cd /home/celphin/scratch/Oxyria/GeneSpace/Oxyria_genomes/genomes/syri

mkdir Sval2_DToL1_out
mv Sval2_DToL1_out.sam Sval2_DToL1_out
cd Sval2_DToL1_out

#---------------
module load StdEnv/2020 python/3.9 scipy-stack/2022a
module load minimap2/2.24
source ~/syri/bin/activate

syri -c Sval2_DToL1_out.sam -r ../Oxy_Sval_h2_short2.fasta -q ../DToL_h1_short.fasta -k -F S

Reading Coords - WARNING - Chromosomes IDs do not match.
Reading Coords - WARNING - Matching them automatically. For each reference genome, most similar query genome will be selected. Check mapids.txt for mapping used.
Reading Coords - WARNING - Reference chromosome Sval_h22 has homologous chromosome in the query genome (DToL_h12). Ensure that same chromosome-strands are being compared in the two genomes, as different strand can result in unexpected errors.
Reading Coords - WARNING - Reference chromosome Sval_h23 has homologous chromosome in the query genome (DToL_h13). Ensure that same chromosome-strands are being compared in the two genomes, as different strand can result in unexpected errors.
Reading Coords - WARNING - Reference chromosome Sval_h24 has homologous chromosome in the query genome (DToL_h14). Ensure that same chromosome-strands are being compared in the two genomes, as different strand can result in unexpected errors.
Reading Coords - WARNING - Reference chromosome Sval_h25 has homologous chromosome in the query genome (DToL_h17). Ensure that same chromosome-strands are being compared in the two genomes, as different strand can result in unexpected errors.
Reading Coords - WARNING - Reference chromosome Sval_h26 has homologous chromosome in the query genome (DToL_h15). Ensure that same chromosome-strands are being compared in the two genomes, as different strand can result in unexpected errors.

deactivate
#----
# Align genomes to each other - remove backwards chromosomes
# Manually rearrange

# DToL h1
# Oxy-5-73842794  HiC_scaffold_4

# split the fasta by chromosome
mkdir chromosome_files; cd chromosome_files
awk '/^>/ { file=substr($1,2) ".fasta" } { print > file }' ../../DToL_h1_short.fasta
find -type f -size -5000000c -delete

# reverse complement the sequence 
seqkit seq -t DNA --reverse --complement DToL_h12.fasta > RDToL_h12.fasta
seqkit seq -t DNA --reverse --complement DToL_h13.fasta > RDToL_h13.fasta
seqkit seq -t DNA --reverse --complement DToL_h14.fasta > RDToL_h14.fasta
seqkit seq -t DNA --reverse --complement DToL_h17.fasta > RDToL_h17.fasta
seqkit seq -t DNA --reverse --complement DToL_h15.fasta > RDToL_h15.fasta

cat DToL_h11.fasta \
RDToL_h12.fasta \
RDToL_h13.fasta \
RDToL_h14.fasta \
RDToL_h15.fasta \
DToL_h16.fasta \
RDToL_h17.fasta \
> ../../DToL_h1_short2.fasta

cd ../..

#--------
# Rerun syri

mv Sval2_DToL1_out2.sam Sval2_DToL1_out
cd Sval2_DToL1_out
mkdir txt
mv *.txt txt

module load StdEnv/2020 python/3.9 scipy-stack/2022a
module load minimap2/2.24
source ~/syri/bin/activate

syri -c Sval2_DToL1_out2.sam -r ../Oxy_Sval_h2_short2.fasta -q ../DToL_h1_short2.fasta -k -F S

# Reading Coords - WARNING - Chromosomes IDs do not match.
# Reading Coords - WARNING - Matching them automatically. For each reference genome, most similar query genome will be selected. Check mapids.txt for mapping used.

#----------
# Output

#Variation_type Count   Length_ref      Length_qry
Syntenicregions        295     464802512       469421242
Inversions      128     38894568        40836682
Translocations  198     3333317 3307397
Duplications(reference)        68      440149  -
Duplications(query)    210     -       1397440
Notaligned(reference) 613     33280212        -
Notaligned(query)     729     -       36561275
SNPs    784588  784588  784588
Insertions      57813   -       4836532
Deletions       58177   4710535 -
Copygains       101     -       732724
Copylosses      109     560150  -
Highlydiverged 8617    272508669       278764065
Tandemrepeats  20      17914   18531


#################################
# DToL H2

cd /home/celphin/scratch/Oxyria/GeneSpace/Oxyria_genomes/genomes/syri

mkdir DToL1_DToL2_out
mv DToL1_DToL2_out.sam DToL1_DToL2_out
cd DToL1_DToL2_out

#---------------
module load StdEnv/2020 python/3.9 scipy-stack/2022a
module load minimap2/2.24
source ~/syri/bin/activate

syri -c DToL1_DToL2_out.sam -r ../DToL_h1_short2.fasta -q ../DToL_h1_short2.fasta -k -F S


# Reading Coords - WARNING - Chromosomes IDs do not match.
# Reading Coords - WARNING - Matching them automatically. For each reference genome, most similar query genome will be selected. Check mapids.txt for mapping used.
# Reading Coords - WARNING - Reference chromosome DToL_h15   (DToL_h27). Filtering out all corresponding alignments.
# Reading Coords - WARNING - Reference chromosome DToL_h12   (DToL_h22). Ensure that same chromosome-strands are being compared in the two genomes, as different strand can result in unexpected errors.
# Reading Coords - WARNING - Reference chromosome DToL_h16   (DToL_h25). Ensure that same chromosome-strands are being compared in the two genomes, as different strand can result in unexpected errors.

deactivate
#----
# Align genomes to each other - remove backwards chromosomes
# Manually rearrange

# DToL h1
# Oxy-5-73842794  HiC_scaffold_4

# split the fasta by chromosome
mkdir chromosome_files; cd chromosome_files
awk '/^>/ { file=substr($1,2) ".fasta" } { print > file }' ../../DToL_h2_short.fasta
find -type f -size -5000000c -delete

# reverse complement the sequence 
seqkit seq -t DNA --reverse --complement DToL_h27.fasta > RDToL_h27.fasta
seqkit seq -t DNA --reverse --complement DToL_h22.fasta > RDToL_h22.fasta
seqkit seq -t DNA --reverse --complement DToL_h25.fasta > RDToL_h25.fasta

cat DToL_h21.fasta \
RDToL_h22.fasta \
DToL_h23.fasta \
DToL_h24.fasta \
RDToL_h25.fasta \
DToL_h26.fasta \
RDToL_h27.fasta \
> ../../DToL_h2_short2.fasta

cd ../..

#--------
# Rerun syri

mv DToL1_DToL2_out2.sam DToL1_DToL2_out
cd DToL1_DToL2_out
mkdir txt
mv *.txt txt

module load StdEnv/2020 python/3.9 scipy-stack/2022a
module load minimap2/2.24
source ~/syri/bin/activate

syri -c DToL1_DToL2_out2.sam -r ../DToL_h1_short2.fasta -q ../DToL_h2_short2.fasta -k -F S

# Reading Coords - WARNING - Chromosomes IDs do not match.
# Reading Coords - WARNING - Matching them automatically. For each reference genome, most similar query genome will be selected. Check mapids.txt for mapping used.

#----------
# Output

#Variation_type Count   Length_ref      Length_qry
Syntenicregions        24      541353141       535591558
Inversions      13      1314951 1370858
Translocations  24      277487  273910
Duplications(reference)        6       31394   -
Duplications(query)    46      -       326021
Notaligned(reference) 51      8352559 -
Notaligned(query)     68      -       3381632
SNPs    54608   54608   54608
Insertions      4667    -       1342035
Deletions       4848    3787557 -
Copygains       5       -       50833
Copylosses      8       266308  -
Highlydiverged 343     15958037        12976761
Tandemrepeats  2       392     244



######################################
#Plotting
# https://github.com/schneebergerlab/plotsr

cd /home/celphin/scratch/Oxyria/GeneSpace/Oxyria_genomes/genomes/syri/

nano genomes.txt
#file	name	tags
Oxy_Elles_Hap1_short5.fasta	Elles1	lw:1.5
Oxy_Elles_Hap2_short2.fasta	Elles2	lw:1.5
Oxy_Sval_h1_short2.fasta	Sval1	lw:1.5
Oxy_Sval_h2_short2.fasta	Sval2	lw:1.5
DToL_h1_short2.fasta	DToL1	lw:1.5
DToL_h2_short2.fasta	DToL2	lw:1.5

plotsr --sr ./Elles1_Elles2_out/syri.out \
       --sr ./Elles2_Sval1_out/syri.out \
       --sr ./Sval1_Sval2_out/syri.out \
       --sr ./Sval2_DToL1_out/syri.out \
       --sr ./DToL1_DToL2_out/syri.out \
       --genomes genomes.txt \
       -o output_plot.png \
       -S 0.5 -W 7 -H 10 -f 8 

#2024-07-26 21:38:43,340 - Plotsr - INFO - Starting



######################################

# Using SyRI to identify genomic rearrangements from whole-genome alignments generated using MUMmer

tmux new-session -s syri
tmux attach-session -t syri

salloc -c40 --time 2:50:00 --mem 191000M --account def-rieseber


cd /home/celphin/scratch/Oxyria/GeneSpace/Oxyria_genomes/genomes/syri
module load StdEnv/2020 python/3.9 scipy-stack/2022a
module load minimap2/2.24
source ~/syri/bin/activate

module load StdEnv/2020
module load python/3.10.2
module load minimap2/2.24
module load scipy-stack/2022a
module load gcc/9.3.0
module load r-bundle-bioconductor/3.12 # includes plotly
module load seqtk/1.3
module load StdEnv/2020
module load mummer/4.0.0beta2


nucmer --maxmatch -c 100 -b 500 -l 50 Oxyria_Main_short.fasta Oxy_Elles_Hap1_short.fasta       # Whole genome alignment. Any other alignment can also be used.
delta-filter -m -i 90 -l 100 out.delta > out.filtered.delta     # Remove small and lower quality alignments
show-coords -THrd out.filtered.delta > out.filtered.coords      # Convert alignment information to a .TSV format as required by SyRI
syri -c out.filtered.coords -d out.filtered.delta -r Oxyria_Main_short.fasta -q Oxy_Elles_Hap1_short.fasta
plotsr syri.out Oxyria_Main_short.fasta Oxy_Elles_Hap1_short.fasta -H 8 -W 5

