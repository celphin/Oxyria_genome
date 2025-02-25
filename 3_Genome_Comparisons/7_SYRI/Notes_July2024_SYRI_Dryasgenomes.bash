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
cd /home/celphin/scratch/Oxyria/GeneSpace/Dryas_genomes/syri

cp /home/celphin/scratch/Annotation/CAP_Snakemake/DoctH0/Ref/DoctH0_Main.fasta ./DoctH0_Main.fasta
cp ~/scratch/Dryas/Dryas_genomes/Dry-int_ragtag_output/ragtag.scaffold.fasta  Dry-int.fasta
cp ~/scratch/Dryas/Dryas_genomes/Dry-drumm_ragtag_output/ragtag.scaffold.fasta  Dry-drumm.fasta
cp ~/scratch/Dryas/Dryas_genomes/Dry-alask_ragtag_output/ragtag.scaffold.fasta  Dry-alask.fasta
cp ~/scratch/Dryas/Dryas_genomes/Dry-octo-H1_ragtag_output/ragtag.scaffold.fasta  Dry-octo-H1.fasta
cp ~/scratch/Dryas/Dryas_genomes/Dry-octo-H2_ragtag_output/ragtag.scaffold.fasta  Dry-ajan.fasta

#-------------------------
# remove short scaffolds
module load seqkit/2.5.1
seqkit seq -m 5000000 DoctH0_Main.fasta > DoctH0_Main-chr.fasta
seqkit seq -m 5000000 Dry-int.fasta > Dry-int-chr.fasta
seqkit seq -m 5000000 Dry-drumm.fasta > Dry-drumm-chr.fasta
seqkit seq -m 5000000 Dry-alask.fasta > Dry-alask-chr.fasta
seqkit seq -m 5000000 Dry-octo-H1.fasta > Dry-octo-H1-chr.fasta
seqkit seq -m 5000000 Dry-ajan.fasta > Dry-ajan-chr.fasta

##############################
module load bioawk/1.0
bioawk -c fastx '{ print $name, length($seq) }' < DoctH0_Main-chr.fasta
bioawk -c fastx '{ print $name, length($seq) }' < Dry-int-chr.fasta
bioawk -c fastx '{ print $name, length($seq) }' < Dry-drumm-chr.fasta
bioawk -c fastx '{ print $name, length($seq) }' < Dry-alask-chr.fasta
bioawk -c fastx '{ print $name, length($seq) }' < Dry-octo-H1-chr.fasta
bioawk -c fastx '{ print $name, length($seq) }' < Dry-ajan-chr.fasta

#################################
# Rename contigs of all genomes

awk '/^>/{print ">Dry_int" ++i; next}{print}' < Dry-int-chr.fasta > Dry-int-renamed.fasta
awk '/^>/{print ">Dry_ajan" ++i; next}{print}' < Dry-ajan-chr.fasta > Dry-ajan-renamed.fasta
awk '/^>/{print ">Dry_drumm" ++i; next}{print}' < Dry-drumm-chr.fasta > Dry-drumm-renamed.fasta
awk '/^>/{print ">Dry_alask" ++i; next}{print}' < Dry-alask-chr.fasta > Dry-alask-renamed.fasta
awk '/^>/{print ">Dry_octo" ++i; next}{print}' < Dry-octo-H1-chr.fasta > Dry-octo-H1-renamed.fasta
awk '/^>/{print ">DoctH0_Main" ++i; next}{print}' < DoctH0_Main-chr.fasta > DoctH0_Main-renamed.fasta

##################################
# Minimap2

tmux new-session -s syri
tmux attach-session -t syri

salloc -c40 --time 2:55:00 --mem 191000M --account def-rieseber

cd /home/celphin/scratch/Oxyria/GeneSpace/Dryas_genomes/syri
module load StdEnv/2020 python/3.9 scipy-stack/2022a
module load minimap2/2.24

minimap2 -ax asm5 -t 40 --eqx Dry-int-renamed.fasta Dry-ajan-renamed.fasta > Int_Ajan_out.sam

minimap2 -ax asm5 -t 40 --eqx Dry-ajan-renamed.fasta Dry-drumm-renamed.fasta > Ajan_Drumm_out.sam

minimap2 -ax asm5 -t 40 --eqx Dry-drumm-renamed.fasta Dry-alask-renamed.fasta > Drumm_Alask_out.sam

minimap2 -ax asm5 -t 40 --eqx DoctH0_Main-renamed.fasta Dry-int-renamed.fasta > Main_Int_out.sam

minimap2 -ax asm5 -t 40 --eqx Dry-int-renamed.fasta Dry-drumm-renamed.fasta > Int_Drumm_out.sam

minimap2 -ax asm5 -t 40 --eqx Dry-int-renamed.fasta Dry-alask-renamed.fasta > Int_Alask_out.sam

minimap2 -ax asm5 -t 40 --eqx Dry-int-renamed.fasta Dry-octo-H1-renamed.fasta > Int_Octo_out.sam

minimap2 -ax asm5 -t 40 --eqx Dry-ajan-renamed.fasta Dry-alask-renamed.fasta > Ajan_Alask_out.sam

minimap2 -ax asm5 -t 40 --eqx Dry-ajan-renamed.fasta Dry-octo-H1-renamed.fasta > Ajan_Octo_out.sam

minimap2 -ax asm5 -t 40 --eqx Dry-drumm-renamed.fasta Dry-octo-H1-renamed.fasta > Drumm_Octo_out.sam

minimap2 -ax asm5 -t 40 --eqx Dry-drumm-renamed.fasta DoctH0_Main-renamed.fasta > Drumm_Main_out.sam

minimap2 -ax asm5 -t 40 --eqx Dry-alask-renamed.fasta DoctH0_Main-renamed.fasta > Alask_Main_out.sam

minimap2 -ax asm5 -t 40 --eqx Dry-octo-H1-renamed.fasta DoctH0_Main-renamed.fasta > Octo_Main_out.sam



minimap2 -ax asm5 -t 40 --eqx Dry-alask-renamed.fasta Dry-octo-H1-renamed.fasta > Alask_Octo_out.sam

minimap2 -ax asm5 -t 40 --eqx Dry-octo-H1-renamed.fasta DoctH0_Main-renamed.fasta > Octo_Main_out.sam
####################################
# Run SyRI
# Elles Hap 1
cd /home/celphin/scratch/Oxyria/GeneSpace/Dryas_genomes/syri

tmux new-session -s syri1
tmux attach-session -t syri1

module load StdEnv/2020 python/3.9 scipy-stack/2022a
module load minimap2/2.24
source ~/syri/bin/activate

# Create directories and run syri for each combination
mkdir Int_Ajan_out
mv Int_Ajan_out.sam Int_Ajan_out/
cd Int_Ajan_out
syri -c Int_Ajan_out.sam -r ../Dry-int-renamed.fasta -q ../Dry-ajan-renamed.fasta -k -F S
cd ..

mkdir Int_Drumm_out
mv Int_Drumm_out.sam Int_Drumm_out/
cd Int_Drumm_out
syri -c Int_Drumm_out.sam -r ../Dry-int-renamed.fasta -q ../Dry-drumm-renamed.fasta -k -F S
cd ..

mkdir Int_Alask_out
mv Int_Alask_out.sam Int_Alask_out/
cd Int_Alask_out
syri -c Int_Alask_out.sam -r ../Dry-int-renamed.fasta -q ../Dry-alask-renamed.fasta -k -F S
cd ..

mkdir Int_Octo_out
mv Int_Octo_out.sam Int_Octo_out/
cd Int_Octo_out
syri -c Int_Octo_out.sam -r ../Dry-int-renamed.fasta -q ../Dry-octo-H1-renamed.fasta -k -F S
cd ..

mkdir Main_Int_out
mv Main_Int_out.sam Main_Int_out/
cd Main_Int_out
syri -c Main_Int_out.sam -r ../DoctH0_Main-renamed.fasta -q ../Dry-int-renamed.fasta -k -F S
cd ..

mkdir Ajan_Drumm_out
mv Ajan_Drumm_out.sam Ajan_Drumm_out/
cd Ajan_Drumm_out
syri -c Ajan_Drumm_out.sam -r ../Dry-ajan-renamed.fasta -q ../Dry-drumm-renamed.fasta -k -F S
cd ..

#---------------
mkdir Ajan_Alask_out
mv Ajan_Alask_out.sam Ajan_Alask_out/
cd Ajan_Alask_out
syri -c Ajan_Alask_out.sam -r ../Dry-ajan-renamed.fasta -q ../Dry-alask-renamed.fasta -k -F S
cd ..

mkdir Ajan_Octo_out
mv Ajan_Octo_out.sam Ajan_Octo_out/
cd Ajan_Octo_out
syri -c Ajan_Octo_out.sam -r ../Dry-ajan-renamed.fasta -q ../Dry-octo-H1-renamed.fasta -k -F S
cd ..

mkdir Drumm_Alask_out
mv Drumm_Alask_out.sam Drumm_Alask_out/
cd Drumm_Alask_out
syri -c Drumm_Alask_out.sam -r ../Dry-drumm-renamed.fasta -q ../Dry-alask-renamed.fasta -k -F S
cd ..

mkdir Drumm_Octo_out
mv Drumm_Octo_out.sam Drumm_Octo_out/
cd Drumm_Octo_out
syri -c Drumm_Octo_out.sam -r ../Dry-drumm-renamed.fasta -q ../Dry-octo-H1-renamed.fasta -k -F S
cd ..

mkdir Drumm_Main_out
mv Drumm_Main_out.sam Drumm_Main_out/
cd Drumm_Main_out
syri -c Drumm_Main_out.sam -r ../Dry-drumm-renamed.fasta -q ../DoctH0_Main-renamed.fasta -k -F S
cd ..

mkdir Alask_Main_out
mv Alask_Main_out.sam Alask_Main_out/
cd Alask_Main_out
syri -c Alask_Main_out.sam -r ../Dry-alask-renamed.fasta -q ../DoctH0_Main-renamed.fasta -k -F S
cd ..

mkdir Alask_Octo_out
mv Alask_Octo_out.sam Alask_Octo_out/
cd Alask_Octo_out
syri -c Alask_Octo_out.sam -r ../Dry-alask-renamed.fasta -q ../Dry-octo-H1-renamed.fasta -k -F S
cd ..

mkdir Octo_Main_out
mv Octo_Main_out.sam Octo_Main_out/
cd Octo_Main_out
syri -c Octo_Main_out.sam -r ../Dry-octo-H1-renamed.fasta -q ../DoctH0_Main-renamed.fasta -k -F S
cd ..




#----------
# Output
more syri.summary

# #Variation_type Count   Length_ref      Length_qry
# Syntenicregions        374     451340179       441587762
# Inversions      265     27528133        26560023
# Translocations  195     7177247 6553530
# Duplications(reference)        39      383625  -
# Duplications(query)    278     -       1726194
# Notaligned(reference) 855     53073235        -
# Notaligned(query)     1055    -       34770628
# SNPs    310180  310180  310180
# Insertions      22643   -       1482208
# Deletions       22100   1522759 -
# Copygains       39      -       200282
# Copylosses      50      179949  -
# Highlydiverged 3824    119572015       108246914
# Tandemrepeats  11      4643    3655

#----------------------------
# Join output files

# Create an output file
output_file="syri_summary_combined.txt"
echo -n > $output_file  # Clear the file if it exists

# Process each directory
for dir in Int_Ajan_out Ajan_Drumm_out Drumm_Alask_out Alask_Octo_out Octo_Main_out Main_Int_out Int_Drumm_out Int_Alask_out Int_Octo_out Ajan_Alask_out Ajan_Octo_out Drumm_Octo_out Drumm_Main_out Alask_Main_out Octo_Main_out; do
    summary_file="${dir}/syri.summary"
    
    # Check if the summary file exists
    if [[ -f $summary_file ]]; then
        # Extract the genome name
        genome_name="${dir%%_out}"
        
        # Process the summary file
        {
            awk -v name="$genome_name" '/^#Variation_type/ {print name "\t" $0; next} /^[^#]/ {print name "\t" $1,$2,$3,$4,$5}' $summary_file
            awk -v name="$genome_name" '/^#Variation_type/ {print name "\t" $0; next} /^[^#]/ {print name "\t" $1,$2,$3,$4,$5}' $summary_file
        } >> $output_file
    else
        echo "Warning: $summary_file not found."
    fi
done

# edit manually to remove spaces in row names

######################################
#Plotting
# https://github.com/schneebergerlab/plotsr

cd /home/celphin/scratch/Oxyria/GeneSpace/Oxyria_genomes/genomes/syri/

nano genomes.txt
#file	name	tags
Dry-int-renamed.fasta	Dint	lw:1.5
Dry-ajan-renamed.fasta	Dajan	lw:1.5
Dry-alask-renamed.fasta	Dalask	lw:1.5
DoctH0_Main-renamed.fasta	DoctoH0	lw:1.5
Dry-octo-H1-renamed.fasta	DoctoH1	lw:1.5
Dry-drumm-renamed.fasta	Ddrumm	lw:1.5


plotsr --sr ./Int_Ajan_out/syri.out \
       --sr ./Ajan_Alask_out/syri.out \
       --sr ./Alask_Main_out/syri.out \
       --genomes genomes.txt \
       -o output_plot.png \
       -S 0.5 -W 7 -H 10 -f 8 

#2024-07-26 21:38:43,340 - Plotsr - INFO - Starting

# ImportError: Cannot read annotations for genome: Dry-octo-H1-renamed.fasta. Make sure that structural annotations for all genomes are provided in the same order as genomes. Exiting.


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

