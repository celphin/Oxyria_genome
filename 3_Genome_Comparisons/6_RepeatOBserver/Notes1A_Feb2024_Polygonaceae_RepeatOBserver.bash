##########################
# Running RepeatOBserver on Polygonaceae genomes on beluga 
# Feb 2024
############################

cd /home/celphin/scratch/Oxyria/RepeatOBserver/

# download the script
wget https://raw.githubusercontent.com/celphin/RepeatOBserverV1/main/Setup_Run_Repeats.sh
chmod +x Setup_Run_Repeats.sh
dos2unix Setup_Run_Repeats.sh

#-----------------------------
# copy over the genomes into one folder

cp F_escelentum_H1.fasta /home/celphin/scratch/Oxyria/RepeatOBserver/
# start the script
cat << EOF > Auto_Fagoesc_H1.sh
#!/bin/bash
#SBATCH --account=def-cronk
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=191000M

module load StdEnv/2020
module load seqkit/2.3.1
module load emboss/6.6.0
module load r/4.1.2

srun Setup_Run_Repeats.sh -i Fagoesc -f F_escelentum_H1.fasta -h H1 -c 19 -m 191000M -g FALSE

EOF

sbatch Auto_Fagoesc_H1.sh

#------------------------- 

cp F_escelentum_H2.fasta /home/celphin/scratch/Oxyria/RepeatOBserver/
# start the script
cat << EOF > Auto_Fagoesc_H2.sh
#!/bin/bash
#SBATCH --account=def-cronk
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=191000M

module load StdEnv/2020
module load seqkit/2.3.1
module load emboss/6.6.0
module load r/4.1.2

srun Setup_Run_Repeats.sh -i Fagoesc -f F_escelentum_H2.fasta -h H2 -c 19 -m 191000M -g FALSE

EOF

sbatch Auto_Fagoesc_H2.sh

#------------------------- 
cp Fagopyrum_tataricum_Main.fasta /home/celphin/scratch/Oxyria/RepeatOBserver/
# start the script
cat << EOF > Auto_Fagotat_H0.sh
#!/bin/bash
#SBATCH --account=def-cronk
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=191000M

module load StdEnv/2020
module load seqkit/2.3.1
module load emboss/6.6.0
module load r/4.1.2

srun Setup_Run_Repeats.sh -i Fagotat -f Fagopyrum_tataricum_Main.fasta -h H0 -c 19 -m 191000M -g FALSE

EOF

sbatch Auto_Fagotat_H0.sh

#------------------------- 
cp F_tataricum_H1.fasta /home/celphin/scratch/Oxyria/RepeatOBserver/
# start the script
cat << EOF > Auto_Fagotat_H1.sh
#!/bin/bash
#SBATCH --account=def-cronk
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=191000M

module load StdEnv/2020
module load seqkit/2.3.1
module load emboss/6.6.0
module load r/4.1.2

srun Setup_Run_Repeats.sh -i Fagotat -f F_tataricum_H1.fasta -h H1 -c 19 -m 191000M -g FALSE

EOF

sbatch Auto_Fagotat_H1.sh

#------------------------- 
cp F_tataricum_H2.fasta /home/celphin/scratch/Oxyria/RepeatOBserver/

# start the script
cat << EOF > Auto_Fagotat_H2.sh
#!/bin/bash
#SBATCH --account=def-cronk
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=191000M

module load StdEnv/2020
module load seqkit/2.3.1
module load emboss/6.6.0
module load r/4.1.2

srun Setup_Run_Repeats.sh -i Fagotat -f F_tataricum_H2.fasta -h H2 -c 19 -m 191000M -g FALSE

EOF

sbatch Auto_Fagotat_H2.sh

#------------------------- 
cp Oxyria_Main.fasta /home/celphin/scratch/Oxyria/RepeatOBserver/
# start the script
cat << EOF > Auto_Oxydig_H0.sh
#!/bin/bash
#SBATCH --account=def-cronk
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=191000M

module load StdEnv/2020
module load seqkit/2.3.1
module load emboss/6.6.0
module load r/4.1.2

srun Setup_Run_Repeats.sh -i Oxydig -f Oxyria_Main.fasta -h H0 -c 19 -m 191000M -g FALSE

EOF

sbatch Auto_Oxydig_H0.sh

#------------------------- 
cp Oxyria_ragtag.fasta /home/celphin/scratch/Oxyria/RepeatOBserver/

# start the script
cat << EOF > Auto_Oxydig_H1.sh
#!/bin/bash
#SBATCH --account=def-cronk
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=191000M

module load StdEnv/2020
module load seqkit/2.3.1
module load emboss/6.6.0
module load r/4.1.2

srun Setup_Run_Repeats.sh -i Oxydig -f Oxyria_ragtag.fasta -h H1 -c 19 -m 191000M -g FALSE

EOF

sbatch Auto_Oxydig_H1.sh

#-------------------------
cp Polavi_Main.fasta /home/celphin/scratch/Oxyria/RepeatOBserver/

# start the script
cat << EOF > Auto_Polavi_H0.sh
#!/bin/bash
#SBATCH --account=def-cronk
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=191000M

module load StdEnv/2020
module load seqkit/2.3.1
module load emboss/6.6.0
module load r/4.1.2

srun Setup_Run_Repeats.sh -i Polavi -f Polavi_Main.fasta -h H0 -c 19 -m 191000M -g FALSE

EOF

sbatch Auto_Polavi_H0.sh

#------------------------- 
cp R_nobile.fasta /home/celphin/scratch/Oxyria/RepeatOBserver/

# start the script
cat << EOF > Auto_Rhunob_H0.sh
#!/bin/bash
#SBATCH --account=def-cronk
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=191000M

module load StdEnv/2020
module load seqkit/2.3.1
module load emboss/6.6.0
module load r/4.1.2

srun Setup_Run_Repeats.sh -i Rhunob -f R_nobile.fasta -h H0 -c 19 -m 191000M -g FALSE

EOF

sbatch Auto_Rhunob_H0.sh

#-------------------------
cp R_tangaticum.fasta /home/celphin/scratch/Oxyria/RepeatOBserver/

# start the script
cat << EOF > Auto_Rhutan_H0.sh
#!/bin/bash
#SBATCH --account=def-cronk
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=191000M

module load StdEnv/2020
module load seqkit/2.3.1
module load emboss/6.6.0
module load r/4.1.2

srun Setup_Run_Repeats.sh -i Rhutan -f R_tangaticum.fasta -h H0 -c 19 -m 191000M -g FALSE

EOF

sbatch Auto_Rhutan_H0.sh

#------------------------- 

cp R_hastalus_linkage.fasta /home/celphin/scratch/Oxyria/RepeatOBserver/

# start the script
cat << EOF > Auto_RumHaus_H0.sh
#!/bin/bash
#SBATCH --account=def-cronk
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=191000M

module load StdEnv/2020
module load seqkit/2.3.1
module load emboss/6.6.0
module load r/4.1.2

srun Setup_Run_Repeats.sh -i RumHaus_H0 -f R_hastalus_linkage.fasta -h H0 -c 19 -m 191000M -g FALSE

EOF

sbatch Auto_RumHaus_H0.sh

#------------------------- 

sq

          JOBID     USER      ACCOUNT           NAME  ST  TIME_LEFT NODES CPUS TRES_PER_N MIN_MEM NODELIST (REASON)
       44910124  celphin def-cronk_cp Auto_Oxydig_H0  PD   10:00:00     1   20        N/A 191000M  (Priority)
       44910135  celphin def-cronk_cp Auto_Oxydig_H1  PD   10:00:00     1   20        N/A 191000M  (Priority)
       44910136  celphin def-cronk_cp Auto_Polavi_H0  PD   10:00:00     1   20        N/A 191000M  (Priority)
       44910137  celphin def-cronk_cp Auto_Rhunob_H0  PD   10:00:00     1   20        N/A 191000M  (Priority)
       44910138  celphin def-cronk_cp Auto_Rhutan_H0  PD   10:00:00     1   20        N/A 191000M  (Priority)
       44909990  celphin def-cronk_cp Auto_Fagoesc_H  PD   10:00:00     1   20        N/A 191000M  (Priority)
       44909991  celphin def-cronk_cp Auto_Fagoesc_H  PD   10:00:00     1   20        N/A 191000M  (Priority)
       44909992  celphin def-cronk_cp Auto_Fagotat_H  PD   10:00:00     1   20        N/A 191000M  (Priority)
       44909993  celphin def-cronk_cp Auto_Fagotat_H  PD   10:00:00     1   20        N/A 191000M  (Priority)
       44909994  celphin def-cronk_cp Auto_Fagotat_H  PD   10:00:00     1   20        N/A 191000M  (Priority)
       44926857  celphin def-cronk_cp Auto_RumHaus_H   R    9:59:17     1   20        N/A 191000M bc11606 (None)


#################################################
# Genomes info

# Species	Chromosome count	N50 (Mbp)	Genome size (Mbp)
#-----------------------------------------------------
#1 Fagopyrum escelentum(H1)	x=8	156(9.8 in paper)	1263 (in paper 1230mbp)
#2 Fagopyrum escelentum(H2)	x=8	150(12.4 in paper)	1190mbp
#3 Fagopyrum tataricum (NCBI)	x=8	53.8	506 mbp
#4 Fagopyrum tataricum (H1)	x=8	56 (50 according to paper)	525 mpb (In paper 453.7 mbp) 
#5 Fagopyrum tataricum (H2)	x=8	56 (30 according to paper)	484 mpb (in paper:446.2 mbp)
#6 Oxyria digyna (our version)	x=7	76	590
#7 Oxyria digyna (version in NCBI)	x=7	36.87586	"	561.2 Mb"
#8 Rheum nobile	x=11	137.450686	1503(in paper 1570 Mb)
#9 Rheum tangaticum	x=?	273(in paper 7.16)	2710(in paper 2740 Mbp)
#10 Polygunum aviculare	x=10	33.8	352

#11 Rumex hastalus	x=5 (2n=10)	158	1647 mbp


#--------------------------------------
# copy to computer
# https://stackoverflow.com/questions/15121337/recursively-use-scp-but-excluding-some-folders

scp -r celphin@beluga.computecanada.ca:/home/celphin/scratch/Oxyria/RepeatOBserver/output_chromosomes/*/Summary_output/DNAwalks* .
scp -r celphin@beluga.computecanada.ca:/home/celphin/scratch/Oxyria/RepeatOBserver/output_chromosomes/*/Summary_output/histograms* .
scp -r celphin@beluga.computecanada.ca:/home/celphin/scratch/Oxyria/RepeatOBserver/output_chromosomes/*/Summary_output/Shannon_div* .

scp -r celphin@beluga.computecanada.ca:/home/celphin/scratch/Oxyria/RepeatOBserver/output_chromosomes/*/Summary_output/spectra/spectra_total_merged* .

scp -r celphin@beluga.computecanada.ca:/home/celphin/scratch/Oxyria/RepeatOBserver/output_chromosomes/*/Summary_output/spectra/spectra_parts_15-35* .
scp -r celphin@beluga.computecanada.ca:/home/celphin/scratch/Oxyria/RepeatOBserver/output_chromosomes/*/Summary_output/spectra/spectra_parts_35-2000* .
scp -r celphin@beluga.computecanada.ca:/home/celphin/scratch/Oxyria/RepeatOBserver/output_chromosomes/*/Summary_output/spectra/spectra_parts_2-8* .

###################################
# Run summary plots of chromosomes

cp -v -u -r /home/celphin/scratch/Oxyria/RepeatOBserver/output_chromosomes/*/Summary_output/pdfs/*.pdf /home/celphin/scratch/Oxyria/RepeatOBserver/Summary_pdfs

cp -v -u -r /home/celphin/scratch/Oxyria/RepeatOBserver/output_chromosomes/*/Summary_output/pdfs/Shannon_div/*Shannon_500.pdf /home/celphin/scratch/Oxyria/RepeatOBserver/Summary_Shannon

cp -v -u -r /home/celphin/scratch/Oxyria/RepeatOBserver/output_chromosomes/*/Summary_output/Shannon_div/Shannon_div_2.5Mbp/*.png /home/celphin/scratch/Oxyria/RepeatOBserver/Summary_Shannon

cp -v -u -r /home/celphin/scratch/Oxyria/RepeatOBserver/output_chromosomes/*/Summary_output/spectra/spectra_total_merged/*TRUE.png /home/celphin/scratch/Oxyria/RepeatOBserver/Summary_spectra


#####################################
# look at repeats

# To run this script EMBOSS6.6.0 is required. 
module load nixpkgs/16.09  intel/2018.3  emboss/6.6.0

# download the script
wget https://raw.githubusercontent.com/celphin/RepeatOBserverV1/main/repeat_seq_finder.sh
chmod +x repeat_seq_finder.sh
dos2unix repeat_seq_finder.sh

# uncomment seqkit load
nano repeat_seq_finder.sh

# direct the script to your chromosome files: "/home/celphin/scratch/repeats/auto_script/input_chromosomes/NerLuet_H0-AT/chromosome_files"
# Choose the Chromosome fasta file that you want to look at the repeat in: "NerLuet_H0-AT_Chr1part01.fasta"
# Define the bp range: 13200000 13500000
# Define the repeat length you are looking for: 4 (will work for repeats up to 600bp long)
# Define the number of bp per line in your fasta file: 60 (default for RepeatOBserverV1)

./repeat_seq_finder.sh "/home/celphin/scratch/repeats/auto_script/input_chromosomes/NerLuet_H0-AT/chromosome_files" \
"NerLuet_H0-AT_Chr1part01.fasta" 10000000 10500000 90 60 500


