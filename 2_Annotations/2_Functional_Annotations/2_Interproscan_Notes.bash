########################################################
# Running Interproscan on Annotations
#####################################

# January 2024

#1. Combine protein files from maker for:
    #Oxyria digyna
    #Oxyria digyna NCBI version
    #Fagopyrum tataricum
    #Polygonum aviculare
#2. Replace character (*) with X
#3. Run interproscan:
        #Oxyria 
        #Fagopyrum tataricum
        #Polygonum aviculare
        #Rheum tanguticum
        #Rheum nobile
        #Fagopyrum escelentum
        #Fagopurum tataricum (another version)
#Requires: run_interproscan.sh
#Note: needs to be run in cedar (connection online look ups doesn't work in compute nodes of graham)
########################################################
#In Cedar:
cd ~/scratch/Oxyria_interproscan
mkdir Interproscan
cd Interproscan
##########################################################
#Combine protein files:
final_fasta="Oxyria_proteins_AED0.6.fasta"
for i in {1..8}; do
    current_file="Oxyria_Proteins/predicted_proteins_MAKER_ORF_Filtered_Oxyria.scaffold_${i}.AED_0.6.fasta"
    cat "$current_file" >> "$final_fasta"
done
grep ">" $final_fasta |wc -l

final_fasta="Fagopyrum_proteins_AED0.6.fasta"
for i in {1..9}; do
    current_file="Fagopyrum_Proteins/predicted_proteins_MAKER_ORF_Filtered_Fagopyrum.scaffold_${i}.AED_0.6.fasta"
    cat "$current_file" >> "$final_fasta"
done
grep ">" $final_fasta |wc -l

final_fasta="Polavi_proteins_AED0.6.fasta"
for i in {1..11}; do
    current_file="Polavi_Proteins/predicted_proteins_MAKER_ORF_Filtered_Polavi.scaffold_${i}.AED_0.6.fasta"
    cat "$current_file" >> "$final_fasta"
done
grep ">" $final_fasta |wc -l

final_fasta="Oxyria_NCBI_proteins_AED0.6.fasta"
for i in {1..8}; do
    current_file="predicted_proteins_MAKER_ORF_Filtered_Oxyria_NCBI.scaffold_${i}.AED_0.6.fasta"
    cat "$current_file" >> "$final_fasta"
done
grep ">" $final_fasta |wc -l

######################################################
#Make sure no "*" characters:
#Oxyria
awk '{ gsub(/\*/, "X"); print }' Oxyria_proteins_AED0.6.fasta > Oxyria_proteins_AED0.6_interproscan_input.fasta
awk '{ gsub(/\*/, "X"); print }' Fagopyrum_proteins_AED0.6.fasta > Fagopyrum_proteins_AED0.6_interproscan_input.fasta
awk '{ gsub(/\*/, "X"); print }' Polavi_proteins_AED0.6.fasta > Polavi_proteins_AED0.6_interproscan_input.fasta
awk '{ gsub(/\*/, "X"); print }' R_nobile_proteins.fasta > R_nobile_interproscan_input.fasta
awk '{ gsub(/\*/, "X"); print }' Oxyria_NCBI_proteins_AED0.6.fasta > Oxyria_NCBI_proteins_AED0.6_interproscan_input.fasta

####################################################
#Run Interproscan
sbatch run_interproscan.sh Oxyria_proteins_AED0.6_interproscan_input.fasta Oxyria_AED0.6_interproscan.tsv
sbatch run_interproscan.sh Fagopyrum_proteins_AED0.6_interproscan_input.fasta Fagopyrum_AED0.6_interproscan.tsv
sbatch run_interproscan.sh Polavi_proteins_AED0.6_interproscan_input.fasta Polavi_AED0.6_interproscan.tsv
sbatch run_interproscan.sh Rheum_tangaticum_proteins.fasta Rheum_tangaticum_interproscan_out.tsv
sbatch run_interproscan.sh Oxyria_NCBI_proteins_AED0.6_interproscan_input.fasta Oxyria_NCBI_interproscan_out.tsv
sbatch run_interproscan.sh F_tataricum_H1.proteins.fasta F_tataricum_H1_interproscan_out.tsv 
sbatch run_interproscan.sh F_escelentum_H1.proteins.fasta F_escelentum_H1_interproscan_out.tsv
sbatch run_interproscan.sh R_nobile_interproscan_input.fasta R_nobile_interproscan_output.tsv

####################################################
# July 2024

# Run Interproscan on:

#1. Combine protein files from maker for:
	# Dryas_octopetala
	# Draba_nivalis
	# Arabis_alpina
	# Oxyria_diyna
	# Rheum nobile
	# Cochlearia_groenlandica

#2. Replace character (*) with X

#3. Run interproscan:
		# Dryas_octopetala
		# Draba_nivalis
		# Arabis_alpina
		# Oxyria_diyna H1
		# Rheum nobile
		# Cochlearia_groenlandica

#Requires: run_interproscan.sh
#Note: needs to be run in cedar (connection online look ups doesn't work in compute nodes of graham)

########################################################

#In Cedar: need internet on compute nodes
cd ~/scratch/Oxyria
mkdir Interproscan
cd Interproscan

# https://interproscan-docs.readthedocs.io/en/latest/HowToRun.html

nano run_interproscan.sh
#!/bin/bash
#SBATCH --account=rrg-rieseber-ac_cpu
#SBATCH --job-name=interproscan
#SBATCH --time=1-12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=0

#Required parameters (in order):
    #inputfile name
    #outputfile name
#to run sh run_interproscan.sh inputfile.fasta outputfile.tsv

module load StdEnv/2020
module load interproscan/5.64-96.0

input_file=$1
output_file=$2
srun interproscan.sh -i "$input_file"  -f tsv -cpu 46 -o "$output_file" -dp  --goterms

##########################################################
# Copy over protein files - copy Beluga to Cedar
cp /home/celphin/scratch/Oxyria/GeneSpace/Total_genomes/peptide/*.fa ~/scratch/Oxyria/Interproscan/

gunzip cg.h1.aa.gz
mv cg.h1.aa Cochlearia_groenlandica.fa

######################################################
#Make sure no "*" characters:
#Oxyria
awk '{ gsub(/\*/, "X"); print }' Arabis_alpina.fa > Arabis_alpina_interproscan_input.fasta
awk '{ gsub(/\*/, "X"); print }' Dryas_octopetala.fa > Dryas_octopetala_interproscan_input.fasta
awk '{ gsub(/\*/, "X"); print }' Oxyria_digyna_H1.fa > Oxyria_digyna_H1_interproscan_input.fasta
awk '{ gsub(/\*/, "X"); print }' Rheum_nobile_H0.fa > Rheum_nobile_H0_interproscan_input.fasta
awk '{ gsub(/\*/, "X"); print }' Draba_nivalis.fa > Draba_nivalis_interproscan_input.fasta
awk '{ gsub(/\*/, "X"); print }' Cochlearia_groenlandica.fa > Cochlearia_groenlandica_interproscan_input.fasta

####################################################
#Run Interproscan
sbatch run_interproscan.sh Oxyria_digyna_H1_interproscan_input.fasta Oxyria_digyna_H1_interproscan_output.tsv
sbatch run_interproscan.sh Dryas_octopetala_interproscan_input.fasta Dryas_octopetala_interproscan_output.tsv
sbatch run_interproscan.sh Arabis_alpina_interproscan_input.fasta Arabis_alpina_interproscan_output.tsv
sbatch run_interproscan.sh Rheum_nobile_H0_interproscan_input.fasta Rheum_nobile_H0_interproscan_output.tsv
sbatch run_interproscan.sh Draba_nivalis_interproscan_input.fasta Draba_nivalis_interproscan_output.tsv
sbatch run_interproscan.sh Cochlearia_groenlandica_interproscan_input.fasta Cochlearia_groenlandica_interproscan_output.tsv

# Update the allocations
scontrol update job=37103206 account=rrg-rieseber-ac_cpu
scontrol update job=37103215 account=rrg-rieseber-ac_cpu
scontrol update job=37103217 account=rrg-rieseber-ac_cpu
scontrol update job=37103219 account=rrg-rieseber-ac_cpu
scontrol update job=37103221 account=rrg-rieseber-ac_cpu

sq
         # JOBID     USER      ACCOUNT           NAME  ST  TIME_LEFT NODES CPUS TRES_PER_N MIN_MEM NODELIST (REASON)
       # 37103206  celphin rrg-rieseber   interproscan   R 1-10:31:54     1   48        N/A 187.50G cdr1977 (None)
       # 37103215  celphin rrg-rieseber   interproscan   R 1-10:31:54     1   48        N/A 187.50G cdr1984 (None)
       # 37103217  celphin rrg-rieseber   interproscan   R 1-11:45:05     1   48        N/A 187.50G cdr1959 (None)
       # 37395871  celphin rrg-rieseber   interproscan  PD 1-12:00:00     1   48        N/A       0  (None)
       # 37103219  celphin rrg-rieseber   interproscan  PD 1-12:00:00     1   48        N/A 187.50G  (Priority)
       # 37103221  celphin rrg-rieseber   interproscan  PD 1-12:00:00     1   48        N/A 187.50G  (Priority)
