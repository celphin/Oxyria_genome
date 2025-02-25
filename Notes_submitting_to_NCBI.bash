###############################
# Submitting the Oxyria digyna haplotypes to NCBI
# Oct 2024
# Instructions: https://www.ncbi.nlm.nih.gov/genbank/diploid_haps/
###############################

#---------------------
# Make a genome submission
# https://submit.ncbi.nlm.nih.gov/subs/genome/

# Login with ORCID ID or another account

# preloading files - same as uploading large individual GBS or WGS files to the sequnece read archive
# https://www.ncbi.nlm.nih.gov/genbank/preloadfiles/

# Request a preload folder once logged in
# Download the key file from the website above. Please do not share this key file. 
# Do not include this information or key file on a public page. This information is regularly updated as per our security policies.
# transfer the aspera key to compute canada through globus or scp 
more /home/celphin/scratch/ aspera.openssh

#-----------
# Download and install the Aspera Connect software from IBM
# https://docs.alliancecan.ca/wiki/Transferring_files_with_Aspera_Connect/ascp
# Find updated download link here: https://www.ibm.com/products/aspera/downloads#connect
cd ~
wget https://d3gcli72yxqn2z.cloudfront.net/downloads/connect/latest/bin/ibm-aspera-connect_4.2.12.780_linux_x86_64.tar.gz
tar -zxvf ibm-aspera-connect_4.2.12.780_linux_x86_64.tar.gz

bash ibm-aspera-connect_4.2.12.780_linux_x86_64.sh

chmod u+x ~/.aspera/connect/plugins/*/*.so ~/.aspera/connect/lib/*
setrpaths.sh --path $HOME/.aspera 
export PATH=~/.aspera/connect/bin:$PATH

#---------------------
# copy genomes and any other files we want to upload to a specific  folder on the server
mkdir /home/celphin/scratch/Oxyria_genome_upload/
cd /home/celphin/scratch/Oxyria_genome_upload/

cp /home/celphin/scratch/Oxyria/Genome_comparisons_Oxyria/Oxy_draft_assembly_Hap1.FINAL.fasta /home/celphin/scratch/Oxyria_genome_upload/
cp /home/celphin/scratch/Oxyria/Genome_comparisons_Oxyria/Oxy_draft_assembly_Hap2.FINAL.fasta /home/celphin/scratch/Oxyria_genome_upload/


###########################
# FASTA file formatting (see more details below)
# https://www.ncbi.nlm.nih.gov/genbank/diploid_haps/#submit
# https://www.ncbi.nlm.nih.gov/genbank/genomesubmit/#fasta

# Chromosome labels 
grep ">" /home/celphin/scratch/Oxyria/Genome_comparisons_Oxyria/Oxy_draft_assembly_Hap1.FINAL.fasta | head
# >HiC_scaffold_1
# >HiC_scaffold_2
# >HiC_scaffold_3
# >HiC_scaffold_4
# >HiC_scaffold_5
# >HiC_scaffold_6
# >HiC_scaffold_7
# >HiC_scaffold_8
# >HiC_scaffold_9
# >HiC_scaffold_10

# Want to label chromosomes as:
# >OxydigH1_Chr1 [organism=Oxyria digyna] [location=chromosome] [chromosome=1]
# >OxydigH1_Chr2 [organism=Oxyria digyna] [location=chromosome] [chromosome=2]
# >OxydigH1_Chr3 [organism=Oxyria digyna] [location=chromosome] [chromosome=3]
# >OxydigH1_Chr4 [organism=Oxyria digyna] [location=chromosome] [chromosome=4]
# >OxydigH1_Chr5 [organism=Oxyria digyna] [location=chromosome] [chromosome=5]
# >OxydigH1_Chr6 [organism=Oxyria digyna] [location=chromosome] [chromosome=6]
# >OxydigH1_Chr7 [organism=Oxyria digyna] [location=chromosome] [chromosome=7]
# >OxydigH1_scaffold_8 [organism=Oxyria digyna] [location=chromosome]
# >OxydigH1_scaffold_9 [organism=Oxyria digyna] [location=chromosome]
# >OxydigH1_scaffold_10 [organism=Oxyria digyna] [location=chromosome]

# >OxydigH2_Chr1 [organism=Oxyria digyna] [location=chromosome] [chromosome=1]
# >OxydigH2_Chr2 [organism=Oxyria digyna] [location=chromosome] [chromosome=2]
# >OxydigH2_Chr3 [organism=Oxyria digyna] [location=chromosome] [chromosome=3]
# >OxydigH2_Chr4 [organism=Oxyria digyna] [location=chromosome] [chromosome=4]
# >OxydigH2_Chr5 [organism=Oxyria digyna] [location=chromosome] [chromosome=5]
# >OxydigH2_Chr6 [organism=Oxyria digyna] [location=chromosome] [chromosome=6]
# >OxydigH2_Chr7 [organism=Oxyria digyna] [location=chromosome] [chromosome=7]
# >OxydigH2_scaffold_8 [organism=Oxyria digyna] [location=chromosome]
# >OxydigH2_scaffold_9 [organism=Oxyria digyna] [location=chromosome]
# >OxydigH2_scaffold_10 [organism=Oxyria digyna] [location=chromosome]

# see current headers
grep ">" Oxy_draft_assembly_Hap1.FINAL.fasta |head
grep ">" Oxy_draft_assembly_Hap2.FINAL.fasta |head

# rename chromosomes
sed -i 's/HiC_scaffold_1$/OxydigH1_Chr1 \[organism=Oxyria digyna\] \[location\=chromosome\] \[chromosome\=1\]/g' Oxy_draft_assembly_Hap1.FINAL.fasta
sed -i 's/HiC_scaffold_2$/OxydigH1_Chr2 \[organism=Oxyria digyna\] \[location\=chromosome\] \[chromosome\=2\]/g' Oxy_draft_assembly_Hap1.FINAL.fasta
sed -i 's/HiC_scaffold_3$/OxydigH1_Chr3 \[organism=Oxyria digyna\] \[location\=chromosome\] \[chromosome\=3\]/g' Oxy_draft_assembly_Hap1.FINAL.fasta
sed -i 's/HiC_scaffold_4$/OxydigH1_Chr4 \[organism=Oxyria digyna\] \[location\=chromosome\] \[chromosome\=4\]/g' Oxy_draft_assembly_Hap1.FINAL.fasta
sed -i 's/HiC_scaffold_5$/OxydigH1_Chr5 \[organism=Oxyria digyna\] \[location\=chromosome\] \[chromosome\=5\]/g' Oxy_draft_assembly_Hap1.FINAL.fasta
sed -i 's/HiC_scaffold_6$/OxydigH1_Chr6 \[organism=Oxyria digyna\] \[location\=chromosome\] \[chromosome\=6\]/g' Oxy_draft_assembly_Hap1.FINAL.fasta
sed -i 's/HiC_scaffold_7$/OxydigH1_Chr7 \[organism=Oxyria digyna\] \[location\=chromosome\] \[chromosome\=7\]/g' Oxy_draft_assembly_Hap1.FINAL.fasta

sed -i 's/HiC_scaffold_1$/OxydigH2_Chr1 \[organism=Oxyria digyna\] \[location\=chromosome\] \[chromosome\=1\]/g' Oxy_draft_assembly_Hap2.FINAL.fasta
sed -i 's/HiC_scaffold_2$/OxydigH2_Chr2 \[organism=Oxyria digyna\] \[location\=chromosome\] \[chromosome\=2\]/g' Oxy_draft_assembly_Hap2.FINAL.fasta
sed -i 's/HiC_scaffold_3$/OxydigH2_Chr3 \[organism=Oxyria digyna\] \[location\=chromosome\] \[chromosome\=3\]/g' Oxy_draft_assembly_Hap2.FINAL.fasta
sed -i 's/HiC_scaffold_4$/OxydigH2_Chr4 \[organism=Oxyria digyna\] \[location\=chromosome\] \[chromosome\=4\]/g' Oxy_draft_assembly_Hap2.FINAL.fasta
sed -i 's/HiC_scaffold_5$/OxydigH2_Chr5 \[organism=Oxyria digyna\] \[location\=chromosome\] \[chromosome\=5\]/g' Oxy_draft_assembly_Hap2.FINAL.fasta
sed -i 's/HiC_scaffold_6$/OxydigH2_Chr6 \[organism=Oxyria digyna\] \[location\=chromosome\] \[chromosome\=6\]/g' Oxy_draft_assembly_Hap2.FINAL.fasta
sed -i 's/HiC_scaffold_7$/OxydigH2_Chr7 \[organism=Oxyria digyna\] \[location\=chromosome\] \[chromosome\=7\]/g' Oxy_draft_assembly_Hap2.FINAL.fasta

# rename scaffolds
sed -i 's/HiC_scaffold_\(.*$\)/OxydigH1_scaffold_\1 \[organism=Oxyria digyna\] \[location=chromosome\]/g' Oxy_draft_assembly_Hap1.FINAL.fasta
sed -i 's/HiC_scaffold_\(.*$\)/OxydigH2_scaffold_\1 \[organism=Oxyria digyna\] \[location=chromosome\]/g' Oxy_draft_assembly_Hap2.FINAL.fasta

# check again
grep ">" Oxy_draft_assembly_Hap1.FINAL.fasta |head
grep ">" Oxy_draft_assembly_Hap2.FINAL.fasta |head

grep "location" Oxy_draft_assembly_Hap1.FINAL.fasta
grep "location" Oxy_draft_assembly_Hap2.FINAL.fasta

#------------------------------------
# make upload session to not be interupted
tmux new-session -s Genome
tmux attach-session -t Genome

# run aspera
# You may use the following command to upload files via Aspera Command-Line:
# ascp -i <path/to/key_file> -QT -l100m -k1 -d <path/to/folder/containing files> subasp@upload.ncbi.nlm.nih.gov:uploads/cassandra.elphinstone_alumni.ubc.ca_dWEjZQuq
export PATH=~/.aspera/connect/bin:$PATH

ascp -i /home/celphin/scratch/aspera.openssh -QT -l100m -k1 -d /home/celphin/scratch/Oxyria_genome_upload subasp@upload.ncbi.nlm.nih.gov:uploads/cassandra.elphinstone_alumni.ubc.ca_dWEjZQuq/Oxyria_genome_upload


############################
# Once files are uploading 

# Start new submission

# Select option 1 (non-wgs genomes)

# Template files can be found here: https://submit.ncbi.nlm.nih.gov/templates/
# or in the attached github directory

# Fill out a Biosample form for each individual or multiple individuals
# includes lat long information, culivar names, etc.

# Fill out a Genome submission form for each indiviudal and its haplotypes
# includes sequencing technologies and coverage and filenames you are uploading

# Connect the preloaded folder to your submission

# Submit the genome

###################################
# After submission remove contamination
# Removing contamination

# SUBID     	BioProject	BioSample	Organism
# --------------------------------------------------------
# SUB14768376	PRJNA1175703	SAMN44369891	Oxyria digyna

# [] We ran your sequences through our Contamination Screen. The screen found 
# contigs that need to be trimmed and/or excluded. The results are in the 
# Contamination.txt file posted in your submission on the WGS submission portal 
# https://submit.ncbi.nlm.nih.gov/subs/genome/. More details about the 
# contamination screening process are available at https://github.com/ncbi/fcs

# GenBank staff will automatically remove contaminants that are found to be 
# the entire sequence or at the end of a sequence, and will post the reports 
# and edited fasta file to the submission portal. Note that internal contamination 
# will not be automatically removed since the sequence may be misassembled and 
# therefore should be split at the contamination and resubmitted as separate sequences.
# In addition, we do not automatically remove mitochondrial sequences in 
# eukaryotic submissions. 

# If you selected the submission portal option "Do not automatically trim or 
# remove sequences identified as contamination" then you will need 
# to adjust the sequences appropriately and then resubmit your sequences. 
# After you remove the contamination, trim any Ns at the ends of the sequence 
# and remove any sequences that are shorter than 200 nt and not part of a 
# multi-component scaffold.

# WARNING: If we do not hear from you by $(add.days,14), your 
# submission will be deleted from the processing queue.

# Note that mismatches between the name of the adaptor/primer identified in the screen 
# and the sequencing technology used to generate the sequencing data should not be used 
# to discount the validity of the screen results as the adaptors/primers of many 
# different sequencing platforms share sequence similarity.


# Screened 886 sequences, 566,594,759 bp.
# 1 sequence to exclude

# Exclude:
# Sequence name, length, apparent source
# OxydigH2_scaffold_14	3894583	prok:firmicutes

wget https://submit.ncbi.nlm.nih.gov/api/2.0/files/of9nuamy/remainingcontamination_oxy_draft_assembly_hap1_final.txt/?format=inline

mv 'index.html?format=inline' Oxyria_H1_contamination

#-----------------------------------------------
# SUBID     	BioProject	BioSample	Organism
# --------------------------------------------------------
# SUB14768376	PRJNA1175704	SAMN44369891	Oxyria digyna

# [] We ran your sequences through our Contamination Screen. The screen found 
# contigs that need to be trimmed and/or excluded. The results are in the 
# Contamination.txt file posted in your submission on the WGS submission portal 
# https://submit.ncbi.nlm.nih.gov/subs/genome/. More details about the 
# contamination screening process are available at https://github.com/ncbi/fcs

# GenBank staff will automatically remove contaminants that are found to be 
# the entire sequence or at the end of a sequence, and will post the reports 
# and edited fasta file to the submission portal. Note that internal contamination 
# will not be automatically removed since the sequence may be misassembled and 
# therefore should be split at the contamination and resubmitted as separate sequences.
# In addition, we do not automatically remove mitochondrial sequences in 
# eukaryotic submissions. 

# If you selected the submission portal option "Do not automatically trim or 
# remove sequences identified as contamination" then you will need 
# to adjust the sequences appropriately and then resubmit your sequences. 
# After you remove the contamination, trim any Ns at the ends of the sequence 
# and remove any sequences that are shorter than 200 nt and not part of a 
# multi-component scaffold.

# WARNING: If we do not hear from you by $(add.days,14), your 
# submission will be deleted from the processing queue.

# Note that mismatches between the name of the adaptor/primer identified in the screen 
# and the sequencing technology used to generate the sequencing data should not be used 
# to discount the validity of the screen results as the adaptors/primers of many 
# different sequencing platforms share sequence similarity.


# Screened 2,649 sequences, 582,616,305 bp.
# 1 sequence to exclude

# Exclude:
# Sequence name, length, apparent source
# OxydigH1_scaffold_126	3894583	prok:firmicutes

wget https://submit.ncbi.nlm.nih.gov/api/2.0/files/2hvx4tfy/contamination_oxy_draft_assembly_hap2_final.txt/?format=inline
mv 'index.html?format=inline' Oxyria_H2_contamination.txt


#####################################
# Remove the scaffolds from H1 and H2
# https://jasonbiology.tokyo/2022/01/03/fasta_remove_scaffold_en/comment-page-1/

module load StdEnv/2023 seqkit/2.5.1

# make a copy
cp Oxy_draft_assembly_Hap1.FINAL.fasta Oxy_draft_assembly_Hap1.FINAL_total.fasta
cp Oxy_draft_assembly_Hap2.FINAL.fasta Oxy_draft_assembly_Hap2.FINAL_total.fasta

# check curent list 
grep ">" Oxy_draft_assembly_Hap1.FINAL_total.fasta | head -n 50
grep ">" Oxy_draft_assembly_Hap2.FINAL_total.fasta | head -n 150

# remove OxydigH1_scaffold_126
seqkit grep -nrvip "^OxydigH1_scaffold_126 \[organism\=Oxyria digyna\] \[location\=chromosome\]$" Oxy_draft_assembly_Hap1.FINAL_total.fasta  > Oxy_draft_assembly_Hap1.FINAL.fasta

# remove OxydigH2_scaffold_14
seqkit grep -nrvip "^OxydigH2_scaffold_14 \[organism\=Oxyria digyna\] \[location\=chromosome\]$" Oxy_draft_assembly_Hap2.FINAL_total.fasta  > Oxy_draft_assembly_Hap2.FINAL.fasta

# check after list 
grep ">" Oxy_draft_assembly_Hap1.FINAL.fasta | head -n 150
grep ">" Oxy_draft_assembly_Hap2.FINAL.fasta | head -n 50

mkdir upload
mv Oxy_draft_assembly_Hap*.FINAL.fasta upload/

#---------------
# reupload

# make upload session to not be interupted
tmux new-session -s Genome
tmux attach-session -t Genome

cd /home/celphin/scratch/Oxyria_genome_upload/upload
# move the final versions of the genomes here

# run aspera
# You may use the following command to upload files via Aspera Command-Line:
# ascp -i <path/to/key_file> -QT -l100m -k1 -d <path/to/folder/containing files> subasp@upload.ncbi.nlm.nih.gov:uploads/cassandra.elphinstone_alumni.ubc.ca_dWEjZQuq
export PATH=~/.aspera/connect/bin:$PATH

ascp -i /home/celphin/scratch/aspera.openssh -QT -l100m -k1 -d /home/celphin/scratch/Oxyria_genome_upload/upload subasp@upload.ncbi.nlm.nih.gov:uploads/cassandra.elphinstone_alumni.ubc.ca_dWEjZQuq/Oxyria_genome_no_contam2


###############################
# FASTA file formatting Details:
# https://www.ncbi.nlm.nih.gov/genbank/diploid_haps/#submit
# https://www.ncbi.nlm.nih.gov/genbank/genomesubmit/#fasta

# How to submit

# [1] In the simple case fasta files can be uploaded because there is no annotation and having 
# the same linkage_evidence for all the gaps is acceptable. However, if annotation or different kinds of gaps are included, then you will need to use the command line program table2asn to create an ASN (.sqn) file, as explained in the Genome Submission Guide.

# [2] If any of the sequences belong to chromosomes or organelles and the "Batch" or 
# "Pseudohaplotypes of a diploid/polyploid assembly" submission option is used, then that
# assignment information must be included in the definition lines of the fasta sequences,
# as described at 'IMPORTANT: Additional requirements for batch submissions'. Specifically:

    # Unlocalized organelle sequence, use [location=xxx], eg:

        # [location=mitochondrion]

        # [location=chloroplast]

    # The complete circular organelle sequence, then add the topology and completeness, eg:
        # [location=mitochondrion] [completeness=complete] [topology=circular]

    # Unlocalized sequence that belongs to a chromosome, eg chromosome 2:
        # [chromosome=2]

    # The sequence represents the chromosome, even if gaps may be present, then add the location:
        # [location=chromosome] [chromosome=2]

# [3] If the files are very big, you may want to upload them before you begin your genome 
# submission, as described at https://www.ncbi.nlm.nih.gov/genbank/preloadfiles/. FYI, this 
# is the same process that exists for preloading files for SRA submissions or any genome submission.

# [4] To submit the haplotypes of an individual, start a new submission in the Genomes Submission Portal.

#########################################
# Annotations - did not test
# .sqn files

# These are generally required only when the submitter wants to include annotation.
# Annotation is optional for GenBank genome submissions.

# https://www.ncbi.nlm.nih.gov/genbank/table2asn/

# Prepare a .sqn file for submission using table2asn. table2asn reads a template file along with the fasta sequence and annotation table files, and outputs an ASN (.sqn) file for submission to GenBank. Follow these three steps:

# 1) Prepare data files
# Prepare fasta files as above, with one file per genome.

# Prepare these additional files:
    # a template file with submitter and publication information.
    # annotation files. These correspond to and have the same basenames as the .fsa files. There are two different 
	# file formats:
        # 5-column feature table files that have the suffix .tbl. Be sure to read the annotation requirements in 
		# the appropriate annotation guidelines:
            # Prokaryotic Annotation Guidelines
            # Eukaryotic Annotation Guidelines
        # GFF files in GenBank-specific format that have the suffix .gff. Be sure to read the instructions at 
		# Genome Annotation with GFF or GTF files.
    # Genome-Assembly-Data Structured Comment. This information can be provided during the genome submission, 
	# but if many genomes are being submitted it could be simpler to include this in the .sqn file itself. To 
	# do that, use the Genome-Assembly-Data Structured Comment Template to create the file and then have it 
	# included with -w genasm.cmt in the tbl2asn commandline, below.
    # quality scores of the sequences. These files correspond to and have the same basenames as the .fsa 
	# files, but have the suffix .qvl. The quality scores are optional.

# 2) Run table2asn
# A. Annotation is in GenBank-specific GFF files: follow the instructions for GFF files.
# B. Annotation is in .tbl files: follow these instructions. Note that a few of the arguments in table2asn have 
# changed relative to tbl2asn, eg -indir instead of -p. The table2asn page provides more details. Here are 
# the instructions for creating annotated genome files when the annotation is in .tbl files:

# Sample command line when the sequences are contigs (overlapping reads with no Ns representing gaps) is
table2asn -indir path_to_files -t template -M n -Z

# If the sequences contain Ns that represent gaps, then run the appropriate table2asn command line with the -l 
# and -gaps-min arguments, as described in the Gapped Genome Submission page. The command line for the most common 
# situation (runs of 10 or more Ns represent a gap, and there are no gaps of completely unknown size, and the 
# evidence for linkage across the gaps is "paired-ends") is:
table2asn -indir path_to_files -t template -M n -Z -gaps-min 10 -l paired-ends

# For either case you can include the source information in the definition line of each contig, as described in the 
# fasta defline components section, above. Alternatively, the organism and strain (or breed or isolate) can be included 
# with -j in the table2asn command line. The additional source qualifiers will be obtained from the registered BioSample.
 # However, chromosome, plasmid & organelle assignment information must be included in the fasta definition lines. 
 # In addition, if the submission is an annnotated prokaryotic genome, then include the genetic code with -j in the 
 # commandline, for example:
table2asn -indir path_to_files -t template -M n -Z -j "[organism=Clostridium difficile ABDC] [strain=ABDC] [gcode=11]"