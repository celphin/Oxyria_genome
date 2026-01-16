#########################
# Prep orthogroups data for CODEML or RELAX
# Jan 2026
############################

# follow instructions here
# https://github.com/siribi/CODEML_WORKFLOW
# https://github.com/siribi/RELAX-workflow

# other tutorial: https://bioinformaticsworkbook.org/dataAnalysis/ComparativeGenomics/Finding_Positive_Selection_With_Codeml.html#gsc.tab=0

###############################
# To run either program we need:

# Make a giant fasta file containing all proteins from every species
# Make a giant fasta file containing all transcripts from every species

# Alignments of genes (MAFFT) and then codon based alignment (pal2nal)

# A set of files with gene identifiers made from  the Orthogroups.txt file

# Extract the newick tree and gene trees
	# https://timetree.org/ 
	# OR use STAG tree

# Make CODEML control files

######################

# Load modules
module load cufflinks/2.2.1 
module load cdbfasta/2017-03-16
module load orthofinder/2.5.2 
module load diamond/2.0.4 
module load clustalo/1.2.4
module load pal2nal.v14 
module load paml/4.9h

# After running orthofinder
# Make a giant fasta file containing all proteins from every species
cat 01_CleanProteins/*fasta >AllCleanProteins.fasta
ml cdbfasta; cdbfasta AllCleanProteins.fasta

# Make a giant fasta file containing all transcripts from every species
for f in *_Transcripts.fasta; do awk '{print $1}' $f >Clean${f};done
ml cdbfasta; cdbfasta AllCleanTranscripts.fasta

# Extract Extract protein and transcript names by looping through the tree files and formatting to get gene names.
for f in 01_CleanProteins/Results_Feb09/Orthologues_Feb09/Gene_Trees/*tree.txt; do sed 's/(//g' $f |sed 's/,/\n/g' |sed 's/:.*//g' |sort|uniq|sed 's/Proteins_/\t/1'|cut -f 2  |cdbyank AllCleanProteins.fasta.cidx >${f}_Proteins.fasta;done
for f in 01_CleanProteins/Results_Feb09/Orthologues_Feb09/Gene_Trees/*; do sed 's/(//g' $f |sed 's/,/\n/g' |sed 's/:.*//g' |sort|uniq|sed 's/Proteins_/\t/1'|cut -f 2  |cdbyank AllCleanTranscripts.fasta.cidx >${f}_Transcripts.fasta;done

# Generate protein alignments with clustalo
# OR MAFFT

20_TestDataset/01_CleanProteins/Results_Feb09/Orthologues_Feb09/Gene_Trees/

ml clustalo/1.2.4-iy2mf6y
for f in *_Proteins.fasta; do echo "ml clustalo;clustalo -i "$f" -o "${f%.*}"alignment";done >generateAlignment.sh


# Run pal2nal – convert alignments into paml format
# The suggestion is to use pal2nal with “-nogap”, though codeml has a way to deal with these.
# However, if you do not have sequence in your pal2nal output, likely you had gaps
# along the entire sequence across all of the species analyzed.

20_TestDataset/01_CleanProteins/Results_Feb09/Orthologues_Feb09/Gene_Trees/

wget http://www.bork.embl.de/pal2nal/distribution/pal2nal.v14.tar.gz
tar -zxvf pal2nal.v14.tar.gz

for f in *alignment; \
do echo "pal2nal.v14/pal2nal.pl  "$f" "${f%_*}"_Transcripts.fasta \
-output paml -nogap >"${f%_*}"_pal2nal" ; done >pal2nal.sh














