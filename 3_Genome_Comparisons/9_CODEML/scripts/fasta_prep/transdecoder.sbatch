#!/bin/bash
# Job name: TransDecoder
#SBATCH --job-name=TransDecoder
#
# Project:
#SBATCH --account=nn8023k
#
# Wall clock limit:
#SBATCH --time=6:00:00
#
# Number of tasks/processes to run:
#SBATCH --ntasks=1
#
# Number of CPUs (this is what I learned from Saga support if I want to use the -t flag):
#SBATCH --cpus-per-task=1
#
# Max memory usage per core (MB):
#SBATCH --mem-per-cpu=16G
#
## Recommended safety settings:
set -o errexit # Make bash exit on any error
set -o nounset # Treat unset variables as errors
#
# Set up job environment:
module purge
#
#These are the previous required modules:
#load required modules
#
#Trying new way of loading all modules in Saga: 
module load TransDecoder/5.5.0-intel-2018b-Perl-5.28.0
#
#Listing modules for easier error-check
module list

## Copy inputfiles to the workdirectory
## Note: can run with precomputed blastdbs by uncommenting the blastrestultsdir copy command below, then you don't need to cp the pep files to scratch
#cp -r pep_files $SCRATCH

## Make sure the results are copied back to the submit directory:
#chkfile "pep_files/Results*" 

## Now do something useful:
cd /cluster/home/siribi/work/EVOTREE/transdecoder/

##Step 1: Extract the long open reading frames:
TransDecoder.LongOrfs -t Poptri.SHORT.cds.fasta

##Step 2 (optional): Identify ORFs with homology to known protein via blast or pfam searches

##Step 3: Predict the likely coding regions (possible to add strand info)
TransDecoder.Predict -t Poptri.SHORT.cds.fasta
