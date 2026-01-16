#!/bin/bash
# Job name: Orthofinder
#SBATCH --job-name=OrthoFinder
#
# Project:
#SBATCH --account=nn8023k
#
# Wall clock limit:
#SBATCH --time=12:00:00
#
# Number of tasks/processes to run:
#SBATCH --ntasks=1
#
# Number of CPUs (this is what I learned from Saga support if I want to use the -t flag):
#SBATCH --cpus-per-task=40
#
# Max memory usage per core (MB):
#SBATCH --mem-per-cpu=4G
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
#All modules should be loaded by Orthofinder
#
#Trying new way of loading all modules in Saga: 
module load OrthoFinder/2.5.1-foss-2020a-Python-3.8.2
#
#Listing modules for easier error-check
module list

## Copy inputfiles to the workdirectory
## Note: can run with precomputed blastdbs by uncommenting the blastrestultsdir copy command below, then you don't need to cp the pep files to scratch
#cp -r pep_files $SCRATCH

## Make sure the results are copied back to the submit directory:
#chkfile "pep_files/Results*" 

## Now do something useful:
cd /cluster/home/siribi/work/EVOTREE/

# Run Orthofinder, standard mode
orthofinder -b /cluster/home/siribi/work/EVOTREE/TRANS_All40/OrthoFinder/Results_Jan03/WorkingDirectory/ -f /cluster/home/siribi/work/EVOTREE/TRANS_All40/NewSpeciesFasta/ -t 24 

## To run if you have pre-computed blastdbs
#orthofinder -b $BLASTRESULTSDIR -t 8 -a 8
