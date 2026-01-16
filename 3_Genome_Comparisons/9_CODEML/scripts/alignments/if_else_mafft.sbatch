#!/bin/bash
# Job name: Prepare for distmat array job
#SBATCH --job-name=if_mafft
#
# Project:
#SBATCH --account=nn8023k
#
# Wall clock limit:
#SBATCH --time=48:00:00
#
# Number of tasks/processes to run:
#SBATCH --ntasks=1
#
# Number of CPUs (this is what I learned from Saga support if I want to use the -t flag):
#SBATCH --cpus-per-task=4
#
# Max memory usage per core (MB):
#SBATCH --mem-per-cpu=5G
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
module load MAFFT/7.427-intel-2018b-with-extensions
#
#Listing modules for easier error-check
module list
#
## Copy inputfiles to the workdirectory
#cp -r pep_files $SCRATCH
#
## Now do something useful:
cd /cluster/home/siribi/work/EVOTREE/CODEML/alignments/All40/
FALIST=`ls -1 *.fa`
i=1
for FAFILE in $FALIST
        do
        FILEBASE=${FAFILE%.fa}
        echo $FILEBASE
        if mafft --thread -1 --auto $FILEBASE.fa > $FILEBASE.aligned.fa | grep -m 1 "Illegal"; then
        mv $FILEBASE.aligned.fa $FILEBASE.ERROR.fa        
        continue
        fi
        done
