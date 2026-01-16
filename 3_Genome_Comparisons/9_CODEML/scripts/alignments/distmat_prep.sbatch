#!/bin/bash
# Job name: Prepare for distmat array job
#SBATCH --job-name=prep_distmat
#
# Project:
#SBATCH --account=XXXXXXX
#
# Wall clock limit:
#SBATCH --time=6:00:00
#
# Number of tasks/processes to run:
#SBATCH --ntasks=1
#
# Number of cores:
#SBATCH --cpus-per-task=1
#
# Max memory usage per core (MB):
#SBATCH --mem-per-cpu=5G
#
## Set up job environment:
## Recommended safety settings:
set -o errexit # Make bash exit on any error
set -o nounset # Treat unset variables as errors
#
#purge and load modules
module purge
module load Python/2.7.16-GCCcore-8.3.0
module load MAFFT/7.453-GCC-8.3.0-with-extensions  
#
#Listing modules for easier error-check
module list
#
#Define parent directory
PARENT=/cluster/home/siribi/work/CODEML_WORKFLOW_TEST/

#Define script directory
SCRIPTS=/cluster/home/siribi/evotree_project/TEST_CODEML_WORKFLOW/

## Enter PARENT directory:
cd $PARENT
cd orthofasta
$SCRIPTS/if_else_mafft.sh #Possible to skip this step if you already have aligned sequences
python $SCRIPTS/sort_gene_of_interest_v8.py 
mkdir sorted_aligned
cd sorted_aligned
rsync -av ../*.sorted.aligned.fa .
sleep 10m
$SCRIPTS/make_100files_directories.sh
cp $SCRIPTS/distance_* .
mkdir distmat
echo MISSION COMPLETED! Please start array run of distmat.
