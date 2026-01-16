#!/bin/bash
# Job name: Prepare for distmat array job
#SBATCH --job-name=Ros_prep_guidance
#
# Project:
#SBATCH --account=nn8023k
#
# Wall clock limit:
#SBATCH --time=72:00:00
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
#module load Python/2.7.15-intel-2018b
module load Python/2.7.15-foss-2018b
#
#Listing modules for easier error-check
module list
#
#Define parent directory
PARENT=/cluster/home/siribi/work/EVOTREE/CODEML/alignments/ROSALES/
#
#Define script directory
SCRIPTS=/cluster/home/siribi/evotree_project/ALIGNMENTS/
#
#Define distmat directory
DISTMAT=/cluster/home/siribi/work/EVOTREE/CODEML/alignments/ROSALES/orthofasta/sorted_aligned/distmat/
#
#REMEMBER TO:
#Change the name of the grep_and_transpose_*.sh script below
#Change "path", "paths", "fdna" and "pahti" in  GUIDANCE_prep_ROSALES.py

##Process distmat output files, make single-copy index files and make files for guidance:
cd $PARENT
mkdir single_copy_index_files
cd $DISTMAT
$SCRIPTS/file_processing_distmat_2.sh
$SCRIPTS/grep_and_transpose_ROSALES.sh
cd $PARENT
cd single_copy_index_files
rsync -av $DISTMAT/*.transposed.txt .
sleep 15m
python $SCRIPTS/GUIDANCE_prep_ROSALES.py
echo SECOND MISSION COMPLETED! Please filter alignments with guidance.
