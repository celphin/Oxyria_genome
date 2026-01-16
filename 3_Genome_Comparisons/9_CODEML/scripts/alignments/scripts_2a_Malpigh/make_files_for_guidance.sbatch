#!/bin/bash
# Job name: Prepare for distmat array job
#SBATCH --job-name=prep_guidance
#
# Project:
#SBATCH --account=nn9370k
#
# Wall clock limit:
#SBATCH --time=48:00:00
#
# Max memory usage per core (MB):
#SBATCH --mem-per-cpu=5G
#
# Number of cores:
#SBATCH --cpus-per-task=1
#
## Set up job environment:
source /cluster/bin/jobsetup
module purge   # clear any inherited modules
set -o errexit # exit on errors
unset LC_CTYPE

#load modules
module load python2/2.7.10.gnu

#Define parent directory
PARENT=/usit/abel/u1/siribi/nobackup/COLD/PIPELINE_RUN_DECEMBER2018/Building_alignments_for_CODEML/draba_MCOGs/

#Define script directory
SCRIPTS=/usit/abel/u1/siribi/nobackup/COLD/PIPELINE_RUN_DECEMBER2018/Building_alignments_for_CODEML/draba_MCOGs/scripts_2a_Draba/

#Define distmat directory
DISTMAT=/usit/abel/u1/siribi/nobackup/COLD/PIPELINE_RUN_DECEMBER2018/Building_alignments_for_CODEML/draba_MCOGs/orthofasta/sorted_aligned/distmat/

#REMEMBER TO:
#Change the name of the grep_and_transpose_*.sh script below
#Change "path", "paths", "fdna" and "pahti" in  GUIDANCE_prep.py

##Process distmat output files, make single-copy index files and make files for guidance:
cd $PARENT
mkdir single_copy_index_files
cd $DISTMAT
$SCRIPTS/file_processing_distmat_2.sh 
$SCRIPTS/grep_and_transpose_DRABA.sh
cd $PARENT
cd single_copy_index_files
rsync -av ../orthofasta/sorted_aligned/distmat/*.transposed.txt .
sleep 10m
python $SCRIPTS/GUIDANCE_prep.py
echo SECOND MISSION COMPLETED! Please filter alignments with guidance.
