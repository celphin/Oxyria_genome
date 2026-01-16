#!/bin/bash
# Job name: Prepare for distmat array job
#SBATCH --job-name=prep_distmat
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
module load mafft

#Define parent directory
PARENT=/usit/abel/u1/siribi/nobackup/COLD/PIPELINE_RUN_DECEMBER2018/Building_alignments_for_CODEML/draba_MCOGs/

#Define script directory
SCRIPTS=/usit/abel/u1/siribi/nobackup/COLD/PIPELINE_RUN_DECEMBER2018/Building_alignments_for_CODEML/draba_MCOGs/scripts_2a_Draba/

## Enter PARENT directory:
cd $PARENT
## Run the python scripts and mafft.
##remember to modify paths and fasta names in the two python scripts. Target directory in MAFFTprep_allOrthogroups.py should be called "orthofasta".  
python $SCRIPTS/MAFFTprep_allOrthogroups.py
cd orthofasta 
$SCRIPTS/run_mafft.sh
python $SCRIPTS/sort_gene_of_interest_v8.py 
mkdir sorted_aligned
cd sorted_aligned
rsync -av ../*.sorted.aligned.fasta .
sleep 10m
$SCRIPTS/make_directories.sh
cp $SCRIPTS/distance_* .
echo MISSION COMPLETED! Please start array run of distmat.
