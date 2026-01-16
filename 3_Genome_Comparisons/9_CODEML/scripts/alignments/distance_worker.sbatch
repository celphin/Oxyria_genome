#!/bin/bash
# Project:EVOTREE
## Run like this: sbatch --array=1-87 distance_worker.sbatch  
#SBATCH --account=XXXXXXX
# Wall clock limit:
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
## Recommended safety settings:
set -o errexit # Make bash exit on any error
set -o nounset # Treat unset variables as errors
module load EMBOSS/6.6.0-intel-2018b
module list

# Define some stuff
DISTMAT_OUT=/cluster/home/siribi/work/EVOTREE/CODEML/alignments/All40/All40_aligned/sorted_aligned/distmat/

# Copy directories to scratch
cp -r $SUBMITDIR/dir.$SLURM_ARRAY_TASK_ID $SCRATCH

cd $SCRATCH

#run distmat:
cd dir.$SLURM_ARRAY_TASK_ID
FALIST=`ls -1 *.sorted.aligned.fa`
i=1
for FAFILE in $FALIST
        do
        FILEBASE=${FAFILE%.sorted.aligned.fa}
        echo $FILEBASE
        distmat $FILEBASE.sorted.aligned.fa -protmethod 2 -outfile $FILEBASE.distmat          
        cp $FILEBASE.distmat $DISTMAT_OUT
	done
