#!/bin/bash
# Project:COLD 
#SBATCH --account=nn9370k
# Wall clock limit:
#SBATCH --time=6:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=1
source /cluster/bin/jobsetup
module purge   # clear any inherited modules
set -o errexit # exit on errors

module load emboss

# Define some stuff
DISTMAT_OUT=/usit/abel/u1/siribi/nobackup/COLD/PIPELINE_RUN_DECEMBER2018/Building_alignments_for_CODEML/draba_MCOGs/orthofasta/sorted_aligned/distmat/

# Copy directories to scratch
cp -r $SUBMITDIR/dir.$TASK_ID $SCRATCH

cd $SCRATCH

#run distmat:
cd dir.$TASK_ID
FALIST=`ls -1 *.sorted.aligned.fasta`
i=1
for FAFILE in $FALIST
        do
        FILEBASE=${FAFILE%.sorted.aligned.fasta}
        echo $FILEBASE
        distmat $FILEBASE.sorted.aligned.fasta -protmethod 2 -outfile $FILEBASE.distmat          
        cp $FILEBASE.distmat $DISTMAT_OUT
	done
