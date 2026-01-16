#!/bin/bash
# Project:EVOTREE 
## How to run: sbatch --array=1-340 worker_guidance_fi.sbatch
#SBATCH --account=nn8023k
#SBATCH --job-name=guidance
# Wall clock limit:
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=5G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3

## Recommended safety settings:
set -o errexit # Make bash exit on any error
set -o nounset # Treat unset variables as errors
#
module purge
module load Anaconda3/2019.03
set +eu
source activate codeml
set -eu

#Load modules not installed in environment
module load Ruby/2.7.2-GCCcore-10.2.0
module load BioPerl/1.7.8-GCCcore-10.2.0

# It is also recommended to to list loaded modules, for easier debugging:
module list

#make directories
#mkdir FASTAS_Guidance_Edits 
#mkdir MSA_Scores_Guidance
#mkdir Seq_Scores_Guidance

# Define directory 
DIR=/cluster/home/siribi/work/EVOTREE/CODEML/alignments/MALPIGH/poptrem_MCOGs_MALPIGH/single_copy_orthofasta/

# Run Guidance and copy files back
cd $DIR
cd dir.$SLURM_ARRAY_TASK_ID
FALIST=`ls -1 *.fasta`
i=1
for FAFILE in $FALIST
        do
        FILEBASE=$(echo $FAFILE | cut -f1 -d.)
        echo $FILEBASE
        if perl /cluster/home/siribi/project/program_dirs/guidance.v2.02/www/Guidance/guidance.pl --seqFile $FILEBASE.fasta --seqType codon --msaProgram PRANK --bootstraps 10 --outDir $DIR/dir.$SLURM_ARRAY_TASK_ID/$FILEBASE.guidance.out --outOrder as_input --proc_num 3 | grep -m 1 "ERROR"; then
		continue
	fi			
        cd $FILEBASE.guidance.out/
        mv MSA.PRANK.Without_low_SP_Col.With_Names $FILEBASE.guidance.edit.fasta
        cp $FILEBASE.guidance.edit.fasta ../../FASTAS_Guidance_Edits 
        mv MSA.PRANK.Guidance2_msa.scr $FILEBASE.MSA.PRANK.Guidance2_msa.scr 
        cp $FILEBASE.MSA.PRANK.Guidance2_msa.scr ../../MSA_Scores_Guidance
        mv MSA.PRANK.Guidance2_res_pair_seq.scr_with_Names $FILEBASE.MSA.PRANK.Guidance2_res_pair_seq.scr_with_Names
        cp $FILEBASE.MSA.PRANK.Guidance2_res_pair_seq.scr_with_Names ../../Seq_Scores_Guidance
        cd ..
        done
