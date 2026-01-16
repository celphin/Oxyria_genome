#!/bin/bash
# Project:EVOTREE
## How to run: sbatch --array=1-340 worker_codeml_alt.sbatch
#SBATCH --account=nn8023k
# Wall clock limit:
#SBATCH --time=2:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
## Recommended safety settings:
set -o errexit # Make bash exit on any error
set -o nounset # Treat unset variables as errors

#Purge and list modules
module purge
#module load Anaconda3/2019.03
#set +eu
#source activate codeml
#set -eu

# It is also recommended to to list loaded modules, for easier debugging:
module list

# Define some stuff
MLC_OUT=/cluster/home/siribi/work/EVOTREE/CODEML/alignments/MALPIGH/CODEML/CODEML_runs/mlc_files_ALT_omega_0.2_kappa_1/

# If needed copy up to scratch
#cp -r $SUBMITDIR/dir.$TASK_ID $SCRATCH

##cd to right directory (either scratch or submit directory) 
cd /cluster/home/siribi/work/EVOTREE/CODEML/alignments/MALPIGH/CODEML/CODEML_runs/alternative_model_omega_0.2_kappa_1/
#cd $SCRATCH

#run codeml:
cd dir.$SLURM_ARRAY_TASK_ID
FALIST=`ls -1 *.phy`
i=1
for FAFILE in $FALIST
        do
        FILEBASE=$(echo $FAFILE | cut -f1 -d.)
        echo $FILEBASE
        mkdir $FILEBASE.alt.files
	if codeml $FILEBASE.alt.ctl | grep -m 1 "Giving up."; then 
		mv r* $FILEBASE.alt.files
        	mv 2* $FILEBASE.alt.files
        	mv 4* $FILEBASE.alt.files
		continue
	fi
	        mv r* $FILEBASE.alt.files
                mv 2* $FILEBASE.alt.files
                mv 4* $FILEBASE.alt.files
                mv lnf $FILEBASE.alt.files
        cp $FILEBASE.alt.mlc $MLC_OUT
	done
