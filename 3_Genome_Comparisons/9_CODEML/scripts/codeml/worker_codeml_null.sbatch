#!/bin/bash
# Project:COLD 
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

# It is also recommended to to list loaded modules, for easier debugging:
module list

# Define some stuff
MLC_OUT=/cluster/home/siribi/work/EVOTREE/CODEML/alignments/MALPIGH/CODEML/CODEML_runs/mlc_files_NULL_omega_1_kappa_1/

# If needed copy up to scratch
#cp -r $SUBMITDIR/dir.$TASK_ID $SCRATCH

## CD to right directory (scratch or submit)
#cd $SCRATCH
cd /cluster/home/siribi/work/EVOTREE/CODEML/alignments/MALPIGH/CODEML/CODEML_runs/null_model_omega_1_kappa_1/

#run codeml:
cd dir.$SLURM_ARRAY_TASK_ID
FALIST=`ls -1 *.phy`
i=1
for FAFILE in $FALIST
        do
        FILEBASE=$(echo $FAFILE | cut -f1 -d.)
        echo $FILEBASE
        mkdir $FILEBASE.fixed.files
	if codeml $FILEBASE.fixed.ctl | grep -m 1 "Giving up."; then
		mv r* $FILEBASE.fixed.files
                mv 2* $FILEBASE.fixed.files
                mv 4* $FILEBASE.fixed.files 
		continue
	fi
	mv r* $FILEBASE.fixed.files
        mv 2* $FILEBASE.fixed.files
        mv 4* $FILEBASE.fixed.files
        mv lnf $FILEBASE.fixed.files
        cp $FILEBASE.fixed.mlc $MLC_OUT
	done
