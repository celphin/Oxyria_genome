#!/bin/bash
# Project:COLD 
#SBATCH --account=nn9370k
# Wall clock limit:
#SBATCH --time=5:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=1
source /cluster/bin/jobsetup
module purge   # clear any inherited modules
set -o errexit # exit on errors

# Define some stuff
MLC_OUT=/usit/abel/u1/siribi/work/CODEML/Cochlearia_groenlandica/CODEML_runs/mlc_files_NULL_omega_1_kappa_1/

# If needed copy up to scratch
#cp -r $SUBMITDIR/dir.$TASK_ID $SCRATCH

## CD to right directory (scratch or submit)
#cd $SCRATCH
cd /usit/abel/u1/siribi/work/CODEML/Cochlearia_groenlandica/CODEML_runs/null_model_omega_1_kappa_1/

#run codeml:
cd dir.$TASK_ID
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
