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
MLC_OUT=/usit/abel/u1/siribi/work/CODEML/Cochlearia_groenlandica/CODEML_runs/mlc_files_ALT_omega_0.2_kappa_1/

# If needed copy up to scratch
#cp -r $SUBMITDIR/dir.$TASK_ID $SCRATCH

##cd to right directory (either scratch or submit directory) 
cd /usit/abel/u1/siribi/work/CODEML/Cochlearia_groenlandica/CODEML_runs/alternative_model_omega_0.2_kappa_1/
#cd $SCRATCH

#run codeml:
cd dir.$TASK_ID
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
