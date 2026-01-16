#!/bin/bash
# Job name: Prepare for distmat array job
#SBATCH --job-name=prep_codeml
#
# Project:
#SBATCH --account=nn8023k
#
# Wall clock limit:
#SBATCH --time=10:00:00
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
#module load XXX
#
#Listing modules for easier error-check
module list

#Define parent directory
PARENT=/cluster/home/siribi/work/EVOTREE/CODEML/alignments/MALPIGH/CODEML/

#Define script directory
SCRIPTS=/cluster/home/siribi/evotree_project/ALIGNMENTS/scripts_2a_Malpigh/

## Enter PARENT directory:
cd $PARENT
#mkdir codeml_input_files
#cd codeml_input_files
#rsync -av ../FASTAS_Guidance_Edits/*.fasta .
#sleep 10m
#$SCRIPTS/make_phylip_hybrids.sh
#cd ..
mkdir CODEML_runs
cd CODEML_runs
mkdir alternative_model_omega_0.2_kappa_1
mkdir alternative_model_omega_1_kappa_2
mkdir alternative_model_omega_2_kappa_3
mkdir alternative_model_omega_4_kappa_2
mkdir null_model_omega_1_kappa_1
mkdir null_model_omega_1_kappa_2
mkdir null_model_omega_1_kappa_3
mkdir mlc_files_ALT_omega_0.2_kappa_1
mkdir mlc_files_ALT_omega_1_kappa_2
mkdir mlc_files_ALT_omega_2_kappa_3 
mkdir mlc_files_ALT_omega_4_kappa_2
mkdir mlc_files_NULL_omega_1_kappa_1
mkdir mlc_files_NULL_omega_1_kappa_2
mkdir mlc_files_NULL_omega_1_kappa_3
echo directories are up!

#FIRST DIRECTORY
echo Preparing OMEGA_0.2 x KAPPA_1
cd alternative_model_omega_0.2_kappa_1
rsync -av ../../codeml_input_files/*.phy .
sleep 10m
cp $SCRIPTS/ALT_omega_0.2_kappa_1/BULKRUN.Malpigh.alt.ctl . 
$SCRIPTS/makeCTLfiles.sh
rm BULKRUN.Malpigh.alt.ctl
$SCRIPTS/make_directories.sh
cp $SCRIPTS/Bulk_run_Malpigh.tree .
for d in */; do cp Bulk_run_Malpigh.tree "$d"; done
cp $SCRIPTS/worker_codeml_alt.sbatch .
echo Alternative model OMEGA_0.2 x KAPPA_1 ready to go!
cd ..
echo Moving on...

#NEXT DIRECTORY
echo Preparing OMEGA_1 x KAPPA_2
cd alternative_model_omega_1_kappa_2
rsync -av ../../codeml_input_files/*.phy .
sleep 10m
cp $SCRIPTS/ALT_omega_1_kappa_2/BULKRUN.Malpigh.alt.ctl .
$SCRIPTS/makeCTLfiles.sh
rm BULKRUN.Malpigh.alt.ctl
$SCRIPTS/make_directories.sh
cp $SCRIPTS/Bulk_run_Malpigh.tree .
for d in */; do cp Bulk_run_Malpigh.tree "$d"; done
cp $SCRIPTS/worker_codeml_alt.sbatch .
echo Alternative model OMEGA_1 x KAPPA_2 ready to go!
cd ..
echo moving on...

#NEXT DIRECTORY
echo Preparing OMEGA_2 x KAPPA_3
cd alternative_model_omega_2_kappa_3
rsync -av ../../codeml_input_files/*.phy .
sleep 10m
cp $SCRIPTS/ALT_omega_2_kappa_3/BULKRUN.Malpigh.alt.ctl .
$SCRIPTS/makeCTLfiles.sh
rm BULKRUN.Malpigh.alt.ctl
$SCRIPTS/make_directories.sh
cp $SCRIPTS/Bulk_run_Malpigh.tree .
for d in */; do cp Bulk_run_Malpigh.tree "$d"; done
cp $SCRIPTS/worker_codeml_alt.sbatch .
echo Alternative model OMEGA_2 x KAPPA_3 ready to go!
cd ..
echo moving on...

#NEXT DIRECTORY
echo Preparing OMEGA_4 x KAPPA_2
cd alternative_model_omega_4_kappa_2
rsync -av ../../codeml_input_files/*.phy .
sleep 10m
cp $SCRIPTS/ALT_omega_4_kappa_2/BULKRUN.Malpigh.alt.ctl .
$SCRIPTS/makeCTLfiles.sh
rm BULKRUN.Malpigh.alt.ctl
$SCRIPTS/make_directories.sh
cp $SCRIPTS/Bulk_run_Malpigh.tree .
for d in */; do cp Bulk_run_Malpigh.tree "$d"; done
cp $SCRIPTS/worker_codeml_alt.sbatch .
echo Alternative model OMEGA_4 x KAPPA_2 ready to go!
cd ..
echo moving on...

echo Here comes the null models!
#FIRST NULL MODEL DIRECTORY
echo Preparing OMEGA_1 x KAPPA_1
cd null_model_omega_1_kappa_1
rsync -av ../../codeml_input_files/*.phy .
sleep 10m
cp $SCRIPTS/NULL_omega_1_kappa_1/BULKRUN.Malpigh.fixed.ctl .
$SCRIPTS/makeCTLfiles_null.sh
rm BULKRUN.Malpigh.fixed.ctl
$SCRIPTS/make_directories.sh
cp $SCRIPTS/Bulk_run_Malpigh.tree .
for d in */; do cp Bulk_run_Malpigh.tree "$d"; done
cp $SCRIPTS/worker_codeml_null.sbatch .
echo Null model OMEGA_1 x KAPPA_1 ready to go!
cd ..
echo moving on...

#NEXT DIRECTORY
echo Preparing OMEGA_1 x KAPPA_2
cd null_model_omega_1_kappa_2
rsync -av ../../codeml_input_files/*.phy .
sleep 10m
cp $SCRIPTS/NULL_omega_1_kappa_2/BULKRUN.Malpigh.fixed.ctl .
$SCRIPTS/makeCTLfiles_null.sh
rm BULKRUN.Malpigh.fixed.ctl
$SCRIPTS/make_directories.sh
cp $SCRIPTS/Bulk_run_Malpigh.tree .
for d in */; do cp Bulk_run_Malpigh.tree "$d"; done
cp $SCRIPTS/worker_codeml_null.sbatch .
echo Null model OMEGA_1 x KAPPA_2 ready to go!
cd ..
echo moving on...

#NEXT DIRECTORY
echo Preparing OMEGA_1 x KAPPA_3
cd null_model_omega_1_kappa_3
rsync -av ../../codeml_input_files/*.phy .
sleep 10m
cp $SCRIPTS/NULL_omega_1_kappa_3/BULKRUN.Malpigh.fixed.ctl .
$SCRIPTS/makeCTLfiles_null.sh
rm BULKRUN.Malpigh.fixed.ctl
$SCRIPTS/make_directories.sh
cp $SCRIPTS/Bulk_run_Malpigh.tree .
for d in */; do cp Bulk_run_Malpigh.tree "$d"; done
cp $SCRIPTS/worker_codeml_null.sbatch .
echo Null model OMEGA_1 x KAPPA_3 ready to go!
cd ..
echo Takk for meg.


