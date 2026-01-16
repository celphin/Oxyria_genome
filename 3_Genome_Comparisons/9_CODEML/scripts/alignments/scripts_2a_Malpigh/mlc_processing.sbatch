#!/bin/bash
# Job name: Prepare for distmat array job
#SBATCH --job-name=mlc_proc
#
# Project:
#SBATCH --account=nn8023k
#
# Wall clock limit:
#SBATCH --time=24:00:00
#
# Number of tasks/processes to run:
#SBATCH --ntasks=1
#
## Number of cores:
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
#
#Define parent directory
PARENT=/cluster/home/siribi/work/EVOTREE/CODEML/alignments/MALPIGH/CODEML/Linum_CODEML_runs/

#Define script directory
SCRIPTS=/cluster/home/siribi/evotree_project/ALIGNMENTS/scripts_2a_Malpigh/

#Define directory where highest LnL will be determined
HIGHEST_LNL=/cluster/home/siribi/work/EVOTREE/CODEML/alignments/MALPIGH/CODEML/Linum_CODEML_runs/HIGHEST_lnL/
R_DIR=/cluster/home/siribi/work/EVOTREE/CODEML/alignments/MALPIGH/CODEML/Linum_CODEML_runs/R_processing/
SITE_CLASSES=/cluster/home/siribi/work/EVOTREE/CODEML/alignments/MALPIGH/CODEML/Linum_CODEML_runs/Site_Class_Tables/

## Enter PARENT directory:
cd $PARENT
mkdir HIGHEST_lnL
mkdir R_processing
mkdir Site_Class_Tables

cd mlc_files_ALT_omega_0.2_kappa_1
mkdir USE_THESE
$SCRIPTS/check_sequence_number.sh
grep -l "9" *.NoSeq.txt >> Files_with_all_sequences.txt
sleep 10m
sed -i 's/.NoSeq.txt/.alt.mlc/' ./Files_with_all_sequences.txt
sleep 5m
cat Files_with_all_sequences.txt | xargs mv -t USE_THESE/
sleep 5m
cd USE_THESE/
$SCRIPTS/ONE_mlc_file.sh 
sleep 10m
mv SORTED_lnL.file SORTED_ALT_omega_0.2_kappa_1.mlc
cp SORTED_ALT_omega_0.2_kappa_1.mlc $HIGHEST_LNL
$SCRIPTS/make_site_class_tables.sh
$SCRIPTS/add_extention_to_site_class_tables_0.2x1.sh
cp *.site.classes.txt $SITE_CLASSES/
sleep 10m 
cd ../..
echo mlc files from 0.2x1 finished! 

cd mlc_files_ALT_omega_1_kappa_2
mkdir USE_THESE
$SCRIPTS/check_sequence_number.sh
grep -l "9" *.NoSeq.txt >> Files_with_all_sequences.txt
sleep 10m
sed -i 's/.NoSeq.txt/.alt.mlc/' ./Files_with_all_sequences.txt
sleep 5m
cat Files_with_all_sequences.txt | xargs mv -t USE_THESE/
sleep 5m
cd USE_THESE/
$SCRIPTS/ONE_mlc_file.sh
sleep 10m
mv SORTED_lnL.file SORTED_ALT_omega_1_kappa_2.mlc
cp SORTED_ALT_omega_1_kappa_2.mlc $HIGHEST_LNL
$SCRIPTS/make_site_class_tables.sh
$SCRIPTS/add_extention_to_site_class_tables_1x2.sh
cp *.site.classes.txt $SITE_CLASSES/
sleep 10m
cd ../..
echo mlc files from 1x2 finished!  

cd mlc_files_ALT_omega_2_kappa_3
mkdir USE_THESE
$SCRIPTS/check_sequence_number.sh
grep -l "9" *.NoSeq.txt >> Files_with_all_sequences.txt
sleep 10m
sed -i 's/.NoSeq.txt/.alt.mlc/' ./Files_with_all_sequences.txt
sleep 5m
cat Files_with_all_sequences.txt | xargs mv -t USE_THESE/
sleep 5m
cd USE_THESE/
$SCRIPTS/ONE_mlc_file.sh
sleep 10m
mv SORTED_lnL.file SORTED_ALT_omega_2_kappa_3.mlc
cp SORTED_ALT_omega_2_kappa_3.mlc $HIGHEST_LNL
$SCRIPTS/make_site_class_tables.sh
$SCRIPTS/add_extention_to_site_class_tables_2x3.sh
cp *.site.classes.txt $SITE_CLASSES/
sleep 10m 
cd ../..
echo mlc files from 2x3 finished! 

cd mlc_files_ALT_omega_4_kappa_2
mkdir USE_THESE
$SCRIPTS/check_sequence_number.sh
grep -l "9" *.NoSeq.txt >> Files_with_all_sequences.txt
sleep 10m
sed -i 's/.NoSeq.txt/.alt.mlc/' ./Files_with_all_sequences.txt
sleep 5m
cat Files_with_all_sequences.txt | xargs mv -t USE_THESE/
sleep 5m
cd USE_THESE/
$SCRIPTS/ONE_mlc_file.sh
sleep 10m
mv SORTED_lnL.file SORTED_ALT_omega_4_kappa_2.mlc
cp SORTED_ALT_omega_4_kappa_2.mlc $HIGHEST_LNL
$SCRIPTS/make_site_class_tables.sh
$SCRIPTS/add_extention_to_site_class_tables_4x2.sh
cp *.site.classes.txt $SITE_CLASSES/
sleep 10m 
cd ../..
echo mlc files from 2x4 finished! 

echo Here comes the NULL models! 
cd mlc_files_NULL_omega_1_kappa_1
mkdir USE_THESE
$SCRIPTS/check_sequence_number_NULL.sh
grep -l "9" *.NoSeq.txt >> Files_with_all_sequences.txt
sleep 10m
sed -i 's/.NoSeq.txt/.fixed.mlc/' ./Files_with_all_sequences.txt
sleep 5m
cat Files_with_all_sequences.txt | xargs mv -t USE_THESE/
sleep 5m
cd USE_THESE/
$SCRIPTS/ONE_mlc_file.sh
sleep 10m
mv SORTED_lnL.file SORTED_NULL_omega_1_kappa_1.mlc
cp SORTED_NULL_omega_1_kappa_1.mlc $HIGHEST_LNL
cd ../..
echo mlc files from 1x1 finished! 

cd mlc_files_NULL_omega_1_kappa_2
mkdir USE_THESE
$SCRIPTS/check_sequence_number_NULL.sh
grep -l "9" *.NoSeq.txt >> Files_with_all_sequences.txt
sleep 10m
sed -i 's/.NoSeq.txt/.fixed.mlc/' ./Files_with_all_sequences.txt
sleep 5m
cat Files_with_all_sequences.txt | xargs mv -t USE_THESE/
sleep 5m
cd USE_THESE/
$SCRIPTS/ONE_mlc_file.sh
sleep 10m
mv SORTED_lnL.file SORTED_NULL_omega_1_kappa_2.mlc
cp SORTED_NULL_omega_1_kappa_2.mlc $HIGHEST_LNL
cd ../..
echo mlc files from 1x2 finished! 

cd mlc_files_NULL_omega_1_kappa_3
mkdir USE_THESE
$SCRIPTS/check_sequence_number_NULL.sh
grep -l "9" *.NoSeq.txt >> Files_with_all_sequences.txt
sleep 10m
sed -i 's/.NoSeq.txt/.fixed.mlc/' ./Files_with_all_sequences.txt
sleep 5m
cat Files_with_all_sequences.txt | xargs mv -t USE_THESE/
sleep 5m
cd USE_THESE/
$SCRIPTS/ONE_mlc_file.sh
sleep 10m
mv SORTED_lnL.file SORTED_NULL_omega_1_kappa_3.mlc
cp SORTED_NULL_omega_1_kappa_3.mlc $HIGHEST_LNL
cd ../..
echo mlc files from 1x3 finished! 

echo Moving on to calculate highest lnL!
cd $HIGHEST_LNL
#Making one table with all lnLs from alternative models: 
join -j 1 -o 1.1,1.2,2.2 SORTED_ALT_omega_0.2_kappa_1.mlc SORTED_ALT_omega_1_kappa_2.mlc > nearly_ALL_lnL_ALT.mlc
join -j 1 -o 1.1,1.2,1.3,2.2 nearly_ALL_lnL_ALT.mlc SORTED_ALT_omega_2_kappa_3.mlc > nearly_nearly_ALL_lnL_ALT.mlc
join -j 1 -o 1.1,1.2,1.3,1.4,2.2 nearly_nearly_ALL_lnL_ALT.mlc SORTED_ALT_omega_4_kappa_2.mlc > ALL_lnL_ALT.mlc

#removing "nan"
sed '/nan/d' ALL_lnL_ALT.mlc > noNAN_ALL_lnL_ALT.mlc 

##removing "OG-number"
awk '{print $2,$3,$4,$5}' noNAN_ALL_lnL_ALT.mlc > noOG_noNAN_ALL_ALT_lnL

##sorting each line hortizontally
awk ' {split( $0, a, " " ); asort( a ); for( i = 1; i <= length(a); i++ ) printf( "%s ", a[i] ); printf( "\n" ); }' noOG_noNAN_ALL_ALT_lnL > sorted_noOG_noNAN_ALT.file

##make a file containing only OG-number:
awk '{print $1}' noNAN_ALL_lnL_ALT.mlc > OG_numbers_ALT 

##paste sorted lnL-values in OG-number file:
paste -d' ' OG_numbers_ALT sorted_noOG_noNAN_ALT.file > sorted_with_OG_numbers_ALT.mlc

##Choose desired columns (OG-number and highest lnL):
awk '{print $1,$5}' sorted_with_OG_numbers_ALT.mlc > HIGHEST_lnL_ALT.mlc

#####SAME FOR NULL MODEL:###########################################################
join -j 1 -o 1.1,1.2,2.2 SORTED_NULL_omega_1_kappa_1.mlc SORTED_NULL_omega_1_kappa_2.mlc > nearly_ALL_lnL_NULL.mlc
join -j 1 -o 1.1,1.2,1.3,2.2 nearly_ALL_lnL_NULL.mlc SORTED_NULL_omega_1_kappa_3.mlc > ALL_lnL_NULL.mlc

#removing "nan"
sed '/nan/d' ALL_lnL_NULL.mlc > noNAN_ALL_lnL_NULL.mlc

##removing "OG-number"
awk '{print $2,$3,$4}' noNAN_ALL_lnL_NULL.mlc > noOG_noNAN_ALL_NULL_lnL

##sorting each line hortizontally
awk ' {split( $0, a, " " ); asort( a ); for( i = 1; i <= length(a); i++ ) printf( "%s ", a[i] ); printf( "\n" ); }' noOG_noNAN_ALL_NULL_lnL > sorted_noOG_noNAN_NULL.file

##make a file containing only OG-number:
awk '{print $1}' noNAN_ALL_lnL_NULL.mlc > OG_numbers_NULL

##paste sorted lnL-values in OG-number file:
paste -d' ' OG_numbers_NULL sorted_noOG_noNAN_NULL.file > sorted_with_OG_numbers_NULL.mlc

##Choose desired columns (OG-number and highest lnL):
awk '{print $1,$4}' sorted_with_OG_numbers_NULL.mlc > HIGHEST_lnL_NULL.mlc

#########JOIN THE TWO TABLES##########################################
join -j 1 -o 1.1,1.2,2.2 HIGHEST_lnL_ALT.mlc HIGHEST_lnL_NULL.mlc > BOTH.mlc
echo BOTH.mlc is finished.
echo Remember that ALT lnL is in column 2 and NULL lnL is in column 3.

cp BOTH.mlc $R_DIR
cp BOTH.mlc $SITE_CLASSES
cd $SITE_CLASSES 
mkdir index_files
cd index_files
mv ../BOTH.mlc .
awk '{print $1, $2}' BOTH.mlc > ALT_lnL.mlc
awk '{for(i=2;i<=NF;i++){printf "%s ", $i >> $1".txt"};printf "\n" >> $1".txt"; close($1".txt")}' ALT_lnL.mlc
sleep 10m
echo You can now grep for lnLs to obtain site class tables. Run make_list_of_site_tables.sh
 
cd $R_DIR
echo Getting ready to run R! So exciting...
module load R/4.0.0-foss-2020a
R CMD BATCH $SCRIPTS/LogLikelihoodTest.R

echo Your results are ready. Takk for meg.
