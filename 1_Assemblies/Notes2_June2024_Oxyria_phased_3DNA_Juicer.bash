# May 5 2022
#-----------------------
# Extra PacBio data download

wget PacBio_data.zip https://sbcshare.sequencing.ubc.ca/s/axxN4S3jXWwNaAq/download
unzip download

#----------------------
#Quality check

#https://www.bioinformatics.babraham.ac.uk/projects/fastqc/INSTALL.txt

cd /home/celphin/scratch/Oxyria/PacBio_April22

module load StdEnv/2020
module load fastqc/0.11.9

fastqc *fastq.gz -o ./fastQC_reports


#--------------------
# total PacBio data download

cd /home/celphin/projects/def-henryg/celphin/Oxyria/PacBio/PacBio_Aug2022

tmux new-session -s download
tmux attach-session -t download

read -p "Login: " login && read -p "Password: " -s password && echo -n "j_username=$login&j_password=$password" > .auth.txt && chmod 600 .auth.txt && wget -O - "https://ces.genomequebec.com/nanuqMPS/readsetList?projectId=21515&tech=Sequel" --no-cookies --no-check-certificate --post-file .auth.txt | wget --no-cookies --no-check-certificate --post-file .auth.txt -ci -; rm -f .auth.txt

# FINISHED --2022-08-15 04:02:12--
# Total wall clock time: 12h 3m 15s
# Downloaded: 8 files, 533G in 11h 47m 42s (12.8 MB/s)

md5sum -c readSets.md5

#-----------------------------------
# bam to fastq
module load StdEnv/2020
module load samtools/1.15.1

cd /home/celphin/projects/def-henryg/celphin/Oxyria/PacBio/PacBio_Aug2022/
samtools bam2fq SAMPLE.bam > SAMPLE.fastq

samtools bam2fq Sequel.RunS218_S2.003.Oxyria_C.subreads.bam > Oxyria_subreads0003.fastq
samtools bam2fq Sequel.RunS218_S2.004.Oxyria_C.subreads.bam > Oxyria_subreads0004.fastq 

# [M::bam2fq_mainloop] discarded 0 singletons
# [M::bam2fq_mainloop] processed 21389047 reads
# [M::bam2fq_mainloop] discarded 0 singletons
# [M::bam2fq_mainloop] processed 23466654 reads

tar -cvzf HiFi-subreads.tar.gz Oxyria_subreads0003.fastq Oxyria_subreads0004.fastq


#-------------------------------------
# Hifiasm
# https://hifiasm.readthedocs.io/en/latest/hic-assembly.html 
# https://github.com/chhylp123/hifiasm

tmux new-session -s Oxy
tmux attach-session -t Oxy

cd /home/celphin/projects/def-henryg/celphin/Oxyria/

# Install hifiasm (requiring g++ and zlib)
git clone https://github.com/chhylp123/hifiasm
cd hifiasm && make

# Hi-C phasing with paired-end short reads in two FASTQ files
hifiasm 
-o Oxyria1_Sept4.asm 
-t32 
--h1 NS.1871.003.IDT_i7_217---IDT_i5_217.Oxyria_Hi-C_R1.fastq.gz 
--h2 NS.1871.003.IDT_i7_217---IDT_i5_217.Oxyria_Hi-C_R2.fastq.gz
HiFi-reads.tar.gz

#--------------------------------

# https://github.com/chhylp123/hifiasm/issues/11 
# should be CCS reads
# maybe format is not right
# test with April reads??

#--------------------------------
# Run ccs reads not subreads

# need ccs reads
tmux new-session -s bam
tmux attach-session -t bam

cd /home/celphin/projects/def-henryg/celphin/Oxyria/PacBio/PacBio_Aug2022

module load StdEnv/2020
module load samtools/1.15.1

samtools bam2fq Sequel.RunS218_S2.003.Oxyria_C.ccs.bam > Oxyria_0003_ccs.fastq 
[M::bam2fq_mainloop] discarded 0 singletons
[M::bam2fq_mainloop] processed 1138973 reads

samtools bam2fq Sequel.RunS218_S2.004.Oxyria_C.ccs.bam > Oxyria_0004_ccs.fastq 
[M::bam2fq_mainloop] discarded 0 singletons
[M::bam2fq_mainloop] processed 1237599 reads

# zip up data
salloc -c1 --time 2:50:00 --mem 1510G --account def-henryg
tar -cvzf HiFi-ccs_reads.tar.gz Oxyria_0003_ccs.fastq Oxyria_0004_ccs.fastq 


# try again with ccs reads
salloc -c32 --time 23:50:00 --mem 1510G --account def-henryg

/home/celphin/projects/def-henryg/celphin/Oxyria/hifiasm/hifiasm
 -o /home/celphin/projects/def-henryg/celphin/Oxyria/Oxyria1_Sept6.asm
 -t32 --h1 /home/celphin/projects/def-henryg/celphin/Oxyria/Hi-C/NS.1871.003.IDT_i7_217---IDT_i5_217.Oxyria_Hi-C_R1.fastq.gz
 --h2 /home/celphin/projects/def-henryg/celphin/Oxyria/Hi-C/NS.1871.003.IDT_i7_217---IDT_i5_217.Oxyria_Hi-C_R2.fastq.gz
 /home/celphin/projects/def-henryg/celphin/Oxyria/PacBio/PacBio_Aug2022/HiFi-ccs_reads.tar.gz

[M::ha_analyze_count] left: none
[M::ha_analyze_count] right: count[56] = 7429776
[M::ha_ft_gen] peak_hom: 56; peak_het: 28
[M::ha_ct_shrink::979.499*7.87] ==> counted 5093711 distinct minimizer k-mers
[M::ha_ft_gen::982.930*7.85@20.872GB] ==> filtered out 5093711 k-mers occurring 280 or more times
[M::ha_opt_update_cov] updated max_n_chain to 280
[M::yak_count] collected 878997595 minimizers
[M::ha_pt_gen::1514.717*7.90] ==> counted 38661657 distinct minimizer k-mers
[M::ha_pt_gen] count[4095] = 0 (for sanity check)
[M::ha_analyze_count] lowest: count[10] = 26177
[M::ha_analyze_count] highest: count[27] = 630362
[M::ha_hist_line]     1: ****************************************************************************************************> 18966805
[M::ha_hist_line]     2: ****************************************************************************************************> 1124580
[M::ha_hist_line]     3: ************************************************ 304805
[M::ha_hist_line]     4: ********************* 133135
[M::ha_hist_line]     5: ************ 74255
[M::ha_hist_line]     6: ******** 47920
[M::ha_hist_line]     7: ****** 35530
[M::ha_hist_line]     8: ***** 29171
[M::ha_hist_line]     9: **** 26359
[M::ha_hist_line]    10: **** 26177
[M::ha_hist_line]    11: ***** 28906
[M::ha_hist_line]    12: ***** 33250
[M::ha_hist_line]    13: ****** 39507
[M::ha_hist_line]    14: ******** 50732
[M::ha_hist_line]    15: *********** 67151
[M::ha_hist_line]    16: ************** 89142
[M::ha_hist_line]    17: ******************* 121318
[M::ha_hist_line]    18: ************************** 161233
[M::ha_hist_line]    19: ********************************* 209539                                                                                    [156/395]
[M::ha_hist_line]    20: ****************************************** 263922
[M::ha_hist_line]    21: **************************************************** 329234
[M::ha_hist_line]    22: *************************************************************** 395978
[M::ha_hist_line]    23: ************************************************************************** 466064
[M::ha_hist_line]    24: ************************************************************************************ 527763
[M::ha_hist_line]    25: ******************************************************************************************* 576581
[M::ha_hist_line]    26: ************************************************************************************************* 613578
[M::ha_hist_line]    27: **************************************************************************************************** 630362
[M::ha_hist_line]    28: **************************************************************************************************** 629509
[M::ha_hist_line]    29: ************************************************************************************************ 606825
[M::ha_hist_line]    30: ******************************************************************************************* 570854
[M::ha_hist_line]    31: *********************************************************************************** 522042
[M::ha_hist_line]    32: ************************************************************************** 465909
[M::ha_hist_line]    33: **************************************************************** 406468
[M::ha_hist_line]    34: ******************************************************* 346852
[M::ha_hist_line]    35: ********************************************** 291324
[M::ha_hist_line]    36: ************************************** 240437
[M::ha_hist_line]    37: ******************************** 199549
[M::ha_hist_line]    38: *************************** 168181
[M::ha_hist_line]    39: *********************** 144095
[M::ha_hist_line]    40: ********************* 129572
[M::ha_hist_line]    41: ******************** 123532
[M::ha_hist_line]    42: ******************** 125919
[M::ha_hist_line]    43: ********************* 135373
[M::ha_hist_line]    44: *********************** 147966
[M::ha_hist_line]    45: ************************** 165869
[M::ha_hist_line]    46: ****************************** 186248
[M::ha_hist_line]    47: ********************************* 209911
[M::ha_hist_line]    48: ************************************* 232500
[M::ha_hist_line]    49: ***************************************** 256126
[M::ha_hist_line]    50: ******************************************** 276333
[M::ha_hist_line]    51: *********************************************** 295641
[M::ha_hist_line]    52: ************************************************* 311695
[M::ha_hist_line]    53: *************************************************** 323590
[M::ha_hist_line]    54: **************************************************** 329718
[M::ha_hist_line]    55: ***************************************************** 332978
[M::ha_hist_line]    56: ***************************************************** 331651
[M::ha_hist_line]    57: **************************************************** 328722
[M::ha_hist_line]    58: ************************************************** 315036
[M::ha_hist_line]    59: *********************************************** 298332
[M::ha_hist_line]    60: ********************************************* 282808
[M::ha_hist_line]    61: ****************************************** 261853
[M::ha_hist_line]    62: ************************************** 240872
[M::ha_hist_line]    63: *********************************** 219614
[M::ha_hist_line]    64: ******************************* 198130
[M::ha_hist_line]    65: **************************** 175598
[M::ha_hist_line]    66: ************************ 153721
[M::ha_hist_line]    67: ********************* 132674
[M::ha_hist_line]    68: ****************** 114684
[M::ha_hist_line]    69: **************** 99156
[M::ha_hist_line]    70: ************* 84544
[M::ha_hist_line]    71: ************ 72509
[M::ha_hist_line]    72: ********** 62353
[M::ha_hist_line]    73: ********* 54281
[M::ha_hist_line]    74: ******** 47588
[M::ha_hist_line]    75: ******* 42134
[M::ha_hist_line]    76: ****** 37711
[M::ha_hist_line]    77: ***** 34569
[M::ha_hist_line]    78: ***** 32293
[M::ha_hist_line]    79: ***** 30532
[M::ha_hist_line]    80: ***** 29084
[M::ha_hist_line]    81: **** 28160
[M::ha_hist_line]    82: **** 27338
[M::ha_hist_line]    83: **** 26909
[M::ha_hist_line]    84: **** 26319
[M::ha_hist_line]    85: **** 25717
[M::ha_hist_line]    86: **** 25739
[M::ha_hist_line]    87: **** 24374
[M::ha_hist_line]    88: **** 24130
[M::ha_hist_line]    89: **** 23094
[M::ha_hist_line]    90: **** 22365
[M::ha_hist_line]    91: *** 21566
[M::ha_hist_line]    92: *** 20933
[M::ha_hist_line]    93: *** 20154
[M::ha_hist_line]    99: *** 17130
[M::ha_hist_line]   100: *** 16586
[M::ha_hist_line]   101: *** 16288
[M::ha_hist_line]   102: *** 15896
[M::ha_hist_line]   103: ** 15711
[M::ha_hist_line]   104: *** 15771
[M::ha_hist_line]   105: ** 15444
[M::ha_hist_line]   106: ** 15572
[M::ha_hist_line]   107: *** 15797
[M::ha_hist_line]   108: ** 15597
[M::ha_hist_line]   109: ** 15155
[M::ha_hist_line]   110: ** 15383
[M::ha_hist_line]   111: ** 15439
[M::ha_hist_line]   112: ** 15124
[M::ha_hist_line]   113: ** 14983
[M::ha_hist_line]   114: ** 14922
[M::ha_hist_line]   115: ** 14795
[M::ha_hist_line]   116: ** 14731
[M::ha_hist_line]   117: ** 14533
[M::ha_hist_line]   118: ** 13704
[M::ha_hist_line]   119: ** 13550
[M::ha_hist_line]   120: ** 13060
[M::ha_hist_line]   121: ** 12985
[M::ha_hist_line]   122: ** 12542
[M::ha_hist_line]   123: ** 12083
[M::ha_hist_line]   124: ** 11720
[M::ha_hist_line]   125: ** 11500
[M::ha_hist_line]   126: ** 10983
[M::ha_hist_line]   127: ** 10604
[M::ha_hist_line]   128: ** 10132
[M::ha_hist_line]   129: ** 9837
[M::ha_hist_line]   130: ** 9638
[M::ha_hist_line]   131: * 9147
[M::ha_hist_line]   132: * 9131
[M::ha_hist_line]   133: * 8903
[M::ha_hist_line]   134: * 8560
[M::ha_hist_line]   135: * 8440
[M::ha_hist_line]   136: * 8225
[M::ha_hist_line]   137: * 7989



[M::ha_hist_line]   176: * 4685
[M::ha_hist_line]   177: * 4749
[M::ha_hist_line]   178: * 4676
[M::ha_hist_line]   179: * 4629
[M::ha_hist_line]   180: * 4536
[M::ha_hist_line]   181: * 4613
[M::ha_hist_line]   182: * 4435
[M::ha_hist_line]   183: * 4313
[M::ha_hist_line]   184: * 4193
[M::ha_hist_line]   185: * 4206
[M::ha_hist_line]   186: * 4149
[M::ha_hist_line]   187: * 4105
[M::ha_hist_line]   188: * 4018
[M::ha_hist_line]   189: * 3976
[M::ha_hist_line]   190: * 3793
[M::ha_hist_line]   191: * 3825
[M::ha_hist_line]   192: * 3845
[M::ha_hist_line]   193: * 3631
[M::ha_hist_line]   194: * 3637
[M::ha_hist_line]   195: * 3643
[M::ha_hist_line]   196: * 3548
[M::ha_hist_line]   197: * 3386
[M::ha_hist_line]   198: * 3469
[M::ha_hist_line]   199: * 3409
[M::ha_hist_line]   200: * 3353
[M::ha_hist_line]   201: * 3301
[M::ha_hist_line]   202: * 3182
[M::ha_hist_line]   203: * 3177
[M::ha_hist_line]   204: * 3226
[M::ha_hist_line]   205: * 3241
[M::ha_hist_line]   206: * 3174
[M::ha_hist_line]  rest: ************************** 162545
[M::ha_analyze_count] left: none
[M::ha_analyze_count] right: count[55] = 332978
[M::ha_pt_gen] peak_hom: 55; peak_het: 27
[M::ha_ct_shrink::1514.835*7.90] ==> counted 19694852 distinct minimizer k-mers
[M::ha_pt_gen::] counting in normal mode
[M::yak_count] collected 878997595 minimizers
[M::ha_pt_gen::1711.608*9.32] ==> indexed 860030790 positions, counted 19694852 distinct minimizer k-mers



[M::ha_hist_line]  rest: **************************** 175767
[M::ha_analyze_count] left: none
[M::ha_analyze_count] right: count[56] = 331532
[M::ha_pt_gen] peak_hom: 56; peak_het: 28
[M::ha_ct_shrink::7463.273*26.53] ==> counted 18079381 distinct minimizer k-mers
[M::ha_pt_gen::] counting in normal mode
[M::yak_count] collected 871849600 minimizers
[M::ha_pt_gen::7655.420*26.38] ==> indexed 870601322 positions, counted 18079381 distinct minimizer k-mers



Writing /home/celphin/projects/def-henryg/celphin/Oxyria/Oxyria1_Sept6.asm.hic.p_ctg.gfa to disk...
[M::build_unitig_index::44.750] ==> Counting
[M::build_unitig_index::8.382] ==> Memory allocating
[M::build_unitig_index::63.636] ==> Filling pos
[M::build_unitig_index::0.277] ==> Sorting pos
[M::build_unitig_index::117.051] ==> HiC index has been built
[M::write_hc_pt_index] Index has been written.
[M::alignment_worker_pipeline::231.830] ==> Qualification
[M::dedup_hits::1.196] ==> Dedup
[M::mc_solve_core::0.129] ==> Partition
[M::dedup_hits::0.658] ==> Dedup
[M::stat] # misjoined unitigs: 7 (N50: 9858488); # corrected unitigs: 14 (N50: 8339145)
[M::adjust_weight_kv_u_trans_advance::0.697]
[M::mb_solve_core::1.460] ==> Partition
[M::mc_solve_core::1.966] ==> Partition
[M::adjust_weight_kv_u_trans_advance::4.504]
[M::mb_solve_core::1.432] ==> Partition
[M::mc_solve_core::1.929] ==> Partition
[M::adjust_weight_kv_u_trans_advance::4.492]
[M::mb_solve_core::1.421] ==> Partition
[M::mc_solve_core::1.902] ==> Partition
[M::stat] # heterozygous bases: 1064182092; # homozygous bases: 101045007
[M::reduce_hamming_error::0.695] # inserted edges: 3992, # fixed bubbles: 7
[M::adjust_utg_by_trio] primary contig coverage range: [47, infinity]
Writing /home/celphin/projects/def-henryg/celphin/Oxyria/Oxyria1_Sept6.asm.hic.hap1.p_ctg.gfa to disk...
[M::adjust_utg_by_trio] primary contig coverage range: [47, infinity]
Writing /home/celphin/projects/def-henryg/celphin/Oxyria/Oxyria1_Sept6.asm.hic.hap2.p_ctg.gfa to disk...
Inconsistency threshold for low-quality regions in BED files: 70%
[M::main] Version: 0.16.1-r375
[M::main] CMD: /home/celphin/projects/def-henryg/celphin/Oxyria/hifiasm/hifiasm -o /home/celphin/projects/def-henryg/celphin/Oxyria/Oxyria1_Sept6.asm -t32 --h1 /home/celphin/projects/def-henryg/celphin/Oxyria/Hi-C/NS.1871.003.IDT_i7_217---IDT_i5_217.Oxyria_Hi-C_R1.fastq.gz --h2 /home/celphin/projects/def-henryg/celphin/Oxyria/Hi-C/NS.1871.003.IDT_i7_217---IDT_i5_217.Oxyria_Hi-C_R2.fastq.gz /home/celphin/projects/def-henryg/celphin/Oxyria/PacBio/PacBio_Aug2022/HiFi-ccs_reads.tar.gz
[M::main] Real time: 20581.846 sec; CPU: 571748.638 sec; Peak RSS: 67.930 GB

#---------------------------------
# get fasta from gfa
# https://hifiasm.readthedocs.io/en/latest/interpreting-output.html

cd /home/celphin/projects/rpp-rieseber/celphin/Oxyria/Oxyria_haplotigs_Sept6/

awk '/^S/{print ">"$2;print $3}' Oxyria1_Sept6.asm.hic.p_ctg.gfa > Oxyria1_Sept6_haplotigs.p_ctg.fa

# Hap 1
awk '/^S/{print ">"$2;print $3}' Oxyria1_Sept6.asm.hic.hap1.p_ctg.gfa > Oxyria1_Sept6_haplotigs.Hap1.fa

# Hap 2
awk '/^S/{print ">"$2;print $3}' Oxyria1_Sept6.asm.hic.hap2.p_ctg.gfa > Oxyria1_Sept6_haplotigs.Hap2.fa

# -rw-r-----  1 celphin rrg-rieseber-ac   582631719 Jun 18 13:16 Oxyria1_Sept6_haplotigs.Hap1.fa
# -rw-r-----  1 celphin rrg-rieseber-ac   566600387 Jun 18 13:18 Oxyria1_Sept6_haplotigs.Hap2.fa
# -rw-r--r--  1 celphin rrg-rieseber-ac   590386466 Nov 10  2023 Oxyria1_Sept6_haplotigs.p_ctg.fa

#--------------------------
# check fasta N50 etc.

# Assembly stats

module load StdEnv/2020 
module load gcc/9.3.0
module load r-bundle-bioconductor/3.14

cd /home/celphin/assembly-stats

./assembly-stats /home/celphin/projects/rpp-rieseber/celphin/Oxyria/Oxyria_haplotigs_Sept6/Oxyria1_Sept6_haplotigs.p_ctg.fa

stats for /home/celphin/projects/rpp-rieseber/celphin/Oxyria/Oxyria_haplotigs_Sept6/Oxyria1_Sept6_haplotigs.p_ctg.fa
sum = 590 373 804, n = 974, ave = 606 133.27, largest = 100 326 276
N50 = 56 351 117, n = 4
N60 = 39 295 833, n = 6
N70 = 34 734 330, n = 7
N80 = 28 174 138, n = 9
N90 = 10 901 797, n = 13
N100 = 14349, n = 974
N_count = 0
Gaps = 0

2n = 14 and 42

Genome size: 800-1000 Mbp
#----------------------------
cd /home/celphin/assembly-stats

./assembly-stats /home/celphin/projects/rpp-rieseber/celphin/Oxyria/Oxyria_haplotigs_Sept6/Oxyria1_Sept6_haplotigs.Hap1.fa

stats for /home/celphin/projects/rpp-rieseber/celphin/Oxyria/Oxyria_haplotigs_Sept6/Oxyria1_Sept6_haplotigs.Hap1.fa
sum = 582 616 305, n = 1101, ave = 529170.12, largest = 76049245
N50 = 38 327 256, n = 6
N60 = 26 598 077, n = 8
N70 = 17086624, n = 11
N80 = 10595219, n = 15
N90 = 2919708, n = 27
N100 = 12798, n = 1101
N_count = 0
Gaps = 0


./assembly-stats /home/celphin/projects/rpp-rieseber/celphin/Oxyria/Oxyria_haplotigs_Sept6/Oxyria1_Sept6_haplotigs.Hap2.fa

stats for /home/celphin/projects/rpp-rieseber/celphin/Oxyria/Oxyria_haplotigs_Sept6/Oxyria1_Sept6_haplotigs.Hap2.fa
sum = 566 594 759, n = 402, ave = 1409439.70, largest = 49230237
N50 = 27 807 003, n = 9
N60 = 26766067, n = 11
N70 = 22934860, n = 13
N80 = 15449834, n = 16
N90 = 5566760, n = 23
N100 = 14349, n = 402
N_count = 0
Gaps = 0


###########################
# send HiC data to Christine

scp -rv /home/celphin/projects/def-henryg/celphin/Oxyria/Hi-C/NS.1871.003.IDT_i7_217---IDT_i5_217.Oxyria_Hi-C_R*.fastq.gz* aphill@kebnekaise.hpc2n.umu.se:/proj/nobackup/snic2021-22-543/oxyria_ref/

folder: aphill@kebnekaise.hpc2n.umu.se:/proj/nobackup/snic2021-22-543/oxyria_ref
Password: 0xyr1a_ref$_144


#################################################################
# Hi-C mapping and scaffolding - run on Cedar
# https://github.com/rieseberglab/haplotype_aware_scaffolding
########################################################

# move data needed over

#try copying input files to main directory
cd /scratch/celphin/Hi-C_scaffold/

scp -v celphin@cedar.computecanada.ca:/home/celphin/projects/rpp-rieseber/celphin/Oxyria/Oxyria_haplotigs_Sept6/Oxyria1_Sept6_haplotigs.Hap1.fa \
/home/celphin/scratch/celphin/Hi-C_scaffold/hic_scaffolding/Oxy_draft_assembly_Hap1.fa

scp -v celphin@cedar.computecanada.ca:/home/celphin/projects/rpp-rieseber/celphin/Oxyria/Oxyria_haplotigs_Sept6/Oxyria1_Sept6_haplotigs.Hap2.fa \
/home/celphin/scratch/Hi-C_scaffold/hic_scaffolding/Oxy_draft_assembly_Hap2.fa

scp -v celphin@cedar.computecanada.ca:/home/celphin/projects/rpp-rieseber/celphin/Oxyria/hic_scaffolding/Oxy_HiC_R1.fastq.gz .
scp -v celphin@cedar.computecanada.ca:/home/celphin/projects/rpp-rieseber/celphin/Oxyria/hic_scaffolding/Oxy_HiC_R2.fastq.gz .

#-------------------------
cd /home/celphin/scratch/Hi-C_scaffold

# follow notes here for mapping:
# https://github.com/rieseberglab/haplotype_aware_scaffolding

##########################################
# install 3D DNA and Juicer

# Juicer
# https://github.com/aidenlab/juicer/wiki
#-----------------------------
# Questions about code on CC please contact: cassandra.elphinstone@shaw.ca
#-----------------------
# Fasta file: /home/celphin/projects/def-henryg/celphin/Oxyria/Oxy_dig_1.scaffolds_FINAL.fasta
# HiC files: /home/celphin/projects/def-henryg/celphin/Oxyria/Hi-C_raw_data/NS.1871.003.IDT_i7_217---IDT_i5_217.Oxyria_Hi-C_R*.fastq.gz

#-------------------
# Dependancies
    # For alignment and creation of the Hi-C pairs file merged_nodups.txt:
        # GNU CoreUtils
        # Burrows-Wheeler Aligner (BWA)
    # For .hic file creation and Juicer tools analysis:
        # Java 1.7 or 1.8 JDK. (Alternative link for Ubuntu/LinuxMint). Minimum system requirements for running Java can be found at http://java.com/en/download/help/sysreq.xml
        # Latest Juicer Tools jar
    # For peak calling:
        # CUDA and an NVIDIA GPU
        # The native libraries included with Juicer are compiled for CUDA 7. Other versions of CUDA can be used, but you will need to download the respective native libraries from JCuda.
        # For best performance, use a dedicated GPU. You may also be able to obtain access to GPU clusters through Amazon Web Services or a local research institution.

#You should have a Juicer directory containing scripts/, references/, and optionally restriction_sites/, and a different working directory containing fastq/. 
#You should download the Juicer tools jar and install it in your scripts/ directory.

#---------------------------
# Dependancies
module load StdEnv/2020
module load nixpkgs/16.09 
module load bwa/0.7.15
module load java/1.8.0_192 
module load gcc/4.8.5
module load cuda/7.5.18

#------------------------
# Setup
# https://github.com/theaidenlab/juicer/wiki/Running-Juicer-on-a-cluster

mkdir Juicer; cd Juicer

# Get Juicer
git clone https://github.com/theaidenlab/juicer.git

cd juicer
# symbolic link to scripts
ln -s SLURM/scripts/ scripts
cd scripts
# Download Juicer tools jar
wget https://github.com/aidenlab/Juicebox/releases/download/v2.18.00/juicer_tools.2.18.00.jar
ln -s juicer_tools.2.18.00.jar juicer_tools.jar

cd /home/celphin/scratch/Hi-C_scaffold/Juicer/juicer
mkdir references/
mkdir restriction_sites
mkdir chromosome_size
mkdir fastq

#-----------------------------
# replace jucier scripts with new CC specific juicer scripts 
git clone https://github.com/rieseberglab/Juicer_Compute_Canada_scripts.git

# copy original scripts and replace with scripts for Cedar SLURM specifically
mv ./scripts/jucier.sh ./scripts/jucier-copy.sh
mv ./scripts/generate_site_positions.py ./scripts/generate_site_positions-copy.py
mv ./scripts/juicer_arrowhead.sh ./scripts/juicer_arrowhead-copy.sh
mv ./scripts/juicer_postprocessing.sh ./scripts/juicer_postprocessing-copy.sh
mv ./scripts/mega.sh ./scripts/mega-copy.sh

#----------------------------
# set default SLURM submission account 
# https://docs.alliancecan.ca/wiki/Running_jobs#Accounts_and_projects
export SLURM_ACCOUNT=def-rieseber
export SBATCH_ACCOUNT=$SLURM_ACCOUNT
export SALLOC_ACCOUNT=$SLURM_ACCOUNT

#----------------
# get 3d-dna
cd /home/celphin/scratch/Hi-C_scaffold/Juicer/juicer
git clone https://github.com/aidenlab/3d-dna.git

########################
git clone https://github.com/rieseberglab/haplotype_aware_scaffolding.git


#-------------------------
# Fix scripts for working on compute canada and slurm
cd /home/celphin/scratch/Hi-C_scaffold/haplotype_aware_scaffolding
dos2unix *.sh

cd /home/celphin/scratch/Hi-C_scaffold/
chmod -R +770 *

nano run_step1_scaffolding.sh
#!/bin/bash

bin=/home/celphin/scratch/Hi-C_scaffold/haplotype_aware_scaffolding 

ID_HAP1=$1 
ID_HAP2=$2 
SAMPLE=$3 
ASSEMBLY_hap1=$4 
ASSEMBLY_hap2=$5 
JUCIER_DIR=$6 
RESULT_DIR=$7 
HIC_R1=$8 #R1.fastq.gz
HIC_R2=$9 #R2.fastq.gz

[ -d ${RESULT_DIR}/${SAMPLE} ] || mkdir -p ${RESULT_DIR}/${SAMPLE}

RESULT_DIR_ID=${RESULT_DIR}/${SAMPLE}

#hap1
bash ${bin}/scaffolding.sh $bin $ID_HAP1 $ASSEMBLY_hap1 $JUCIER_DIR $HIC_R1 $HIC_R2 1 \
$RESULT_DIR_ID &

#hap2
bash ${bin}/scaffolding.sh $bin $ID_HAP2 $ASSEMBLY_hap2 $JUCIER_DIR $HIC_R1 $HIC_R2 1 \
$RESULT_DIR_ID &

#mix hap
[ -d ${RESULT_DIR}/${SAMPLE}/sequences ] || mkdir -p ${RESULT_DIR}/${SAMPLE}/sequences

bash ${bin}/cat_fasta.sh $ASSEMBLY_hap1 $ASSEMBLY_hap2 ${RESULT_DIR}/${SAMPLE}/sequences

bash ${bin}/scaffolding.sh $bin ${ID_HAP1}_${ID_HAP2} \
${RESULT_DIR}/${SAMPLE}/sequences/$(basename ${ASSEMBLY_hap1%%.fasta})_$(basename ${ASSEMBLY_hap2%%.fasta}).fasta \
$JUCIER_DIR $HIC_R1 $HIC_R2 10 $RESULT_DIR_ID

#mini map alignment
[ -d ${RESULT_DIR}/${SAMPLE}/minimap ] || mkdir -p ${RESULT_DIR}/${SAMPLE}/minimap

sbatch ${bin}/minimap.sh $ASSEMBLY_hap1 $ASSEMBLY_hap2 ${RESULT_DIR}/${SAMPLE}/minimap


#--------------------------------------
nano scaffolding.sh
#!/bin/bash
bin=$1 
ID=$2 
ASSEMBLY=$3 
JUCIER_DIR=$4 
HIC_R1=$5 
HIC_R2=$6 
MAPQ=$7 
RESULT_DIR_ID=$8

#preparation for first round scaffolding: indexes the assembly, 
#finds the restriction enzyme locations, calculate the contig sizes

sbatch ${bin}/JUICER_PREP.sh $ASSEMBLY $ID $JUCIER_DIR \
${RESULT_DIR_ID} $HIC_R1 $HIC_R2 ${MAPQ}


#-------------------------------------
nano JUICER_PREP.sh
#!/bin/bash
#SBATCH --account=def-rieseber
#SBATCH --time=0-1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=192000m

#eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"
#conda activate ASSEM_SCAFF2

module load StdEnv/2020
module load nixpkgs/16.09 
module load bwa/0.7.15
module load java/1.8.0_192 
module load gcc/4.8.5
module load cuda/7.5.18
module load python/3.8.2

ASSEMBLY=$1   #assembly name 
ID=$2     #scaffold name in working directory 
JUCIER_DIR=$3 #/lustre04/scratch/mjahani/opt/juicer
RESULT_DIR_ID=$4
HIC_R1=$5
HIC_R2=$6
MAPQ=$7

cp -R -u -p $ASSEMBLY ${JUCIER_DIR}/references/

[ -d ${RESULT_DIR_ID}/LOG ] || mkdir -p ${RESULT_DIR_ID}/LOG

bwa index ${JUCIER_DIR}/references/${ASSEMBLY##*/} 2> \
${RESULT_DIR_ID}/LOG/${ID}_BWA_INDEX.log

cd ${JUCIER_DIR}/restriction_sites

python ${JUCIER_DIR}/misc/generate_site_positions.py MboI \
$ID ${JUCIER_DIR}/references/${ASSEMBLY##*/} 2> \
${RESULT_DIR_ID}/LOG/${ID}_generate_site_positions.log

awk 'BEGIN{OFS="\t"}{print $1, $NF}' ${ID}_MboI.txt \
>${JUCIER_DIR}/chromosome_size/${ID}.chrom.sizes 2> \
${RESULT}/LOG/${ID}_chromosomesize.log

#Script to map HiC data on assembly for 3D DNA scaffolding
cd /home/celphin/scratch/Hi-C_scaffold
sbatch /home/celphin/scratch/Hi-C_scaffold/haplotype_aware_scaffolding/JUICER.sh \
$ASSEMBLY $ID $JUCIER_DIR $HIC_R1 $HIC_R2 ${RESULT_DIR_ID} ${MAPQ}

#-----------------------------------
nano JUICER.sh

#!/bin/bash
#SBATCH --account=def-rieseber
#SBATCH --time=0-1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=192000m

#eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"
#conda activate ASSEM_SCAFF2
module load StdEnv/2020
module load nixpkgs/16.09 
module load bwa/0.7.15
module load java/1.8.0_192 
module load gcc/4.8.5
module load cuda/7.5.18

ASSEMBLY=$1   #assembly name #Kalea.hic.hap1.p_ctg.fasta
ID=$2     #scaffoldinf name in working directory 
JUCIER_DIR=$3 #/lustre04/scratch/mjahani/opt/juicer
HIC_R1=$4 
HIC_R2=$5
RESULT_DIR_ID=$6
MAPQ=$7

[ -d ${RESULT_DIR_ID}/LOG ] || mkdir -p ${RESULT_DIR_ID}/LOG

[ -d ${JUCIER_DIR}/scaffolding/${ID}/fastq ] || \
mkdir -p ${JUCIER_DIR}/scaffolding/${ID}/fastq

cp -R -u -p $HIC_R1 ${JUCIER_DIR}/scaffolding/${ID}/fastq

cp -R -u -p $HIC_R2 ${JUCIER_DIR}/scaffolding/${ID}/fastq

cd ${JUCIER_DIR}/scaffolding/${ID}

${JUCIER_DIR}/scripts/juicer.sh \
    -g $ID \
    -s MboI \
    -t 46 \
    -z ${JUCIER_DIR}/references/${ASSEMBLY##*/} \
    -y ${JUCIER_DIR}/restriction_sites/${ID}_MboI.txt \
    -p ${JUCIER_DIR}/chromosome_size/${ID}.chrom.sizes \
    -D $JUCIER_DIR 2> ${RESULT_DIR_ID}/LOG/${ID}_JUICER.log

rm -rf ${JUCIER_DIR}/scaffolding/${ID}/fastq

#First round of 3D DNA scaffolding
cd /home/celphin/scratch/Hi-C_scaffold
sbatch /home/celphin/scratch/Hi-C_scaffold/haplotype_aware_scaffolding/3DDNA.sh \
$ASSEMBLY $ID $JUCIER_DIR ${MAPQ} ${RESULT_DIR_ID}

#------------------------------------
nano 3DDNA.sh
#!/bin/bash
#SBATCH --account=def-rieseber
#SBATCH --time=0-5:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=192000m

#eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"
#conda activate ASSEM_SCAFF2

module load StdEnv/2020
module load nixpkgs/16.09 
module load bwa/0.7.15
module load java/1.8.0_192 
module load gcc/4.8.5
module load cuda/7.5.18

ASSEMBLY=$1   #assembly name #Kalea.hic.hap1.p_ctg.fasta
ID=$2     #scaffold
JUCIER_DIR=$3 #/lustre04/scratch/mjahani/opt/juicer
MAPQ=$4
RESULT_DIR_ID=$5

[ -d ${RESULT_DIR_ID}/LOG ] || mkdir -p ${RESULT_DIR_ID}/LOG

[ -d ${JUCIER_DIR}/scaffolding/${ID}/3D_DNA ] || \
mkdir -p ${JUCIER_DIR}/scaffolding/${ID}/3D_DNA

[ -d ${RESULT_DIR_ID}/HiC_map ] || mkdir -p ${RESULT_DIR_ID}/HiC_map

cd ${JUCIER_DIR}/scaffolding/${ID}/3D_DNA

export LANGUAGE=en_US.UTF-8
export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8
export LC_CTYPE=en_US.UTF-8

bash ${JUCIER_DIR}/3d-dna/run-asm-pipeline.sh -i 5000 -r 0 -q ${MAPQ} \
--early-exit \
--editor-repeat-coverage 5 \
--editor-coarse-resolution 100000 \
--editor-coarse-region 500000 \
${JUCIER_DIR}/references/${ASSEMBLY##*/} \
${JUCIER_DIR}/scaffolding/${ID}/aligned/merged_nodups.txt 2> \
${RESULT_DIR_ID}/LOG/${ID}_3DDNA.log

cp ${JUCIER_DIR}/scaffolding/${ID}/3D_DNA/$(basename ${ASSEMBLY%%fasta})0.assembly \
${RESULT_DIR_ID}/HiC_map

cp ${JUCIER_DIR}/scaffolding/${ID}/3D_DNA/$(basename ${ASSEMBLY%%fasta})0.hic \
${RESULT_DIR_ID}/HiC_map


#-------------------------------

#cat_fasta.sh 
#--no change

#------------------------------
nano minimap.sh
#!/bin/bash
#SBATCH --account=def-rieseber
#SBATCH --time=0-5:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=192000m

#eval "$(command conda 'shell.bash' 'hook' 2>/dev/null)"
#conda activate scaff_env
module load StdEnv/2023 minimap2/2.28

ASSEMBLY1=$1 #/DATA/home/mjahani/opt/juicer/scaffolding/AGA10_hap1/3D_DNA/AGA10.hic.hap1.p_ctg.final.fasta
ASSEMBLY2=$2 #/DATA/home/mjahani/opt/juicer/scaffolding/AGA10_hap2/3D_DNA/AGA10.hic.hap2.p_ctg.final.fasta
SAVE_DIR=$3  #/DATA/home/mjahani/curation_AGA10

minimap2 -x asm10 -t 46 $ASSEMBLY1 $ASSEMBLY2 \
>${SAVE_DIR}/$(basename "${ASSEMBLY1%%.fasta}")_$(basename "${ASSEMBLY2%%.fasta}").paf

######################################
# to run 
cd ..
dos2unix ./haplotype_aware_scaffolding/*.sh

module load StdEnv/2020
module load nixpkgs/16.09 
module load bwa/0.7.15
module load java/1.8.0_192 
module load gcc/4.8.5
module load cuda/7.5.18

cd /home/celphin/scratch/Hi-C_scaffold

bash ./haplotype_aware_scaffolding/run_step1_scaffolding.sh \
Oxydig_HAP1 \ 
Oxydig_HAP2 \ 
Oxydig \
Oxy_draft_assembly_Hap1.fa \
Oxy_draft_assembly_Hap2.fa \
/home/celphin/scratch/Hi-C_scaffold/Juicer/juicer \
/home/celphin/scratch/Hi-C_scaffold/output \
Oxy_HiC_R1.fastq.gz \
Oxy_HiC_R2.fastq.gz

# ID_HAP1 = ID for haplotype one, example: MK_ultra_hap1
# ID_HAP2 = ID for haplotype two, example: MK_ultra_hap2
# SAMPLE = ID for the genotype, example: MK_ultra
# ASSEMBLY_hap1 = genome assembly of haplotype 1 in fasta format, example: MK_ultra.hic.hap1.p_ctg.fasta
# ASSEMBLY_hap2 = genome assembly of haplotype 2 in fasta format, example: MK_ultra.hic.hap2.p_ctg.fasta
# JUCIER_DIR = Path to where Juicer was installed
# RESULT_DIR = path for saving results
# HIC_R1 = HiC read `R1` in fastq format
# HIC_R2 = HiC read `R2` in fastq format


#-----------------------------------
# Errors

mkdir: cannot create directory '//sequences': Permission denied
/home/celphin/scratch/Hi-C_scaffold/haplotype_aware_scaffolding/cat_fasta.sh: line 10: /sequences_.fasta: Permission denied
basename: missing operand
Try 'basename --help' for more information.
basename: missing operand
Try 'basename --help' for more information.


########################################
# try running directly 

bin=/home/celphin/scratch/Hi-C_scaffold/haplotype_aware_scaffolding 
ID_HAP1=Oxydig_HAP1
ID_HAP2=Oxydig_HAP2
SAMPLE=Oxydig
ASSEMBLY_hap1=Oxy_draft_assembly_Hap1.fa 
ASSEMBLY_hap2=Oxy_draft_assembly_Hap2.fa 
JUCIER_DIR=/home/celphin/scratch/Hi-C_scaffold/Juicer/juicer 
RESULT_DIR=/home/celphin/scratch/Hi-C_scaffold/output
HIC_R1=Oxy_HiC_R1.fastq.gz 
HIC_R2=Oxy_HiC_R2.fastq.gz 

[ -d ${RESULT_DIR}/${SAMPLE} ] || mkdir -p ${RESULT_DIR}/${SAMPLE}

RESULT_DIR_ID=${RESULT_DIR}/${SAMPLE}

#hap1
bash ${bin}/scaffolding.sh $bin $ID_HAP1 $ASSEMBLY_hap1 $JUCIER_DIR $HIC_R1 $HIC_R2 1 \
$RESULT_DIR_ID

#hap2hap
bash ${bin}/scaffolding.sh $bin $ID_HAP2 $ASSEMBLY_hap2 $JUCIER_DIR $HIC_R1 $HIC_R2 1 \
$RESULT_DIR_ID

#----------------------
# Issues

# make slurm scripts as dependancies
# https://stackoverflow.com/questions/65672312/chain-multiple-slurm-jobs-with-dependency
# https://bioinformaticsworkbook.org/Appendix/HPC/SLURM/submitting-dependency-jobs-using-slurm.html#gsc.tab=0

#!/bin/bash

ID=$(sbatch --parsable $1)
shift 
for script in "$@"; do
  ID=$(sbatch --parsable --dependency=afterok:${ID} $script)
done

# for now try just adding sbatch call to end of last script

#----------------------
/scratch/celphin/Hi-C_scaffold/Juicer/juicer/3d-dna/run-asm-pipeline.sh -i 5000 -r 0 -q 1 --early-ex
it --editor-repeat-coverage 5 --editor-coarse-resolution 100000 --editor-coarse-region 500000 /home/
celphin/scratch/Hi-C_scaffold/Juicer/juicer/references/Oxy_draft_assembly_Hap1.fa /home/celphin/scra
tch/Hi-C_scaffold/Juicer/juicer/scaffolding/Oxydig_HAP1/aligned/merged_nodups.txt

version: 180922
 -i|--input flag was triggered, filtering draft contigs/scaffolds smaller than 5000.
 -r|--rounds flag was triggered, will run 0 round(s) of misjoin correction.
 -q|--mapq flag was triggered, scaffolding using reads with at least 1 mapping quality.
 -e|--early-exit flag was triggered, will do early exit.
 --editor-repeat-coverage flag was triggered, threshold repeat coverage parameter set to 5.
 --editor-coarse-resolution flag was triggered, misjoin editor coarse matrix resolution set to 100000.
 --editor-coarse-region flag was triggered, misjoin editor coarse resolution depletion region size set to 500000.
 

# perl: warning: Setting locale failed.
# perl: warning: Please check that your locale settings:
        # LANGUAGE = (unset),
        # LC_ALL = (unset),
        # LANG = "C.UTF-8"
    # are supported and installed on your system.
# perl: warning: Falling back to the standard locale ("C").
# perl: warning: Setting locale failed.
# perl: warning: Please check that your locale settings:
        # LANGUAGE = (unset),
        # LC_ALL = (unset),
        # LANG = "C.UTF-8"
    # are supported and installed on your system.
# perl: warning: Falling back to the standard locale ("C").
# Not sure how to parse your input: files not listed or not found at expected locations. Exiting!


# try adding to the 3D DNA script
export LANGUAGE=en_US.UTF-8
export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8
export LC_CTYPE=en_US.UTF-8

#-------------------------------
# cp: cannot stat '/home/celphin/scratch/Hi-C_scaffold/Juicer/juicer/scaffolding/Oxydig_HAP1/3D_DNA/Ox
# y_draft_assembly_Hap1.fa0.assembly': No such file or directory
# cp: cannot stat '/home/celphin/scratch/Hi-C_scaffold/Juicer/juicer/scaffolding/Oxydig_HAP1/3D_DNA/Ox
# y_draft_assembly_Hap1.fa0.hic': No such file or directory

module load StdEnv/2020
module load nixpkgs/16.09 
module load bwa/0.7.15
module load java/1.8.0_192 
module load gcc/4.8.5
module load cuda/7.5.18
cd /home/celphin/scratch/Hi-C_scaffold/Juicer/juicer/scaffolding/Oxydig_HAP1

/scratch/celphin/Hi-C_scaffold/Juicer/juicer/3d-dna/run-asm-pipeline.sh -i 5000 \
-r 0 -q 1 --early-exit --editor-repeat-coverage 5 --editor-coarse-resolution 100000 \
--editor-coarse-region 500000 \
/home/celphin/scratch/Hi-C_scaffold/Juicer/juicer/references/Oxy_draft_assembly_Hap1.fa \
/home/celphin/scratch/Hi-C_scaffold/Juicer/juicer/scaffolding/Oxydig_HAP1/aligned/merged_nodups.txt

# missing
more /home/celphin/scratch/Hi-C_scaffold/Juicer/juicer/scaffolding/Oxydig_HAP1/aligned/merged_nodups.txt
# need to adjust the number of cpu to 46 in JUICER.sh

#--------------------------
# Issue running Juicer
# try running directly not in the script

tmux new-session -s Juicer
tmux attach-session -t Juicer

cd /home/celphin/scratch/Hi-C_scaffold

#eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"
#conda activate ASSEM_SCAFF2
module load StdEnv/2023 bwa/0.7.17 java/1.8.0_292 samtools/1.18 nixpkgs/16.09 gcc/5.4.0 cuda/8.0.44 

#ID=Oxydig_HAP1
ID=Oxydig_HAP2
#ASSEMBLY=Oxy_draft_assembly_Hap1.fa 
ASSEMBLY=Oxy_draft_assembly_Hap2.fa 
JUCIER_DIR=/home/celphin/scratch/Hi-C_scaffold/Juicer/juicer 
RESULT_DIR_ID=/home/celphin/scratch/Hi-C_scaffold/output/Oxydig
HIC_R1=Oxy_HiC_R1.fastq.gz 
HIC_R2=Oxy_HiC_R2.fastq.gz 


[ -d ${RESULT_DIR_ID}/LOG ] || mkdir -p ${RESULT_DIR_ID}/LOG

[ -d ${JUCIER_DIR}/scaffolding/${ID}/fastq ] || \
mkdir -p ${JUCIER_DIR}/scaffolding/${ID}/fastq

cp -R -u -p $HIC_R1 ${JUCIER_DIR}/scaffolding/${ID}/fastq

cp -R -u -p $HIC_R2 ${JUCIER_DIR}/scaffolding/${ID}/fastq

cd ${JUCIER_DIR}/scaffolding/${ID}

${JUCIER_DIR}/scripts/juicer.sh \
    -g $ID \
    -s MboI \
    -t 46 \
    -z ${JUCIER_DIR}/references/${ASSEMBLY##*/} \
    -y ${JUCIER_DIR}/restriction_sites/${ID}_MboI.txt \
    -p ${JUCIER_DIR}/chromosome_size/${ID}.chrom.sizes \
    -D $JUCIER_DIR 2> ${RESULT_DIR_ID}/LOG/${ID}_JUICER.log

# to check logs
cd /home/celphin/scratch/Hi-C_scaffold/Juicer/juicer/scaffolding/Oxydig_HAP2/debug/mergesort-34712045.out



#-------------------- 
# Error
#https://github.com/aidenlab/juicer/issues/9

in queue default to genome /home/celphin/scratch/Hi-C_scaffold/Juicer/juicer/references/Oxy_draft_a
ssembly_Hap2.fa with no fragment delimited maps.
(-: Created /home/celphin/scratch/Hi-C_scaffold/Juicer/juicer/scaffolding/Oxydig_HAP2/splits and /ho
me/celphin/scratch/Hi-C_scaffold/Juicer/juicer/scaffolding/Oxydig_HAP2/aligned.


***! Failure during chimera handling of /home/celphin/scratch/Hi-C_scaffold/Juicer/juicer/scaffolding/Oxydig_HAP1/splits/Oxy_HiC.fastq.gz

# need newer version of awk
# can replace with gawk
# i've replaced awk by gawk in juicer.sh line 423,


grep -n "***! Failure during chimera handling of" juicer.sh
#886:   echo "***! Failure during chimera handling of $name${ext}"

less -N juicer.sh

nano +423 juicer.sh

# try newer StdEnv/2023 in module loads at the start 
# check other differences
# mostly the same but for slurm calls
# https://www.diffchecker.com/text-compare/

# old
load_bwa="module load StdEnv/2020 bwa/0.7.17"
load_java="module load nixpkgs/16.09  java/1.8.0_192"
load_gpu="module load nixpkgs/16.09 gcc/4.8.5 cuda/7.5.18"
load_samtools="module load StdEnv/2020 samtools/1.16.1"

module load StdEnv/2023 bwa java samtools

load_bwa="module load StdEnv/2023 bwa/0.7.17"
load_java="module load StdEnv/2023  java/1.8.0_292"
load_gpu="module load nixpkgs/16.09 gcc/5.4.0 cuda/8.0.44"
load_samtools="module load StdEnv/2020  gcc/9.3.0 samtools/1.13"
# try juicer again with Hap2

#-------------------------

# worked?
rm -rf ${JUCIER_DIR}/scaffolding/${ID}/fastq

# Launch the 3D DNA script

#First round of 3D DNA scaffolding
cd /home/celphin/scratch/Hi-C_scaffold
sbatch /home/celphin/scratch/Hi-C_scaffold/haplotype_aware_scaffolding/3DDNA.sh \
$ASSEMBLY $ID $JUCIER_DIR ${MAPQ} ${RESULT_DIR_ID}



#---------------------
# wait here for scaffolding to finish

#mix hap
[ -d ${RESULT_DIR}/${SAMPLE}/sequences ] || mkdir -p ${RESULT_DIR}/${SAMPLE}/sequences

bash ${bin}/cat_fasta.sh $ASSEMBLY_hap1 $ASSEMBLY_hap2 ${RESULT_DIR}/${SAMPLE}/sequences

bash ${bin}/scaffolding.sh $bin ${ID_HAP1}_${ID_HAP2} \
${RESULT_DIR}/${SAMPLE}/sequences/$(basename ${ASSEMBLY_hap1%%.fasta})_$(basename ${ASSEMBLY_hap2%%.fasta}).fasta \
$JUCIER_DIR $HIC_R1 $HIC_R2 10 $RESULT_DIR_ID

#mini map alignment
[ -d ${RESULT_DIR}/${SAMPLE}/minimap ] || mkdir -p ${RESULT_DIR}/${SAMPLE}/minimap

sbatch ${bin}/minimap.sh $ASSEMBLY_hap1 $ASSEMBLY_hap2 ${RESULT_DIR}/${SAMPLE}/minimap

##############################
# try Kaede's code for Juicer
# https://github.com/ericgonzalezs/ASSEMBLIES/blob/main/Juicer 

## Juicer and 3D-DNA pipelne for genome scaffolding ## 
#1. Run juicer to produce .hic and .assembly
#2. Open the .hic file in Juicebox (can use cloud version for bigger files)
#3. 3d-dna first step keeps all contigs, 2nd and 3rd step polishes - removes contigs to trash (on debris. 
#4. For placing cnntigs in the right place, look for lines of red/white. If there is a thick white line, move contigs. 

#do this to run it as a pipeline based on Eric's script (https://github.com/ericgonzalezs/ASSEMBLIES/blob/main/Juicer)
#haplotype 1 
ln -s ../Juicer/juicer/CPU/ scripts
cd scripts/common
wget https://hicfiles.tc4ga.com/public/juicer/juicer_tools.1.9.9_jcuda.0.8.jar
ln -s juicer_tools.1.9.9_jcuda.0.8.jar  juicer_tools.jar
cd ../..
mkdir fastq
cd fastq
mv ../../Oxy_HiC_R1.fastq.gz Oxy_HiC_R1.fastq.gz
mv ../../Oxy_HiC_R2.fastq.gz Oxy_HiC_R2.fastq.gz
cd ../
mkdir references
cd references/
mv ../Oxy_draft_assembly_Hap1.fa   Oxy_draft_assembly_Hap1.fa
mv ../Oxy_draft_assembly_Hap2.fa   Oxy_draft_assembly_Hap2.fa

#-----------------------------
tmux new-session -s bwa
tmux attach-session -t bwa

cd /home/celphin/scratch/Hi-C_scaffold/Juicer_run/references
salloc -c1 --time 2:55:00 --mem 120000m --account def-rieseber

module load bwa
bwa index Oxy_draft_assembly_Hap1.fa #Run this in a job
bwa index Oxy_draft_assembly_Hap2.fa #Run this in a job

#--------------------------------
cd ..
mkdir restriction_sites 
cd restriction_sites
wget https://raw.githubusercontent.com/aidenlab/juicer/main/misc/generate_site_positions.py

#---- here download the juicer/misc/generate_site_positions.py and edit accordingly
  filenames = {
    'hg19': '/seq/references/Homo_sapiens_assembly19.fasta',
    'mm9' : '/seq/references/Mus_musculus_assembly9.fasta',
    'mm10': '/seq/references/Mus_musculus_assembly10.fasta',
    'hg18': '/seq/references/Homo_sapiens_assembly18.fasta',
    'Oxy_Hap1': '../references/Oxy_draft_assembly_Hap1.fa', 
    'Oxy_Hap2': '../references/Oxy_draft_assembly_Hap2.fa',#here you put your contig/assembly and its path
  }

tmux new-session -s python
tmux attach-session -t python

/home/celphin/scratch/Hi-C_scaffold/Juicer_run/restriction_sites
salloc -c1 --time 2:55:00 --mem 120000m --account def-rieseber

module load python
python generate_site_positions.py DpnII Oxy_Hap1 #Run this in a job
python generate_site_positions.py DpnII Oxy_Hap2 #Run this in a job

#--------------------------------
#generate a file Chromosome_sizes.sh with this:
for i in $(ls *_DpnII.txt)
do
name=$(echo $i | cut -d "." -f 1 )
awk 'BEGIN{OFS="\t"}{print $1, $NF}'  $i > "$name"".chrom.sizes"
done

#-----------------------
# tree looks like

tree
.
├── fastq
│   ├── Oxy_HiC_R1.fastq.gz
│   └── Oxy_HiC_R2.fastq.gz
├── juicer_run_Oxyria_Hap1.sh
├── juicer_run_Oxyria_Hap2.sh
├── references
│   ├── Oxy_draft_assembly_Hap1.fa
│   ├── Oxy_draft_assembly_Hap1.fa.amb
│   ├── Oxy_draft_assembly_Hap1.fa.ann
│   ├── Oxy_draft_assembly_Hap1.fa.bwt
│   ├── Oxy_draft_assembly_Hap1.fa.pac
│   ├── Oxy_draft_assembly_Hap1.fa.sa
│   ├── Oxy_draft_assembly_Hap2.fa
│   ├── Oxy_draft_assembly_Hap2.fa.amb
│   ├── Oxy_draft_assembly_Hap2.fa.ann
│   ├── Oxy_draft_assembly_Hap2.fa.bwt
│   ├── Oxy_draft_assembly_Hap2.fa.pac
│   └── Oxy_draft_assembly_Hap2.fa.sa
├── restriction_sites
│   ├── Oxy_Hap1_DpnII.chrom.sizes
│   ├── Oxy_Hap1_DpnII.txt
│   ├── Oxy_Hap2_DpnII.chrom.sizes
│   ├── Oxy_Hap2_DpnII.txt
│   └── generate_site_positions.py
└── scripts -> ../Juicer/juicer/CPU/

#----------------------------
# run Juicer 

#run juicier with out GPUS and more cores
nano juicer_run_Oxyria_Hap1.sh

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --account=def-rieseber
#SBATCH --time=3-0
#SBATCH --ntasks-per-node=32
#SBATCH --mem=149G
module load StdEnv/2020 bwa/0.7.17 java/17.0.2 samtools/1.15.1
cd /home/celphin/scratch/Hi-C_scaffold/Juicer_run
bash scripts/juicer.sh -D $PWD -g Oxy_Hap1 -s DpnII -p restriction_sites/Oxy_Hap1_DpnII.chrom.sizes -y restriction_sites/Oxy_Hap1_DpnII.txt -z references/Oxy_draft_assembly_Hap1.fa -t 32 -S early

sbatch juicer_run_Oxyria_Hap1.sh
(-:  Align of /home/celphin/scratch/Hi-C_scaffold/Juicer_run/splits/Oxy_HiC.fastq.gz.sam done successfully

#------------------
nano juicer_run_Oxyria_Hap2.sh

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --account=def-rieseber
#SBATCH --time=3-0
#SBATCH --ntasks-per-node=32
#SBATCH --mem=149G
module load StdEnv/2020 bwa/0.7.17 java/17.0.2 samtools/1.15.1
cd /home/celphin/scratch/Hi-C_scaffold/Juicer_run
bash scripts/juicer.sh -D $PWD -g Oxy_Hap2 -s DpnII -p restriction_sites/Oxy_Hap2_DpnII.chrom.sizes -y restriction_sites/Oxy_Hap2_DpnII.txt -z references/Oxy_draft_assembly_Hap2.fa -t 32 -S early

sbatch juicer_run_Oxyria_Hap2.sh

# worried that the files overlap between haplotypes

##################################
# try running again in separate folders
cd /home/celphin/scratch/Hi-C_scaffold/Juicer_run
cp -rv fastq ../Oxy_Hap1_juicer
cp -rv references/ ../Oxy_Hap1_juicer
cp -rv restriction_sites ../Oxy_Hap1_juicer
cp -rv scripts ../Oxy_Hap1_juicer
cd ..
cp -rv Oxy_Hap1_juicer Oxy_Hap2_juicer

#-------------------
cd /home/celphin/scratch/Hi-C_scaffold/Oxy_Hap1_juicer
nano juicer_run_Oxyria_Hap1.sh

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --account=def-rieseber
#SBATCH --time=0-22:00:00
#SBATCH --ntasks-per-node=32
#SBATCH --mem=149G
module load StdEnv/2020 bwa/0.7.17 java/17.0.2 samtools/1.15.1
cd /home/celphin/scratch/Hi-C_scaffold/Oxy_Hap1_juicer
bash scripts/juicer.sh -D $PWD -g Oxy_Hap1 -s DpnII -p restriction_sites/Oxy_Hap1_DpnII.chrom.sizes -y restriction_sites/Oxy_Hap1_DpnII.txt -z references/Oxy_draft_assembly_Hap1.fa -t 32 -S early

sbatch juicer_run_Oxyria_Hap1.sh

#------------------
cd /home/celphin/scratch/Hi-C_scaffold/Oxy_Hap2_juicer
nano juicer_run_Oxyria_Hap2.sh

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --account=def-rieseber
#SBATCH --time=0-22:00:00
#SBATCH --ntasks-per-node=32
#SBATCH --mem=149G
module load StdEnv/2020 bwa/0.7.17 java/17.0.2 samtools/1.15.1
cd /home/celphin/scratch/Hi-C_scaffold/Oxy_Hap2_juicer
bash scripts/juicer.sh -D $PWD -g Oxy_Hap2 -s DpnII -p restriction_sites/Oxy_Hap2_DpnII.chrom.sizes -y restriction_sites/Oxy_Hap2_DpnII.txt -z references/Oxy_draft_assembly_Hap2.fa -t 32 -S early

sbatch juicer_run_Oxyria_Hap2.sh




########################################
# try 3D DNA
# https://github.com/ericgonzalezs/ASSEMBLIES/blob/main/3D-DNA.sh

cd /home/celphin/scratch/Hi-C_scaffold/
git clone https://github.com/aidenlab/3d-dna.git
cd 3d-dna
chmod -R 770 * # give execute permissions


#create python env
module load StdEnv/2020 python/3.11.2

virtualenv 3ddna

source 3ddna/bin/activate

pip install scipy numpy matplotlib #libraries required for 3d-dna 

deactivate

#----------------------
cd /home/celphin/scratch/Hi-C_scaffold/3d-dna
# copy over files as symbolic links
ln -s /home/celphin/scratch/Hi-C_scaffold/Oxy_Hap1_juicer/references/Oxy_draft_assembly_Hap1.fa
ln -s /home/celphin/scratch/Hi-C_scaffold/Oxy_Hap1_juicer/references/Oxy_draft_assembly_Hap2.fa

ln -s /home/celphin/scratch/Hi-C_scaffold/Oxy_Hap1_juicer/aligned/merged_nodups_Oxy_Hap1.txt 
ln -s /home/celphin/scratch/Hi-C_scaffold/Oxy_Hap2_juicer/aligned/merged_nodups_Oxy_Hap2.txt

#----------------
# run 3DDNA

nano 3DDNA_Oxy_Hap1.sh

#!/bin/bash
#SBATCH --account=def-rieseber
#SBATCH --time=1-22:00:00
#SBATCH --cpus-per-task=15
#SBATCH --mem=100G

module load StdEnv/2020 python/3.10.2 java/17.0.2 lastz/1.04.03
export PATH="/home/celphin/scratch/Hi-C_scaffold/3d-dna:$PATH"
source /home/celphin/scratch/Hi-C_scaffold/3d-dna/3ddna/bin/activate

cd /home/celphin/scratch/Hi-C_scaffold/3d-dna/Oxy_hap1
/home/celphin/scratch/Hi-C_scaffold/3d-dna/run-asm-pipeline.sh -r 2 Oxy_draft_assembly_Hap1.fa merged_nodups_Oxy_Hap1.txt

deactivate

sbatch 3DDNA_Oxy_Hap1.sh

#----------------------------
nano 3DDNA_Oxy_Hap2.sh

#!/bin/bash
#SBATCH --account=def-rieseber
#SBATCH --time=1-22:00:00
#SBATCH --cpus-per-task=15
#SBATCH --mem=100G

module load StdEnv/2020 python/3.10.2 java/17.0.2 lastz/1.04.03
export PATH="/home/celphin/scratch/Hi-C_scaffold/3d-dna:$PATH"
source /home/celphin/scratch/Hi-C_scaffold/3d-dna/3ddna/bin/activate

cd /home/celphin/scratch/Hi-C_scaffold/3d-dna/Oxy_hap2
/home/celphin/scratch/Hi-C_scaffold/3d-dna/run-asm-pipeline.sh -r 2 Oxy_draft_assembly_Hap2.fa merged_nodups_Oxy_Hap2.txt

deactivate

sbatch 3DDNA_Oxy_Hap2.sh

#------------------------------
# taking a long time
Finished writing norms
:| Warning: No input for label1 was provided. Default for label1 is ":::fragment_"
:| Warning: No input for label2 was provided. Default for label2 is ":::debris"

###########################
# Finalize output

nano 3DDNA_Oxy_Hap1_final.sh

#!/bin/bash
#SBATCH --account=def-rieseber
#SBATCH --time=0-22:00:00
#SBATCH --cpus-per-task=15
#SBATCH --mem=100G

module load StdEnv/2020 python/3.10.2 java/17.0.2 lastz/1.04.03
export PATH="/home/celphin/scratch/Hi-C_scaffold/3d-dna:$PATH"
source /home/celphin/scratch/Hi-C_scaffold/3d-dna/3ddna/bin/activate

run-asm-pipeline.sh -r 0 --stage finalize Oxy_draft_assembly_Hap1.fa merged_nodups_Oxy_Hap1.txt 

deactivate

sbatch 3DDNA_Oxy_Hap1_final.sh

#-----------------------------
nano 3DDNA_Oxy_Hap2_final.sh

#!/bin/bash
#SBATCH --account=def-rieseber
#SBATCH --time=1-22:00:00
#SBATCH --cpus-per-task=15
#SBATCH --mem=100G

module load StdEnv/2020 python/3.10.2 java/17.0.2 lastz/1.04.03
export PATH="/home/celphin/scratch/Hi-C_scaffold/3d-dna:$PATH"
source /home/celphin/scratch/Hi-C_scaffold/3d-dna/3ddna/bin/activate

run-asm-pipeline.sh -r 0 --stage seal Oxy_draft_assembly_Hap2.fa merged_nodups_Oxy_Hap2.txt 

deactivate

sbatch 3DDNA_Oxy_Hap2_final.sh

#######################################
# check assembly stats

module load StdEnv/2020 
module load gcc/9.3.0
module load r-bundle-bioconductor/3.14
cd ~/assembly-stats/build

./assembly-stats /home/celphin/scratch/Hi-C_scaffold/3d-dna/Oxy_hap1/Oxy_draft_assembly_Hap1.FINAL.fasta

stats for /home/celphin/scratch/Hi-C_scaffold/3d-dna/Oxy_hap1/Oxy_draft_assembly_Hap1.FINAL.fasta
sum = 582 998 305, n = 1886, ave = 309118.93, largest = 79 385 870
N50 = 74 109 000, n = 4
N60 = 70 648 606, n = 5
N70 = 69 949 313, n = 6
N80 = 67 411 532, n = 7

N90 = 100188, n = 65
N100 = 1000, n = 1886
N_count = 382000
Gaps = 764

#----------------------------
./assembly-stats /home/celphin/scratch/Hi-C_scaffold/3d-dna/Oxy_hap2/Oxy_draft_assembly_Hap2.FINAL.fasta

stats for /home/celphin/scratch/Hi-C_scaffold/3d-dna/Oxy_hap2/Oxy_draft_assembly_Hap2.FINAL.fasta
sum = 566 728 759, n = 615, ave = 921510.18, largest = 86 137 797
N50 = 74 731 636, n = 4
N60 = 72 720 770, n = 5
N70 = 72 112 163, n = 6
N80 = 72 112 163, n = 6
N90 = 70 815 000, n = 7

N100 = 1000, n = 615
N_count = 134000
Gaps = 268

#################################################