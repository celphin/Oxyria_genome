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

# Setup
# https://github.com/theaidenlab/juicer/wiki/Running-Juicer-on-a-cluster
cd /home/celphin/projects/def-henryg/celphin/Oxyria/
mkdir Juicer
cd Juicer
mkdir references/

# copy over reference fasta
cp -v /home/celphin/projects/def-henryg/celphin/Oxyria/Oxy_dig_1.scaffolds_FINAL.fasta ./references/

# Get Juicer
git clone https://github.com/theaidenlab/juicer.git

# symbolic link to scripts
ln -s juicer/SLURM/scripts/ scripts
cd scripts
# Download Juicer tools jar
wget https://github.com/aidenlab/Juicebox/releases/download/v2.18.00/juicer_tools.2.18.00.jar
ln -s juicer_tools.2.18.00.jar juicer_tools.jar

cd ..
mkdir fastq; cd fastq
# copy over the HiC files into fastq
cp -v /home/celphin/projects/def-henryg/celphin/Oxyria/Hi-C_raw_data/NS.1871.003.IDT_i7_217---IDT_i5_217.Oxyria_Hi-C_R*.fastq.gz .
cd ../..

#---------------------
# Add genome to Juicer
# https://github.com/aidenlab/juicer/wiki/Usage

tmux new-session -s Juicer # Ctrl+B D to dettach
tmux attach-session -t Juicer

# get allocation 
salloc -c1 --time 2:50:00 --mem 120000m --account def-rieseber

cd /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer

module load StdEnv/2020
module load nixpkgs/16.09 
module load bwa/0.7.15
module load java/1.8.0_192 
module load gcc/4.8.5
module load cuda/7.5.18

# index fasta 
bwa index ./references/Oxy_dig_1.scaffolds_FINAL.fasta
# [main] Real time: 725.084 sec; CPU: 711.967 sec

#------------------------
# add restriction sites to reference

# copy over python script
cp -v /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer/juicer/miscgenerate_site_positions.py /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer/scripts

module load StdEnv/2020
module load python
python ./scripts/generate_site_positions.py DpnII Oxyria '/project/6064374/celphin/Oxyria/Juicer/references/Oxy_dig_1.scaffolds_FINAL.fasta'
mkdir restriction_sites
mv Oxyria_DpnII.txt restriction_sites

#---------------------
# get chromosome sizes list
awk 'BEGIN{OFS="\t"}{print $1, $NF}' ./restriction_sites/Oxyria_DpnII.txt > ./restriction_sites/Oxyria.chrom.sizes

##############################
# try running

module load StdEnv/2020
module load nixpkgs/16.09 
module load bwa/0.7.15
module load java/1.8.0_192 
module load gcc/4.8.5
module load cuda/7.5.18
cd /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer

# run juicer
./scripts/juicer.sh -z ./references/Oxy_dig_1.scaffolds_FINAL.fasta -p ./restriction_sites/Oxyria.chrom.sizes -y ./restriction_sites/Oxyria_DpnII.txt

sbatch: error: Batch job submission failed: Unspecified error
sbatch: error: ----------------------------------------
sbatch: error: You are associated with multiple _cpu allocations...
sbatch: error: Please specify one of the following accounts to submit this job:
sbatch: error:   RAS default accounts: def-cronk, def-henryg, def-rieseber,
sbatch: error:           RAC accounts: rpp-rieseber,
sbatch: error: Compute-Burst accounts:
sbatch: error:         Other accounts:
sbatch: error: Use the parameter --account=desired_account when submitting your job
sbatch: error: ----------------------------------------

# set default SLURM submission account 
# https://docs.alliancecan.ca/wiki/Running_jobs#Accounts_and_projects
export SLURM_ACCOUNT=def-rieseber
export SBATCH_ACCOUNT=$SLURM_ACCOUNT
export SALLOC_ACCOUNT=$SLURM_ACCOUNT

#------------------------
# try again
./scripts/juicer.sh -z ./references/Oxy_dig_1.scaffolds_FINAL.fasta -p ./restriction_sites/Oxyria.chrom.sizes -y ./restriction_sites/Oxyria_DpnII.txt 

sbatch: error: ----------------------------------------
sbatch: error: The specified partition does not exist, or the submitted job cannot fit in it...
sbatch: error: Please specify a different partition, or simply submit the job without the --partition option,
sbatch: error: the scheduler will redirect it to the most suitable partition automatically
sbatch: error: ----------------------------------------
sbatch: error: Batch job submission failed: Unspecified error

# https://docs.alliancecan.ca/wiki/Running_jobs#Do_not_specify_a_partition
--partition=default

#-----------------------
# try again setting the partitions to default
./scripts/juicer.sh -z ./references/Oxy_dig_1.scaffolds_FINAL.fasta -p ./restriction_sites/Oxyria.chrom.sizes -y ./restriction_sites/Oxyria_DpnII.txt -q default -l default

# Works
(-: Looking for fastq files...fastq files exist
(-: Aligning files matching /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer/fastq/*_R*.fastq*
 in queue default to genome ./references/Oxy_dig_1.scaffolds_FINAL.fasta with no fragment delimited maps.
(-: Created /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer/splits and /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer/aligned.
(-: Starting job to launch other jobs once splitting is complete
sbatch: NOTE: Your memory request of 2048M was likely submitted as 2G. Please note that Slurm interprets memory requests denominated in G as multiples of 1024M, not 1000M.
sbatch: NOTE: Your memory request of 65536M was likely submitted as 64G. Please note that Slurm interprets memory requests denominated in G as multiples of 1024M, not 1000M.
sbatch: NOTE: Your memory request of 2048M was likely submitted as 2G. Please note that Slurm interprets memory requests denominated in G as multiples of 1024M, not 1000M.
sbatch: NOTE: Your memory request of 8192M was likely submitted as 8G. Please note that Slurm interprets memory requests denominated in G as multiples of 1024M, not 1000M.
sbatch: NOTE: Your memory request of 1024M was likely submitted as 1G. Please note that Slurm interprets memory requests denominated in G as multiples of 1024M, not 1000M.
sbatch: NOTE: Your memory request of 4096M was likely submitted as 4G. Please note that Slurm interprets memory requests denominated in G as multiples of 1024M, not 1000M.
sbatch: NOTE: Your memory request of 8192M was likely submitted as 8G. Please note that Slurm interprets memory requests denominated in G as multiples of 1024M, not 1000M.
sbatch: NOTE: Your memory request of 2048M was likely submitted as 2G. Please note that Slurm interprets memory requests denominated in G as multiples of 1024M, not 1000M.
(-: Finished adding all jobs... Now is a good time to get that cup of coffee... Last job id 60717442


#-------------------------
sq
          JOBID     USER      ACCOUNT           NAME  ST  TIME_LEFT NODES CPUS TRES_PER_N MIN_MEM NODELIST (REASON)
       60710509  celphin def-rieseber    interactive   R    1:08:50     1   32        N/A 120000M cdr1308 (None)
       60717187  celphin def-rieseber a1677103668_cm  PD       2:00     1    1        N/A    256M  (Priority)
       60717229  celphin def-rieseber a1677103668_NS  PD 2-00:00:00     1    1        N/A      5G  (Priority)
       60717261  celphin def-rieseber a1677103668_al  PD 2-00:00:00     1    8        N/A  62.50G  (Nodes required for job are DOWN, DRAINED or reserved for jobs in higher priority partitions)
       60717295  celphin def-rieseber a1677103668_me  PD 5-00:00:00     1    1        N/A     10G  (Dependency)
       60717305  celphin def-rieseber a1677103668_me  PD 5-00:00:00     1    1        N/A     10G  (Dependency)
       60717312  celphin def-rieseber a1677103668_me  PD 5-00:00:00     1    8        N/A      2G  (Dependency)
       60717324  celphin def-rieseber a1677103668_ch  PD 2-00:00:00     1    1        N/A    256M  (Dependency)
       60717331  celphin def-rieseber a1677103668_fr  PD 7-00:00:00     1    8        N/A     64G  (Dependency)
       60717334  celphin def-rieseber a1677103668_de  PD      10:00     1    1        N/A    256M  (JobHeldUser)
       60717342  celphin def-rieseber a1677103668_de  PD 5-00:00:00     1    1        N/A      2G  (Dependency)
       60717347  celphin def-rieseber a1677103668_po  PD    1:40:00     1    1        N/A    256M  (Dependency)
       60717352  celphin def-rieseber a1677103668_du  PD 2-00:00:00     1    1        N/A      8G  (Dependency)
       60717355  celphin def-rieseber a1677103668_me  PD 2-00:00:00     1    8        N/A     10G  (Dependency)
       60717361  celphin def-rieseber a1677103668_me  PD 2-00:00:00     1    8        N/A     10G  (Dependency)
       60717369  celphin def-rieseber a1677103668_pr  PD 2-00:00:00     1    8        N/A      1G  (Dependency)
       60717387  celphin def-rieseber a1677103668_ba  PD 2-00:00:00     1    8        N/A     10G  (Dependency)
       60717399  celphin def-rieseber a1677103668_st  PD 5-00:00:00     1    1        N/A     25G  (Dependency)
       60717404  celphin def-rieseber a1677103668_st  PD 5-00:00:00     1    1        N/A     25G  (Dependency)
       60717428  celphin def-rieseber a1677103668_hi  PD 5-00:00:00     1    1        N/A    150G  (Dependency)
       60717435  celphin def-rieseber a1677103668_hi  PD 5-00:00:00     1    1        N/A    150G  (Dependency)
       60717438  celphin def-rieseber a1677103668_hi  PD 2-00:00:00     1    1        N/A      4G  (Dependency)
       60717440  celphin def-rieseber a1677103668_ar  PD 2-00:00:00     1    1        N/A      8G  (Dependency)
       60717442  celphin def-rieseber a1677103668_pr  PD   20:00:00     1    1        N/A      2G  (Dependency)


#---------------------------
# module load error

/var/spool/slurmd/job60717261/slurm_script: spack: command not found
/var/spool/slurmd/job60717261/slurm_script: spack: command not found
mem: invalid option -- '5'

# https://groups.google.com/g/3d-genomics/c/GxihVHgCqD4?pli=1 

nano ./scripts/juicer.sh

# change start to:
# from:
#---- 
juicer_version="2.0"
## Set the following variables to work with your system
## path additionals, make sure paths are correct for your system
## use cluster load commands
load_bwa="module load  StdEnv/2020 bwa/0.7.17"
load_java="module load nixpkgs/16.09  java/1.8.0_192"
load_gpu="module load nixpkgs/16.09 gcc/4.8.5 cuda/7.5.18"
load_samtools="module load StdEnv/2020 samtools/1.16.1"

# Juicer directory, contains scripts/, references/, and restriction_sites/
# can also be set in options via -D
juiceDir="/home/celphin/projects/def-henryg/celphin/Oxyria/Juicer" ### RICE
# default queue, can also be set in options via -q
queue="default"
queue_time="23:50:00"
# default long queue, can also be set in options via -l
long_queue="default"
long_queue_time="2-00:00:00"
#---- 
# to start sizes
# size to split fastqs. adjust to match your needs. 4000000=1M reads per split

#------------------
# also change 
#SBATCH --mem=$alloc_mem
# to
#SBATCH --mem=${alloc_mem}m

#---------------------
# also line 370

if [ -z "$threads" ]
then
	threads=8
	sortthreads=8
	threadstring="-t $threads"
	sthreadstring="-@ $threads"
else
	threadstring="-t $threads"
	sthreadstring="-@ $threads"
	sortthreads=$threads
fi

#-----------------
# remove lines 415

if [ $isBCM -eq 1 ] || [ $isRice -eq 1 ]
then
     alloc_mem=50000
fi

#--------------------
# remove 
./scripts/juicer.sh: line 484: [: -ne: unary operator expected
./scripts/juicer.sh: line 950: [: -eq: unary operator expected
./scripts/juicer.sh: line 956: [: -eq: unary operator expected


    484 if [ $isVoltron -ne 1 ]
    485 then
    486     if [ -z $splitme ]
    487     then
    488         fastqsize=$(ls -lgGL  ${fastqdir} | awk '{sum+=$3}END{print sum}')
    489         if [ "$fastqsize" -gt "2592410750" ]
    490         then
    491             splitme=1
    492         fi
    493     fi
    494 fi

#----------------------------
# remove

950     if [ $isVoltron -eq 1 ]
    951     then
    952         sbatch_time="#SBATCH -t 10080"
    953     else
    954         sbatch_time="#SBATCH -t 1440"
    955     fi
    956     if [ $isBCM -eq 1 ]
    957     then
    958         sbatch_cpu_alloc="#SBATCH -c 1"
    959         sbatch_mem_alloc="#SBATCH --mem=80G"
    960     else
    961         sbatch_cpu_alloc="#SBATCH -c 8"
    962         sbatch_mem_alloc="#SBATCH --mem=64G"
    963     fi

# replace with
sbatch_cpu_alloc="#SBATCH -c 8"
sbatch_mem_alloc="#SBATCH --mem=80G"

#--------------------------
# update threads to 32 instead of 8
# change back  - cannot do more than 8?


#------------------------
# try again with added -D to main directory

module load StdEnv/2020
module load nixpkgs/16.09 
module load bwa/0.7.15
module load java/1.8.0_192 
module load gcc/4.8.5
module load cuda/7.5.18
cd /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer

./scripts/juicer.sh \
-D /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer/ \
-z ./references/Oxy_dig_1.scaffolds_FINAL.fasta \
-p ./restriction_sites/Oxyria.chrom.sizes \
-y ./restriction_sites/Oxyria_DpnII.txt


(-: Looking for fastq files...fastq files exist
(-: Aligning files matching /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer/fastq/*_R*.fastq*
 in queue default to genome ./references/Oxy_dig_1.scaffolds_FINAL.fasta with no fragment delimited maps.
(-: Created /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer/splits and /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer/aligned.
(-: Starting job to launch other jobs once splitting is complete
sbatch: NOTE: Your memory request of 2048M was likely submitted as 2G. Please note that Slurm interprets memory requests denominated in G as multiples of 1024M, not 1000M.
sbatch: NOTE: Your memory request of 2048M was likely submitted as 2G. Please note that Slurm interprets memory requests denominated in G as multiples of 1024M, not 1000M.
sbatch: NOTE: Your memory request of 8192M was likely submitted as 8G. Please note that Slurm interprets memory requests denominated in G as multiples of 1024M, not 1000M.
sbatch: NOTE: Your memory request of 1024M was likely submitted as 1G. Please note that Slurm interprets memory requests denominated in G as multiples of 1024M, not 1000M.
sbatch: NOTE: Your memory request of 8192M was likely submitted as 8G. Please note that Slurm interprets memory requests denominated in G as multiples of 1024M, not 1000M.
sbatch: NOTE: Your memory request of 2048M was likely submitted as 2G. Please note that Slurm interprets memory requests denominated in G as multiples of 1024M, not 1000M.
(-: Finished adding all jobs... Now is a good time to get that cup of coffee... Last job id 60732560


sq
          JOBID     USER      ACCOUNT           NAME  ST  TIME_LEFT NODES CPUS TRES_PER_N MIN_MEM NODELIST (REASON)
       60722892  celphin def-rieseber    interactive   R    1:03:01     1   32        N/A 120000M cdr1594 (None)
       60732536  celphin def-rieseber a1677113991_cm  PD       2:00     1    1        N/A    256M  (Priority)
       60732537  celphin def-rieseber a1677113991_NS  PD   23:50:00     1    1        N/A      5G  (Priority)
       60732538  celphin def-rieseber a1677113991_al  PD   23:50:00     1    8        N/A  62.50G  (Priority)
       60732539  celphin def-rieseber a1677113991_me  PD 2-00:00:00     1    1        N/A     10G  (Dependency)
       60732540  celphin def-rieseber a1677113991_me  PD 2-00:00:00     1    1        N/A     10G  (Dependency)
       60732541  celphin def-rieseber a1677113991_me  PD 2-00:00:00     1    8        N/A      2G  (Dependency)
       60732542  celphin def-rieseber a1677113991_ch  PD   23:50:00     1    1        N/A    256M  (Dependency)
       60732543  celphin def-rieseber a1677113991_fr  PD    1:00:00     1    8        N/A     80G  (Dependency)
       60732544  celphin def-rieseber a1677113991_de  PD      10:00     1    1        N/A    256M  (JobHeldUser)
       60732545  celphin def-rieseber a1677113991_de  PD 2-00:00:00     1    1        N/A      2G  (Dependency)
       60732546  celphin def-rieseber a1677113991_po  PD    1:40:00     1    1        N/A    256M  (Dependency)
       60732547  celphin def-rieseber a1677113991_du  PD   23:50:00     1    1        N/A      8G  (Dependency)
       60732548  celphin def-rieseber a1677113991_me  PD   23:50:00     1    8        N/A     10G  (Dependency)
       60732549  celphin def-rieseber a1677113991_me  PD   23:50:00     1    8        N/A     10G  (Dependency)
       60732552  celphin def-rieseber a1677113991_pr  PD   23:50:00     1    8        N/A      1G  (Dependency)
       60732553  celphin def-rieseber a1677113991_ba  PD   23:50:00     1    8        N/A     10G  (Dependency)
       60732554  celphin def-rieseber a1677113991_st  PD 2-00:00:00     1    1        N/A     25G  (Dependency)
       60732556  celphin def-rieseber a1677113991_st  PD 2-00:00:00     1    1        N/A     25G  (Dependency)
       60732557  celphin def-rieseber a1677113991_hi  PD 2-00:00:00     1    1        N/A    150G  (Dependency)
       60732558  celphin def-rieseber a1677113991_hi  PD 2-00:00:00     1    1        N/A    150G  (Dependency)
       60732559  celphin def-rieseber a1677113991_ar  PD   23:50:00     1    1        N/A      8G  (Dependency)
       60732560  celphin def-rieseber a1677113991_pr  PD   20:00:00     1    1        N/A      2G  (Dependency)



#######################################

# add samtools modules to juicer.sh
# try to resume
# Relaunch via the same script. Type juicer.sh [options] -S stage where "stage" is one of chimeric, merge, dedup, final, postproc, or early. 
"chimeric" is when alignment is done but not chimera read handling. See below for more information. 
"merge" is for when alignment has finished but merged_sort hasn't been created; 
"dedup" is for when merged_sort is there but not merged_nodups (this will relaunch all dedup jobs); 
"final" is for when merged_nodups is there and you want the stats and hic files; 
"postproc" is for when you have the hic files and just want feature annotations; 
"early" is for early exit, before hic file creation. You can also assign early exit via the "-e" flag and start at one of the other stages. 
If your jobs failed at the alignment stage, run relaunch_prep.sh and then run juicer.sh.

module load StdEnv/2020
module load nixpkgs/16.09 
module load bwa/0.7.15
module load java/1.8.0_192 
module load gcc/4.8.5
module load cuda/7.5.18
cd /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer

export SLURM_ACCOUNT=def-rieseber
export SBATCH_ACCOUNT=$SLURM_ACCOUNT
export SALLOC_ACCOUNT=$SLURM_ACCOUNT

./scripts/juicer.sh \
-D /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer/ \
-z ./references/Oxy_dig_1.scaffolds_FINAL.fasta \
-p ./restriction_sites/Oxyria.chrom.sizes \
-y ./restriction_sites/Oxyria_DpnII.txt \
-S chimeric

(-: Aligning files matching /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer/fastq/*_R*.fastq*
 in queue default to genome ./references/Oxy_dig_1.scaffolds_FINAL.fasta with no fragment delimited maps.
---  Using already created files in /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer/splits

(-: Starting job to launch other jobs once splitting is complete
sbatch: NOTE: Your memory request of 2048M was likely submitted as 2G. Please note that Slurm interprets memory requests denominated in G as multiples of 102
4M, not 1000M.
sbatch: NOTE: Your memory request of 2048M was likely submitted as 2G. Please note that Slurm interprets memory requests denominated in G as multiples of 102
4M, not 1000M.
sbatch: NOTE: Your memory request of 8192M was likely submitted as 8G. Please note that Slurm interprets memory requests denominated in G as multiples of 102
4M, not 1000M.
sbatch: NOTE: Your memory request of 1024M was likely submitted as 1G. Please note that Slurm interprets memory requests denominated in G as multiples of 102
4M, not 1000M.
sbatch: NOTE: Your memory request of 8192M was likely submitted as 8G. Please note that Slurm interprets memory requests denominated in G as multiples of 102
4M, not 1000M.
sbatch: NOTE: Your memory request of 2048M was likely submitted as 2G. Please note that Slurm interprets memory requests denominated in G as multiples of 102
4M, not 1000M.
(-: Finished adding all jobs... Now is a good time to get that cup of coffee... Last job id 60765337


          JOBID     USER      ACCOUNT           NAME  ST  TIME_LEFT NODES CPUS TRES_PER_N MIN_MEM NODELIST (REASON)
       60765152  celphin def-rieseber a1677137803_me  PD 2-00:00:00     1    1        N/A     10G  (Dependency)
       60765208  celphin def-rieseber a1677137803_me  PD 2-00:00:00     1    1        N/A     10G  (Dependency)
       60765218  celphin def-rieseber a1677137803_me  PD 2-00:00:00     1    8        N/A      2G  (Dependency)
       60765224  celphin def-rieseber a1677137803_ch  PD   23:50:00     1    1        N/A    256M  (Dependency)
       60765230  celphin def-rieseber a1677137803_fr  PD    1:00:00     1    8        N/A     80G  (Dependency)
       60765257  celphin def-rieseber a1677137803_de  PD      10:00     1    1        N/A    256M  (JobHeldUser)
       60765275  celphin def-rieseber a1677137803_de  PD 2-00:00:00     1    1        N/A      2G  (Dependency)
       60765281  celphin def-rieseber a1677137803_po  PD    1:40:00     1    1        N/A    256M  (Dependency)
       60765284  celphin def-rieseber a1677137803_du  PD   23:50:00     1    1        N/A      8G  (Dependency)
       60765305  celphin def-rieseber a1677137803_me  PD   23:50:00     1    8        N/A     10G  (Dependency)
       60765311  celphin def-rieseber a1677137803_me  PD   23:50:00     1    8        N/A     10G  (Dependency)
       60765317  celphin def-rieseber a1677137803_pr  PD   23:50:00     1    8        N/A      1G  (Dependency)
       60765321  celphin def-rieseber a1677137803_ba  PD   23:50:00     1    8        N/A     10G  (Dependency)
       60765324  celphin def-rieseber a1677137803_st  PD 2-00:00:00     1    1        N/A     25G  (Dependency)
       60765327  celphin def-rieseber a1677137803_st  PD 2-00:00:00     1    1        N/A     25G  (Dependency)
       60765328  celphin def-rieseber a1677137803_hi  PD 2-00:00:00     1    1        N/A    150G  (Dependency)
       60765330  celphin def-rieseber a1677137803_hi  PD 2-00:00:00     1    1        N/A    150G  (Dependency)
       60765332  celphin def-rieseber a1677137803_ar  PD   23:50:00     1    1        N/A      8G  (Dependency)
       60765337  celphin def-rieseber a1677137803_pr  PD   20:00:00     1    1        N/A      2G  (Dependency)
       60765146  celphin def-rieseber a1677137803_NS  CG   23:49:48     1    1        N/A      5G cdr360 (None)

# did not find the correct files

########################
# try again from beginning

module load StdEnv/2020
module load nixpkgs/16.09 
module load bwa/0.7.15
module load java/1.8.0_192 
module load gcc/4.8.5
module load cuda/7.5.18
cd /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer

./scripts/juicer.sh \
-D /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer/ \
-z ./references/Oxy_dig_1.scaffolds_FINAL.fasta \
-p ./restriction_sites/Oxyria.chrom.sizes \
-y ./restriction_sites/Oxyria_DpnII.txt


(-: Looking for fastq files...fastq files exist
(-: Aligning files matching /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer/fastq/*_R*.fastq*
 in queue default to genome ./references/Oxy_dig_1.scaffolds_FINAL.fasta with no fragment delimited maps.
(-: Created /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer/splits and /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer/aligned.
(-: Starting job to launch other jobs once splitting is complete
sbatch: NOTE: Your memory request of 2048M was likely submitted as 2G. Please note that Slurm interprets memory requests denominated in G as multiples of 1024M, not 1000M.
sbatch: NOTE: Your memory request of 2048M was likely submitted as 2G. Please note that Slurm interprets memory requests denominated in G as multiples of 1024M, not 1000M.
sbatch: NOTE: Your memory request of 8192M was likely submitted as 8G. Please note that Slurm interprets memory requests denominated in G as multiples of 1024M, not 1000M.
sbatch: NOTE: Your memory request of 1024M was likely submitted as 1G. Please note that Slurm interprets memory requests denominated in G as multiples of 1024M, not 1000M.
sbatch: NOTE: Your memory request of 8192M was likely submitted as 8G. Please note that Slurm interprets memory requests denominated in G as multiples of 1024M, not 1000M.
sbatch: NOTE: Your memory request of 2048M was likely submitted as 2G. Please note that Slurm interprets memory requests denominated in G as multiples of 1024M, not 1000M.
(-: Finished adding all jobs... Now is a good time to get that cup of coffee... Last job id 60768160


          JOBID     USER      ACCOUNT           NAME  ST  TIME_LEFT NODES CPUS TRES_PER_N MIN_MEM NODELIST (REASON)
       60768139  celphin def-rieseber a1677140213_al   R   23:48:48     1    8        N/A  62.50G cdr570 (None)
       60768140  celphin def-rieseber a1677140213_me  PD 2-00:00:00     1    1        N/A     10G  (Dependency)
       60768141  celphin def-rieseber a1677140213_me  PD 2-00:00:00     1    1        N/A     10G  (Dependency)
       60768142  celphin def-rieseber a1677140213_me  PD 2-00:00:00     1    8        N/A      2G  (Dependency)
       60768143  celphin def-rieseber a1677140213_ch  PD   23:50:00     1    1        N/A    256M  (Dependency)
       60768144  celphin def-rieseber a1677140213_fr  PD    1:00:00     1    8        N/A     80G  (Dependency)
       60768145  celphin def-rieseber a1677140213_de  PD      10:00     1    1        N/A    256M  (JobHeldUser)
       60768146  celphin def-rieseber a1677140213_de  PD 2-00:00:00     1    1        N/A      2G  (Dependency)
       60768147  celphin def-rieseber a1677140213_po  PD    1:40:00     1    1        N/A    256M  (Dependency)
       60768148  celphin def-rieseber a1677140213_du  PD   23:50:00     1    1        N/A      8G  (Dependency)
       60768149  celphin def-rieseber a1677140213_me  PD   23:50:00     1    8        N/A     10G  (Dependency)
       60768151  celphin def-rieseber a1677140213_me  PD   23:50:00     1    8        N/A     10G  (Dependency)
       60768152  celphin def-rieseber a1677140213_pr  PD   23:50:00     1    8        N/A      1G  (Dependency)
       60768153  celphin def-rieseber a1677140213_ba  PD   23:50:00     1    8        N/A     10G  (Dependency)
       60768154  celphin def-rieseber a1677140213_st  PD 2-00:00:00     1    1        N/A     25G  (Dependency)
       60768155  celphin def-rieseber a1677140213_st  PD 2-00:00:00     1    1        N/A     25G  (Dependency)
       60768156  celphin def-rieseber a1677140213_hi  PD 2-00:00:00     1    1        N/A    150G  (Dependency)
       60768157  celphin def-rieseber a1677140213_hi  PD 2-00:00:00     1    1        N/A    150G  (Dependency)
       60768158  celphin def-rieseber a1677140213_ar  PD   23:50:00     1    1        N/A      8G  (Dependency)
       60768160  celphin def-rieseber a1677140213_pr  PD   20:00:00     1    1        N/A      2G  (Dependency)


total 47964371
drwxr-s---  2 celphin def-henryg        4096 Feb 23 05:20 .
drwxr-s--- 10 celphin def-henryg        4096 Feb 23 00:16 ..
-rw-r-----  1 celphin def-henryg          10 Feb 23 00:17 NS.1871.003.IDT_i7_217---IDT_i5_217.Oxyria_Hi-C.fastq_linecount.txt
-rw-r-----  1 celphin def-henryg          58 Feb 23 05:14 NS.1871.003.IDT_i7_217---IDT_i5_217.Oxyria_Hi-C.fastq_norm.txt.res.txt
-rw-r-----  1 celphin def-henryg           0 Feb 23 05:50 NS.1871.003.IDT_i7_217---IDT_i5_217.Oxyria_Hi-C.fastq.sam
-rw-r-----  1 celphin def-henryg 69473425327 Feb 23 05:14 NS.1871.003.IDT_i7_217---IDT_i5_217.Oxyria_Hi-C.fastq.sam2
-rw-r-----  1 celphin def-henryg 69473424263 Feb 23 05:45 NS.1871.003.IDT_i7_217---IDT_i5_217.Oxyria_Hi-C.fastq.sam3
lrwxrwxrwx  1 celphin def-henryg         118 Feb 23 00:16 NS.1871.003.IDT_i7_217---IDT_i5_217.Oxyria_Hi-C_R1.fastq -> /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer/fastq/NS.1871.003.IDT_i7_217---IDT_i5_217.Oxyria_Hi-C_R1.fastq
lrwxrwxrwx  1 celphin def-henryg         118 Feb 23 00:16 NS.1871.003.IDT_i7_217---IDT_i5_217.Oxyria_Hi-C_R2.fastq -> /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer/fastq/NS.1871.003.IDT_i7_217---IDT_i5_217.Oxyria_Hi-C_R2.fastq


***! Failure during chimera handling of /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer/splits/NS.1871.003.IDT_i7_217---IDT_i5_217.Oxyria_Hi-C.fastq

Command terminated by signal 9
43.91user 20.10system 0:26.67elapsed 239%CPU (0avgtext+0avgdata 16777908maxresident)k
573976inputs+0outputs (147major+4193748minor)pagefaults 0swaps
slurmstepd: error: Detected 1 oom-kill event(s) in StepId=60768142.batch. Some of your processes may have been killed by the cgroup out-of-memory handler.


#----------------------------
# edit memory to be more (120G over 8 threads)

    864 #!/bin/bash -l
    865 #SBATCH -p $long_queue
    866 #SBATCH -o $debugdir/mergesort-%j.out
    867 #SBATCH -e $debugdir/mergesort-%j.err
    868 #SBATCH --mem=120G
    869 #SBATCH -t $long_queue_time
    870 #SBATCH -c $sortthreads
    871 #SBATCH --ntasks=1
    872 #SBATCH -d $dependalign
    873 #SBATCH -J "${groupname}_mergesort_${jname}"
    874 #SBATCH --threads-per-core=1
    875 $userstring
    876 ${load_samtools}
    877 #we should probably set the -m based on memory / num of threads
    878 if time samtools sort -t cb -n -O SAM -@ $sortthreads -l 0 -m 10G $name$ext.sam3 >  ${name}${ext}.sam
    879 then
    880    rm -f $name$ext.sam2 $name$ext.sam3
    881    touch $touchfile
    882 else
    883    echo "***! Failure during chimera handling of $name${ext}"
    884    touch $errorfile
    885    exit 1
    886 fi


########################
# try again 

module load StdEnv/2020
module load nixpkgs/16.09 
module load bwa/0.7.15
module load java/1.8.0_192 
module load gcc/4.8.5
module load cuda/7.5.18
cd /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer


./scripts/juicer.sh \
-D /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer/ \
-z ./references/Oxy_dig_1.scaffolds_FINAL.fasta \
-p ./restriction_sites/Oxyria.chrom.sizes \
-y ./restriction_sites/Oxyria_DpnII.txt \
-S chimeric

(-: Aligning files matching /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer/fastq/*_R*.fastq*
 in queue default to genome ./references/Oxy_dig_1.scaffolds_FINAL.fasta with no fragment delimited maps.
---  Using already created files in /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer/splits

(-: Starting job to launch other jobs once splitting is complete
sbatch: NOTE: Your memory request of 2048M was likely submitted as 2G. Please note that Slurm interprets memory requests denominated in G as multiples of 1024M, not 1000M.
sbatch: NOTE: Your memory request of 8192M was likely submitted as 8G. Please note that Slurm interprets memory requests denominated in G as multiples of 1024M, not 1000M.
sbatch: NOTE: Your memory request of 1024M was likely submitted as 1G. Please note that Slurm interprets memory requests denominated in G as multiples of 1024M, not 1000M.
sbatch: NOTE: Your memory request of 8192M was likely submitted as 8G. Please note that Slurm interprets memory requests denominated in G as multiples of 1024M, not 1000M.
sbatch: NOTE: Your memory request of 2048M was likely submitted as 2G. Please note that Slurm interprets memory requests denominated in G as multiples of 1024M, not 1000M.
(-: Finished adding all jobs... Now is a good time to get that cup of coffee... Last job id 60845146

          JOBID     USER      ACCOUNT           NAME  ST  TIME_LEFT NODES CPUS TRES_PER_N MIN_MEM NODELIST (REASON)
       60845001  celphin def-rieseber a1677188694_NS  PD   23:50:00     1    1        N/A      5G  (Priority)
       60845006  celphin def-rieseber a1677188694_me  PD 2-00:00:00     1    1        N/A     10G  (Dependency)
       60845014  celphin def-rieseber a1677188694_me  PD 2-00:00:00     1    1        N/A     10G  (Dependency)
       60845016  celphin def-rieseber a1677188694_me  PD 2-00:00:00     1    8        N/A    120G  (Dependency)
       60845021  celphin def-rieseber a1677188694_ch  PD   23:50:00     1    1        N/A    256M  (Dependency)
       60845036  celphin def-rieseber a1677188694_fr  PD    1:00:00     1    8        N/A     80G  (Dependency)
       60845050  celphin def-rieseber a1677188694_de  PD      10:00     1    1        N/A    256M  (JobHeldUser)
       60845058  celphin def-rieseber a1677188694_de  PD 2-00:00:00     1    1        N/A      2G  (Dependency)
       60845067  celphin def-rieseber a1677188694_po  PD    1:40:00     1    1        N/A    256M  (Dependency)
       60845071  celphin def-rieseber a1677188694_du  PD   23:50:00     1    1        N/A      8G  (Dependency)
       60845092  celphin def-rieseber a1677188694_me  PD   23:50:00     1    8        N/A     10G  (Dependency)
       60845104  celphin def-rieseber a1677188694_me  PD   23:50:00     1    8        N/A     10G  (Dependency)
       60845109  celphin def-rieseber a1677188694_pr  PD   23:50:00     1    8        N/A      1G  (Dependency)
       60845114  celphin def-rieseber a1677188694_ba  PD   23:50:00     1    8        N/A     10G  (Dependency)
       60845120  celphin def-rieseber a1677188694_st  PD 2-00:00:00     1    1        N/A     25G  (Dependency)
       60845127  celphin def-rieseber a1677188694_st  PD 2-00:00:00     1    1        N/A     25G  (Dependency)
       60845129  celphin def-rieseber a1677188694_hi  PD 2-00:00:00     1    1        N/A    150G  (Dependency)
       60845137  celphin def-rieseber a1677188694_hi  PD 2-00:00:00     1    1        N/A    150G  (Dependency)
       60845143  celphin def-rieseber a1677188694_ar  PD   23:50:00     1    1        N/A      8G  (Dependency)
       60845146  celphin def-rieseber a1677188694_pr  PD   20:00:00     1    1        N/A      2G  (Dependency)

# blanked sam files again??

#########################################
# up memory and lower time for other parts

##################################################
# try again from start??

module load StdEnv/2020
module load nixpkgs/16.09 
module load bwa/0.7.15
module load java/1.8.0_192 
module load gcc/4.8.5
module load cuda/7.5.18

cd /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer

./scripts/juicer.sh \
-D /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer \
-z ./references/Oxy_dig_1.scaffolds_FINAL.fasta \
-p ./restriction_sites/Oxyria.chrom.sizes \
-y ./restriction_sites/Oxyria_DpnII.txt


(-: Looking for fastq files...fastq files exist
(-: Aligning files matching /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer/fastq/*_R*.fastq*
 in queue default to genome ./references/Oxy_dig_1.scaffolds_FINAL.fasta with no fragment delimited maps.
(-: Created /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer/splits and /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer/aligned.
(-: Starting job to launch other jobs once splitting is complete
(-: Finished adding all jobs... Now is a good time to get that cup of coffee... Last job id 60850608


          JOBID     USER      ACCOUNT           NAME  ST  TIME_LEFT NODES CPUS TRES_PER_N MIN_MEM NODELIST (REASON)
       60850524  celphin def-rieseber a1677193290_cm  PD       2:00     1    1        N/A    256M  (Priority)
       60850530  celphin def-rieseber a1677193290_NS  PD   12:50:00     1    1        N/A      5G  (Priority)
       60850535  celphin def-rieseber a1677193290_al  PD   12:50:00     1    8        N/A  62.50G  (Priority)
       60850542  celphin def-rieseber a1677193290_me  PD 2-00:00:00     1    1        N/A     10G  (Dependency)
       60850546  celphin def-rieseber a1677193290_me  PD 2-00:00:00     1    1        N/A     10G  (Dependency)
       60850552  celphin def-rieseber a1677193290_me  PD 2-00:00:00     1    8        N/A    120G  (Dependency)
       60850559  celphin def-rieseber a1677193290_ch  PD   12:50:00     1    1        N/A    256M  (Dependency)
       ***60850563  celphin def-rieseber a1677193290_fr  PD    1:00:00     1    8        N/A    120G  (Dependency)
       60850566  celphin def-rieseber a1677193290_de  PD      10:00     1    1        N/A    256M  (JobHeldUser)
       60850570  celphin def-rieseber a1677193290_de  PD 2-00:00:00     1    1        N/A    120G  (Dependency)
       60850571  celphin def-rieseber a1677193290_po  PD    5:00:00     1    1        N/A    256M  (Dependency)
       60850592  celphin def-rieseber a1677193290_du  PD   12:50:00     1    1        N/A    120G  (Dependency)
       60850593  celphin def-rieseber a1677193290_me  PD   12:50:00     1    8        N/A     20G  (Dependency)
       60850594  celphin def-rieseber a1677193290_me  PD   12:50:00     1    8        N/A     20G  (Dependency)
       60850595  celphin def-rieseber a1677193290_pr  PD   12:50:00     1    8        N/A     20G  (Dependency)
       60850602  celphin def-rieseber a1677193290_ba  PD   12:50:00     1    8        N/A     20G  (Dependency)
       60850603  celphin def-rieseber a1677193290_st  PD 2-00:00:00     1    1        N/A    120G  (Dependency)
       60850604  celphin def-rieseber a1677193290_st  PD 2-00:00:00     1    1        N/A    120G  (Dependency)
       60850605  celphin def-rieseber a1677193290_hi  PD 2-00:00:00     1    1        N/A    150G  (Dependency)
       60850606  celphin def-rieseber a1677193290_hi  PD 2-00:00:00     1    1        N/A    150G  (Dependency)
       60850607  celphin def-rieseber a1677193290_ar  PD   12:50:00     1    1        N/A     20G  (Dependency)
       60850608  celphin def-rieseber a1677193290_pr  PD    2:00:00     1    1        N/A     20G  (Dependency)

#################################
# worked further
[celphin@cedar1 aligned]$ ls
inter_30_hists.m  inter_30.txt  inter_hists.m  inter.txt  merged1.txt  merged30.txt  merged_dedup.bam  tmp

# problems

JOBID     STATE    NAME                        DEPENDENCY                                                   NODELIST(REASON)
60850571  RUNNING  a1677193290_post_dedup      (null)                                                       cdr621
60850592  PENDING  a1677193290_dupcheck        afterok:60850571(unfulfilled)                                (Dependency)
60850593  PENDING  a1677193290_merged1         afterok:60850592(unfulfilled)                                (Dependency)
60850594  PENDING  a1677193290_merged30        afterok:60850592(unfulfilled)                                (Dependency)
60850595  PENDING  a1677193290_prestats        afterok:60850592(unfulfilled)                                (Dependency)
60850602  PENDING  a1677193290_bamrm           afterok:60850593(unfulfilled),afterok:60850594(unfulfilled)  (Dependency)
60850603  PENDING  a1677193290_stats           afterok:60850593(unfulfilled),afterok:60850595(unfulfilled)  (Dependency)
60850604  PENDING  a1677193290_stats30         afterok:60850594(unfulfilled),afterok:60850602(unfulfilled)  (Dependency)
60850605  PENDING  a1677193290_hic             afterok:60850603(unfulfilled)                                (Dependency)
60850606  PENDING  a1677193290_hic30           afterok:60850604(unfulfilled)                                (Dependency)
60850607  PENDING  a1677193290_arrowhead_wrap  afterok:60850605(unfulfilled),afterok:60850606(unfulfilled)  (Dependency)
60850608  PENDING  a1677193290_prep_done       afterok:60850607(unfulfilled)                                (Dependency)

more fragmerge-60850563.err
/var/spool/slurmd/job60850563/slurm_script: samtools: command not found

[E::hts_open_format] Failed to open file "/home/celphin/projects/def-henryg/celphin/Oxyria/Juicer/aligned/merged_dedup.sam" : No such file or directory
samtools view: failed to open "/home/celphin/projects/def-henryg/celphin/Oxyria/Juicer/aligned/merged_dedup.sam" for reading: No such file or directory

awk: /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer//scripts/split_rmdups_sam.awk:36: fatal: cannot open file `/home/celphin/projects/def-henryg/cel
phin/Oxyria/Juicer/aligned/merged_sort.sam' for reading (No such file or directory)

#-------------------
# problem: other scripts need samtools

adjust_insert_size.awk  dups_sam.awk                                                                          juicer_postprocessing.sh  relaunch_prep.sh
check.sh                fragment_4dnpairs.pl                                                                  juicer.sh                 sam_to_mnd.sh
chimeric_sam.awk        generate_site_positions.py                                                            juicer_tools              sam_to_pre.awk
cleanup.sh              GSE63525_GM12878_primary_replicate_HiCCUPS_looplist_with_motifs_unique_localized.txt  juicer_tools.2.18.00.jar  split_rmdups_sam.awk
conversion.sh           index_by_chr.awk                                                                      juicer_tools.jar          stats_sub.awk
countligations.sh       juicer_arrowhead.sh                                                                   lib64
diploid.pl              juicer-copy.sh                                                                        makemega_addstats.awk
diploid_split.awk       juicer_hiccups.sh                                                                     mega.sh

#----------------------
more juicer_arrowhead.sh

# Aiden Lab specific check
isRice=$(hostname | awk '{if ($1~/rice/){print 1}else {print 0}}')
isBCM=$(hostname | awk '{if ($1~/bcm/){print 1}else {print 0}}')
isVoltron=0
# Set default appropriately
if [ $isRice -eq 1 ]
then
    juicer_tools_path="/projects/ea14/juicer/scripts/juicer_tools"
elif [ $isBCM -eq 1 ]
then
    juicer_tools_path="/storage/aiden/juicer/scripts/juicer_tools"
else
    isVoltron=1
    juicer_tools_path="/gpfs0/juicer2/scripts/juicer_tools"
fi

# change to 
juicer_tools_path="/home/celphin/projects/def-henryg/celphin/Oxyria/Juicer/scripts/juicer_tools"

#----------------------------
grep isVoltron *

juicer_hiccups.sh:isVoltron=0
juicer_hiccups.sh:    isVoltron=1

nano juicer_hiccups.sh

# change to
juiceDir="/home/celphin/projects/def-henryg/celphin/Oxyria/Juicer"
juicer_tools_path="/home/celphin/projects/def-henryg/celphin/Oxyria/Juicer/scripts/juicer_tools"
bed_file_dir="/home/celphin/projects/def-henryg/celphin/Oxyria/Juicer/references/motif"

#---------------------
juicer_postprocessing.sh:isVoltron=0
juicer_postprocessing.sh:    isVoltron=1

nano juicer_postprocessing.sh
# change to
juiceDir="/home/celphin/projects/def-henryg/celphin/Oxyria/Juicer"
juicer_tools_path="/home/celphin/projects/def-henryg/celphin/Oxyria/Juicer/scripts/juicer_tools"
bed_file_dir="/home/celphin/projects/def-henryg/celphin/Oxyria/Juicer/references/motif"

#------------------------
mega.sh:isVoltron=0
mega.sh:    isVoltron=1
mega.sh:    if [ $isRice -eq 1 ] || [ $isVoltron -eq 1 ]

nano mega.sh

load_bwa="module load  StdEnv/2020 bwa/0.7.17"
load_java="module load nixpkgs/16.09  java/1.8.0_192"
load_gpu="module load nixpkgs/16.09 gcc/4.8.5 cuda/7.5.18"
load_samtools="module load StdEnv/2020 samtools/1.16.1"
# Juicer directory, contains scripts/, references/, and restriction_sites/
# can also be set in options via -D
juiceDir="/home/celphin/projects/def-henryg/celphin/Oxyria/Juicer" ### RICE
# default queue, can also be set in options via -q
queue=default
# default long queue, can also be set in options via -l
long_queue=default
long_queue_time=2-00:00:00


    myPath=/bin:$PATH
    load_bwa="module load BioBuilds/2015.04"
    load_java="module load Java/8.0.3.22"
    load_gpu="module load gcccuda/2016a;module load CUDA/8.0.54;"
    # Juicer directory, contains scripts/, references/, and restriction_sites/
    # can also be set in options via -D
    juiceDir=""
    # default queue, can also be set in options via -q
    queue="commons"
    # default long queue, can also be set in options via -l
    long_queue="commons"
    long_queue_time="1440"
	
	
######################
# samtools needed to load line 966 in juicer.sh
${load_samtools}

#########################################
# try again from further along

module load StdEnv/2020
module load nixpkgs/16.09 
module load bwa/0.7.15
module load java/1.8.0_192 
module load gcc/4.8.5
module load cuda/7.5.18

cd /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer

./scripts/juicer.sh \
-D /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer \
-z ./references/Oxy_dig_1.scaffolds_FINAL.fasta \
-p ./restriction_sites/Oxyria.chrom.sizes \
-y ./restriction_sites/Oxyria_DpnII.txt \
-S merge

sq
          JOBID     USER      ACCOUNT           NAME  ST  TIME_LEFT NODES CPUS TRES_PER_N MIN_MEM NODELIST (REASON)
       60881379  celphin def-rieseber a1677216310_cm  PD       2:00     1    1        N/A    256M  (Priority)
       60881381  celphin def-rieseber a1677216310_fr  PD    1:00:00     1    8        N/A    120G  (Priority)
       60881383  celphin def-rieseber a1677216310_de  PD      10:00     1    1        N/A    256M  (JobHeldUser)
       60881384  celphin def-rieseber a1677216310_de  PD 2-00:00:00     1    1        N/A    120G  (Dependency)
       60881385  celphin def-rieseber a1677216310_po  PD    5:00:00     1    1        N/A    256M  (Dependency)
       60881387  celphin def-rieseber a1677216310_du  PD   12:50:00     1    1        N/A    120G  (Dependency)
       60881388  celphin def-rieseber a1677216310_me  PD   12:50:00     1    8        N/A     20G  (Dependency)
       60881390  celphin def-rieseber a1677216310_me  PD   12:50:00     1    8        N/A     20G  (Dependency)
       60881391  celphin def-rieseber a1677216310_pr  PD   12:50:00     1    8        N/A     20G  (Dependency)
       60881392  celphin def-rieseber a1677216310_ba  PD   12:50:00     1    8        N/A     20G  (Dependency)
       60881393  celphin def-rieseber a1677216310_st  PD 2-00:00:00     1    1        N/A    120G  (Dependency)
       60881394  celphin def-rieseber a1677216310_st  PD 2-00:00:00     1    1        N/A    120G  (Dependency)
       60881395  celphin def-rieseber a1677216310_hi  PD 2-00:00:00     1    1        N/A    150G  (Dependency)
       60881396  celphin def-rieseber a1677216310_hi  PD 2-00:00:00     1    1        N/A    150G  (Dependency)
       60881397  celphin def-rieseber a1677216310_ar  PD   12:50:00     1    1        N/A     20G  (Dependency)
       60881398  celphin def-rieseber a1677216310_pr  PD    2:00:00     1    1        N/A     20G  (Dependency)

more fincln-60881398.out

Fri Feb 24 00:43:23 PST 2023
(-: Pipeline successfully completed (-:
Run cleanup.sh to remove the splits directory
Check /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer/aligned for results
Fri Feb 24 00:43:24 PST 2023

#################################################
# output

cd /home/celphin/projects/rpp-rieseber/celphin/Cannabis/Juicer_Mojtaba/

HiC_XandY_ragtag.scaffold.0.assembly  HiC_XandY_ragtag.scaffold.0.hic  inter.hic

###########################
# need merged_nodups.txt still so try rerunning again

module load StdEnv/2020
module load nixpkgs/16.09 
module load bwa/0.7.15
module load java/1.8.0_192 
module load gcc/4.8.5
module load cuda/7.5.18

cd /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer

./scripts/juicer.sh \
-D /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer \
-z ./references/Oxy_dig_1.scaffolds_FINAL.fasta \
-p ./restriction_sites/Oxyria.chrom.sizes \
-y ./restriction_sites/Oxyria_DpnII.txt


(-: Looking for fastq files...fastq files exist
(-: Aligning files matching /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer/fastq/*_R*.fastq*
 in queue default to genome ./references/Oxy_dig_1.scaffolds_FINAL.fasta with no fragment delimited maps.
(-: Created /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer/splits and /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer/aligned.
(-: Starting job to launch other jobs once splitting is complete
(-: Finished adding all jobs... Now is a good time to get that cup of coffee... Last job id 60897757


          JOBID     USER      ACCOUNT           NAME  ST  TIME_LEFT NODES CPUS TRES_PER_N MIN_MEM NODELIST (REASON)
       60897699  celphin def-rieseber a1677231328_cm  PD       2:00     1    1        N/A    256M  (Priority)
       60897701  celphin def-rieseber a1677231328_NS  PD   12:50:00     1    1        N/A      5G  (Priority)
       60897703  celphin def-rieseber a1677231328_al  PD   12:50:00     1    8        N/A  62.50G  (Priority)
       60897706  celphin def-rieseber a1677231328_me  PD 2-00:00:00     1    1        N/A     10G  (Dependency)
       60897708  celphin def-rieseber a1677231328_me  PD 2-00:00:00     1    1        N/A     10G  (Dependency)
       60897711  celphin def-rieseber a1677231328_me  PD 2-00:00:00     1    8        N/A    120G  (Dependency)
       60897713  celphin def-rieseber a1677231328_ch  PD   12:50:00     1    1        N/A    256M  (Dependency)
       60897716  celphin def-rieseber a1677231328_fr  PD    1:00:00     1    8        N/A    120G  (Dependency)
       60897718  celphin def-rieseber a1677231328_de  PD      10:00     1    1        N/A    256M  (JobHeldUser)
       60897722  celphin def-rieseber a1677231328_de  PD 2-00:00:00     1    1        N/A    120G  (Dependency)
       60897726  celphin def-rieseber a1677231328_po  PD    5:00:00     1    1        N/A    256M  (Dependency)
       60897729  celphin def-rieseber a1677231328_du  PD   12:50:00     1    1        N/A    120G  (Dependency)
       60897731  celphin def-rieseber a1677231328_me  PD   12:50:00     1    8        N/A     20G  (Dependency)
       60897734  celphin def-rieseber a1677231328_me  PD   12:50:00     1    8        N/A     20G  (Dependency)
       60897738  celphin def-rieseber a1677231328_pr  PD   12:50:00     1    8        N/A     20G  (Dependency)
       60897741  celphin def-rieseber a1677231328_ba  PD   12:50:00     1    8        N/A     20G  (Dependency)
       60897743  celphin def-rieseber a1677231328_st  PD 2-00:00:00     1    1        N/A    120G  (Dependency)
       60897747  celphin def-rieseber a1677231328_st  PD 2-00:00:00     1    1        N/A    120G  (Dependency)
       60897749  celphin def-rieseber a1677231328_hi  PD 2-00:00:00     1    1        N/A    150G  (Dependency)
       60897751  celphin def-rieseber a1677231328_hi  PD 2-00:00:00     1    1        N/A    150G  (Dependency)
       60897753  celphin def-rieseber a1677231328_ar  PD   12:50:00     1    1        N/A     20G  (Dependency)
       60897757  celphin def-rieseber a1677231328_pr  PD    2:00:00     1    1        N/A     20G  (Dependency)

#----------------
# output is the same

#------------------
grep merged_nodups.txt *
juicer-copy.sh: samtools view $sthreadstring -O SAM -F 1024 $outputdir/merged_dedup.*am | awk -v mnd=1 -f ${juiceDir}/scripts/sam_to_pre.awk > ${outputdir}/merged_nodups.txt
juicer.sh:      samtools view $sthreadstring -O SAM -F 1024 $outputdir/merged_dedup.*am | awk -v mnd=1 -f ${juiceDir}/scripts/sam_to_pre.awk > ${outputdir}/merged_nodups.txt
grep: lib64: Is a directory
sam_to_mnd.sh:# How to create old style merged_nodups.txt for 3D-DNA


more sam_to_mnd.sh
# How to create old style merged_nodups.txt for 3D-DNA
# File should be dedupped SAM
# NB: you should likely put the directory in from of the sam_to_pre.awk file
samtools view -O SAM -F 1024 $1 | awk -v mnd=1 -f sam_to_pre.awk


####################
# make merged_nodups.txt from bam 

tmux attach-session -t Juicer

salloc -c10 --time 2:50:00 --mem 120000m --account def-rieseber

module load StdEnv/2020 samtools/1.16.1

samtools view -@ 8 -O SAM -F 1024 /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer/aligned/merged_dedup.*am | \
awk -v mnd=1 -f /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer/scripts/sam_to_pre.awk > \
/home/celphin/projects/def-henryg/celphin/Oxyria/Juicer/aligned/merged_nodups.txt

###################################
# need to run through 3D DNA to get movable scaffolds .assembly file ...
# https://aidenlab.org/assembly/manual_180322.pdf 

# downloading 3D DNA
# https://github.com/aidenlab/3d-dna

#--------------------------
# Prerequisites

    # LastZ (version 1.03.73 released 20150708) – for diploid mode only
    # Java version >=1.7
    # Bash >=4
    # GNU Awk >=4.0.2
    # GNU coreutils sort >=8.11
    # Python >=2.7 - for chromosome number-aware splitter module only
    # scipy numpy matplotlib - for chromosome number-aware splitter module only

module load StdEnv/2020 lastz/1.04.03 java/1.8.0_192 python/3.10.2 scipy-stack/2022a

#--------------
# git clone
cd /home/celphin/projects/def-henryg/celphin/Oxyria/3D_DNA
git clone https://github.com/aidenlab/3d-dna.git

chmod -R 770 /home/celphin/projects/def-henryg/celphin/Oxyria/3D_DNA/3d-dna/*

#-----------------------
# get .assembly file from  merged_nodups.txt and fasta
tmux attach-session -t Juicer

cd /home/celphin/projects/def-henryg/celphin/Oxyria/3D_DNA
awk -f ./3d-dna/utils/generate-assembly-file-from-fasta.awk /home/celphin/projects/def-henryg/celphin/Oxyria/Oxy_dig_1.scaffolds_FINAL.fasta > Oxyria_draft.assembly

salloc -c10 --time 2:50:00 --mem 120000m --account def-rieseber
module load StdEnv/2020 lastz/1.04.03 java/1.8.0_192 python/3.10.2 scipy-stack/2022a
cd /home/celphin/projects/def-henryg/celphin/Oxyria/3D_DNA
/home/celphin/projects/def-henryg/celphin/Oxyria/3D_DNA/3d-dna/run-asm-pipeline-post-review.sh \
-r /home/celphin/projects/def-henryg/celphin/Oxyria/3D_DNA/Oxyria_draft.assembly \
/home/celphin/projects/def-henryg/celphin/Oxyria/Oxy_dig_1.scaffolds_FINAL.fasta \
/home/celphin/projects/def-henryg/celphin/Oxyria/Juicer/aligned/merged_nodups.txt

 -r|--review flag was triggered, treating file /home/celphin/projects/def-henryg/celphin/Oxyria/3D_DNA/Oxyria_draft.assembly as a JB4A review file for draft fasta in arguments.
###############
Finilizing output:
:| Warning: No input for label1 was provided. Default for label1 is ":::fragment_"
:| Warning: No input for label2 was provided. Default for label2 is ":::debris"
:| Warning: No input for label1 was provided. Default for label1 is ":::fragment_"
:| Warning: No input for label2 was provided. Default for label2 is ":::debris"
:| Warning: No input for label1 was provided. Default for label1 is ":::fragment_"
:| Warning: No input for label2 was provided. Default for label2 is ":::debris"
:| Warning: No input for label1 was provided. Default for label1 is ":::fragment_"
:| Warning: No input for label2 was provided. Default for label2 is ":::debris"
:| Warning: No input for label1 was provided. Default for label1 is ":::fragment_"
:| Warning: No input for label2 was provided. Default for label2 is ":::debris"
:| Warning: No input for label1 was provided. Default for label1 is ":::fragment_"
:| Warning: No input for label2 was provided. Default for label2 is ":::debris"
:| Warning: No input for label1 was provided. Default for label1 is ":::fragment_"
:| Warning: No input for label2 was provided. Default for label2 is ":::debris"
:| Warning: No input for label1 was provided. Default for label1 is ":::fragment_"
:| Warning: No input for label2 was provided. Default for label2 is ":::debris"
:| Warning: No input for label1 was provided. Default for label1 is ":::fragment_"
:| Warning: No input for label2 was provided. Default for label2 is ":::debris"
:| Warning: No input for label1 was provided. Default for label1 is ":::fragment_"
:| Warning: No input for label2 was provided. Default for label2 is ":::debris"
:| Warning: No input for label1 was provided. Default for label1 is ":::fragment_"
:| Warning: No input for label2 was provided. Default for label2 is ":::debris"
:| Warning: No input for label1 was provided. Default for label1 is ":::fragment_"
:| Warning: No input for label2 was provided. Default for label2 is ":::debris"
:| Warning: No input for label1 was provided. Default for label1 is ":::fragment_"
:| Warning: No input for label2 was provided. Default for label2 is ":::debris"
:) -p flag was triggered. Running with GNU Parallel support parameter set to true.
:) -q flag was triggered, starting calculations for 1 threshold mapping quality
:) -i flag was triggered, building mapq without
:) -c flag was triggered, will remove temporary files after completion
...Remapping contact data from the original contig set to assembly

...Building track files
...Building the hic file
Picked up JAVA_TOOL_OPTIONS: -Xmx2g
Feb 24, 2023 3:36:10 PM java.util.prefs.FileSystemPreferences$1 run
INFO: Created user preferences directory.
Not including fragment map
Start preprocess
Writing header
Writing body
..
Writing footer

Finished preprocess
HiC file version: 8

Calculating norms for zoom BP_2500000
Calculating norms for zoom BP_1000000
Calculating norms for zoom BP_500000
Calculating norms for zoom BP_250000
Calculating norms for zoom BP_100000
Calculating norms for zoom BP_50000
Calculating norms for zoom BP_25000
Calculating norms for zoom BP_10000


#########################
# open in Juicebox 
# https://aidenlab.org/juicebox/

# use desktop version to edit manually
# https://www.youtube.com/watch?v=Nj7RhQZHM18
# need .assembly and .hic files

#################################
# after editing 
# convert back to fasta: juicebox_assembly_converter.py
# copy review.assembly into Cedar
# https://github.com/phasegenomics/juicebox_scripts

git clone https://github.com/phasegenomics/juicebox_scripts.git
chmod -R 770 /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer/juicebox_scripts/*

tmux new-session -s Juicer 
tmux attach-session -t Juicer

salloc -c1 --time 2:50:00 --mem 120000m --account def-rieseber
module load StdEnv/2020 lastz/1.04.03 java/1.8.0_192 python/3.10.2 scipy-stack/2022a

cd /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer/juicebox_scripts/

python /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer/juicebox_scripts/juicebox_scripts/juicebox_assembly_converter.py \
-a /home/celphin/projects/def-henryg/celphin/Oxyria/Juicer/juicebox_scripts/data/Oxy_dig_1.scaffolds_FINAL.final.review.assembly  \
-f /home/celphin/projects/def-henryg/celphin/Oxyria/Oxy_dig_1.scaffolds_FINAL.fasta 


########################
module load StdEnv/2020 
module load gcc/9.3.0
module load r-bundle-bioconductor/3.14

cd /home/celphin/projects/def-cronk/celphin/Plectritis_analysis/Assembly/GenomeQC/assembly-stats

./assembly-stats /home/celphin/scratch/repeats/Oxyria/Oxy_dig_1.scaffolds_FINAL.final.review.fasta

stats for /home/celphin/scratch/repeats/Oxyria/Oxy_dig_1.scaffolds_FINAL.final.review.fasta
sum = 590427020, n = 1560, ave = 378478.86, largest = 88146246
N50 = 76036264, n = 4
N60 = 73842794, n = 5
N70 = 73289246, n = 6
N80 = 70136240, n = 7
N90 = 70136240, n = 7
N100 = 7, n = 1560
N_count = 87158
Gaps = 250

###############################
# BUSCO

module load StdEnv/2020 gcc python augustus hmmer blast+ metaeuk prodigal r
source /project/6019339/celphin/Cannabis/busco_env/

cd /home/celphin/projects/def-henryg/celphin/Oxyria/BUSCO/

# Final assembly
busco --offline --in /home/celphin/scratch/repeats/Oxyria/Oxy_dig_1.scaffolds_FINAL.final.review.fasta \
--out Final_BUSCO_eudicots --lineage_dataset eudicots_odb10 --mode genome --cpu 32 \
--download_path /project/6019339/celphin/Cannabis/BUSCO_downloads/BUSCO_downloads/

 Results:        C:93.4%[S:85.3%,D:8.1%],F:1.8%,M:4.8%,n:2326

2023-02-26 21:07:25 INFO:

        --------------------------------------------------
        |Results from dataset eudicots_odb10              |
        --------------------------------------------------
        |C:93.4%[S:85.3%,D:8.1%],F:1.8%,M:4.8%,n:2326     |
        |2172   Complete BUSCOs (C)                       |
        |1983   Complete and single-copy BUSCOs (S)       |
        |189    Complete and duplicated BUSCOs (D)        |
        |43     Fragmented BUSCOs (F)                     |
        |111    Missing BUSCOs (M)                        |
        |2326   Total BUSCO groups searched               |
        --------------------------------------------------


# HiC scaffolds
busco --offline --in /home/celphin/projects/def-henryg/celphin/Oxyria/Oxy_dig_1.scaffolds_FINAL.fasta \
--out HiC_Scaffolds_BUSCO_eudicots --lineage_dataset eudicots_odb10 --mode genome --cpu 32 \
--download_path /project/6019339/celphin/Cannabis/BUSCO_downloads/BUSCO_downloads/

Results:        C:93.4%[S:85.2%,D:8.2%],F:1.8%,M:4.8%,n:2326

        --------------------------------------------------
        |Results from dataset eudicots_odb10              |
        --------------------------------------------------
        |C:93.4%[S:85.2%,D:8.2%],F:1.8%,M:4.8%,n:2326     |
        |2171   Complete BUSCOs (C)                       |
        |1981   Complete and single-copy BUSCOs (S)       |
        |190    Complete and duplicated BUSCOs (D)        |
        |43     Fragmented BUSCOs (F)                     |
        |112    Missing BUSCOs (M)                        |
        |2326   Total BUSCO groups searched               |
        --------------------------------------------------

# longstitch
busco --offline --in /home/celphin/projects/def-cronk/celphin/Plectritis_analysis/Assembly/Polish/LongStitch/Oxyria1_Sept6_haplotigs_oxy3.fa \
--out Longstitch_BUSCO_eudicots --lineage_dataset eudicots_odb10 --mode genome --cpu 32 \
--download_path /project/6019339/celphin/Cannabis/BUSCO_downloads/BUSCO_downloads/

 Results:        C:93.4%[S:85.2%,D:8.2%],F:1.8%,M:4.8%,n:2326

2023-02-26 22:34:51 INFO:

        --------------------------------------------------
        |Results from dataset eudicots_odb10              |
        --------------------------------------------------
        |C:93.4%[S:85.2%,D:8.2%],F:1.8%,M:4.8%,n:2326     |
        |2172   Complete BUSCOs (C)                       |
        |1981   Complete and single-copy BUSCOs (S)       |
        |191    Complete and duplicated BUSCOs (D)        |
        |43     Fragmented BUSCOs (F)                     |
        |111    Missing BUSCOs (M)                        |
        |2326   Total BUSCO groups searched               |
        --------------------------------------------------

# Contig assembly
busco --offline --in /home/celphin/projects/def-henryg/celphin/Oxyria/Oxyria1_Sept6_haplotigs.p_ctg.fa \
--out Contig_BUSCO_eudicots --lineage_dataset eudicots_odb10 --mode genome --cpu 32 \
--download_path /project/6019339/celphin/Cannabis/BUSCO_downloads/BUSCO_downloads/

Results:        C:93.3%[S:85.1%,D:8.2%],F:1.8%,M:4.9%,n:2326


        --------------------------------------------------
        |Results from dataset eudicots_odb10              |
        --------------------------------------------------
        |C:93.3%[S:85.1%,D:8.2%],F:1.8%,M:4.9%,n:2326     |
        |2171   Complete BUSCOs (C)                       |
        |1980   Complete and single-copy BUSCOs (S)       |
        |191    Complete and duplicated BUSCOs (D)        |
        |43     Fragmented BUSCOs (F)                     |
        |112    Missing BUSCOs (M)                        |
        |2326   Total BUSCO groups searched               |
        --------------------------------------------------

# try another busco with fewer genes
cd /home/celphin/projects/def-henryg/celphin/Oxyria/BUSCO/

busco --offline --in /home/celphin/projects/def-henryg/celphin/Oxyria/Oxyria1_Sept6_haplotigs.p_ctg.fa \
--out Contig_BUSCO_viridiplantae --lineage_dataset viridiplantae_odb10 --mode genome --cpu 20 \
--download_path /project/6019339/celphin/Cannabis/BUSCO_downloads/BUSCO_downloads/



################################
# quality check of initial reads

cd /project/6064374/celphin/Oxyria/PacBio_raw_reads/PacBio_Aug2022/

module load StdEnv/2020
module load fastqc/0.11.9

fastqc *fastq -o ./fastQC_reports

# done copied to google drive


##############################################
# try annotation with Jose's pipeline

# https://github.com/megahitokiri/CAP_Snakemake
