#############################################################
#Install/modify snakemake + dependencies + configs for Pipeline:
#August 2023
#Follow conda set up notes first
    #Repeat Masker download
    #Setting up conda environments
    #Snakemake set up
    #Pipeline download + modifications
    #Notes on running pipeline + some trouble shooting
#############################################################
#Getting Maker: 
# on Compute Canada  - simple
# email a screen clipping of your registration to support@tech.alliancecan.ca

module load StdEnv/2020 gcc/9.3.0 openmpi/4.0.3 maker/3.01.03

# test 
maker -h

####################################################
# Git clone needed repos
tmux new-session -s Annotation
tmux attach-session -t Annotation

cd ~/scratch/
mkdir Annotation; cd Annotation

# clone the needed repos
git clone https://github.com/megahitokiri/CAP_Snakemake.git
git clone https://github.com/billzt/gff3sort.git
git clone https://github.com/oushujun/EDTA.git

#######################################################
# copy over new files from: https://github.com/rieseberglab/Genome_assemblies_annotations/tree/main/Oxyria_Assembly_Annotations/2_Annotations/1_Annotation_Maker_Pipeline/MS_Annotation_Files_Oxyria

cd ~/scratch/Annotation/MS_Annotation_Files_Oxyria/
cp * ../CAP_Snakemake/
 
##########################################################
#Conda environment set up 
~/miniconda2/bin/conda info --envs
cd ~/scratch/Annotation/EDTA

#----------------
# activate base environment
source ~/miniconda2/bin/activate

# update conda
conda update -n base -c defaults conda
# Your installed version is: 2.17

# install EDTA
conda env create -f EDTA.yml
Your installed version is: 2.17
#install EDTA (CHECK FILE NAME OF .yml file)    
    #Trouble shooting note: 
        #depending on updates might require miniconda redownload
        #Works with Original EDTA.yml file, available in this repo 
conda env create -n EDTA -f  EDTA.yml
#to remove and reinstall: conda remove --name EDTA --all

########################

# you can open this environment here
source ~/miniconda2/bin/activate EDTA
# to escape
# conda deactivate
#-------------------------
#create  - maybe not needed - since not working
# just module load StdEnv/2020 gcc/9.3.0 python augustus hmmer blast+ metaeuk prodigal r 
~/miniconda2/bin/conda create --name BUSCO  
# open BUSCO
source ~/miniconda2/bin/activate BUSCO
# install BUSCO
# https://anaconda.org/bioconda/busco
conda install -c bioconda busco
# fails to properly install (?)

# to escape
# conda deactivate
#---
# to remove an enivronment
# source ~/miniconda2/bin/activate
# conda remove -n EDTA --all
###############################

#Set up snakemake profile:
module load StdEnv/2020
module load python/3.10.2
module load scipy-stack/2022a
module load gcc/9.3.0
module load cuda/11.7
module load hdf5/1.12.1
module load parasail/2.5

pip install cookiecutter --user

module load  StdEnv/2020 intel/2020.1.217 metagenome-atlas/2.5.0
cookiecutter gh:khanlab/cc-slurm -o ~/.config/snakemake -f

You've downloaded /home/celphin/.cookiecutters/cc-slurm before. Is it okay to delete and re-download it? [yes]: yes
profile_name [cc-slurm]:
account [ctb-akhanf]: def-rieseber
gpu_type [t4]:
default_time [60]: 420
default_mem_mb [4000]: 12000
default_gpus [0]:
max_time [4320]:
max_mem_mb [64000]:
max_gpus [1]:
max_threads [16]: 32
immediate_submit [True]:
use_singularity [True]:
use_envmodules [True]:
verbose [False]:
singularity_prefix [/project/6050199/akhanf/singularity/snakemake_containers/]: 'singularity'
Success! You can run the profile with:'

# edit one part of profile manually

# edit singularity-prefix: 'singularity'
# comment out: notebook-listen: '$(hostname -f):8888'

nano ~/.config/snakemake/cc-slurm/config.yaml

cluster: "~/.config/snakemake/cc-slurm/slurm-submit.py {dependencies}"
#cluster-status: "~/.config/snakemake/cc-slurm/slurm-status.py"  #not fully tested yet.. better to rely on filesystem checks..
jobscript: "slurm-jobscript.sh"
jobs: 500
verbose: False
use-singularity: True
use-envmodules: True
singularity-args: ' -e '
immediate-submit: True
notemp: True #notemp set to the value of immediate_submit (since both have to be true together)
max-jobs-per-second: 2
#max-status-checks-per-second: 0.2
#restart-times: 1
default-resources: "mem_mb=12000" #this is for ensuring grouped-jobs have enough memory (ie. sum of memory request)
local-cores: 1

singularity-prefix: 'singularity'
#notebook-listen: '$(hostname -f):8888'

rerun-incomplete: True
keep-going: True
show-failed-logs: True


#####################################
#Repeat Masker

# Install and compile RepeatMasker
# Follow instructions here: https://www.repeatmasker.org/RepeatMasker/
# In desired folder:
mkdir ~/scratch/Annotation/RepeatMasker
cd ~/scratch/Annotation/RepeatMasker

#1 Install +compile RMBlast: https://www.repeatmasker.org/rmblast/
wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.14.0/ncbi-blast-2.14.0+-src.tar.gz
wget https://www.repeatmasker.org/rmblast/isb-2.14.0+-rmblast.patch.gz

tar zxvf ncbi-blast-2.14.0+-src.tar.gz
gunzip isb-2.14.0+-rmblast.patch.gz

cd ncbi-blast-2.14.0+-src
patch -p1 < ../isb-2.14.0+-rmblast.patch

cd c++
# be sure to adjust the prefix to the location of your installation
./configure --with-mt \
--without-debug \
--without-krb5 \
--without-openssl \
--with-projects=scripts/projects/rmblastn/project.lst \
--prefix=/home/celphin/scratch/Annotation/RepeatMasker/rmblast

make
#make -j # can be used to parallelize the build on multiprocessor systems, e.g. make -j2 to dedicate two cores to the build process.
make install

#2 Get TRF: https://github.com/Benson-Genomics-Lab/TRF/blob/master/README.md#instructions-for-compiling
cd ~/scratch/Annotation/RepeatMasker
git clone https://github.com/Benson-Genomics-Lab/TRF.git
cd TRF
mkdir build; cd build
../configure
make

#3 Configure repeat masker
cd ~/scratch/Annotation/RepeatMasker
wget https://www.repeatmasker.org/RepeatMasker/RepeatMasker-4.1.5.tar.gz
gunzip RepeatMasker-4.1.5.tar.gz
tar xvf RepeatMasker-4.1.5.tar

module load StdEnv/2020
module load python/3.10.2
module load scipy-stack/2022a
module load gcc/9.3.0
module load cuda/11.7
module load hdf5/1.12.1
module load parasail/2.5
module load  h5py/3.6.0

cd ~/scratch/Annotation/RepeatMasker/RepeatMasker
perl ./configure

# gives prompts below
# needs to be the full path not using the ~/scratch
The full path including the name for the TRF program.
TRF_PRGM:

/home/celphin/scratch/Annotation/RepeatMasker/TRF/build/src/trf
/home/msandler/scratch/Annotation/RepeatMasker/TRF/build/src/trf

Add a Search Engine:
   1. Crossmatch: [ Un-configured ]
   2. RMBlast: [ Un-configured ]
   3. HMMER3.1 & DFAM: [ Un-configured ]
   4. ABBlast: [ Un-configured ]

   5. Done


Enter Selection:2

The path to the installation of the RMBLAST sequence alignment program.
RMBLAST_DIR: 

/home/celphin/scratch/Annotation/RepeatMasker/rmblast/bin/


Do you want RMBlast to be your default
search engine for Repeatmasker? (Y/N)  [ Y ]: Y

Add a Search Engine:
   1. Crossmatch: [ Un-configured ]
   2. RMBlast: [ Configured, Default ]
   3. HMMER3.1 & DFAM: [ Un-configured ]
   4. ABBlast: [ Un-configured ]

   5. Done


Enter Selection:2

The path to the installation of the RMBLAST sequence alignment program.
RMBLAST_DIR [/home/celphin/scratch/Annotation/RepeatMasker/rmblast/bin/]:

/home/celphin/scratch/Annotation/RepeatMasker/rmblast/bin/

Select default - Yes
#-------------------------
Add a Search Engine:
   1. Crossmatch: [ Un-configured ]
   2. RMBlast: [ Configured, Default ]
   3. HMMER3.1 & DFAM: [ Un-configured ]
   4. ABBlast: [ Un-configured ]

   5. Done
Enter Selection: 5
Building FASTA version of RepeatMasker.lib .....


#----------------------------------
#Update the maker rule in the snakefile, to match wherever you downloaded Repeat Masker:
cd ~/scratch/Annotation/CAP_Snakemake/

grep "/home/msandler/scratch/Oxyria/" Snakefile_Oxyria_Rheum
                sed -i 's|RepeatMasker=/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Compiler/gcc9/repeatmasker/4.1.1/RepeatMasker|RepeatMasker=/home/msandler/scratch/Oxyria/Annotation/RepeatMasker/RepeatMasker|g'  Chr{wildcards.Chrs}/maker_exe.ctl

sed -i 's|\/home\/msandler\/scratch\/Oxyria\/|\/home\/celphin\/scratch\/|g' Snakefile_Oxyria_Rheum

######################################################################
#Load R libraries
module load  StdEnv/2020 r/4.1.0
export R_LIBS_USER=/home/celphin/R/x86_64-pc-linux-gnu-library/4.1.0/

R

# run one line at a time - will have some prompts
install.packages("dplyr")

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.14")

BiocManager::install("GenomicRanges")
BiocManager::install("Biostrings")
BiocManager::install("optparse")
BiocManager::install("GenomicFeatures")
BiocManager::install("Biostrings")
BiocManager::install("ORFik")
BiocManager::install("BSgenome")
BiocManager::install("rtracklayer")

#-----------
# check
library(GenomicRanges)
library(Biostrings)
library(optparse)
library(GenomicFeatures)
library(Biostrings)
library(ORFik)
library(BSgenome)
library(rtracklayer)

######################################################################
#Running + Some trouble shooting notes

module load StdEnv/2020 intel/2020.1.217 metagenome-atlas/2.5.0 
snakemake --snakefile Snakefile_XX --profile cc-slurm 

#If maker crashes with fix-nucleotides (see log), request some memory and run within the chromosome
    #once maker -fix-nucleotides 
    #once maker -> as soon as running can start 
    #current time too little may require restart
    salloc -c1 --time 4:00:00 --mem 187000m --account def-rieseber
    module load StdEnv/2020 gcc/9.3.0 openmpi/4.0.3 maker/3.01.03
    maker -fix-nucleotides
    #test
    maker 
