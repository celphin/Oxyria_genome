##########################
# Install OrthoFinder and GeneSpace on Beluga
# Feb 2024
##################################################
# OrthoFinder

# Simlar papers
# https://www.nature.com/articles/s42003-023-05044-1#Sec9
# https://onlinelibrary.wiley.com/doi/full/10.1111/1755-0998.13280

# Orthofinder
# https://davidemms.github.io/
# https://github.com/davidemms/OrthoFinder

# Genespace Paper: https://elifesciences.org/articles/78526 
# https://github.com/jtlovell/GENESPACE

################################
# Installation

# needs  OrthoFinder, MCScanX and R

cd /home/celphin/scratch/Oxyria/GeneSpace

#---------------------------
# OrthoFinder Installation
# https://github.com/davidemms/OrthoFinder#installing-orthofinder-on-linux

# No python
wget https://github.com/davidemms/OrthoFinder/releases/download/2.5.5/OrthoFinder.tar.gz
tar xzf OrthoFinder.tar.gz

module load StdEnv/2020 python/3.11.5 scipy-stack/2021a

# Test
./OrthoFinder/orthofinder -h
./OrthoFinder/orthofinder -f /OrthoFinder/ExampleData/

#[88567] Error loading Python lib '/tmp/_MEIVMHG2R/libpython3.10.so.1.0': dlopen: /lib64/libm.so.6: version `GLIBC_2.35' not found (required by /tmp/_MEIVMHG2R/libpython3.10.so.1.0)

#---------------------------
# with python 

wget https://github.com/davidemms/OrthoFinder/releases/download/2.5.5/OrthoFinder_source.tar.gz
tar xzf OrthoFinder_source.tar.gz

module load StdEnv/2020 python/3.11.5 scipy-stack/2021a

python OrthoFinder_source/orthofinder.py -h

# Results:
    # /project/6003374/Dryas_shared_data/Oxyria/GeneSpace/OrthoFinder_source/ExampleData/OrthoFinder/Results_Jan18/

#---------------------------------
# to run
module load StdEnv/2020 python/3.11.5 scipy-stack/2021a
cd /home/celphin/scratch/Oxyria/GeneSpace/
python OrthoFinder_source/orthofinder.py -f .../Protein_files/

#--------------------
#MCScanX
# https://github.com/wyp1125/MCScanX

git clone https://github.com/wyp1125/MCScanX.git

module load StdEnv/2023 java/17.0.6
cd MCScanX
make

#-----------------------------
# load R
module load StdEnv/2020 r/4.2.2 glpk/5.0

R

# install library
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("jtlovell/GENESPACE")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("Biostrings", "rtracklayer"))

library(GENESPACE)
# failed dependencies ‘XML’, ‘restfulr’ are not available for package ‘rtracklayer’

install.packages("restfulr")
rm -r /home/celphin/R/x86_64-pc-linux-gnu-library/4.2/00LOCK-XML
install.packages("XML")

# try again 
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("Biostrings", "rtracklayer"))

BiocManager::install("Biostrings",force = TRUE)

install.packages("igraph") # needs glpk

devtools::install_github("jtlovell/GENESPACE")
library(GENESPACE)

# GENESPACE v1.3.1: synteny and orthology constrained comparative genomics

###################################
# after initial start

# Could not find a valid path to the orthofinder program from R. To run  orthofinder, ensure that the 
# orthofinder program is in the         $PATH, then call the following from the shell: 
 orthofinder -f /home/celphin/scratch/Oxyria/GeneSpace/tmp -t 16 -a 1 -X \
 -o /home/celphin/scratch/Oxyria/GeneSpace/orthofinderError 
# in run_orthofinder(gsParam = gsParam, verbose = TRUE) :

# Once OrthoFinder has been run, re-call run_genespace

# check PATH
echo $PATH

# This sets your PATH variable to the existing PATH plus what you add to the end. 
# Check that it has been added 
# (Caveat: it presist only in the current session of the terminal):
export PATH=$PATH:/home/celphin/scratch/Oxyria/GeneSpace/OrthoFinder_source/orthofinder.py

# check PATH again
echo $PATH

# to make permanent add export PATH=$PATH:/path/to/my/program to:
nano ~/.bashrc
export PATH=$PATH:/home/celphin/scratch/Oxyria/GeneSpace/OrthoFinder/orthofinder
#alias orthofinder='python /home/celphin/scratch/Oxyria/GeneSpace/OrthoFinder_source/orthofinder.py'

# try running from bash?
 orthofinder -f /home/celphin/scratch/Oxyria/GeneSpace/tmp -t 16 -a 1 -X \
 -o /home/celphin/scratch/Oxyria/GeneSpace/orthofinder 
# works

################################