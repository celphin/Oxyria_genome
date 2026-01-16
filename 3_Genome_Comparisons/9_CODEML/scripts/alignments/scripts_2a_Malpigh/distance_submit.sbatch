#!/bin/sh
#SBATCH --account=nn9370k
#SBATCH --time=160:0:0
#SBATCH --mem-per-cpu=200M
source /cluster/bin/jobsetup
module purge   # clear any inherited modules
set -o errexit # exit on errors

arrayrun 001-160 distance_worker.sbatch
