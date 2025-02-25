#!/bin/bash
#SBATCH --account=rpp-rieseber
#SBATCH --job-name=interproscan
#SBATCH --time=1-12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=120000m

#Required parameters (in order): 
    #inputfile name
    #outputfile name
#to run sh run_interproscan.sh inputfile.fasta outputfile.tsv

module load StdEnv/2020
module load interproscan/5.64-96.0

input_file=$1
output_file=$2
srun interproscan.sh -i "$input_file"  -f tsv -cpu 30 -o "$output_file" -dp  --goterms
