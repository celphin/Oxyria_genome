#!/usr/bin/env bash

eval "$(conda shell.bash hook)"

conda activate snakemake

snakemake --cluster "sbatch -t {cluster.time} --mem {cluster.mem} -N {cluster.nodes} --ntasks-per-node {cluster.ntask} --cpus-per-task {cluster.cpus} -o {cluster.output} -e {cluster.error}" --cluster-config cluster_config.yaml --jobs 48 --rerun-incomplete -R MAKER3

conda deactivate
