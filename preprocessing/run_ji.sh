#!/bin/bash
#SBATCH --job-name=run_ji
#SBATCH --output=/scratch/users/tbencomo/logs/run_ji.out
#SBATCH --nodes=1
#SBATCH --time=03-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=300
#SBATCH --qos=long
#SBATCH --mail-type=END
#SBATCH --mail-user=tbencomo@stanford.edu


set -e
cd /scratch/groups/carilee/cs191w/preprocessing
snakemake -s preproc_ji.smk targets --cluster-config cluster.json -j 499 \
    --use-conda --use-singularity \
    --cluster 'sbatch -p {cluster.partition} -t {cluster.time} --mem {cluster.mem} -c {cluster.ncpus} -o {cluster.out}'
