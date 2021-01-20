#!/bin/bash
#SBATCH --job-name=run_snake
#SBATCH --output=/scratch/users/tbencomo/logs/run_snake.out
#SBATCH --nodes=1
#SBATCH --time=03-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=300
#SBATCH --qos=long
#SBATCH --mail-type=END
#SBATCH --mail-user=tbencomo@stanford.edu


set -e
cd /home/users/tbencomo/cs191w/preprocessing
snakemake --cluster-config cluster.json -j 499 \
    --use-conda --use-singularity \
    --cluster 'sbatch -p {cluster.partition} -t {cluster.time} --mem {cluster.mem} -c {cluster.ncpus} -o {cluster.out}'
