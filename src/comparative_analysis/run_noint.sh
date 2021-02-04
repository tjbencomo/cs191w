#!/bin/bash
#SBATCH --job-name=run_noint
#SBATCH --output=/scratch/users/tbencomo/logs/run_noint.out
#SBATCH --nodes=1
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=64000
#SBATCH --mail-type=END
#SBATCH --mail-user=tbencomo@stanford.edu

set -euo pipefail
ml R/4.0.2
cd /scratch/groups/carilee/cs191w/
Rscript src/comparative_analysis/no_integration.R
