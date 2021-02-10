#!/bin/bash
#SBATCH --job-name=run_harmony
#SBATCH --output=/scratch/users/tbencomo/logs/run_harmony.out
#SBATCH --nodes=1
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=20000
#SBATCH --mail-type=END
#SBATCH --mail-user=tbencomo@stanford.edu

set -euo pipefail
ml R/4.0.2
cd /scratch/groups/carilee/cs191w/
Rscript src/comparative_analysis/harmony_integration.R
