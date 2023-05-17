#!/bin/sh
#SBATCH --job-name=msn_cc
#SBATCH --output=res_msn_cc.txt
#
#SBATCH --ntasks=1
#SBATCH --time=50:00
#SBATCH --mem-per-cpu=2000
#
#SBATCH --array=1-250
#
#SBATCH --mail-user=laure.brenot@ulb.be
#SBATCH --mail-type=ALL
#
source /home/ulb/gtime/lbrenot/.bashrc
conda activate msnoise_env
srun msnoise compute_cc

