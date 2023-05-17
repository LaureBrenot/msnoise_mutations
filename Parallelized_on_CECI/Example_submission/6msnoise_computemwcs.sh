#!/bin/sh
#SBATCH --job-name=msn_mwcs
#SBATCH --output=res_msn_mwcs.txt
#
#SBATCH --ntasks=1
#SBATCH --time=06:20:00
#SBATCH --mem-per-cpu=2000
#
#SBATCH --array=1-250
#
#SBATCH --mail-user=laure.brenot@ulb.be
#SBATCH --mail-type=ALL
#
source /home/ulb/gtime/lbrenot/.bashrc
conda activate msnoise_env
#srun msnoise -t 4 compute_mwcs
srun msnoise -t 4 daily_compute_mwcs
