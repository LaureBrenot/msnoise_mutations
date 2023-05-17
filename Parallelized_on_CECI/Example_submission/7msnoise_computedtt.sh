#!/bin/sh
#SBATCH --job-name=msn_dtt
#SBATCH --output=res_msn_dtt.txt
#
#SBATCH --ntasks=1
#SBATCH --time=60:00
#SBATCH --mem-per-cpu=2000
#
#SBATCH --array=1-250
#
#SBATCH --mail-user=laure.brenot@ulb.be
#SBATCH --mail-type=ALL
#
source /home/ulb/gtime/lbrenot/.bashrc
conda activate msnoise_env
#srun msnoise compute_dtt
srun msnoise daily_compute_dtt
