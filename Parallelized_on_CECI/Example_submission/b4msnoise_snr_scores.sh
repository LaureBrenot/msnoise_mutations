#!/bin/sh
#SBATCH --job-name=msn_snr_score
#SBATCH --output=res_msn_snr_score.txt
#
#SBATCH --ntasks=1
#SBATCH --time=40:00:00
#SBATCH --mem-per-cpu=10000
#
#SBATCH --mail-user=laure.brenot@ulb.be
#SBATCH --mail-type=ALL
#
source /home/ulb/gtime/lbrenot/.bashrc
conda activate msnoise_env
srun msnoise snr_score

