#!/bin/sh
#SBATCH --job-name=msn_scan
#SBATCH --output=res_msn_scan.txt
#
#SBATCH --ntasks=1
#SBATCH --time=30:00
#SBATCH --mem-per-cpu=500
#
#SBATCH --cpus-per-task=8
#
#SBATCH --mail-user=laure.brenot@ulb.be
#SBATCH --mail-type=ALL
#
source /home/ulb/gtime/lbrenot/.bashrc 
#cd ./Msnoise_para/Aso_msn
conda activate msnoise_env
srun msnoise -t 8 scan_archive --init
srun msnoise db execute "update stations set used_location_codes='--'"
