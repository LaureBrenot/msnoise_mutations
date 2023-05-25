#!/bin/sh
#SBATCH --job-name=msn_stack
#SBATCH --output=res_msn_stack.txt
#
#SBATCH --ntasks=1
#SBATCH --time=50:00
#SBATCH --mem-per-cpu=2000
#
#SBATCH --array=1-250
#
#SBATCH --mail-user=@
#SBATCH --mail-type=ALL
#
source /home/ulb/gtime/lbrenot/.bashrc
conda activate msnoise_env
#srun msnoise stack -r
#srun msnoise reset STACK
srun msnoise stack -m

