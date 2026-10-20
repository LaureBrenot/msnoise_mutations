#!/bin/bash
#SBATCH --job-name=msn_wct
#SBATCH --output=logs/wct_%A_%a.out
#SBATCH --array=1-40%20          # 40 workers, at most 20 running at once
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4        # 4 msnoise processes per worker
#SBATCH --mem-per-cpu=3000
#SBATCH --time=12:00:00
#SBATCH --signal=B:USR1@900      # warn the script 15 min before the time limit

# ---- environment ----------------------------------------------------------
# One thread per process: numpy/scipy/pycwt would otherwise each spawn as many
# threads as the node has cores and fight for CPUs already used by others.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1

source ~/miniconda3/etc/profile.d/conda.sh
conda activate msnoise_stable
cd /globalscratch/users/c/c/ccaudron        # MSNoise project folder (db.ini)

# ---- stop starting new rounds when the time limit approaches ---------------
STOP=0
trap 'echo "Time limit close, finishing current round"; STOP=1' USR1

# ---- check the database is reachable from this compute node ----------------
python - <<'EOF' || { echo "Cannot reach the MSNoise database from $(hostname)"; exit 1; }
from msnoise.api import connect
from sqlalchemy import text
db = connect(); db.execute(text("select 1")); db.close()
EOF

has_jobs() {
python - <<'EOF'
import sys
from msnoise.api import connect, is_dtt_next_job
db = connect()
ok = is_dtt_next_job(db, flag='T', jobtype='WCT')
db.close()
sys.exit(0 if ok else 1)
EOF
}

# ---- work loop -------------------------------------------------------------
# Each msnoise process exits after 1000 jobs or 3 empty batches, so keep
# relaunching while WCT jobs remain to be done.
while [ "$STOP" -eq 0 ] && has_jobs; do
    msnoise -t "$SLURM_CPUS_PER_TASK" -d 2 cc dvv compute_wct -b 5 &
    wait $!        # 'wait' lets the USR1 trap fire while msnoise runs
    wait $! 2>/dev/null
done

echo "Worker $SLURM_ARRAY_TASK_ID done on $(hostname)"
