#!/bin/bash
#SBATCH --job-name=msn_wct
#SBATCH --output=logs/wct_%A_%a.out
#SBATCH --array=1-40%20
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=3000
#SBATCH --time=12:00:00
#SBATCH --signal=B:USR1@900
#
# MSNoise WCT on Slurm, single file.  Submit with:
#
#     mkdir -p logs; sbatch msnoise_wct.sh
#
# 40 workers (at most 20 running at once), each running 4 msnoise processes.
# Array task 1 also queues a small cleanup job that runs once every worker
# has ended (done, failed or killed at the time limit): it puts WCT jobs left
# "I"n progress back to "T"odo, so a later resubmission picks them up.
# USR1 is sent 15 min before the time limit so workers stop starting new work.

PROJECT=/globalscratch/users/c/c/ccaudron     # MSNoise project folder (db.ini)

source ~/miniconda3/etc/profile.d/conda.sh
conda activate msnoise_stable
cd "$PROJECT" || exit 1

# One thread per process: avoids oversubscribing CPUs shared with other users
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1

# ---------------------------------------------------------------------------
# Cleanup role (queued by task 1, see below)
# ---------------------------------------------------------------------------
if [ "$WCT_ROLE" = "cleanup" ]; then
    msnoise reset WCT --no-reset-e
    msnoise info -j
    exit 0
fi

# ---------------------------------------------------------------------------
# Task 1 queues the cleanup job, dependent on the whole array
# ---------------------------------------------------------------------------
if [ "$SLURM_ARRAY_TASK_ID" = "1" ]; then
    sbatch --dependency=afterany:"$SLURM_ARRAY_JOB_ID" \
           --array=1 --job-name=msn_wct_reset \
           --output=logs/wct_reset_%j.out \
           --cpus-per-task=1 --mem-per-cpu=1000 --time=0:10:00 \
           --export=ALL,WCT_ROLE=cleanup \
           "$PROJECT/msnoise_wct.sh" \
      && echo "Cleanup job queued after array $SLURM_ARRAY_JOB_ID"
fi

# ---------------------------------------------------------------------------
# Worker
# ---------------------------------------------------------------------------
STOP=0
trap 'echo "Time limit close, finishing current round"; STOP=1' USR1

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

# Each msnoise process exits after 1000 jobs or 3 empty batches:
# relaunch while WCT jobs remain to be done.
while [ "$STOP" -eq 0 ] && has_jobs; do
    msnoise -t "$SLURM_CPUS_PER_TASK" -d 2 cc dvv compute_wct -b 5 &
    wait $!
    wait $! 2>/dev/null
done

echo "Worker $SLURM_ARRAY_TASK_ID done on $(hostname)"
