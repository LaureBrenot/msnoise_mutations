#!/bin/bash
# Submit the WCT workers, then a small cleanup job that runs once all workers
# have ended (finished, failed or hit the time limit). The cleanup puts jobs
# left "I"n progress by killed workers back to "T"odo, so a resubmission
# picks them up. Jobs flagged "E"rror are left alone (they failed for real).
#
# Usage:  conda activate msnoise_stable; bash submit_wct.sh

PROJECT=/globalscratch/users/c/c/ccaudron
mkdir -p "$PROJECT/logs"
cd "$PROJECT"

jid=$(sbatch --parsable wct_array.sh)
echo "Submitted WCT workers: $jid"

cid=$(sbatch --parsable --dependency=afterany:$jid \
      --job-name=msn_wct_reset --output=logs/wct_reset_%j.out \
      --ntasks=1 --cpus-per-task=1 --mem=1000 --time=0:10:00 \
      --wrap "cd $PROJECT && msnoise reset WCT --no-reset-e && msnoise info -j")
echo "Submitted cleanup: $cid (runs after $jid)"
