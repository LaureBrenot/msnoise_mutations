# Move from msnoise_c back to stable msnoise while keeping the CCs you already computed
**Why**: This moves the project back to the prior version of msnoise **without recomputing the cross-correlations**. The msnoise_c (msnoise_current) project (database + OUTPUT folder) is only read, never modified, so you can always go back.
However, if cc are only covering a few years and not on an extended network, this may be shorter than just starting from scratch in the DB with a cloned environment that has the msnoise folder overwritten with the prior one.
**What carries over**: CCs (daily stacks + keep_all windows), config values, filters, stations, data_availability, CC job status.
**What is recomputed by prior msnoise**: REF/MOV stacks, MWCS, DTT, WCT, stretching. These are cheap compared to the CCs.

> **You don't need to know msnoise_c to do this safely.** See the section *"How to be sure it worked"* further down: (1) nothing can be lost, (2) three questions to ask before you start, (3) a `verify` command that checks everything automatically and writes a PASS/FAIL report you can hand over.

The `change_without_kill.md` copied tables one to one. That no longer works here: msnoise_c changed the schema (config sets, a `lineages` table, `workflow_steps`, no `filters` table, relative data paths) and the file layout (NetCDF under `OUTPUT/preprocess_1/cc_1/filter_1/_output/...`). The script `migrate_c_to_stable.py` does the translation.

| | msnoise_c | stable msnoise |
|---|---|---|
| daily CC | `OUTPUT/preprocess_1/cc_1/filter_1/_output/daily/ZZ/S1_S2/DAY.nc` | `STACKS/01/001_DAYS/ZZ/S1_S2/DAY.MSEED` |
| keep_all windows | `OUTPUT/preprocess_1/cc_1/filter_1/_output/all/ZZ/S1_S2/DAY.nc` | `CROSS_CORRELATIONS/01/S1/S2/ZZ/DAY.h5` |
| filters | config category `filter`, sets 1..N | `filters` table, ref 1..N |
| config | `(name, category, set_number)` | `(name, value)` |
| CC jobs | `jobtype = cc_1` + lineage | `jobtype = CC` |
| DA path | relative to the data source URI | absolute |

---

## Step 0: Stop msnoise_c and back up its DB

### Before you start: three questions
Run `inspect` (Step 3) and check the output with these questions:
1. **Is msnoise_c stopped?** Nobody should be computing while you migrate (`msnoise info -j` should show no `I` jobs).
2. **Which CC setup is the "real" one?** If `inspect` lists several `cc` or `preprocess` sets, the stable version can only keep one. Get the number (`--cc-set N`, `--preprocess-set N`).
```
conda activate msnoise_c
python migrate_c_to_stable.py inspect --src /path/to/msnoise_c_project > inspect_output.txt
cat inspect_output.txt
```
The output looks like this:
```
Config sets found:
  cc                 set=1   21 params  maxlag=60 cc_sampling_rate=20 corr_duration=1800 keep_all=Y keep_days=Y comps=ZZ
  filter             set=1    6 params  freqmin=0.1 freqmax=1.0 used=Y
  filter             set=2    6 params  freqmin=1.0 freqmax=4.0 used=Y
  preprocess         set=1    9 params
  ...
Jobs per step / flag:
  cc_1                 D        9
  preprocess_1         D        9
  stack_1              T       18
Data sources:
  ref=1 name=local uri=/path/to/SDS structure=SDS
CC output folders on disk:
  OUTPUT/preprocess_1/cc_1/filter_1/_output/daily  -> 9 files
  OUTPUT/preprocess_1/cc_1/filter_1/_output/all  -> 9 files
  ...
```
4. **Are these the right filters?** `inspect` lists every `filter` set with its frequency band. Each one becomes stable filter N (folder `0N`).

With those three answers, nothing else about msnoise_c needs deciding by you.\

! Use appropriate IP ! 
check with
```
python -c "import pickle; print(pickle.load(open('db.ini','rb')))"
```
You'll get something like:
```
[2, '10.44.2.1:5020', 'static_c', 'root', 'noise', '']
```
That list is tech, host:port, database, user, password and table prefix.

Make sure no msnoise_c job is running (check with `msnoise info -j`), then:
```
./mariadb/bin/mysqldump -h 10.44.2.1 -P 5020 -u root -p static_c > static_c_backup.sql
```
(replace `static_c` with your msnoise_c database name)

## Step 1: Create the empty database for the stable version and the environement
```./mariadb/bin/mysql```
```
CREATE DATABASE IF NOT EXISTS static_stable;
```
Create or clone an environment and check it has no msnoise already
```
conda create -n msnoise_stable --clone msnoise_clean
conda activate msnoise_stable
python -c "import msnoise; print(msnoise.__file__)"
```
Replace the unzip msnoise.zip folder in the site-packages folder to replace the existing msnoise folder (=msnoise_c). Usually path is miniconda/env/env_name/lib/python3.smth/site-package/msnoise

## Step 2: Create the stable project folder and initialise it with the STABLE msnoise
Use the conda env that has the stable `msnoise.zip` installed (e.g. `msnoise_clean`), **not** the msnoise_c env.
```
conda activate msnoise_clean
mkdir /globalscratch/ulb/gtime/lbrenot/msnoise_stable
cd /globalscratch/ulb/gtime/lbrenot/msnoise_stable
msnoise db init
```
Answer: mysql, host `10.44.2.1:5050`, database `static_stable`, user `root`, password, empty prefix.
This creates the stable tables and the default config. Do **not** add filters or stations; the script does that.

Copy `migrate_c_to_stable.py` into this folder.

## Step 3: Inspect what is in the msnoise_c project (read-only)
```
python migrate_c_to_stable.py inspect --src /path/to/msnoise_c_project
```
It lists the config sets, lineages, job counts per step, and the CC folders on disk with file counts.
Check:
- how many `cc` / `preprocess` sets you have. Stable msnoise supports only one CC configuration, so pick one with `--cc-set N` (and `--preprocess-set N` if you have several).
- that the `filter` sets are the ones you want. Each `filter_N` becomes filter ref `N`, folder `0N`.

## Step 4: Fill the stable DB
```
python migrate_c_to_stable.py db \
    --src /path/to/msnoise_c_project \
    --dst . \
    --cc-set 1 \
    --classic-ref /globalscratch/ulb/gtime/lbrenot/msnoise__static
```
`--classic-ref` points to your **old stable project** (the one with the `static` DB). Your custom `s08compute_wct.py` reads `coda_safety_factor` and `ref_type`, which `msnoise db init` does not create. The script copies them from there. Without it, it inserts `coda_safety_factor=1.2` and `ref_type=static` and prints a warning.

What it does:
- **config**: each stable parameter takes its value from the msnoise_c global config first, then preprocess, cc, stack, refstack, mwcs, mwcs_dtt, wavelet, wavelet_dtt, etc. Renamed parameters are mapped (`dtt_maxdtt → dtt_maxdt`, `wct_codacycles → dtt_codacycles`, `wct_min_nonzero → dvv_min_nonzero`). `data_folder`, `data_structure`, `network` and `channels` come from the msnoise_c data source. `output_folder` is set to `CROSS_CORRELATIONS`.
- **filters**: from `filter_N` sets. `mwcs_low/high/wlen/step` come from the mwcs set linked to that filter.
- **stations**: copied.
- **data_availability**: paths made absolute, flag set to `A`, so `new_jobs` does not recreate CC jobs for data already processed.
- **jobs**: `cc_1` jobs become `CC` jobs with the same flag. Each done CC also gets a `STACK` job set to `T`.

It prints which parameters it copied and which stayed at stable defaults. Check with:
```
msnoise info
```
To redo this step: `--force`.

## Step 5: Convert the CC files
Check the disk space first. The keep_all `.h5` files take about as much space as the `_output/all` NetCDF files:
```
du -sh /path/to/msnoise_c_project/OUTPUT/preprocess_1/cc_1
```
Do a dry run:
```
python migrate_c_to_stable.py files --src /path/to/msnoise_c_project --dst . --cc-set 1 --dry-run
```
Then run it (preferably in a Slurm job, it can take a while):
```
python migrate_c_to_stable.py files --src /path/to/msnoise_c_project --dst . --cc-set 1 --threads 8
```
- It is **resumable**: if the job is killed, rerun the same command. Files already converted are skipped.
- It stops before converting anything if the lag axis of the files does not match `maxlag` / `cc_sampling_rate` in the stable DB.
- If disk space is tight, use `--keep-days-only` in **both** steps 4 and 5. This only converts the daily stacks and sets `keep_all=N`, so stacking works on daily CCs.

## Step 5b: Verify (automatic)
```
python migrate_c_to_stable.py verify --src /path/to/msnoise_c_project --dst . --cc-set 1 --threads 8
```
(add `--keep-days-only` if you used it before). It must end with `RESULT : PASS`. The full report goes to `migration_check.txt`. See below for what it checks.

## Step 6: Continue with stable msnoise
```
msnoise cc stack -r
msnoise reset STACK
msnoise cc stack -m
msnoise cc dvv compute_mwcs      # and/or
msnoise cc dvv compute_wct
msnoise cc dvv compute_dtt
```
For new data, use the normal routine: `msnoise scan_archive`, `msnoise new_jobs`, `msnoise cc compute_cc`. Only new or modified files create CC jobs.

---

## How to be sure it worked (when you have never used msnoise_c)

### 1. Nothing can be lost
- The script only **reads** the msnoise_c database and `OUTPUT` folder. It never writes to them.
- Step 0 also gives you a SQL backup of the msnoise_c DB.
- If anything looks wrong, drop `static_stable`, delete the stable folder, and start again from Step 1. The worst case is lost time, never lost CCs.

### 2. After the migration: the `verify` command (Step 5b)
It compares the stable project with the msnoise_c one and checks:

| check | what it proves |
|---|---|
| config | the parameters that define the CCs (maxlag, sampling rate, corr_duration, whitening, winsorizing, components, stack method, dates, mov_stack, ref period…) are identical in both |
| filters | same filters, same frequency bands |
| stations | same stations, same `used` flags |
| data_availability | same number of rows, and a sample of 50 raw data files exists at the new absolute paths (so new data scans will work) |
| jobs | every day/pair that msnoise_c finished is marked done in stable (so nothing gets recomputed, and nothing is missed) |
| files | every msnoise_c CC file has its stable counterpart, and no stray files from another setup |
| lag axis | the files' lag axis matches the stable maxlag / cc_sampling_rate |
| content | each converted file holds **exactly the same numbers**, window times and lag axis as the msnoise_c file (bit-for-bit comparison, not approximate) |

On a very large project you can check the numbers in a random subset first with `--sample 2000`, then run the full check overnight.

**How I tested `verify`:** I ran it on a copy of my test migration where I broke 8 things on purpose: deleted one CC file, changed one number by 0.000001 in one file, added a stray file, changed maxlag, changed a filter band, disabled a station, broke the data paths, and reset 2 finished jobs. It caught all 8 and returned `FAIL`. On the correct migration it returns `PASS`.

### 4. What to show
- `migration_check.txt` (the PASS report).
- Once the stable version has run stack/MWCS/DTT (or WCT) on a few pairs, compare those dv/v curves with any dv/v plots you already have from msnoise_c. The CCs are byte-identical, so any difference comes from the stable version's post-processing, not from the migration.

---

## Things to know
- **Tested end to end** on MariaDB with a synthetic 3-station, 3-day, 2-filter project: msnoise_c computed the CCs, the script migrated them, and stable msnoise then ran stack -r/-m, MWCS, DTT and WCT on them. `scan_archive` + `new_jobs` recreated 0 jobs.
- **Formats identical**: the converted files have exactly the same layout, number of samples, sampling rate, lag axis and window time labels as files written natively by stable msnoise.
- **Values are close but not identical**: on the same data, CCs computed by msnoise_c and by stable msnoise are 0.994–0.998 correlated (processing differs slightly; msnoise_c also keeps the last 30-min window of the day, 48 vs 47). So:
  - stacking, MWCS and WCT on the migrated CCs is consistent (all CCs come from one version);
  - if you compute **new days** with stable msnoise and append them to msnoise_c CCs, check a few overlap days. Recompute 2–3 days with stable msnoise in a scratch copy and compare dv/v before trusting the join. If it matters, recompute everything with one version.
- WCT on the **first day** fails with "No data for …" in stable msnoise. This also happens when stable msnoise computes the CCs itself: its 1D moving stack is labelled on the next day. It is not caused by the migration.
