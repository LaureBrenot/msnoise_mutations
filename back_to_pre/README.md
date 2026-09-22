# MSNoise older version: install and reuse existing cross-correlations

One of the stable MSNoise version (`msnoise/`, reports itself as **1.6.post1**) and two helper scripts:

| script | use it when |
|---|---|
| `clone_classic.py` | your existing project is a **classic** MSNoise project (tables `config, filters, jobs, stations, data_availability`). This is the usual case. |
| `migrate_c_to_stable.py` | your existing project was made with **msnoise_c** (the new "workflow/lineage" version, with an `OUTPUT\preprocess_1\cc_1\...` folder). |

Both scripts only **read** the original project. It is never modified.

All commands below are for **Windows, Anaconda Prompt**. Replace everything in `<...>`.

Example names used below:
- original project folder: `E:\Kilauea` (contains its `db.ini`)
- original database: `kilauea2`
- new project folder: `E:\Kilauea\stable_dir_LB`
- new database: `kilauea_stable`
- MySQL/MariaDB user / password: `something` / `something`

---

## 1. Install the stable MSNoise

### 1.1 Create a new environment
Use a **new** environment. An environment that already has another msnoise (e.g. the official 1.6.5 or msnoise_c) causes confusing errors.
```
conda create -n msnoise_stable python=3.11
conda activate msnoise_stable
```

### 1.2 Install git and the dependencies
```
conda install -c conda-forge git
conda install -c conda-forge xarray netcdf4 pytables pandas scipy numpy obspy matplotlib sqlalchemy sqlalchemy-utils pymysql click logbook flask flask-admin flask-wtf markdown "setuptools<81"
pip install pycwt
```

### 1.3 Remove any other msnoise, then install this one
```
pip uninstall msnoise
```
Repeat until it says `not installed`. If `conda list msnoise` shows one, run `conda remove msnoise` too. Then:
```
pip install --no-deps "git+https://github.com/LaureBrenot/msnoise_mutations.git#subdirectory=back_to_pre"
```

### 1.4 Check the installation
Run this **from a folder that does not contain an `msnoise` folder** (e.g. `cd %USERPROFILE%`):
```
pip show msnoise
python -c "import msnoise; print(msnoise.__file__)"
python -c "import xarray, netCDF4, tables, obspy, pandas, sqlalchemy, sqlalchemy_utils, pymysql, click, logbook, flask, flask_admin, pycwt; print('all OK')"
msnoise --help
```
Expected:
- `Version: 1.6.post1` (if it says `1.6.5`, the official msnoise is still installed: go back to 1.3)
- a path inside `...\envs\msnoise_stable\Lib\site-packages\msnoise\`
- `all OK`
- the msnoise command list

A `UserWarning: pkg_resources is deprecated` message is normal and can be ignored.

To update later: add `--force-reinstall` to the pip install command in 1.3.

---

## 2. Identify the existing project

### 2.1 Read its db.ini
```
cd <ORIGINAL_FOLDER>
python -c "import pickle; print(pickle.load(open('db.ini','rb')))"
```
Output: `[tech, host, database, user, password, prefix]`, e.g. `[2, 'localhost', 'kilauea2', 'something', 'something', '']`
- `1` = SQLite (the database is a file in the folder)
- `2` = MySQL/MariaDB (a database server, here on this computer since host is `localhost`)

### 2.2 Look at its tables and jobs (MySQL example)
```
python -c "import pymysql; c=pymysql.connect(host='localhost', user='<USER>', password='<PASSWORD>', database='<OLD_DB>'); cur=c.cursor(); cur.execute('SHOW TABLES'); print(cur.fetchall()); cur.execute('SELECT jobtype, flag, COUNT(*) FROM jobs GROUP BY jobtype, flag'); print(cur.fetchall())"
```
- Tables `config, data_availability, filters, jobs, stations` and jobtypes `CC, STACK, MWCS, WCT...` → **classic project → section 3**
- Tables like `lineages, workflow_steps, data_sources` and jobtypes `cc_1, stack_1...` → **msnoise_c project → section 5**

Job flags: `T` = to do, `I` = in progress (or interrupted), `D` = done, `E` = error.

---

## 3. Copy a classic project (keep the CCs, new package, new folder)

### 3.0 Stop and back up
Make sure nothing is running on the original project. Then:
```
mysqldump -u <USER> -p <OLD_DB> > E:\<OLD_DB>_backup.sql
```
If `mysqldump` is not found, use its full path, e.g. `"C:\Program Files\MariaDB 11.x\bin\mysqldump.exe"`.

### 3.1 Create the new database
```
python -c "import pymysql; c=pymysql.connect(host='localhost', user='<USER>', password='<PASSWORD>'); c.cursor().execute('CREATE DATABASE IF NOT EXISTS <NEW_DB>'); print('created')"
```
If you get `Access denied`, run the same command with `user='root'` and the MySQL root password, then give the user access:
```
python -c "import pymysql; c=pymysql.connect(host='localhost', user='root', password='<ROOTPASSWORD>'); c.cursor().execute(\"GRANT ALL ON <NEW_DB>.* TO '<USER>'@'localhost'\"); print('ok')"
```
If you get `Can't connect to MySQL server`, start the MySQL/MariaDB service (Windows "Services").

**Never use the original database name here.**

### 3.2 Create the new project folder and initialise it
```
mkdir <NEW_FOLDER>
cd <NEW_FOLDER>
msnoise db init
```
Answers: `2` (mysql), server `localhost`, database `<NEW_DB>`, user, password, prefix empty (or the same prefix as the original db.ini).
It must end with `Installation Done! - Go to Configuration Step!`. Check with `msnoise info`.

### 3.3 Copy the database content
Put `clone_classic.py` in the new folder (download it from this GitHub folder), then do a dry run:
```
python clone_classic.py --src <ORIGINAL_FOLDER> --dst <NEW_FOLDER>
```
Example output:
```
config            : 77 keys (75 taken from source, 2 new keys added with default: coda_safety_factor, ref_type)
filters           : 1 rows
stations          : 21 rows
data_availability : 345198 rows
jobs              : CC D=911975, CC T=144, STACK T=911975
  output_folder  = CROSS_CORRELATIONS
  keep_all       = N
  ...
```
If it looks right, write it:
```
python clone_classic.py --src <ORIGINAL_FOLDER> --dst <NEW_FOLDER> --go
```
What it does:
- copies config, filters, stations, data_availability;
- copies CC jobs (interrupted `I`/`E` jobs become `T`, so they are recomputed);
- creates one `STACK` job (`T`) per finished CC, so the new package rebuilds the stacks;
- does not copy MWCS/STR/WCT/PSD jobs (the new package recreates them);
- adds config keys the new package needs that the old database lacks, with default values (listed in the output).

`coda_safety_factor` (default 1.2) and `ref_type` (default `static`) are only used by `compute_wct`. Check them with whoever owns the project. To change one:
```
python -c "import pymysql; c=pymysql.connect(host='localhost',user='<USER>',password='<PASSWORD>',database='<NEW_DB>'); c.cursor().execute(\"UPDATE config SET value='moving' WHERE name='ref_type'\"); c.commit(); print('ok')"
```
(`msnoise config set` refuses these two names, because they are not in the default list.)

### 3.4 Find which filter folder holds the CCs
The folder number in `STACKS` is the filter **ref**. Old folders from deleted filters may still be there.
```
python -c "import pymysql; c=pymysql.connect(host='localhost',user='<USER>',password='<PASSWORD>',database='<OLD_DB>'); cur=c.cursor(); cur.execute('SELECT * FROM filters'); print(cur.fetchall())"
```
The first number of each row is the ref (e.g. `3` → folder `03`). Check that folder holds about as many files as the finished CC jobs:
```
dir /s /b <ORIGINAL_FOLDER>\STACKS\03\001_DAYS | find /c ".MSEED"
```

### 3.5 Copy the CC files
Check free disk space first (`dir <ORIGINAL_FOLDER>\STACKS\03\001_DAYS /s`, total at the bottom).

Daily CCs (always), once per filter ref:
```
robocopy <ORIGINAL_FOLDER>\STACKS\03\001_DAYS <NEW_FOLDER>\STACKS\03\001_DAYS /E
```
Per-window CCs: **only if `keep_all = Y`** in the dry run. Copy the `output_folder` folder:
```
robocopy <ORIGINAL_FOLDER>\CROSS_CORRELATIONS <NEW_FOLDER>\CROSS_CORRELATIONS /E
```
If `output_folder` was a full path, point the copy to its own folder: `msnoise config set output_folder=CROSS_CORRELATIONS`.

The other folders (`STACKS\..\REF`, `STACKS2`, `MWCS2`, `DTT2`, `WCT`...) are **not** needed: the new package recomputes them.

### 3.6 Check and run
Always run msnoise **from the new folder**. The `db.ini` in the current folder decides which database is used.
```
cd <NEW_FOLDER>
msnoise info
msnoise info -j
msnoise cc compute_cc
msnoise cc stack -r
msnoise reset STACK
msnoise cc stack -m
msnoise cc dvv compute_mwcs
msnoise cc dvv compute_dtt
msnoise cc dvv compute_wct
```
- `compute_cc` only computes the CC jobs still `T` (interrupted ones and new data); finished CCs are reused.
- With `keep_all = N`, stacking warns that the sampling will be 1 day. That is expected.
- `compute_wct` fails on the very first day with "No data for ..." – known behaviour of this version, not a problem.
- New data later: `msnoise scan_archive`, `msnoise new_jobs`, then the commands above.

---

## 4. Troubleshooting (all seen in practice)

| message | cause | fix |
|---|---|---|
| `No module named 'pkg_resources'` | recent setuptools removed it; an old package still imports it | `pip install "setuptools<81"` |
| `UserWarning: pkg_resources is deprecated` | harmless warning | ignore, or `set PYTHONWARNINGS=ignore::UserWarning` |
| `msnoise.__file__` points to e.g. `E:\Laure\msnoise\__init__.py` | Python found an `msnoise` folder in the current directory | run the check from another folder |
| `pip show msnoise` says `1.6.5` / traceback mentions `scripts\msnoise.py` line 9 `import pkg_resources` | the **official** msnoise is installed, not this one | section 1.3 (uninstall, reinstall) |
| `Cannot find command 'git'` | git not installed | `conda install -c conda-forge git` |
| `No module named 'xarray'` (or other) | dependencies not installed (`--no-deps`) | section 1.2 |
| `Table '<db>.config' doesn't exist` on **any** msnoise command | a `db.ini` in the current folder points to an empty database (every msnoise command reads it at start) | `del db.ini` then `msnoise db init`, or `cd` elsewhere |
| `db init`: "database seems to already exist and is not empty" | a previous attempt left half-created tables | drop and recreate the **new** database (below), then `msnoise db init` again |
| `invalid choice: 'E:\...'` from a script | the step name is missing | `python migrate_c_to_stable.py inspect --src ...` (step name first) |
| path ends up as `E:\Kilauea"` | a Windows path ending with `\` before a quote | write `E:\Kilauea`, not `"E:\Kilauea\"` |
| `migrate_c_to_stable.py`: "no 'category' column: looks like a CLASSIC db" | the project is classic, not msnoise_c | use `clone_classic.py` (section 3) |
| `Access denied` creating a database | user lacks rights | section 3.1 with root |

Reset the **new** database (never the original):
```
python -c "import pymysql; c=pymysql.connect(host='localhost', user='<USER>', password='<PASSWORD>'); cur=c.cursor(); cur.execute('DROP DATABASE <NEW_DB>'); cur.execute('CREATE DATABASE <NEW_DB>'); print('reset')"
```

Useful commands:
```
msnoise info -j                     # job counts per type and flag
msnoise reset STACK                 # I/E jobs -> T
msnoise reset STACK --all           # all jobs -> T
msnoise config get <name>
msnoise config set <name>=<value>
```

---

## 5. If the project is an msnoise_c project

Use `migrate_c_to_stable.py`. It converts the msnoise_c NetCDF CC files and database into this format. Do sections 1, 3.0, 3.1 and 3.2 first, then from the new folder:
```
python migrate_c_to_stable.py inspect --src <C_FOLDER> > inspect_output.txt
python migrate_c_to_stable.py db      --src <C_FOLDER> --dst . --cc-set 1
python migrate_c_to_stable.py files   --src <C_FOLDER> --dst . --cc-set 1 --dry-run
python migrate_c_to_stable.py files   --src <C_FOLDER> --dst . --cc-set 1 --threads 4
python migrate_c_to_stable.py verify  --src <C_FOLDER> --dst . --cc-set 1 --threads 4
```
- `inspect` lists config sets, filters and jobs. If there are several `cc` or `preprocess` sets, choose one with `--cc-set N` / `--preprocess-set N`.
- `files` can be re-run if interrupted; already converted files are skipped.
- `verify` must end with `RESULT : PASS`. The report is written to `migration_check.txt`.
- Then run section 3.6.

CCs computed by msnoise_c and by this version are very close but not identical (correlation ≈ 0.995). Stacking and dv/v on the migrated CCs are consistent. If new days are computed with this version, compare a few overlap days before trusting the join.
