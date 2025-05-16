# Fork database and clone msnoise env
**Why**: try a new WCT implementation without losing your current progress or breaking compatibility with your workflow

```./mariadb/bin/mysql```

## Step 1: Create the new database for your fork
```CREATE DATABASE IF NOT EXISTS static_new_wct;```

## Step 2: Copy the structure and data of each table from "static" to "static_new_wct"

-- Copy config table
```
CREATE TABLE static_new_wct.config LIKE static.config;
INSERT INTO static_new_wct.config SELECT * FROM static.config;
```
-- Copy data_availability table
```
CREATE TABLE static_new_wct.data_availability LIKE static.data_availability;
INSERT INTO static_new_wct.data_availability SELECT * FROM static.data_availability;
```
-- Copy filters table
```
CREATE TABLE static_new_wct.filters LIKE static.filters;
INSERT INTO static_new_wct.filters SELECT * FROM static.filters;
```
-- Copy jobs table
```
CREATE TABLE static_new_wct.jobs LIKE static.jobs;
INSERT INTO static_new_wct.jobs SELECT * FROM static.jobs;
```
-- Copy stations table
```
CREATE TABLE static_new_wct.stations LIKE static.stations;
INSERT INTO static_new_wct.stations SELECT * FROM static.stations;
```
## Step 3: Verify that all data was copied correctly
```
SELECT 'Config table count:', COUNT(*) FROM static.config
UNION ALL
SELECT 'Config table count (new):', COUNT(*) FROM static_new_wct.config
UNION ALL
SELECT 'Data availability count:', COUNT(*) FROM static.data_availability
UNION ALL
SELECT 'Data availability count (new):', COUNT(*) FROM static_new_wct.data_availability
UNION ALL
SELECT 'Filters count:', COUNT(*) FROM static.filters
UNION ALL
SELECT 'Filters count (new):', COUNT(*) FROM static_new_wct.filters
UNION ALL
SELECT 'Jobs count:', COUNT(*) FROM static.jobs
UNION ALL
SELECT 'Jobs count (new):', COUNT(*) FROM static_new_wct.jobs
UNION ALL
SELECT 'Stations count:', COUNT(*) FROM static.stations
UNION ALL
SELECT 'Stations count (new):', COUNT(*) FROM static_new_wct.stations;
```
## Step 4: Create new db.ini file
```mkdir msnoise_new_wct```
```cd msnoise_new_wct```
```vi copy_db.py```

```
import pickle
import sys

# Try to infer the format from your existing db.ini file
try:
    with open('/globalscratch/ulb/gtime/lbrenot/msnoise__static/db.ini', 'rb') as f:
        original_config = pickle.load(f)

    # Print the original format
    print(f"Original config has {len(original_config)} values: {original_config}")

    # Create a new config with the same format but different database
    if len(original_config) == 5:
        new_config = (2, '10.44.2.1:5050', 'static_new_wct', 'root', 'noise')
    else:
        new_config = (2, '10.44.2.1:5050', 'static_new_wct', 'root', 'noise', '')

    # Write the new config
    with open('db.ini', 'wb') as f:
        pickle.dump(new_config, f)

    print(f"Created new db.ini with format matching the original")
except Exception as e:
    print(f"Error: {e}")
    sys.exit(1)
```

check with:
```msnoise info -j```

## Step 5: Clone environment

```conda create -n msnoise_copy_clean --clone msnoise_clean```

