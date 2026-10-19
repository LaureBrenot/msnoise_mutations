#!/usr/bin/env python3
"""
Recover MSNoise settings from an existing STACKS2 folder and create the
MWCS jobs needed to resume processing from there (MWCS -> DTT -> DVV).

Run it from the MSNoise project folder (the one containing db.ini and STACKS2):

    python recover_from_stacks2.py              # inspect only, prints what it found
    python recover_from_stacks2.py --create-jobs  # also inserts MWCS jobs (flag T)
"""
import argparse
import datetime
import glob
import os
from collections import defaultdict

import numpy as np
import pandas as pd
import xarray as xr


def scan(root):
    filters, movstacks, comps, stations = set(), set(), set(), set()
    days = defaultdict(set)
    taxis_info = None
    nfiles = 0

    # STACKS2/<filterid>/<movstack>/<components>/<sta1>_<sta2>.nc
    for fn in sorted(glob.glob(os.path.join(root, "*", "*", "*", "*.nc"))):
        fid, ms, comp, name = os.path.relpath(fn, root).split(os.sep)
        if ms == "REF":
            continue
        s1, s2 = name[:-3].split("_")
        filters.add(int(fid))
        movstacks.add(ms)
        comps.add(comp)
        stations.update([s1, s2])
        pair = "%s:%s" % (s1, s2)
        with xr.open_dataset(fn) as ds:
            t = pd.DatetimeIndex(ds.times.values)
            days[pair].update(t.normalize().strftime("%Y-%m-%d"))
            if taxis_info is None:
                tax = ds.taxis.values
                dt = float(np.median(np.diff(tax)))
                taxis_info = (float(tax.max()), 1.0 / dt, len(tax))
        nfiles += 1

    return dict(filters=sorted(filters), movstacks=sorted(movstacks),
                comps=sorted(comps), stations=sorted(stations),
                days=days, taxis=taxis_info, nfiles=nfiles)


def report(info):
    if not info["nfiles"]:
        raise SystemExit("No STACKS2 files found. Run this from the project folder.")
    alldays = sorted(set().union(*info["days"].values()))
    maxlag, sr, npts = info["taxis"]
    print("Files scanned     : %i" % info["nfiles"])
    print("Filters (refs)    : %s" % info["filters"])
    print("mov_stacks        : %s" % info["movstacks"])
    print("Components        : %s" % ",".join(info["comps"]))
    print("Pairs             : %i" % len(info["days"]))
    print("Stations          : %s" % ", ".join(info["stations"]))
    print("Days              : %s -> %s (%i unique)" % (alldays[0], alldays[-1], len(alldays)))
    print("taxis             : maxlag=%g s, cc_sampling_rate=%g Hz, %i samples" % (maxlag, sr, npts))

    ms = ",".join("('%s','%s')" % tuple(m.split("_")) for m in info["movstacks"])
    cross = [c for c in info["comps"]]
    print("\nSuggested config (check against what you originally used):")
    print("  msnoise config set maxlag=%g" % maxlag)
    print("  msnoise config set cc_sampling_rate=%g" % sr)
    print("  msnoise config set \"mov_stack=(%s)\"" % ms)
    print("  msnoise config set components_to_compute=%s" % ",".join(cross))
    print("  # if some pairs are autocorrelations (sta1 == sta2), set")
    print("  # components_to_compute_single_station accordingly instead")
    print("\nFilters must exist in the DB with refs %s (create them in order)." % info["filters"])


def create_jobs(info, chunk=10000):
    from msnoise.api import connect, massive_insert_job
    from msnoise.msnoise_table_def import Job

    db = connect()
    existing = set(db.query(Job.day, Job.pair).filter(Job.jobtype == "MWCS").all())
    db.close()

    now = datetime.datetime.utcnow()
    jobs = [dict(day=d, pair=p, jobtype="MWCS", flag="T", lastmod=now)
            for p, ds in info["days"].items() for d in sorted(ds)
            if (d, p) not in existing]
    for i in range(0, len(jobs), chunk):
        massive_insert_job(jobs[i:i + chunk])
    print("\nInserted %i MWCS jobs (%i already existed)." % (len(jobs), len(existing)))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default="STACKS2")
    ap.add_argument("--create-jobs", action="store_true")
    args = ap.parse_args()
    info = scan(args.root)
    report(info)
    if args.create_jobs:
        create_jobs(info)
