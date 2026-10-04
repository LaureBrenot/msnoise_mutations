#!/usr/bin/env python
"""
check_ref.py - why does `msnoise cc stack -r` say "Empty dataset"?

Run it from the project folder (the one with db.ini):
    python check_ref.py

It prints the reference window, the filters, and for a few pairs the exact
file names msnoise looks for, whether they exist, and what IS on disk.
Read-only.
"""
import datetime
import glob
import os
import pickle
import sys


def connect(folder="."):
    with open(os.path.join(folder, "db.ini"), "rb") as f:
        v = list(pickle.load(f))
    if len(v) == 5:
        v.append("")
    tech, host, db, user, pw, prefix = v
    if tech == 1:
        import sqlite3
        p = host if os.path.isabs(host) else os.path.join(folder, host)
        return sqlite3.connect(p), prefix
    import pymysql
    h, port = (host.split(":") + ["3306"])[:2] if ":" in str(host) else (host, 3306)
    return pymysql.connect(host=h, port=int(port), user=user, password=pw,
                           database=db), prefix


con, prefix = connect()
p = (prefix + "_") if prefix else ""
cur = con.cursor()
cur.execute("SELECT name, value FROM %sconfig" % p)
cfg = dict(cur.fetchall())
print("project folder :", os.path.abspath("."))
print("ref_begin      :", cfg.get("ref_begin"))
print("ref_end        :", cfg.get("ref_end"))
print("startdate      :", cfg.get("startdate"))
print("enddate        :", cfg.get("enddate"))
print("export_format  :", cfg.get("export_format"))
print("keep_all       :", cfg.get("keep_all"))
print("components     :", cfg.get("components_to_compute"),
      "| single station:", cfg.get("components_to_compute_single_station"))

cur.execute("SELECT ref, low, high, used FROM %sfilters" % p)
filters = cur.fetchall()
print("filters        :", filters)

# ---- reference window, same logic as msnoise ------------------------------
begin, end = cfg.get("ref_begin", ""), cfg.get("ref_end", "")
try:
    if begin.startswith("-"):
        start = datetime.date.today() + datetime.timedelta(days=int(begin))
        stop = datetime.date.today() + datetime.timedelta(days=int(end))
    else:
        start = datetime.datetime.strptime(begin, "%Y-%m-%d").date()
        stop = datetime.datetime.strptime(end, "%Y-%m-%d").date()
    stop = min(stop, datetime.date.today())
    print("REF window     : %s .. %s  (%d days)" % (start, stop, (stop - start).days + 1))
except Exception as e:
    sys.exit("cannot parse ref_begin/ref_end: %s" % e)

# ---- what is on disk ------------------------------------------------------
print("\nSTACKS folder:")
if not os.path.isdir("STACKS"):
    sys.exit("  !! no STACKS folder (or junction) in this project folder")
for ff in sorted(os.listdir("STACKS")):
    sub = os.path.join("STACKS", ff)
    if os.path.isdir(sub):
        print("  STACKS\\%s -> %s" % (ff, sorted(os.listdir(sub))[:5]))

# ---- what msnoise will look for ------------------------------------------
cur.execute("SELECT DISTINCT pair FROM %sjobs WHERE jobtype='CC' AND flag='D' "
            "ORDER BY pair LIMIT 3" % p)
pairs = [r[0] for r in cur.fetchall()]
print("\nchecking a few pairs:")
for fref, low, high, used in filters:
    if not used:
        continue
    for pair in pairs:
        s1, s2 = pair.split(":")
        base = os.path.join("STACKS", "%02i" % int(fref), "001_DAYS")
        for comp in ("ZZ",):
            d = os.path.join(base, comp, "%s_%s" % (s1, s2))
            exists = os.path.isdir(d)
            print("\n  %s  filter %02i  %s" % (pair, int(fref), comp))
            print("    folder          : %s  %s" % (d, "OK" if exists else "MISSING"))
            if not exists:
                # is the folder there under another name?
                alt = glob.glob(os.path.join(base, comp, "*%s*" % s1.split(".")[1]))
                print("    similar folders : %s" % (alt[:3] or "none"))
                continue
            files = sorted(os.listdir(d))
            print("    files on disk   : %d  (first %s, last %s)"
                  % (len(files), files[0], files[-1]))
            want = os.path.join(d, "%s.MSEED" % start)
            print("    wants (1st day) : %s  %s"
                  % (want, "OK" if os.path.isfile(want) else "MISSING"))
            inwin = [f for f in files
                     if str(start) <= f[:10] <= str(stop)]
            print("    files inside REF window: %d" % len(inwin))
            if inwin:
                print("      e.g. %s" % inwin[:3])
