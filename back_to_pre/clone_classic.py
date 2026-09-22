#!/usr/bin/env python
"""
clone_classic.py  -  copy a CLASSIC msnoise project (same DB format) into a new
project folder + new database, keeping the computed CCs, so the copy can be run
with another msnoise package without touching the original.

Copies: config, filters, stations, data_availability, CC jobs.
STACK jobs: one 'T' job per finished CC day/pair (so stacks are recomputed
by the new package). MWCS / DTT / STR / WCT / PSD jobs are NOT copied: they are
recreated by the new package when it stacks.
Adds any config key the new package needs that the old DB lacks (default value).

Usage (from anywhere, env of the NEW package active):
  python clone_classic.py --src E:\\Kilauea --dst E:\\Kilauea\\stable_dir_LB           (dry run)
  python clone_classic.py --src E:\\Kilauea --dst E:\\Kilauea\\stable_dir_LB --go

--dst must already contain a db.ini made by `msnoise db init` (fresh database).
The source database and folder are only READ.
"""
import argparse, csv, os, pickle, sys
from sqlalchemy import create_engine, MetaData, select, func
from sqlalchemy.pool import NullPool

EXTRA = {"coda_safety_factor": "1.2", "ref_type": "static"}


def engine_for(folder):
    with open(os.path.join(folder, "db.ini"), "rb") as f:
        v = list(pickle.load(f))
    if len(v) == 5:
        v.append("")
    tech, host, db, user, pw, prefix = v
    if tech == 1:
        p = host if os.path.isabs(host) else os.path.join(folder, host)
        url = "sqlite:///" + p
    elif tech == 2:
        url = "mysql+pymysql://%s:%s@%s/%s" % (user, pw, host, db)
    else:
        url = "postgresql+psycopg2://%s:%s@%s/%s" % (user, pw, host, db)
    e = create_engine(url, poolclass=NullPool)
    md = MetaData()
    md.reflect(bind=e)
    pre = (prefix + "_") if prefix else ""
    return e, md, pre, "%s (tech %s)" % (db or host, tech)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", required=True, help="original project folder (db.ini)")
    ap.add_argument("--dst", required=True, help="new project folder (db.ini from msnoise db init)")
    ap.add_argument("--go", action="store_true", help="really write (default: dry run)")
    a = ap.parse_args()
    if os.path.abspath(a.src) == os.path.abspath(a.dst):
        sys.exit("--src and --dst are the same folder")

    se, smd, spre, sname = engine_for(a.src)
    de, dmd, dpre, dname = engine_for(a.dst)
    if sname == dname and str(se.url) == str(de.url):
        sys.exit("source and destination point to the SAME database - stop.")
    T = lambda md, pre, n: md.tables.get(pre + n)
    for n in ("config", "filters", "jobs", "stations", "data_availability"):
        if T(smd, spre, n) is None:
            sys.exit("source has no '%s' table" % n)
        if T(dmd, dpre, n) is None:
            sys.exit("destination has no '%s' table - run `msnoise db init` in --dst" % n)
    if "category" in T(smd, spre, "config").c:
        sys.exit("source is an msnoise_c project: use migrate_c_to_stable.py instead")
    with de.connect() as c:
        if c.execute(select(func.count()).select_from(T(dmd, dpre, "jobs"))).scalar():
            sys.exit("destination jobs table is not empty: use a fresh database")

    print("source      :", sname, "  folder", os.path.abspath(a.src))
    print("destination :", dname, "  folder", os.path.abspath(a.dst))
    print("mode        :", "WRITE" if a.go else "DRY RUN (add --go to write)")

    def rows(n):
        with se.connect() as c:
            return [dict(r._mapping) for r in c.execute(select(T(smd, spre, n)))]

    plan = {}
    # config: source values override destination defaults, keep dst-only keys
    src_cfg = {r["name"]: r["value"] for r in rows("config")}
    with de.connect() as c:
        dst_cfg = {r.name: r.value for r in c.execute(select(T(dmd, dpre, "config")))}
    try:
        import msnoise
        dfl = os.path.join(os.path.dirname(msnoise.__file__), "default.csv")
        need = {r["name"]: r["default"] for r in csv.DictReader(open(dfl))}
    except Exception:
        need = {}
    need.update(EXTRA)
    new_cfg = dict(dst_cfg)
    new_cfg.update(src_cfg)
    added = {k: v for k, v in need.items() if k not in new_cfg}
    new_cfg.update(added)
    plan["config"] = [{"name": k, "value": v} for k, v in new_cfg.items()]

    for n in ("filters", "stations", "data_availability"):
        cols = set(T(dmd, dpre, n).c.keys())
        plan[n] = [{k: v for k, v in r.items() if k in cols} for r in rows(n)]

    jobs = rows("jobs")
    jcols = set(T(dmd, dpre, "jobs").c.keys()) - {"ref", "lastmod"}
    cc = [r for r in jobs if r["jobtype"] == "CC"]
    out = []
    for r in cc:
        j = {k: v for k, v in r.items() if k in jcols}
        if j["flag"] in ("I", "E"):
            j["flag"] = "T"            # interrupted -> redo
        out.append(j)
        if j["flag"] == "D":
            out.append({"day": r["day"], "pair": r["pair"], "jobtype": "STACK", "flag": "T"})
    plan["jobs"] = out

    print("\nconfig            : %d keys (%d taken from source, %d new keys added with default: %s)"
          % (len(new_cfg), len(src_cfg), len(added), ", ".join(sorted(added)) or "-"))
    for n in ("filters", "stations", "data_availability"):
        print("%-18s: %d rows" % (n, len(plan[n])))
    cnt = {}
    for j in out:
        cnt[(j["jobtype"], j["flag"])] = cnt.get((j["jobtype"], j["flag"]), 0) + 1
    print("jobs              :", ", ".join("%s %s=%d" % (k[0], k[1], v) for k, v in sorted(cnt.items())))
    for k in ("output_folder", "keep_all", "keep_days", "export_format", "data_folder"):
        print("  %-14s = %s" % (k, new_cfg.get(k)))

    if not a.go:
        print("\nDry run only. Re-run with --go to write.")
        return
    with de.begin() as c:
        for n in ("config", "filters", "stations", "data_availability", "jobs"):
            t = T(dmd, dpre, n)
            c.execute(t.delete())
            data = plan[n]
            for i in range(0, len(data), 5000):
                c.execute(t.insert(), data[i:i + 5000])
    print("\nDone. Database written.")


if __name__ == "__main__":
    main()
