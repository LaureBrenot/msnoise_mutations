#!/usr/bin/env python
"""
migrate_c_to_stable.py
======================

Move an MSNoise project from the *workflow / lineage* version ("msnoise_c",
config sets, OUTPUT/<preprocess_x>/<cc_x>/<filter_x>/_output/...) back to the
*classic* version ("msnoise", single config table, filters table,
STACKS/<ff>/001_DAYS/... and <output_folder>/<ff>/<sta1>/<sta2>/<comp>/<day>.h5),
WITHOUT recomputing the cross-correlations.

The source project (msnoise_c) is only READ, never modified.

Three sub-commands, run in this order:

  inspect  : show what is in the source project (config sets, lineages,
             CC output folders, job counts). Read-only.
  db       : fill a *freshly initialised* classic database
             (created with `msnoise db init` from the classic environment)
             with config, filters, stations, data_availability and jobs.
  verify   : check the result: config, filters, stations, data availability,
             jobs, and that every converted CC file holds exactly the same
             numbers as the msnoise_c file. Writes migration_check.txt.
  files    : convert the CC NetCDF files to the classic on-disk layout
             (daily stacks -> STACKS/<ff>/001_DAYS/<comp>/<s1>_<s2>/<day>.MSEED,
              keep_all windows -> <output_folder>/<ff>/<s1>/<s2>/<comp>/<day>.h5).
             Resumable: existing target files are skipped.

Run it with the python of the CLASSIC msnoise environment (needs sqlalchemy,
pandas, tables, xarray, netCDF4, obspy; pymysql for MySQL/MariaDB).

Example
-------
  python migrate_c_to_stable.py inspect --src /path/msnoise_c_project
  python migrate_c_to_stable.py db      --src /path/msnoise_c_project --dst /path/stable_project --cc-set 1 \
                                        --classic-ref /path/old_classic_project
  python migrate_c_to_stable.py files   --src /path/msnoise_c_project --dst /path/stable_project --cc-set 1 --threads 8
  python migrate_c_to_stable.py verify  --src /path/msnoise_c_project --dst /path/stable_project --cc-set 1 --threads 8
"""
import argparse
import collections
import glob
import os
import pickle
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
from sqlalchemy import MetaData, create_engine, select
from sqlalchemy.pool import NullPool

# --------------------------------------------------------------------------
# DB helpers
# --------------------------------------------------------------------------
IniFile = collections.namedtuple(
    "IniFile", ["tech", "hostname", "database", "username", "password", "prefix"])


def read_ini(project_dir):
    fn = os.path.join(project_dir, "db.ini")
    with open(fn, "rb") as f:
        vals = list(pickle.load(f))
    if len(vals) == 5:
        vals.append("")
    return IniFile(*vals)


def make_engine(project_dir):
    ini = read_ini(project_dir)
    if ini.tech == 1:
        path = ini.hostname
        if not os.path.isabs(path):
            path = os.path.join(project_dir, path)
        return create_engine("sqlite:///%s" % path), ini
    if ini.tech == 2:
        url = "mysql+pymysql://%s:%s@%s/%s" % (
            ini.username, ini.password, ini.hostname, ini.database)
    elif ini.tech == 3:
        url = "postgresql+psycopg2://%s:%s@%s/%s" % (
            ini.username, ini.password, ini.hostname, ini.database)
    else:
        raise ValueError("Unknown tech %s in db.ini" % ini.tech)
    return create_engine(url, poolclass=NullPool), ini


def reflect(engine, ini):
    md = MetaData()
    md.reflect(bind=engine)
    p = ini.prefix or ""

    def t(name):
        # new version: most tables are "<prefix>_<name>", workflow tables
        # are "<prefix><name>" (no underscore). Try both.
        for cand in (p + "_" + name if p else name, p + name, name):
            if cand in md.tables:
                return md.tables[cand]
        return None
    return t


# --------------------------------------------------------------------------
# Source (msnoise_c) readers
# --------------------------------------------------------------------------
class Source:
    def __init__(self, project_dir):
        self.dir = os.path.abspath(project_dir)
        self.engine, self.ini = make_engine(self.dir)
        self.t = reflect(self.engine, self.ini)
        for needed in ("config", "jobs", "stations", "data_availability"):
            if self.t(needed) is None:
                sys.exit("Source DB has no '%s' table - is --src really the "
                         "msnoise_c project folder?" % needed)
        if "category" not in self.t("config").c:
            sys.exit("Source 'config' table has no 'category' column: this "
                     "looks like a CLASSIC db, not an msnoise_c one.")

    def rows(self, table, *where):
        tab = self.t(table)
        q = select(tab)
        for w in where:
            q = q.where(w(tab))
        with self.engine.connect() as c:
            return [dict(r._mapping) for r in c.execute(q)]

    def config(self):
        """{(category, set_number): {name: value}}"""
        # global params are stored with set_number NULL or 1 depending on
        # the version -> always key them as ("global", None)
        out = collections.defaultdict(dict)
        for r in self.rows("config"):
            s = None if r["category"] == "global" else r["set_number"]
            out[(r["category"], s)][r["name"]] = r["value"]
        return out

    def steps(self):
        tab = self.t("workflow_steps")
        if tab is None:
            return {}
        with self.engine.connect() as c:
            return {r.step_id: dict(r._mapping) for r in c.execute(select(tab))}

    def lineages(self):
        tab = self.t("lineages")
        if tab is None:
            return {}
        with self.engine.connect() as c:
            return {r.lineage_id: r.lineage_str for r in c.execute(select(tab))}

    def data_sources(self):
        tab = self.t("data_sources")
        if tab is None:
            return {}
        with self.engine.connect() as c:
            return {r.ref: dict(r._mapping) for r in c.execute(select(tab))}

    def output_folder(self):
        cfg = self.config()
        of = cfg.get(("global", None), {}).get("output_folder", "OUTPUT")
        return of if os.path.isabs(of) else os.path.join(self.dir, of)


def find_cc_roots(src, cc_set, preprocess_set=None):
    """Return list of <.../cc_N> folders on disk holding filter_* outputs."""
    root = src.output_folder()
    hits = []
    for dirpath, dirnames, _ in os.walk(root):
        if os.path.basename(dirpath) == "cc_%d" % cc_set:
            if any(d.startswith("filter_") for d in dirnames):
                hits.append(dirpath)
            dirnames[:] = []  # don't go deeper
        elif "_output" in os.path.basename(dirpath):
            dirnames[:] = []
    if preprocess_set is not None:
        hits = [h for h in hits
                if "preprocess_%d" % preprocess_set in h.split(os.sep)]
    return sorted(hits)


# --------------------------------------------------------------------------
# inspect
# --------------------------------------------------------------------------
def cmd_inspect(a):
    src = Source(a.src)
    print("Source project :", src.dir)
    print("DB             : tech=%s host=%s db=%s prefix=%r" % (
        src.ini.tech, src.ini.hostname, src.ini.database, src.ini.prefix))
    cfg = src.config()
    print("\nConfig sets found:")
    for (cat, s) in sorted(cfg, key=lambda k: (k[0], k[1] or 0)):
        extra = ""
        if cat == "filter":
            c = cfg[(cat, s)]
            extra = "  freqmin=%s freqmax=%s used=%s" % (
                c.get("freqmin"), c.get("freqmax"), c.get("used"))
        if cat == "cc":
            c = cfg[(cat, s)]
            extra = "  maxlag=%s cc_sampling_rate=%s corr_duration=%s keep_all=%s keep_days=%s comps=%s" % (
                c.get("maxlag"), c.get("cc_sampling_rate"), c.get("corr_duration"),
                c.get("keep_all"), c.get("keep_days"), c.get("components_to_compute"))
        print("  %-18s set=%-4s %3d params%s" % (cat, s, len(cfg[(cat, s)]), extra))

    lin = src.lineages()
    if lin:
        print("\nLineages:")
        for k, v in sorted(lin.items()):
            print("  %4d  %s" % (k, v))

    print("\nJobs per step / flag:")
    counts = collections.Counter((r["jobtype"], r["flag"]) for r in src.rows("jobs"))
    for (jt, fl), n in sorted(counts.items()):
        print("  %-20s %s %8d" % (jt, fl, n))

    print("\nData sources:")
    for r in src.data_sources().values():
        print("  ref=%s name=%s uri=%s structure=%s" % (
            r["ref"], r["name"], r["uri"], r["data_structure"]))

    print("\nCC output folders on disk (under %s):" % src.output_folder())
    for cc in sorted({k[1] for k in cfg if k[0] == "cc"}):
        for root in find_cc_roots(src, cc):
            for fdir in sorted(glob.glob(os.path.join(root, "filter_*"))):
                for kind in ("daily", "all"):
                    p = os.path.join(fdir, "_output", kind)
                    if os.path.isdir(p):
                        n = sum(len(f) for _, _, f in os.walk(p))
                        print("  %s  -> %d files" % (os.path.relpath(p, src.dir), n))


# --------------------------------------------------------------------------
# db
# --------------------------------------------------------------------------
# classic name  <-  (msnoise_c category, msnoise_c name)   explicit renames
RENAMES = {
    "dtt_maxdt": ("mwcs_dtt", "dtt_maxdtt"),
    "dtt_codacycles": ("wavelet_dtt", "wct_codacycles"),
    "dvv_min_nonzero": ("wavelet_dtt", "wct_min_nonzero"),
}
# Keys read by the customised classic s08compute_wct.py but absent from the
# classic default.csv
EXTRA_DEFAULTS = {"coda_safety_factor": "1.2", "ref_type": "static"}

# Where to look for a classic parameter, in priority order
CATEGORY_PRIORITY = ["global", "preprocess", "cc", "stack", "refstack",
                     "mwcs", "mwcs_dtt", "mwcs_dtt_dvv", "wavelet",
                     "wavelet_dtt", "wavelet_dtt_dvv", "stretching",
                     "stretching_dvv", "psd", "psd_rms"]


def pick_sets(cfg, a):
    """Choose one set number per category (the ones you want to carry)."""
    chosen = {}
    for cat in CATEGORY_PRIORITY[1:]:
        nums = sorted(s for (c, s) in cfg if c == cat and s is not None)
        if not nums:
            continue
        want = getattr(a, cat.replace("-", "_") + "_set", None)
        chosen[cat] = want if want in nums else nums[0]
    return chosen


def cmd_db(a):
    src = Source(a.src)
    dst_engine, dst_ini = make_engine(a.dst)
    dt = reflect(dst_engine, dst_ini)
    for needed in ("config", "filters", "jobs", "stations", "data_availability"):
        if dt(needed) is None:
            sys.exit("Destination DB has no '%s' table. Run `msnoise db init` "
                     "with the CLASSIC msnoise in %s first." % (needed, a.dst))
    if "category" in dt("config").c:
        sys.exit("Destination looks like an msnoise_c DB, not a classic one.")

    with dst_engine.connect() as c:
        n = c.execute(select(dt("jobs"))).first()
    if n is not None and not a.force:
        sys.exit("Destination 'jobs' table is not empty. Use a fresh database "
                 "(or --force to append).")

    cfg = src.config()
    sets = pick_sets(cfg, a)
    print("Using config sets:", sets)

    # ---------------- config ----------------
    ds_list = list(src.data_sources().values())
    ds = ds_list[0] if ds_list else {}
    if len(ds_list) > 1:
        print("WARNING: %d data sources in source DB; classic msnoise supports "
              "one data_folder. Using '%s'." % (len(ds_list), ds.get("name")))

    def lookup(name):
        if name in RENAMES:
            cat, n = RENAMES[name]
            if cat in sets and n in cfg.get((cat, sets[cat]), {}):
                return cfg[(cat, sets[cat])][n]
        if name in cfg.get(("global", None), {}):
            return cfg[("global", None)][name]
        for cat in CATEGORY_PRIORITY[1:]:
            if cat in sets and name in cfg.get((cat, sets[cat]), {}):
                return cfg[(cat, sets[cat])][name]
        return None

    ds_map = {}
    if ds:
        uri = ds.get("uri") or ""
        if uri.startswith("sds://"):
            uri = uri[len("sds://"):]
        ds_map = {"data_folder": uri,
                  "data_structure": ds.get("data_structure"),
                  "archive_format": ds.get("archive_format") or "",
                  "network": ds.get("network_code"),
                  "channels": ds.get("channels")}

    ctab = dt("config")
    updated, missing = [], []
    with dst_engine.begin() as c:
        classic_names = [r.name for r in c.execute(select(ctab.c.name))]
        for name in classic_names:
            if name in ds_map and ds_map[name] is not None:
                val = ds_map[name]
            elif name == "output_folder":
                val = a.output_folder
            elif name == "export_format":
                val = "MSEED"
            else:
                val = lookup(name)
            if val is None:
                missing.append(name)
                continue
            if a.keep_days_only and name == "keep_all":
                val = "N"
            c.execute(ctab.update().where(ctab.c.name == name).values(value=str(val)))
            updated.append(name)
    print("config: %d params copied, %d left at classic defaults: %s" % (
        len(updated), len(missing), ", ".join(missing)))

    # Extra keys used by your customised classic scripts (e.g. s08compute_wct)
    # that `msnoise db init` does not create.
    extra = {}
    if a.classic_ref:
        ref_engine, ref_ini = make_engine(a.classic_ref)
        rt = reflect(ref_engine, ref_ini)
        with ref_engine.connect() as c:
            ref_cfg = {r.name: r.value for r in c.execute(select(rt("config")))}
        extra.update({k: v for k, v in ref_cfg.items() if k not in classic_names})
        print("classic reference DB: %d extra keys found: %s" % (
            len(extra), ", ".join(sorted(extra))))
    for k, v in EXTRA_DEFAULTS.items():
        if k not in extra and k not in classic_names:
            extra[k] = v
            print("WARNING: '%s' not found anywhere, inserted default %r - "
                  "check it with `msnoise config get %s`" % (k, v, k))
    with dst_engine.begin() as c:
        for k, v in extra.items():
            c.execute(ctab.insert().values(name=k, value=v))

    # ---------------- filters ----------------
    lin = src.lineages()
    ftab = dt("filters")
    filters = []
    for (cat, s), c in sorted(cfg.items(), key=lambda kv: (kv[0][0], kv[0][1] or 0)):
        if cat != "filter" or s is None:
            continue
        # mwcs set downstream of this filter (via lineage strings)
        mw = None
        for l in lin.values():
            parts = l.split("/")
            if "filter_%d" % s in parts:
                m = [p for p in parts if p.startswith("mwcs_") and p[5:].isdigit()]
                if m:
                    mw = int(m[0][5:])
                    break
        if mw is None:
            mw = sets.get("mwcs")
        mc = cfg.get(("mwcs", mw), {})
        low, high = float(c.get("freqmin", 0.1)), float(c.get("freqmax", 1.0))
        used = str(c.get("used", "Y")).upper() in ("Y", "TRUE", "1")
        filters.append(dict(
            ref=s, low=low, high=high,
            mwcs_low=float(mc.get("freqmin", low)),
            mwcs_high=float(mc.get("freqmax", high)),
            mwcs_wlen=float(mc.get("mwcs_wlen", 10.0)),
            mwcs_step=float(mc.get("mwcs_step", 5.0)),
            used=used))
    with dst_engine.begin() as c:
        if a.force:
            c.execute(ftab.delete())
        if filters:
            c.execute(ftab.insert(), filters)
    for f in filters:
        print("filter %02d: %.3f-%.3f Hz  (mwcs %.3f-%.3f, wlen %.1f, step %.1f) used=%s"
              % (f["ref"], f["low"], f["high"], f["mwcs_low"], f["mwcs_high"],
                 f["mwcs_wlen"], f["mwcs_step"], f["used"]))

    # ---------------- stations ----------------
    stab = dt("stations")
    cols = [c for c in ("ref", "net", "sta", "used_location_codes",
                        "used_channel_names", "X", "Y", "altitude",
                        "coordinates", "used") if c in stab.c]
    stations = [{k: r.get(k) for k in cols} for r in src.rows("stations")]
    with dst_engine.begin() as c:
        if a.force:
            c.execute(stab.delete())
        if stations:
            c.execute(stab.insert(), stations)
    print("stations: %d copied" % len(stations))

    # ---------------- data_availability ----------------
    datab = dt("data_availability")
    dcols = [c for c in ("net", "sta", "loc", "chan", "path", "file",
                         "starttime", "endtime", "data_duration",
                         "gaps_duration", "samplerate", "flag") if c in datab.c]
    roots = {k: (v.get("uri") or "").replace("sds://", "")
             for k, v in src.data_sources().items()}
    das = []
    for r in src.rows("data_availability"):
        d = {k: r.get(k) for k in dcols}
        root = roots.get(r.get("data_source_id"), ds_map.get("data_folder", ""))
        if d["path"] and not os.path.isabs(d["path"]) and root:
            d["path"] = os.path.join(root, d["path"])
        # 'A' = already processed: classic `new_jobs` will not recreate CC jobs
        # for these files. Files added later get N/M and will be processed.
        d["flag"] = "A"
        das.append(d)
    with dst_engine.begin() as c:
        if a.force:
            c.execute(datab.delete())
        for i in range(0, len(das), 5000):
            c.execute(datab.insert(), das[i:i + 5000])
    print("data_availability: %d rows copied (flag set to 'A')" % len(das))

    # ---------------- jobs ----------------
    steps = src.steps()
    cc_step_ids = {sid for sid, s in steps.items()
                   if s["category"] == "cc" and s["set_number"] == sets.get("cc")}
    pre = sets.get("preprocess")
    jobs = {}
    for r in src.rows("jobs"):
        if r.get("step_id") not in cc_step_ids:
            continue
        if pre is not None and lin and r.get("lineage_id") in lin:
            if "preprocess_%d" % pre not in lin[r["lineage_id"]].split("/"):
                continue
        key = (r["day"], r["pair"])
        # if duplicates (several lineages), keep "least done" flag
        order = {"T": 0, "I": 1, "F": 2, "D": 3}
        if key not in jobs or order.get(r["flag"], 0) < order.get(jobs[key], 0):
            jobs[key] = r["flag"]
    jtab = dt("jobs")
    out = []
    for (day, pair), flag in jobs.items():
        if flag == "I":   # was running when you stopped -> redo it
            flag = "T"
        out.append(dict(day=day, pair=pair, jobtype="CC", flag=flag))
        if flag == "D":
            out.append(dict(day=day, pair=pair, jobtype="STACK", flag="T"))
    with dst_engine.begin() as c:
        if a.force:
            c.execute(jtab.delete())
        for i in range(0, len(out), 5000):
            c.execute(jtab.insert(), out[i:i + 5000])
    fc = collections.Counter((j["jobtype"], j["flag"]) for j in out)
    print("jobs:", dict(fc))
    print("\nDone. Now run the 'files' step.")


# --------------------------------------------------------------------------
# files
# --------------------------------------------------------------------------
def _convert_daily(src_fn, dst_fn, sr, maxlag, s1, s2, ncorr):
    import xarray as xr
    from obspy import Stream, Trace
    ds = xr.open_dataset(src_fn)
    data = ds["CCF"].values.astype(np.float64)
    ds.close()
    tr = Trace(data=data)
    pair = "%s:%s" % (s1, s2)
    tr.stats["station"] = pair[:11]          # same as classic export_mseed
    tr.stats["sampling_rate"] = sr
    tr.stats["location"] = "%02i" % min(int(ncorr), 99)
    os.makedirs(os.path.dirname(dst_fn), exist_ok=True)
    tmp = dst_fn + ".tmp"
    Stream([tr]).write(tmp, format="MSEED")
    os.replace(tmp, dst_fn)
    return len(data)


def _convert_all(src_fn, dst_fn, taxis):
    import pandas as pd
    import xarray as xr
    ds = xr.open_dataset(src_fn)
    da = ds["CCF"].load()
    ds.close()
    if da.shape[1] != len(taxis):
        raise ValueError("taxis length %d != classic %d" % (da.shape[1], len(taxis)))
    idx = pd.to_datetime(da.coords["times"].values).strftime("%Y-%m-%d %H:%M:%S")
    df = pd.DataFrame(da.values, index=idx, columns=taxis)
    os.makedirs(os.path.dirname(dst_fn), exist_ok=True)
    tmp = dst_fn + ".tmp"
    df.to_hdf(tmp, key="data")
    os.replace(tmp, dst_fn)
    return da.shape[0]


def _task(kind, src_fn, dst_fn, sr, maxlag, s1, s2, taxis):
    try:
        if kind == "daily":
            ncorr = 0
            all_fn = src_fn.replace(os.sep + "daily" + os.sep, os.sep + "all" + os.sep)
            if os.path.isfile(all_fn):
                import xarray as xr
                with xr.open_dataset(all_fn) as d:
                    ncorr = d.sizes.get("times", 0)
            _convert_daily(src_fn, dst_fn, sr, maxlag, s1, s2, ncorr)
        else:
            _convert_all(src_fn, dst_fn, taxis)
        return None
    except Exception as e:  # report, keep going
        return "%s: %s" % (src_fn, e)


def cmd_files(a):
    src = Source(a.src)
    dst_engine, dst_ini = make_engine(a.dst)
    dt = reflect(dst_engine, dst_ini)
    with dst_engine.connect() as c:
        dcfg = {r.name: r.value for r in c.execute(select(dt("config")))}
    maxlag = float(dcfg["maxlag"])
    sr = float(dcfg["cc_sampling_rate"])
    taxis = np.linspace(-maxlag, maxlag, int(2 * maxlag * sr) + 1)
    out_folder = dcfg.get("output_folder", "CROSS_CORRELATIONS")
    if not os.path.isabs(out_folder):
        out_folder = os.path.join(os.path.abspath(a.dst), out_folder)
    stacks = os.path.join(os.path.abspath(a.dst), "STACKS")

    roots = find_cc_roots(src, a.cc_set, a.preprocess_set)
    if not roots:
        sys.exit("No cc_%d folder with filter_* outputs found under %s"
                 % (a.cc_set, src.output_folder()))
    if len(roots) > 1:
        sys.exit("Several cc_%d folders found, choose one with --preprocess-set:\n  %s"
                 % (a.cc_set, "\n  ".join(roots)))
    root = roots[0]
    print("Reading CCs from:", root)
    print("Writing daily  -> %s" % stacks)
    if not a.keep_days_only:
        print("Writing keep_all -> %s" % out_folder)

    tasks, skipped, checked = [], 0, False
    for fdir in sorted(glob.glob(os.path.join(root, "filter_*"))):
        fid = int(os.path.basename(fdir).split("_")[1])
        kinds = ["daily"] if a.keep_days_only else ["daily", "all"]
        for kind in kinds:
            base = os.path.join(fdir, "_output", kind)
            if not os.path.isdir(base):
                continue
            for comp in sorted(os.listdir(base)):
                for pairdir in sorted(glob.glob(os.path.join(base, comp, "*"))):
                    pair = os.path.basename(pairdir)
                    parts = pair.split("_")
                    if len(parts) != 2:
                        print("  ! cannot split pair folder %r, skipped" % pair)
                        continue
                    s1, s2 = parts
                    for fn in sorted(glob.glob(os.path.join(pairdir, "*.nc"))):
                        day = os.path.basename(fn)[:-3]
                        if kind == "daily":
                            dst = os.path.join(stacks, "%02i" % fid, "001_DAYS",
                                               comp, pair, day + ".MSEED")
                        else:
                            dst = os.path.join(out_folder, "%02i" % fid, s1, s2,
                                               comp, day + ".h5")
                        if not checked:
                            import xarray as xr
                            with xr.open_dataset(fn) as d:
                                n = d.sizes["taxis"]
                            if n != len(taxis):
                                sys.exit("Lag axis mismatch: files have %d samples, "
                                         "classic config (maxlag=%s, cc_sampling_rate=%s) "
                                         "expects %d. Fix maxlag / cc_sampling_rate in "
                                         "the classic DB first." % (n, maxlag, sr, len(taxis)))
                            checked = True
                        if os.path.isfile(dst) and not a.overwrite:
                            skipped += 1
                            continue
                        tasks.append((kind, fn, dst, sr, maxlag, s1, s2, taxis))
    print("%d files to convert, %d already done (skipped)" % (len(tasks), skipped))
    if a.dry_run:
        for t in tasks[:10]:
            print("  %s\n   -> %s" % (t[1], t[2]))
        return

    errors = []
    if a.threads > 1:
        with ProcessPoolExecutor(a.threads) as ex:
            futs = [ex.submit(_task, *t) for t in tasks]
            for i, f in enumerate(as_completed(futs), 1):
                r = f.result()
                if r:
                    errors.append(r)
                if i % 1000 == 0:
                    print("  %d / %d" % (i, len(tasks)), flush=True)
    else:
        for i, t in enumerate(tasks, 1):
            r = _task(*t)
            if r:
                errors.append(r)
            if i % 1000 == 0:
                print("  %d / %d" % (i, len(tasks)), flush=True)
    print("Converted %d files, %d errors" % (len(tasks) - len(errors), len(errors)))
    for e in errors[:20]:
        print("  ERROR", e)


# --------------------------------------------------------------------------
# verify
# --------------------------------------------------------------------------
def _num_equal(a, b):
    try:
        return abs(float(a) - float(b)) < 1e-9
    except (TypeError, ValueError):
        return str(a).strip().upper() == str(b).strip().upper()


def _check_file(kind, src_fn, dst_fn):
    """Return None if dst holds exactly the same numbers as src, else a message."""
    import pandas as pd
    import xarray as xr
    try:
        with xr.open_dataset(src_fn) as d:
            ref = d["CCF"].load()
        if kind == "daily":
            from obspy import read
            tr = read(dst_fn, format="MSEED")[0]
            if tr.stats.npts != ref.size:
                return "npts %d != %d" % (tr.stats.npts, ref.size)
            if not np.array_equal(tr.data.astype(np.float64), ref.values.astype(np.float64)):
                return "values differ (max abs diff %.3g)" % np.abs(
                    tr.data - ref.values).max()
            return None
        df = pd.read_hdf(dst_fn, "data")
        if df.shape != ref.shape:
            return "shape %s != %s" % (df.shape, ref.shape)
        if not np.array_equal(df.values, ref.values):
            return "values differ"
        t_src = pd.to_datetime(ref.coords["times"].values)
        t_dst = pd.to_datetime(df.index)
        if not (t_src == t_dst).all():
            return "window times differ"
        if not np.allclose(df.columns.values.astype(float),
                           ref.coords["taxis"].values, atol=1e-6):
            return "lag axis differs"
        return None
    except Exception as e:
        return "cannot read: %s" % e


def cmd_verify(a):
    import random
    src = Source(a.src)
    dst_engine, dst_ini = make_engine(a.dst)
    dt = reflect(dst_engine, dst_ini)
    report = []
    status = {"PASS": 0, "WARN": 0, "FAIL": 0}

    def log(level, msg):
        status[level] += 1
        line = "[%s] %s" % (level, msg)
        report.append(line)
        print(line, flush=True)

    def rows(table):
        with dst_engine.connect() as c:
            return [dict(r._mapping) for r in c.execute(select(dt(table)))]

    cfg = src.config()
    sets = pick_sets(cfg, a)
    dcfg = {r["name"]: r["value"] for r in rows("config")}

    # ---- 1. CC-defining parameters -------------------------------------
    checks = [("global", n) for n in ("startdate", "enddate")] + \
             [("preprocess", n) for n in ("preprocess_lowpass", "preprocess_highpass",
                                          "resampling_method", "remove_response")] + \
             [("cc", n) for n in ("maxlag", "cc_sampling_rate", "corr_duration",
                                  "overlap", "winsorizing", "whitening", "whitening_type",
                                  "cc_normalisation", "stack_method",
                                  "components_to_compute", "cc_type", "keep_days")] + \
             [("stack", "mov_stack"), ("refstack", "ref_begin"), ("refstack", "ref_end")]
    if not a.keep_days_only:
        checks.append(("cc", "keep_all"))
    bad = []
    n_checked = 0
    for cat, name in checks:
        key = (cat, None) if cat == "global" else (cat, sets.get(cat))
        if name not in cfg.get(key, {}) or name not in dcfg:
            continue
        n_checked += 1
        if not _num_equal(cfg[key][name], dcfg[name]):
            bad.append("%s: msnoise_c %s=%r, stable=%r" % (name, cat, cfg[key][name], dcfg[name]))
    if bad:
        for b in bad:
            log("FAIL", "config " + b)
    else:
        log("PASS", "config: %d CC-defining parameters identical" % n_checked)
    pre_sr = cfg.get(("preprocess", sets.get("preprocess")), {}).get("cc_sampling_rate")
    cc_sr = cfg.get(("cc", sets.get("cc")), {}).get("cc_sampling_rate")
    if pre_sr and cc_sr and not _num_equal(pre_sr, cc_sr):
        log("WARN", "msnoise_c preprocess cc_sampling_rate=%s but cc cc_sampling_rate=%s"
            % (pre_sr, cc_sr))

    # ---- 2. filters ------------------------------------------------------
    dfil = {r["ref"]: r for r in rows("filters")}
    fbad = []
    nf = 0
    for (cat, s), c in cfg.items():
        if cat != "filter" or s is None:
            continue
        nf += 1
        f = dfil.get(s)
        if f is None:
            fbad.append("filter_%d missing" % s)
        elif not (_num_equal(f["low"], c.get("freqmin")) and _num_equal(f["high"], c.get("freqmax"))):
            fbad.append("filter_%d: %s-%s Hz vs %s-%s Hz" % (
                s, c.get("freqmin"), c.get("freqmax"), f["low"], f["high"]))
    for b in fbad:
        log("FAIL", b)
    if not fbad:
        log("PASS", "filters: %d filters, same frequency bands" % nf)

    # ---- 3. stations -----------------------------------------------------
    s_src = {(r["net"], r["sta"]): bool(r["used"]) for r in src.rows("stations")}
    s_dst = {(r["net"], r["sta"]): bool(r["used"]) for r in rows("stations")}
    if s_src == s_dst:
        log("PASS", "stations: %d stations, same 'used' flags" % len(s_src))
    else:
        diff_used = sorted(".".join(k) for k in set(s_src) & set(s_dst) if s_src[k] != s_dst[k])
        log("FAIL", "stations differ: only in msnoise_c %s, only in stable %s, "
            "'used' flag differs for %s" % (
                sorted(".".join(k) for k in set(s_src) - set(s_dst)),
                sorted(".".join(k) for k in set(s_dst) - set(s_src)), diff_used))

    # ---- 4. data availability -------------------------------------------
    da = rows("data_availability")
    n_src = len(src.rows("data_availability"))
    if len(da) == n_src:
        log("PASS", "data_availability: %d rows in both" % n_src)
    else:
        log("FAIL", "data_availability: %d rows in msnoise_c, %d in stable" % (n_src, len(da)))
    sample = random.Random(0).sample(da, min(50, len(da)))
    missing = [os.path.join(r["path"], r["file"]) for r in sample
               if not os.path.isfile(os.path.join(r["path"], r["file"]))]
    if missing:
        log("WARN", "data_availability: %d/%d sampled raw files not found on this "
            "machine (e.g. %s) - wrong data_folder, or data not mounted here?"
            % (len(missing), len(sample), missing[0]))
    elif sample:
        log("PASS", "data_availability: %d sampled raw data paths exist on disk" % len(sample))

    # ---- 5. jobs ---------------------------------------------------------
    steps = src.steps()
    cc_ids = {k for k, v in steps.items()
              if v["category"] == "cc" and v["set_number"] == sets.get("cc")}
    src_done = {(r["day"], r["pair"]) for r in src.rows("jobs")
                if r.get("step_id") in cc_ids and r["flag"] == "D"}
    dst_done = {(r["day"], r["pair"]) for r in rows("jobs")
                if r["jobtype"] == "CC" and r["flag"] == "D"}
    if src_done <= dst_done:
        log("PASS", "jobs: all %d finished CC day-pairs are marked done in stable" % len(src_done))
    else:
        log("FAIL", "jobs: %d finished CC day-pairs not marked done in stable"
            % len(src_done - dst_done))

    # ---- 6. CC files ----------------------------------------------------
    maxlag = float(dcfg["maxlag"])
    sr = float(dcfg["cc_sampling_rate"])
    n_lag = int(2 * maxlag * sr) + 1
    out_folder = dcfg.get("output_folder", "CROSS_CORRELATIONS")
    if not os.path.isabs(out_folder):
        out_folder = os.path.join(os.path.abspath(a.dst), out_folder)
    stacks = os.path.join(os.path.abspath(a.dst), "STACKS")
    roots = find_cc_roots(src, a.cc_set, a.preprocess_set)
    if len(roots) != 1:
        log("FAIL", "expected exactly one cc_%d folder, found %d" % (a.cc_set, len(roots)))
        roots = roots[:1]
    pairs = []   # (kind, src, dst)
    per_group = collections.Counter()
    kinds = ["daily"] if a.keep_days_only else ["daily", "all"]
    for root in roots:
        for fdir in sorted(glob.glob(os.path.join(root, "filter_*"))):
            fid = int(os.path.basename(fdir).split("_")[1])
            for kind in kinds:
                base = os.path.join(fdir, "_output", kind)
                for fn in sorted(glob.glob(os.path.join(base, "*", "*", "*.nc"))):
                    comp = fn.split(os.sep)[-3]
                    pair = fn.split(os.sep)[-2]
                    s1, s2 = pair.split("_")
                    day = os.path.basename(fn)[:-3]
                    if kind == "daily":
                        d = os.path.join(stacks, "%02i" % fid, "001_DAYS", comp, pair, day + ".MSEED")
                    else:
                        d = os.path.join(out_folder, "%02i" % fid, s1, s2, comp, day + ".h5")
                    pairs.append((kind, fn, d))
                    per_group[(fid, kind, comp)] += 1
    for (fid, kind, comp), n in sorted(per_group.items()):
        report.append("       filter %02d %-5s %s : %d source files" % (fid, kind, comp, n))

    absent = [p for p in pairs if not os.path.isfile(p[2])]
    if absent:
        log("FAIL", "files: %d of %d CC files were not converted (e.g. %s) - rerun the "
            "'files' step" % (len(absent), len(pairs), absent[0][1]))
    else:
        log("PASS", "files: all %d CC files have their stable counterpart" % len(pairs))

    # orphans: stable files with no msnoise_c source (old runs, wrong cc set...)
    expected = {p[2] for p in pairs}
    found = set(glob.glob(os.path.join(stacks, "*", "001_DAYS", "*", "*", "*.MSEED")))
    if not a.keep_days_only:
        found |= set(glob.glob(os.path.join(out_folder, "*", "*", "*", "*", "*.h5")))
    orphans = found - expected
    if orphans:
        log("WARN", "files: %d files in the stable project do not come from this "
            "msnoise_c cc set (e.g. %s)" % (len(orphans), sorted(orphans)[0]))

    # lag axis
    if pairs:
        import xarray as xr
        with xr.open_dataset(pairs[0][1]) as d:
            n = d.sizes["taxis"]
        if n == n_lag:
            log("PASS", "lag axis: %d samples = 2*maxlag(%g)*cc_sampling_rate(%g)+1"
                % (n, maxlag, sr))
        else:
            log("FAIL", "lag axis: files have %d samples, stable config expects %d" % (n, n_lag))

    # content, exact numbers
    present = [p for p in pairs if os.path.isfile(p[2])]
    if a.sample and a.sample < len(present):
        present = random.Random(1).sample(present, a.sample)
        what = "a random sample of %d" % len(present)
    else:
        what = "all %d" % len(present)
    print("Comparing numbers in %s converted files..." % what, flush=True)
    errs = []
    if a.threads > 1:
        with ProcessPoolExecutor(a.threads) as ex:
            futs = {ex.submit(_check_file, *p): p for p in present}
            for f in as_completed(futs):
                r = f.result()
                if r:
                    errs.append("%s: %s" % (futs[f][2], r))
    else:
        for i, p in enumerate(present, 1):
            r = _check_file(*p)
            if r:
                errs.append("%s: %s" % (p[2], r))
            if i % 2000 == 0:
                print("  %d / %d" % (i, len(present)), flush=True)
    if errs:
        log("FAIL", "content: %d files differ from msnoise_c" % len(errs))
        for e in errs[:20]:
            report.append("       " + e)
    else:
        log("PASS", "content: %s converted files hold exactly the same numbers, "
            "window times and lag axis as the msnoise_c files" % what)

    # ---- summary ---------------------------------------------------------
    verdict = "FAIL" if status["FAIL"] else ("PASS with warnings" if status["WARN"] else "PASS")
    head = ["MSNoise migration check: msnoise_c -> stable",
            "source      : %s" % src.dir,
            "destination : %s" % os.path.abspath(a.dst),
            "cc set      : %d   config sets used: %s" % (a.cc_set, sets),
            "RESULT      : %s  (%d pass, %d warn, %d fail)" % (
                verdict, status["PASS"], status["WARN"], status["FAIL"]),
            ""]
    fn = os.path.join(a.dst, "migration_check.txt")
    with open(fn, "w") as f:
        f.write("\n".join(head + report) + "\n")
    print("\n" + head[4])
    print("Report written to %s" % fn)
    sys.exit(1 if status["FAIL"] else 0)


# --------------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    def common(sp, dst=True):
        sp.add_argument("--src", required=True, help="msnoise_c project folder (contains db.ini)")
        if dst:
            sp.add_argument("--dst", required=True, help="classic project folder (contains db.ini)")
        sp.add_argument("--cc-set", type=int, default=1, help="cc config set to migrate (default 1)")
        sp.add_argument("--preprocess-set", type=int, default=None,
                        help="preprocess set (only needed if you have several)")
        sp.add_argument("--keep-days-only", action="store_true",
                        help="only migrate daily stacks (sets keep_all=N)")

    common(sub.add_parser("inspect"), dst=False)
    sp = sub.add_parser("db")
    common(sp)
    for cat in CATEGORY_PRIORITY[3:]:
        sp.add_argument("--%s-set" % cat.replace("_", "-"), type=int, default=None,
                        dest="%s_set" % cat)
    sp.add_argument("--output-folder", default="CROSS_CORRELATIONS",
                    help="classic output_folder for keep_all files (default CROSS_CORRELATIONS)")
    sp.add_argument("--classic-ref", default=None,
                    help="an OLD classic project folder (with db.ini) to copy "
                         "extra custom config keys from (e.g. coda_safety_factor, ref_type)")
    sp.add_argument("--force", action="store_true",
                    help="wipe filters/stations/DA/jobs in destination and refill")
    sp = sub.add_parser("files")
    common(sp)
    sp.add_argument("--threads", type=int, default=1)
    sp.add_argument("--overwrite", action="store_true")
    sp.add_argument("--dry-run", action="store_true")
    sp = sub.add_parser("verify")
    common(sp)
    sp.add_argument("--threads", type=int, default=1)
    sp.add_argument("--sample", type=int, default=0,
                    help="compare the numbers of only N random files (default: all)")
    a = p.parse_args()
    {"inspect": cmd_inspect, "db": cmd_db, "files": cmd_files,
     "verify": cmd_verify}[a.cmd](a)


if __name__ == "__main__":
    main()
