#!/usr/bin/env python3
"""
Set up an MSNoise project in one run: config values, filters and stations.

Edit the three sections below (CONFIG, FILTERS, STATIONS), then run from the
project folder (the one containing db.ini):

    python setup_msnoise_project.py

It is safe to re-run: existing values, filters and stations are updated in place.
Check the result with:  msnoise info
"""
import sys

from msnoise.api import (connect, get_config, update_config, get_filters,
                         update_filter, update_station, get_station)
from msnoise.default import default

# ---------------------------------------------------------------------------
# 1. CONFIG -- any MSNoise parameter (see `msnoise info` for the full list).
#    Values are given as strings, exactly as you would type them in
#    `msnoise config set name=value`. Comment out what you don't want to touch.
#    MUST match how STACKS2 was built: maxlag, cc_sampling_rate, mov_stack,
#    components (use the values printed by recover_from_stacks2.py).
# ---------------------------------------------------------------------------
CONFIG = {
    # --- must match STACKS2 ---
    "maxlag": "120",
    "cc_sampling_rate": "25",
    "mov_stack": "(('10D','1D'),('1D','1D'))",
    "components_to_compute": "EE,EN,EZ,NN,NZ,ZZ",
    "components_to_compute_single_station": "EE,EN,EZ,NN,NZ,ZZ",

    # --- period / reference ---
    "startdate": "2024-01-01",
    "enddate": "2026-12-31",
    "ref_begin": "2024-01-01",
    "ref_end": "2026-12-31",

    # --- dt/t parameters ---
    "dtt_lag": "dynamic",
    "dtt_v": "1.0",
    "dtt_minlag": "5.0",
    "dtt_width": "30.0",
    "dtt_sides": "both",
    "dtt_mincoh": "0.0",
    "dtt_maxerr": "0.1",
    "dtt_maxdt": "0.1",
}

# ---------------------------------------------------------------------------
# 2. FILTERS -- `ref` must match the folder numbers in STACKS2 (01 -> 1, ...).
#    low/high: whitening band (Hz) used for the CC; mwcs_low/mwcs_high: band
#    of the MWCS linear regression; mwcs_wlen/mwcs_step in seconds.
# ---------------------------------------------------------------------------
FILTERS = [
    dict(ref=1, low=0.1, high=25.0, mwcs_low=0.11, mwcs_high=24.99,
         mwcs_wlen=12.0, mwcs_step=4.0, used=True),
    # dict(ref=2, low=1.0, high=2.0, mwcs_low=1.0, mwcs_high=2.0,
    #      mwcs_wlen=6.0, mwcs_step=2.0, used=True),
]

# ---------------------------------------------------------------------------
# 3. STATIONS -- either fill STATIONS by hand, or set STATIONXML to a file/glob
#    and list only the stations you want in STATION_IDS.
#    loc / chan must reproduce the names in STACKS2 (NET.STA.LOC), e.g. for
#    "BE.UCC." the location code is "" ; for "PF.BON.00" it is "00".
#    coordinates: "DEG" (lon/lat in degrees) or "UTM" (metres).
# ---------------------------------------------------------------------------
STATIONS = [
    # dict(net="PF", sta="BON", loc="00", chan="HH?",
    #      X=55.7, Y=-21.2, altitude=2000.0, coordinates="DEG"),
]

STATIONXML = "Responses/.xml"            # e.g. "inventory/*.xml"
STATION_IDS = []             # e.g. ["PF.BON.00", "PF.FOR.00"] (NET.STA.LOC)
STATIONXML_CHAN = "HH?"      # channel names stored for stations from StationXML

# ---------------------------------------------------------------------------


def apply_config(db):
    print("== Config")
    for name, value in CONFIG.items():
        if name not in default:
            print("  !! unknown parameter '%s', skipped" % name)
            continue
        old = get_config(db, name)
        update_config(db, name, str(value))
        flag = "  " if str(old) == str(value) else " M"
        print("%s %s: %s -> %s" % (flag, name, old, value))


def apply_filters(db):
    print("== Filters")
    for f in sorted(FILTERS, key=lambda f: f["ref"]):
        update_filter(db, f["ref"], f["low"], f["mwcs_low"], f["high"],
                      f["mwcs_high"], f["mwcs_wlen"], f["mwcs_step"], f["used"])
    # New filters get an auto-incremented ref: make sure they ended up right.
    wanted = {f["ref"] for f in FILTERS}
    got = {f.ref for f in get_filters(db, all=True)}
    for f in get_filters(db, all=True):
        print("  %2i: [%.3f:%.3f] mwcs [%.3f:%.3f] wlen=%g step=%g used=%s"
              % (f.ref, f.low, f.high, f.mwcs_low, f.mwcs_high,
                 f.mwcs_wlen, f.mwcs_step, f.used))
    if wanted - got:
        print("  !! filter refs %s were not created with the expected ref. "
              "Create filters in order 1, 2, 3... on an empty filter table."
              % sorted(wanted - got))


def stations_from_stationxml():
    import glob
    from obspy import read_inventory
    inv = None
    for fn in sorted(glob.glob(STATIONXML)):
        inv = read_inventory(fn) if inv is None else inv + read_inventory(fn)
    if inv is None:
        sys.exit("No StationXML found for %s" % STATIONXML)
    out = []
    for sid in STATION_IDS:
        net, sta, loc = sid.split(".")
        sel = inv.select(network=net, station=sta)
        if not len(sel):
            print("  !! %s not found in StationXML, skipped" % sid)
            continue
        s = sel[0][0]
        out.append(dict(net=net, sta=sta, loc=loc, chan=STATIONXML_CHAN,
                        X=s.longitude, Y=s.latitude, altitude=s.elevation,
                        coordinates="DEG"))
    return out


def apply_stations(db):
    print("== Stations")
    stations = list(STATIONS)
    if STATIONXML:
        stations += stations_from_stationxml()

    # Merge location codes when a station appears with several loc codes
    merged = {}
    for s in stations:
        key = (s["net"], s["sta"])
        if key in merged:
            locs = merged[key]["locs"]
            if s["loc"] not in locs:
                locs.append(s["loc"])
        else:
            merged[key] = dict(s, locs=[s["loc"]])

    for (net, sta), s in merged.items():
        update_station(db, net, sta, float(s["X"]), float(s["Y"]),
                       float(s["altitude"]), coordinates=s["coordinates"],
                       used=1)
        st = get_station(db, net, sta)
        st.used_location_codes = ",".join(s["locs"])
        st.used_channel_names = s["chan"]
        db.commit()
        print("  %s.%s  X=%.4f Y=%.4f alt=%.1f %s  locs=%s chans=%s"
              % (net, sta, st.X, st.Y, st.altitude, st.coordinates,
                 st.used_location_codes or "''", st.used_channel_names))


if __name__ == "__main__":
    db = connect()
    apply_config(db)
    apply_filters(db)
    apply_stations(db)
    db.close()
    print("\nDone. Check with: msnoise info")
