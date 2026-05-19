# register_wct_params.py
"""
Run this once to inject WCT custom parameters into the MSNoise config table.
After running, use get_config(db, "coda_velocity") etc. in your WCT code.

Usage:
    python register_wct_params.py
    python register_wct_params.py --coda_velocity 2.5 --dtt_minlag 8.0
"""

import argparse
from msnoise.api import connect, get_config
from msnoise.msnoise_table_def import Config

WCT_DEFAULTS = {
    #"dtt_v":      ("2.0",    "Surface wave velocity in km/s for coda window detection"),
    "coda_safety_factor": ("1.2",    "Safety factor applied to surface wave travel time"),
    #"dtt_minlag":         ("6.0",    "Minimum coda lag time in seconds"),
    "ref_type":           ("moving", "Reference type: 'static' or 'moving'"),
    #"ref_end":            ("auto",   "Reference window end: 'auto' or integer days relative to current date"),
    #"ref_begin":          ("auto",   "Reference window start: 'auto' or integer days before ref_end"),
}

def set_or_update_config(session, name, value):
    """Insert if not exists, update if exists — what update_config lacks."""
    existing = session.query(Config).filter(Config.name == name).first()
    if existing is None:
        session.add(Config(name=name, value=value))
    else:
        existing.value = value
    session.commit()

def register_wct_params(db, overrides=None):
    overrides = overrides or {}
    for key, (default, description) in WCT_DEFAULTS.items():
        existing = get_config(db, key)
        if not existing:  # get_config returns '' (empty string) when not found
            value = overrides.get(key, default)
            set_or_update_config(db, key, value)
            print(f"  [NEW]      {key} = {value}  ({description})")
        else:
            if key in overrides:
                set_or_update_config(db, key, overrides[key])
                print(f"  [UPDATED]  {key} = {overrides[key]}  (was: {existing})")
            else:
                print(f"  [EXISTING] {key} = {existing}  (unchanged)")

def main():
    parser = argparse.ArgumentParser(description="Register WCT parameters in MSNoise config table")
    for key, (default, description) in WCT_DEFAULTS.items():
        parser.add_argument(f"--{key}", default=None, help=f"{description} (default: {default})")
    args = parser.parse_args()
    overrides = {k: v for k, v in vars(args).items() if v is not None}

    db = connect()
    print("Registering WCT parameters in MSNoise config table...")
    register_wct_params(db, overrides)
    db.close()
    print("Done.")

if __name__ == "__main__":
    main()
