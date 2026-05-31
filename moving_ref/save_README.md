s08compute_wct.py replace the one in msnoise package (before Alec version/msnoise2)
The others are local to the work directory.
To run moving reference wct:
1. run add_pram_moving_ref, to modify and add two parameters that where not in the original config.
2. run msnoise cc dvv compute_wct command with the new file in the package
3. once all jobs ran, run mer_wct_day_track to get nc file per station pairs
4. run delta_dvv_integration.py to retrive dv/v from the moving reference (without integration the output of wct are delta-dv/v)


Moving reference parameters:
ref_end / ref_start
│
├── 'auto' / 'auto'  
│   ref_end   = date - stack_size
│   ref_start = implicit (every stack is stack_size long)
│   → always takes most recent stack before ref_end
│   → warns if separation > stack_size + step_days (gap in data)
│   → never skips
│
├── 'auto' / int(N)
│   ref_end   = date - stack_size  
│   ref_start = ref_end - N        (explicit fixed-width window)
│   → takes most recent stack in [ref_start, ref_end]
│   → if empty: fallback to most recent before ref_end + warn
│   → never skips
│
├── int(M) / 'auto'
│   ref_end   = date + M           (warn if M > -stack_size, overlap risk)
│   ref_start = ref_end - stack_size  (width adapts to stack_size)
│   → takes most recent stack in [ref_start, ref_end]
│   → if empty: fallback to most recent before ref_end + warn
│   → never skips
│
└── int(M) / int(N)
    ref_end   = date + M           (warn if M > -stack_size, overlap risk)
    ref_start = ref_end - N        (fully explicit)
    → takes most recent stack in [ref_start, ref_end]
    → if empty: fallback to most recent before ref_end + warn
    → never skips

Specific parameters used:
    "dtt_v":              ("2.0",    "Surface wave velocity in km/s for coda window detection"),
    "coda_safety_factor": ("1.2",    "Safety factor applied to surface wave travel time"),
    "dtt_minlag":         ("6.0",    "Minimum coda lag time in seconds"),
    "ref_type":           ("moving", "Reference type: 'static' or 'moving'"),
    "ref_end":            ("auto",   "Reference window end: 'auto' or integer days relative to current date"),
    "ref_begin":          ("auto",   "Reference window start: 'auto' or integer days before ref_end"),
