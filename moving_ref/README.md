s08compute_wct.py replace the one in msnoise package (before Alec version/msnoise2) The others are local to the work directory. To run moving reference wct:

1. run add_pram_moving_ref, to modify and add two parameters that where not in the original config.
2. run msnoise cc dvv compute_wct command with the new file in the package
3. once all jobs ran, run mer_wct_day_track to get nc file per station pairs
4. run delta_dvv_integration.py to retrive dv/v from the moving reference (without integration the output of wct are delta-dv/v)

# WCT 

A three-stage pipeline for computing, merging, and integrating seismic velocity
changes (dv/v) using WCT within the MSNoise framework.

---

## Pipeline Overview

```
s08compute_wct.py        →    merge_wct_day_track.py    →    delta_dvv_integration.py
(daily .npz per pair)        (merged .nc per pair)           (reconstructed absolute dv/v)
```

---

## Stage 1 — `s08compute_wct.py`: Delta dv/v Computation

### What it does

For each station pair, day, filter, component and moving-stack window, computes
the WCT between a **current** stacked CCF and a **reference** stacked CCF.
The phase delay extracted by the WCT gives:

$$\Delta\frac{\delta v}{v}(t) = \frac{v(t) - v_\mathrm{ref}(t)}{v_\mathrm{ref}(t)}$$

where $v(t)$ is the current stack and $v_\mathrm{ref}(t)$ is the reference stack
averaged over the configured reference window.

Results are saved as daily `.npz` files under:

```
WCT/[filter_id]/[mov_stack]/[component]/[station_pair]/YYYY-MM-DD.npz
```

### Reference types

#### Static reference (`ref_type = static`)

A single pre-computed long-term stack is used as reference for every date.
This is the classical MSNoise approach and produces dv/v relative to a fixed
baseline.

#### Moving reference (`ref_type = moving`)

The reference is recomputed for each date $t$ from stacks falling within a
sliding window defined by two parameters:

| Parameter | Config key | Meaning |
|-----------|------------|---------|
| `ref_end` | `ref_end` | End of the reference window relative to $t$ |
| `ref_begin` | `ref_begin` | Width of the reference window |

Both parameters accept either `'auto'` or an **integer number of days**.

### Moving reference parameter cases

#### Case 1 — `ref_end = auto`, `ref_begin = auto` (default)

```
ref_end   = t - stack_size
ref_start = ref_end - stack_size
```

The reference is the mean of all stacks within the single preceding non-overlapping window. With gap-free regularly-stepped data this window contains exactly one stack, reducing to a simple lag-1 difference. The measured
quantity is the velocity change between two consecutive stacks:

$$\Delta\frac{\delta v}{v}(t) = \frac{\delta v}{v}(t)_\mathrm{curr} - \frac{\delta v}{v}(t-1)_\mathrm{ref}$$

This is the simplest and most common case. With 10-day stacks it gives a
10-day differential velocity measurement. The absolute dv/v is recovered by
cumulative summation (see Stage 3).

#### Case 2 — `ref_end = auto`, `ref_begin = N` (integer days)

```
ref_end   = t - stack_size
ref_start = ref_end - N
```

The reference is the **mean of all stacks within a window of width N days**
ending one stack before the current date. A larger N reduces noise in the
reference at the cost of smoothing rapid velocity changes in the reference
baseline. The integration step uses the same window width to reconstruct the
absolute dv/v consistently.

Useful when individual stacks are noisy and a more stable reference baseline
is needed.

#### Case 3 — `ref_end = M` (integer days), `ref_begin = auto`

```
ref_end   = t + M          (M is negative: window ends M days before t)
ref_start = ref_end - stack_size
```

The reference window is shifted by a fixed lag M relative to the current date.
If M is negative (the typical case), the reference ends further in the past
than the default auto case. This increases the temporal separation between
reference and current stack, making the measurement sensitive to slower
velocity changes while keeping a single-stack reference.

A warning is raised if `ref_end` overlaps with the current stack window.

#### Case 4 — `ref_end = M`, `ref_begin = N` (both integers)

```
ref_end   = t + M
ref_start = ref_end - N
```

Full control over the reference window. The reference is the mean of all
stacks within a window of width N days ending M days relative to the current
date. This generalises all previous cases:

- M controls the **lag** between current and reference states
- N controls the **smoothing width** of the reference

The integration step reconstructs absolute dv/v by anchoring each point on
the mean of already-reconstructed values over the same `[t+M-N, t+M]` window,
which is the correct inverse of the WCT measurement (see Stage 3).

### Adaptive coda window

Station distances are pre-cached from the MSNoise database. For each pair the
coda window start is set to:

```
coda_start = max(distance / velocity * safety_factor, dtt_minlag)
```

For autocorrelations (distance = 0) the minimum lag is used directly.

### Metadata saved per daily file

Each `.npz` file stores, alongside `dvv`, `err`, `coh` and `freqs`:

| Field | Description |
|-------|-------------|
| `curr_window_days` | Width of the current stack in days |
| `ref_window_days` | Width of the reference window (N) in days |
| `step_days` | Step between consecutive stacks in days |
| `ref_lag_days` | Lag from current date to end of reference window (M) in days |
| `actual_separation` | Actual calendar distance to last reference stack (gap monitor) |
| `coda_start` | Coda window start time used for this pair |
| `coda_constraint` | Whether coda start was set by physics or minimum lag |

---

## Stage 2 — `merge_wct_day_track.py`: Merging Daily Files

### What it does

Consolidates all daily `.npz` files for each combination of filter, moving
stack, component and station pair into a single NetCDF file:

```
WCT_MERGED/[filter_id]/[mov_stack]/[component]/[station_pair].nc
```

Each `.nc` file contains three 2-D arrays (`times × freq`): `dvv`, `err`,
`coh`, plus the metadata fields listed above stored as dataset attributes.

### Usage

```bash
python merge_wct_day_track.py [--wct_dir WCT] [--output_dir WCT_MERGED] [--verbose]
```

A progress log (`merge_progress.log`) is maintained in the output directory so
interrupted runs can be resumed without reprocessing completed pairs.

CSV copies of `dvv`, `err` and `coh` are also written under
`WCT_MERGED/csv/...` for quick inspection.

---

## Stage 3 — `delta_dvv_integration.py`: Absolute dv/v Reconstruction

### What it does

Integrates the $\Delta\delta v/v$ time series produced by Stage 1 into an
absolute $\delta v/v$ time series anchored at zero at `start_date`.

### Reconstruction formula

For each time step $t$:

$$\frac{\delta v}{v}(t) = \frac{1}{K} \sum_{k \in [t+M-N,\, t+M]} \frac{\delta v}{v}(t_k) + \Delta\frac{\delta v}{v}(t)$$

where the sum is over already-reconstructed values whose timestamps fall
within the same reference window `[t+M-N, t+M]` that was used during the WCT
step. This is the exact inverse of the WCT measurement and is therefore
self-consistent for all parameter cases.

The window bounds `M` (`ref_lag_days`) and `N` (`ref_window_days`) are read
directly from the dataset attributes written by Stage 2, so no manual
configuration is needed — the integrator automatically adapts to whatever
parameters were used during Stage 1.

### Special cases recovered by the general formula

| ref_end | ref_begin | Reconstruction reduces to |
|---------|-----------|--------------------------|
| auto | auto | Simple cumulative sum: $v(t) = v(t-1) + \Delta\delta v/v(t)$ when step_days = stack_size (one stack per window).|
| auto | N | Mean of N/step past values as anchor |
| M | auto | Single past value lagged by M days |
| M | N | Mean of window `[t-M-N, t-M]` as anchor |

### Gap handling

A gap is detected when the time step between consecutive entries exceeds
`3 × step_days`. At a gap:

- The last reconstructed value before the gap is used as anchor
- If $\Delta\delta v/v(t)$ is NaN at the gap boundary, the anchor is carried
  forward unchanged
- This recovers the net velocity change across the gap, the total difference between the last pre-gap and first post-gap stack, but temporal variations within the gap are not resolved.



Output files mirror the input directory structure under `WCT_INTEGRATED/`.

---

## Configuration reference (MSNoise parameters)

| Parameter | Type | Description |
|-----------|------|-------------|
| `ref_type` | `static` / `moving` | Reference type |
| `ref_end` | `auto` / int (days) | End of reference window relative to current date |
| `ref_begin` | `auto` / int (days) | Width of reference window in days |
| `dtt_minlag` | float (s) | Minimum coda window start time |
| `dtt_v` | float (km/s) | Surface wave velocity for coda window |
| `coda_safety_factor` | float | Safety multiplier on surface wave travel time |
| `dtt_codacycles` | int | Number of coda cycles defining window end |
| `dtt_mincoh` | float | Minimum wavelet coherence for weighting |
| `dtt_maxdt` | float (s) | Maximum allowed time delay |
| `dvv_min_nonzero` | float | Minimum fraction of non-zero weights for valid estimate |
| `wct_norm` | bool | Normalise waveforms before WCT |
| `wavelet_type` | list | Wavelet type and parameter, e.g. `['Morlet', 6]` |
| `wct_ns` | int | Smoothing parameter (scale direction) |
| `wct_nt` | float | Smoothing parameter (time direction) |
| `wct_vpo` | int | Voices per octave (scale resolution) |
| `wct_nptsfreq` | int | Number of frequency points |
