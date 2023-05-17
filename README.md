# msnoise_mutations

# Modification/Addification to MSNoise (https://github.com/ROBelgium/MSNoise)
## Next mutations
Plot dvv pairs

## Zoom in dvv
**Why?**
1 dvv point per day\
**Computation**
Start from compute_cc with config parameter keep_all=Y\
Modified version of s05compute_mwcs.py and s06compute_dtt.py\
**Plot**
Modified version of plots/dvv.py \
**Result** 
Compute and plot dvv with an higher sampling than 1 point per day based on corr_duration config parameter\
**Output** zoom_s05compute_mwcs.py, zoom_s06compute_dtt.py, plots/zoomerrdvv.py\
**To do next**
Compute dvv on a moving window

## Plot error of dvv
**Why?** See when a value reliability \
**Output** plots/errdvv.py

### Plot several dvv curve on the same figure
**Why?:** Compare dvv variation from a year to another. Have a first look on seasonal variations.\
**Output** plots/multidvv.py

### Evaluate the coss-correlation parameters used
**Why?** What parameters use to compute the cross-correlation?\
**Computation:** Use cross-correlation computed with a set of config parameters. Compute SNR of all the CC (=reference) and randomly stacked CC. Compute a ratio between the two SNR considering data availabilty of the station pair.\
**Result:** Give a score for each station pair and set of parameter. The score is the average number of CC stacked to reach 10% of the SNR reference.\ 
**Output** plots/snr_score.py

### Detect drastic change in dvv variations
**Why?** See real change in dvv variations (add an argument to silmutaneous variation)\
**Computation:** Compute all slopes. Determine a stastiticall threshold based on all slopes. Detect slope above the threshold.\
**Results:** Highlight period of dvv outlier changes.\
**Output** plots/drastic_slopes.py
