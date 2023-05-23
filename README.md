# [MSNoise](https://github.com/ROBelgium/MSNoise) mutations
## Next mutations
* Fix plot dvv pairs
* Zoom in dvv: Compute dvv on a moving window

## Zoom in dvv
**Why?**
&nbsp;&nbsp;&nbsp; because 1 dvv point per day\
**Computation**
&nbsp;&nbsp;&nbsp; Start from compute_cc with config parameter ```keep_all=Y```\
Modified version of s05compute_mwcs.py and s06compute_dtt.py\
**Plot**
&nbsp;&nbsp;&nbsp; Modified version of plots/dvv.py \
**Result**
&nbsp;&nbsp;&nbsp; Compute and plot dvv with a higher sampling than 1 point per day based on ```corr_duration``` config parameter\
**Output** &nbsp;&nbsp;&nbsp; zoom_s05compute_mwcs.py, zoom_s06compute_dtt.py, plots/zoomerrdvv.py\
**To do next**
&nbsp;&nbsp;&nbsp; Compute dvv on a moving window

## Detect drastic change in dvv variations
**Why?** &nbsp;&nbsp;&nbsp; See remarkable change in dvv variations (add an argument to simultaneous variation)\
**Computation** &nbsp;&nbsp;&nbsp; Compute all slopes. Determine a statistical threshold based on all slopes. Detect slope above the threshold.\
**Results** &nbsp;&nbsp;&nbsp; Highlight period of dvv outlier changes.\
**Output** &nbsp;&nbsp;&nbsp; plots/drastic_slopes.py

## Evaluate the coss-correlation parameters used
**Why?** &nbsp;&nbsp;&nbsp; What parameters use to compute the cross-correlation?\
**Computation** &nbsp;&nbsp;&nbsp; Use cross-correlation computed with a set of config parameters. Compute SNR of all the CC (=reference) and randomly stacked CC. Compute a ratio between the two SNR considering data availability of the station pair.\
**Result** &nbsp;&nbsp;&nbsp; Give a score for each station pair and set of parameters. The score is the average number of CC stacked to reach 10% of the SNR reference.\
**Output** &nbsp;&nbsp;&nbsp; plots/snr_score.py

## Plot error of dvv
**Why?** &nbsp;&nbsp;&nbsp; See when a value reliability \
**Output** &nbsp;&nbsp;&nbsp; plots/errdvv.py

## Plot several dvv curves on the same figure
**Why?** &nbsp;&nbsp;&nbsp; Compare dvv variation from one year to another. Have a first look at seasonal variations.\
**Output** &nbsp;&nbsp;&nbsp; plots/multidvv.py, local/compare_multidvv.py

## Launch MSNoise on [CECI](https://www.ceci-hpc.be/)
**Why?** &nbsp;&nbsp;&nbsp; Save time\
**Result** &nbsp;&nbsp;&nbsp; Set up and submission scripts\
**Output** &nbsp;&nbsp;&nbsp; Parallelized_on_CECI/* \
**To do next** &nbsp;&nbsp;&nbsp; Further explanations here

![image](https://github.com/LaureBrenot/msnoise_mutations/assets/133853397/4a56b340-278d-4d4d-8d71-1a80b9f04172)

