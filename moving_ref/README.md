s08compute_wct.py replace the one in msnoise package (before Alec version/msnoise2)
The others are local to the work directory.
To run moving reference wct:
1. run add_pram_moving_ref, to modify and add two parameters that where not in the original config.
2. run msnoise cc dvv compute_wct command with the new file in the package
3. once all jobs ran, run mer_wct_day_track to get nc file per station pairs
4. run delta_dvv_integration.py to retrive dv/v from the moving reference (without integration the output of wct are delta-dv/v)
