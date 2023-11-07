# -*- coding: utf-8 -*-
"""
Created on Thu Sep 28 21:49:27 2023

@author: laure
modification based on: https://github.com/LaureBrenot/msnoise_mutations.git

from Lecocq, T., C. Caudron, et F. Brenguier (2014), 
MSNoise, a Python Package for Monitoring Seismic Velocity Changes Using Ambient Seismic Noise, 
Seismological Research Letters, 85(3), 715‑726, doi:10.1785/0220130073.
Doc: http://msnoise.org/doc/master/index.html
Gituhb: https://github.com/ROBelgium/MSNoise.git

NEEDS:
Create an environment of python 3.7
Activate environment
conda install -c conda-forge flask-admin flask-wtf markdown folium pymysql logbook tables
conda install -c conda-forge obspy
pip install msnoise
pip install tables

add in .../site-package/msnoise/msnoise_mutations.py
modify in .../site-package/msnoise/scripts/msnoise.py for msnoise_save.py
add in .../site-package/msnoise/scripts/msnoise.py (the new one)

READY TO WORK
This code plot dv/v with one point for 30min
"""

import os
from msnoise import s000installer, api
from msnoise.msnoise_mutations import process_job_type, zoomerrdvv
import subprocess

def setup_msnoise(data_location, data_output, start, end,  start_ref, end_ref):
    os.chdir(data_output)

    # Initialize MSNoise if it's the first time
    s000installer.main(tech=1)  # Tech 1=SQLite 2=MySQL 3=PostgreSQL', default=None

    db = api.connect()
    api.update_config(db, 'data_folder', data_location)
    api.update_config(db, 'output_folder', os.path.join(data_output, 'CROSS_CORRELATIONS'))
    api.update_config(db, 'startdate', start)
    api.update_config(db, 'enddate', end)
    
    # Populate Station Table and scan archive
    subprocess.run('msnoise populate', check=True, shell=True)
    subprocess.run('msnoise scan_archive --init', check=True, shell=True)
    # !msnoise populate
    # !msnoise scan_archive --init

    # Computation parameters
    api.update_config(db, 'components_to_compute', 'ZZ') # station pairs: cross-correlation
    #api.update_config(db, 'components_to_compute_single_station', 'EN,EZ,NZ') # single station: auto-correlation or cross-components
    api.update_config(db, 'remove_response', 'N') # remove instrument response, give path or in './inventory' folder
    api.update_config(db, 'keep_all', 'Y') # keep cross correlation to zoom in

    # Define filter parameters
    # update_filter(session, ref, low, mwcs_low, high, mwcs_high, rms_threshold, mwcs_wlen, mwcs_step, used), where low and high are the whitening bounds, mwcs_low and mwcs_high the frequency bounds for moving-window cross-spectral analysis, and mwcs_wlen and mwcs_step define the window length and window step used in this same analysis respectively.
    # filter 0.1 to 1Hz of the whitening function, filter 0.15 to 0.95Hz of the linear regression done in MWCS(looks like 5%), 0:not use anymore, 12sec windows to perform MWCS, 2sec step for windows, True activate the filter
    filters = [
        (1, 0.1, 0.15, 1.0, 0.95, 0, 12, 2, True),
        (2, 1.0, 1.05, 2.0, 1.95, 0, 12, 2, True),
        (3, 2.0, 2.05, 4.0, 3.95, 1, 12, 2, True),
        (4, 0.5, 0.55, 5.0, 4.95, 1, 12, 2, True),
        (5, 1.0, 1.05, 20.0, 19.95, 0, 12, 2, True)
    ]

    for params in filters:
        api.update_filter(db, *params)

    # Stack parameters
    api.update_config(db, 'mov_stack', '1') # day length windows, def=5
    api.update_config(db, 'ref_begin', start_ref) # can be relative, define the reference stack,
    api.update_config(db, 'ref_end', end_ref)

    # DTT parameters
    api.update_config(db, 'dtt_minlag', '10') # set minimum lag time, def=5
    api.update_config(db, 'dtt_maxerr', '0.2') # set maximum error, def=0.1
    api.update_config(db, 'dtt_maxdt', '0.2') # set maximum delay time, def=0.1
    api.update_config(db, 'dtt_mincoh', '0.6') # set minimum coherence, def=0.65
    
    # Check for new jobs
    # !msnoise db execute "select sta, count(*) from data_availability group by sta"
    subprocess.run('msnoise new_jobs --init', check=True, shell=True)
    # !msnoise new_jobs --init
    
def modify_msnoise(data_location, data_output, start, end,  start_ref, end_ref):
    db = api.connect()
    api.update_config(db, 'startdate', start)
    api.update_config(db, 'enddate', end)
    
    # Scan the archive
    subprocess.run('msnoise scan_archive --init', check=True, shell=True)
    # !msnoise scan_archive --init

    # Define filter parameters
    # update_filter(session, ref, low, mwcs_low, high, mwcs_high, rms_threshold, mwcs_wlen, mwcs_step, used), where low and high are the whitening bounds, mwcs_low and mwcs_high the frequency bounds for moving-window cross-spectral analysis, and mwcs_wlen and mwcs_step define the window length and window step used in this same analysis respectively.
    # filter 0.1 to 1Hz of the whitening function, filter 0.15 to 0.95Hz of the linear regression done in MWCS(looks like 5%), 0:not use anymore, 12sec windows to perform MWCS, 2sec step for windows, True activate the filter
    filters = [
        (1, 0.1, 0.15, 1.0, 0.95, 0, 12, 2, True),
        (2, 1.0, 1.05, 2.0, 1.95, 0, 12, 2, True),
        (3, 2.0, 2.05, 4.0, 3.95, 1, 12, 2, True),
        (4, 0.5, 0.55, 5.0, 4.95, 1, 12, 2, True),
        (5, 1.0, 1.05, 20.0, 19.95, 0, 12, 2, True)
    ]

    for params in filters:
        api.update_filter(db, *params)

    # Set up stack parameters
    api.update_config(db, 'mov_stack', '1') # day length windows, def=5
    api.update_config(db, 'ref_begin', start_ref) # can be relative, define the reference stack,
    api.update_config(db, 'ref_end', end_ref)

    # Set up DTT parameters
    api.update_config(db, 'dtt_minlag', '10') # set minimum lag time, def=5
    api.update_config(db, 'dtt_maxerr', '0.2') # set maximum error, def=0.1
    api.update_config(db, 'dtt_maxdt', '0.2') # set maximum delay time, def=0.1
    api.update_config(db, 'dtt_mincoh', '0.6') # set minimum coherence, def=0.65
    
    # Check for new jobs
    # !msnoise db execute "select sta, count(*) from data_availability group by sta"
    # !msnoise new_jobs
    subprocess.run('msnoise new_jobs', check=True, shell=True)
    
def compute_msnoise():
    db = api.connect()

    # Execute MSNoise jobs
    process_job_type(db, 'CC', 'msnoise compute_cc')
    process_job_type(db, 'STACK', 'msnoise stack -m')
    process_job_type(db, 'MWCS', 'msnoise compute_zoom_mwcs')
    process_job_type(db, 'DTT', 'msnoise compute_zoom_dtt')


if __name__ == "__main__":
    data_location = "C:/Users/laure/OneDrive - Université Libre de Bruxelles/Documents/AThese/Scripts/Scripts_O1/Compile_msnoise/avo_smart/Okmok_data_test"
    data_output = r"C:\Users\laure\OneDrive - Université Libre de Bruxelles\Documents\AThese\Scripts\Scripts_O1\Compile_msnoise\avo_smart\test1"
    start_date, end_date ='2021-08-05', '2021-08-07'
    start_ref, end_ref = start_date, end_date
    
    try:
        db = api.connect()
        if not api.get_stations(db, all=True):
            # new station added
            setup_msnoise(data_location, data_output, start_date, end_date, start_ref, end_ref) 
        else:
            # when modifying dates or parameters
            modify_msnoise(data_location, data_output, start_date, end_date, start_ref, end_ref) 
    except:
        # first time 
        setup_msnoise(data_location, data_output, start_date, end_date, start_ref, end_ref)

    # compute cross-correlation, stack, mwcs, dt/t
    compute_msnoise()
    # plot network dv/v
    zoomerrdvv(mov_stack=1, show=True)
    # plot dv/v pairs
    #zoomerrdvv(mov_stack=1, pairs=['AV.OKBR._AV.OKCE.', 'AV.OKAK._AV.OKCF.']) 
