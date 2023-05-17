# -*- coding: utf-8 -*-
"""
Created on Sat Dec 17 14:20:21 2022

@author: breno
"""
import matplotlib
from msnoise import api
from msnoise import s000installer
import os
scan 		= True
newjobs 	= False
compute_cc 	= True
stacks 		= True
compute_dtt 	= True

db = api.connect() #connect to database

#api.reset_jobs(db, 'CC', alljobs=True)
#api.reset_jobs(db, 'MWCS', alljobs=True)
#api.reset_jobs(db, 'DTT', alljobs=True)

################################## SCAN ARCHIVE
start 	= '2020-06-22'
end 	= '2020-08-22'
if scan:
    api.update_config(db, 'startdate', start) #set startdate
    api.update_config(db, 'enddate', end) #set enddate
##ready to scan archive

################################# NEWJOBS
#os.system("msnoise db update_loc_chan") # update stations loc and chan from data availability table
#if newjobs:
    #os.system("msnoise db execute 'update stations set used_location_codes='--' '")
    #api.update_config(db, 'used_locations_codes','--')
    #api.update_station(db,'NZ','WIZ',177.1894,-37.5265,40,coordinates='DEG') # enter location of each station
    #api.update_station(db,'NZ','WSRZ',177.178,-37.518,290,coordinates='DEG')

    #api.update_config(db, 'components_to_compute', 'UU') #for station-pair correlation functions
    #api.update_config(db, 'components_to_compute_single_station', 'EN,EZ,NZ')
##ready new_jobs init

################################ COMPUTE CC
if compute_cc:
    api.update_config(db, 'cc_sampling_rate', '20') #set frequency to resample to, default=20
    api.update_config(db, 'preprocess_lowpass', '8') #preprocess lowpass filter frequency, default=8
    api.update_config(db, 'preprocess_highpass', '0.01') #preprocess highpass filter frequency, default=0.01
    api.update_config(db, 'maxlag', '120') #maximum lag of cross-correlation functions, default=120
    api.update_config(db, 'corr_duration', '1800') #length of time slices (seconds) to cross-correlate, default=1800
    api.update_config(db, 'windsorizing', '3') #temporal normalization method (3 = clip 3 times RMS)
    api.update_config(db, 'whitening_type', 'B') #type of whitening (B = whiten to amplitude of 1.0)
    api.update_config(db, 'clip_after_whiten', 'Y') #Boolean to perform temporal norm after spectral whitening
    api.update_config(db, 'stack_method', 'linear') # stack method, linear or pws (phase-weighted)
    api.update_config(db, 'remove_response', 'Y') #remove instrument response
    #api.update_cingig(db, 'response_format', ?)
    api.update_config(db, 'response_path', '/home/ulb/gtime/lbrenot/lscratch/Alaska_volcanoes/Pavlof/data_Pavlof/Response/')
    api.update_config(db, 'keep_all', 'Y')
 
    #update_filter(session, ref, low, mwcs_low, high, mwcs_high, rms_threshold, mwcs_wlen, mwcs_step, used), where low and high are the whitening bounds, mwcs_low and mwcs_high the frequency bounds for moving-window cross-spectral analysis, and mwcs_wlen and mwcs_step define the window length and window step used in this same analysis respectively.
    api.update_filter(db, 1, 0.1, 0.15, 1.0, 0.95, 0, 12, 2, True)
    #filter 0.1 to 1Hz of the whitening function, filter 0.15 to 0.95Hz of the linear regression done in MWCS(looks like 5%), 0:not use anymore, 12sec windows to perform MWCS, 2sec step for windows, True activate the filter
    api.update_filter(db, 2, 1.0, 1.05, 2.0, 1.95, 0, 12, 2, True)
    api.update_filter(db, 3, 1.0, 1.05, 5.0, 4.95, 1, 12, 2, True)
    api.update_filter(db, 4, 5.0, 5.05, 10.0, 9.95, 1, 12, 2, True)
    api.update_filter(db, 5, 0.5, 0.55, 2.0, 1.95, 1, 12, 2, True)

#ready to compute_cc

#################################### STACK
#to record dvv : compare current and reference stacks of ccfs. Parameters to decide of the size stacks
if stacks:
    api.update_config(db, 'mov_stack', '1,5,10') #day length windows
    api.update_config(db, 'ref_begin', start)#define the reference stack
    api.update_config(db, 'ref_end', end)

#os.system("msnoise stack -r") #compute reference stack
#os.system("msnoise reset STACK") #reinitialize the job list to compute
#os.system("msnoise stack -m")#compute the current stack (without the reference)

#ready stack and compute_mwcs

################################### COMPUTE DTT
if compute_dtt:
    api.update_config(db, 'dtt_lag', 'static') #set lag time window to be the same for all CCFs
    api.update_config(db, 'dtt_minlag', '10') #set minimum lag time
    api.update_config(db, 'dtt_width', '30') #set width of window to compute dv/v , in sec
    api.update_config(db, 'dtt_sides', 'both') #choose which side(s) of CCF to compute dv/v
    api.update_config(db, 'dtt_maxerr', '0.2') #set maximum error
    api.update_config(db, 'dtt_maxdt', '0.2') #set maximum delay time
    api.update_config(db, 'dtt_mincoh', '0.6') #set minimum coherence
