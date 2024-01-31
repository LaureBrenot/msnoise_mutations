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
conda install -c conda-forge flask-admin flask-wtf markdown folium pymysql logbook tables pywct
conda install -c conda-forge obspy
pip install msnoise
pip install tables

add in .../site-package/msnoise/msnoise_mutations.py
modify in .../site-package/msnoise/scripts/msnoise.py for msnoise_save.py
add in .../site-package/msnoise/scripts/msnoise.py (the new one)
modify in .../site-package/msnoise/scripts/default.py for default_save.py
add in .../site-package/msnoise/scripts/default.py (the new one)

READY TO WORK
This code compute dv/v with one point for 30min
Plot wavelet
"""
import os
from msnoise.msnoise_mutations import zoomerrdvv, parametrization_msnoise, compute_msnoise, plot_wct

if __name__ == "__main__":
    data_location = "C:/Users/laure/OneDrive - Université Libre de Bruxelles/Documents/AThese/Scripts/Scripts_O1/Compile_msnoise_wavelet/Okmok_test"
    data_output = r"C:\Users\laure\OneDrive - Université Libre de Bruxelles\Documents\AThese\Scripts\Scripts_O1\Compile_msnoise_wavelet\test4"
    start_date, end_date ='2021-08-05', '2021-08-15'
    start_ref, end_ref = start_date, end_date

    config_updates = {'data_folder': data_location,
                      'output_folder': os.path.join(data_output, 'CROSS_CORRELATIONS'),
                      'startdate': start_date,
                      'enddate': end_ref,
                      'components_to_compute': 'ZZ',
                      'remove_response': 'N',
                      'keep_all': 'Y',
                      'mov_stack': '1',
                      'ref_begin': start_ref,
                      'ref_end': end_ref,
                      'dtt_minlag': '10',
                      'dtt_maxerr': '0.2',
                      'dtt_maxdt': '0.2',
                      'dtt_mincoh': '0.6'} # 'components_to_compute_single_station': 'EN,EZ,NZ',
    filters = [
        (1, 0.5, 0.55, 5.0, 4.95, 1, 12, 2, True),
        (2, 1.0, 1.05, 20.0, 19.95, 0, 12, 2, True)
    ]
    
    # parametrization of msnoise (first time or when changing parameters)
    parametrization_msnoise(data_output, config_updates, filters)

    # compute cross-correlation, stack, wct, mwcs, dt/t
    compute_msnoise(['WCT']) #['CC', 'STACK', 'WCT', 'MWCS', 'DTT']

    # plot dv/v pairs with wct computation
    event_list = []
    current_config = {'freqmin':'0.5', 'freqmax':'5.0', 'freqranges':[[0.5, 1.0], [1.0, 2.0], [2.0, 4.0], [0.5, 2.0]],
                  'plot_all_period':False,'startdateplot':start_date, 'enddateplot':end_date, 
                  'same_dvv_scale':False, 'dvv_min': -1.2, 'dvv_max':1,
                  'plot_event':True, 'event_list':event_list, 'seasonal_rm':False, 
                  'save_fig':False}
    wanted_plot = ['dvv', 'coh', 'dvv_curve', 'coherence_curve', 'dvv_coh']
    
    plot_wct(current_config, wanted_plot, pairs=['AV.OKBR._AV.OKCE.', 'AV.OKAK._AV.OKCF.'], mov_stack=1, components='ZZ', filterid=1)


    # plot network dv/v with mwcs computation
    #zoomerrdvv(mov_stack=1, show=True)
    # plot dv/v pairs
    #zoomerrdvv(mov_stack=1, pairs=['AV.OKBR._AV.OKCE.', 'AV.OKAK._AV.OKCF.'])