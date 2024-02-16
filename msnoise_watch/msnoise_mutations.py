# -*- coding: utf-8 -*-
"""
Created on Wed Sep 27 15:45:59 2023

@author: laure
Contain get_ref(), zoom_s05compute_mwcs.py, zoom_s06compute_dtt.py, zoomerrdvv.py
"""

from .api import *
#from .api import get_config, get_extension, get_job_types, get_logger, get_params, get_filters, get_results_all, get_results
from .move2obspy import mwcs

from .s000installer import main as s000installer
from .msnoise_table_def import Job

import logbook
import subprocess

#wct
#from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import pyplot, transforms
from scipy.io import loadmat
import pycwt as wavelet
from pycwt.helpers import find
from scipy.signal import convolve2d
import warnings
import os
from obspy.signal.regression import linear_regression

from matplotlib.dates import DateFormatter
import matplotlib as mpl
import matplotlib.pylab as plt
import numpy as np
import pandas as pd
import os
import itertools
import matplotlib.dates as mdates
from matplotlib.lines import Line2D
import matplotlib.cm as cm
from matplotlib.colors import Normalize


def get_extension(export_format):
    if export_format == "BOTH":
        return ".MSEED"
    elif export_format == "SAC":
        return ".SAC"
    elif export_format == "MSEED":
        return ".MSEED"
    else:
        return ".MSEED"
    
def get_ref(session, station1, station2, filterid, components, params=None):
    """
    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`
    :type station1: str
    :param station1: The name of station 1 (formatted NET.STA)
    :type station2: str
    :param station2: The name of station 2 (formatted NET.STA)
    :type filterid: int
    :param filterid: The ID (ref) of the filter
    :type components: str
    :param components: The name of the components used (ZZ, ZR, ...)
    :type params: dict
    :param params: A dictionnary of MSNoise config parameters as returned by
        :func:`get_params`.
    :rtype: :class:`obspy.trace`
    :return: A Trace object containing the ref
    """
    from obspy import Trace, read
    if not params:
        export_format = get_config(session, 'export_format')
        extension = get_extension(export_format)
    else:
        extension = get_extension(params.export_format)

    ref_name = "%s_%s" % (station1, station2)
    ref_name = ref_name
    rf = os.path.join("\STACKS", "%02i" %
                      filterid, "REF", components,
                      ref_name + extension)
    if not os.path.isfile(rf):
        logging.debug("No REF file named %s, skipping." % rf)
        print(os.sep, rf,'file does not exist')
        return Trace()
    else:
        print('file exists and found')
        print(read(rf)[0])
    return read(rf)[0]

##########################################################################################

def process_job_type(session, job_type, compute_command, reset_command='msnoise reset ', info_command='msnoise info -j'):
    iteration = 0
    if job_type == 'STACK':
        while any('T' in job for job in get_job_types(session, job_type)):
            print(f"After {iteration} iterations, there is work to do for {job_type} Ref")
            print(get_job_types(session, job_type))
            if iteration <= 3:
                subprocess.run(compute_command[:-2]+'-r', check=True, shell=True)#{compute_command[:-2]+'-r'}
                iteration += 1
            else:
                print(f"Tried to compute {job_type} 4 times and there are still jobs to do")
                break
        subprocess.run('msnoise info -j', check=True, shell=True)#!{info_command}
        
    while any('I' in job or 'T' in job for job in get_job_types(session, job_type)):
        print(f"After {iteration} iterations, there is work to do for {job_type}")
        print(get_job_types(session, job_type))
        if iteration <= 3:
            subprocess.run('msnoise reset '+job_type, check=True, shell=True)#{reset_command+job_type}
            subprocess.run(compute_command, check=True, shell=True)#!{compute_command}
            iteration += 1
        else:
            print(f"Tried to compute {job_type} 4 times and there are still jobs to do")
            break
    subprocess.run('msnoise info -j', check=True, shell=True)#{info_command}

def parametrization_msnoise(data_output, config_updates, filters):
    try:
        db = connect()
        if not get_stations(db, all=True): 
            setup = True    # new station added
        else: 
            setup = False   # when modifying dates or parameters
    except: 
        setup = True        # first time
    
    if setup:
        # Initialize MSNoise if it's the first time
        os.chdir(data_output) 
        s000installer(tech=1)  # Tech 1=SQLite 2=MySQL 3=PostgreSQL', default=None
    
    db = connect()
    #update_config(db, **config_updates)
    for key, value in config_updates.items():
        update_config(db, key, value)

    if setup:
        # Populate Station Table
        subprocess.run('msnoise populate', check=True, shell=True)
    # Scan the archive
    subprocess.run('msnoise scan_archive --init', check=True, shell=True)

    for params in filters:
        update_filter(db, *params)

    # Check for new jobs
    if setup:
        subprocess.run('msnoise new_jobs --init', check=True, shell=True)
    else:
        subprocess.run('msnoise new_jobs', check=True, shell=True)
    
def compute_msnoise(job_list):
    db = connect()

    # Execute MSNoise jobs
    if 'CC' in job_list:
        process_job_type(db, 'CC', 'msnoise compute_cc')
    if 'STACK' in job_list:
        process_job_type(db, 'STACK', 'msnoise stack -m')
    if 'WCT' in job_list:
        # update WCT jobs as MWCS
        jobs_mwcs = db.query(Job).filter(Job.jobtype == 'MWCS').filter(Job.flag == 'T')
        for job_mwcs in jobs_mwcs:
            update_job(db, job_mwcs.day, job_mwcs.pair, 'WCT', 'T')
        process_job_type(db, 'WCT', 'msnoise compute_wct')
    if 'MWCS' in job_list:
        process_job_type(db, 'MWCS', 'msnoise compute_mwcs')
    if 'DTT' in job_list:
        process_job_type(db, 'DTT', 'msnoise compute_dtt')
    if 'DVV' in job_list:
        process_job_type(db, 'DVV', 'msnoise compute_dvv')

    if 'zoom_WCT' in job_list:
        # update WCT jobs as MWCS
        jobs_mwcs = db.query(Job).filter(Job.jobtype == 'MWCS').filter(Job.flag == 'T')
        for job_mwcs in jobs_mwcs:
            update_job(db, job_mwcs.day, job_mwcs.pair, 'WCT', 'T')
        process_job_type(db, 'WCT', 'msnoise compute_zoom_wct')
    if 'zoom_MWCS' in job_list:
        process_job_type(db, 'MWCS', 'msnoise compute_zoom_mwcs')
    if 'zoom_DTT' in job_list:
        process_job_type(db, 'DTT', 'msnoise compute_zoom_dtt')
    if 'zoom_DVV' in job_list:
        process_job_type(db, 'DVV', 'msnoise compute_zoom_dvv')
##########################################################################################

def build_movstack_timelist(session, mov_stack):
    """
    Creates a time array for the analysis period with a time step of 1 hour.
    The returned tuple contains a start and an end time, and a list of
    individual times between the two.

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`

    :rtype: tuple
    :returns: (start, end, timelist)
    """
    begin = get_config(session, "startdate")
    end = get_config(session, "enddate")
    freq = mov_stack

    begin_datetime = datetime.datetime.strptime(begin, "%Y-%m-%d")
    begin = begin_datetime.strftime('%Y-%m-%d %H:%M:%S')
    end_datetime = datetime.datetime.strptime(end, "%Y-%m-%d")
    end = end_datetime.strftime('%Y-%m-%d 23:59:59')
    print(begin, end, freq)

    if begin[0] == '-':
        start = datetime.datetime.today() + datetime.timedelta(hours=int(begin))
        end = datetime.datetime.today() + datetime.timedelta(hours=int(end))
    elif begin == "1970-01-01 00:00:00": # TODO this fails when the DA is empty
        start = session.query(DataAvailability).order_by(
            DataAvailability.starttime).first().starttime()
        end = datetime.datetime.strptime(end, '%Y-%m-%d %H:%M:%S')
    else:
        start = datetime.datetime.strptime(begin, '%Y-%m-%d %H:%M:%S')
        end = datetime.datetime.strptime(end, '%Y-%m-%d %H:%M:%S')
    end = min(end, datetime.datetime.now())
    timelist = pd.date_range(start, end, freq=freq).tolist()

    return start, end, timelist


def zoom_mwcs(loglevel="INFO"):
    logger = logbook.Logger(__name__)
    # Reconfigure logger to show the pid number in log records
    logger = get_logger('msnoise.compute_mwcs_child', loglevel,
                        with_pid=True)
    logger.info('*** Starting: Compute MWCS ***')    
    db = connect()
    
    export_format = get_config(db, 'export_format')
    if export_format == "BOTH":
        extension = ".MSEED"
    else:
        extension = "."+export_format
    mov_stack = get_config(db, "mov_stack")
    if mov_stack.count(',') == 0:
        mov_stacks = [int(mov_stack), ]
    else:
        mov_stacks = [int(mi) for mi in mov_stack.split(',')]
    
    goal_sampling_rate = float(get_config(db, "cc_sampling_rate"))
    maxlag = float(get_config(db, "maxlag"))
    params = get_params(db)

    startdate =  datetime.datetime.strptime(params['startdate'],"%Y-%m-%d")
    enddate = datetime.datetime.strptime(params['enddate'], "%Y-%m-%d")
    # print(startdate, enddate)
    # First we reset all DTT jobs to "T"odo if the REF is new for a given pair
    # for station1, station2 in get_station_pairs(db, used=True):
    #     sta1 = "%s.%s" % (station1.net, station1.sta)
    #     sta2 = "%s.%s" % (station2.net, station2.sta)
    #     pair = "%s:%s" % (sta1, sta2)
    #     if is_dtt_next_job(db, jobtype='DTT', ref=pair):
    #         logger.info(
    #             "We will recompute all MWCS based on the new REF for %s" % pair)
    #         reset_dtt_jobs(db, pair)
    #         update_job(db, "REF", pair, jobtype='DTT', flag='D')
    # 
    # logger.debug('Ready to compute')
    # Then we compute the jobs
    outfolders = []
    filters = get_filters(db, all=False)
    time.sleep(np.random.random() * 5)
    while is_dtt_next_job(db, flag='T', jobtype='MWCS'):
        #TODO would it be possible to make the next 8 lines in the API ?
        jobs = get_dtt_next_job(db, flag='T', jobtype='MWCS')
        
        if not len(jobs):
            # edge case, should only occur when is_next returns true, but
            # get_next receives no jobs (heavily parallelised calls).
            time.sleep(np.random.random())
            continue
        pair = jobs[0].pair
        refs, days = zip(*[[job.ref, job.day] for job in jobs])
        date_range = pd.date_range(startdate, enddate, freq='d')

        # print(date_range)
        # print(refs)
        logger.info(
            "There are MWCS jobs for some days to recompute for %s" % pair)
        for f in filters:
            filterid = int(f.ref)
            for components in params.all_components:
                ref_name = pair.replace('.', '_').replace(':', '_')
                rf = os.path.join("STACKS", "%02i" %
                                  filterid, "REF", components,
                                  ref_name + extension)
                if not os.path.isfile(rf):
                    logging.debug("No REF file named %s, skipping." % rf)
                    continue
                ref = read(rf)[0].data
                ref_name1 = pair.replace(':', '_')
                #ref_name = pair.replace('.', '_').replace(':', '_') #
                station1, station2 = pair.split(":")

                #ref = get_ref(db, station1.replace('.', '_'), station2.replace('.', '_'), 
                #              filterid, components, params)

                if not len(ref):
                    #print("error ref")
                    logging.debug("No REF file found for %s.%i.%s, skipping." %
                                   (ref_name, filterid, components))
                    continue
                ref = ref.data
                mov_stacks = [1]
                days2 = []
                for day in date_range:
                    day= datetime.datetime.strptime(str(day)[:10], '%Y-%m-%d')
                    days2.append(day)
                days = days2
                #logger.info(days)

                output_folder = get_config(db, 'output_folder')
                path = os.path.join(output_folder, "%02i" % int(filterid),
                                    station1, station2, components)
                #print(path, days[0].strftime('%Y-%m-%d.h5'))
                for mov_stack in mov_stacks:
                    n, data2 = get_results(db, station1, station2, filterid,
                                          components, days, mov_stack,
                                          format="matrix", params=params)
                    # print(np.shape(data2)) 
                    # print(data2.index)
                    # print(data2.columns)
                    #for j, cur2 in enumerate(data2):
                        #print(cur2)
                    data = get_results_all(db, station1, station2, filterid, components, days)

                    data_idx = data.index
                    #print(data_idx)
                    data=data.to_numpy()
                    for i, cur in enumerate(data):
                        if np.all(np.isnan(cur)):
                            continue
                        logger.debug(
                            'Processing MWCS for: %s.%s.%02i - %s - %02i days' %
                            (ref_name1, components, filterid, data_idx[i], mov_stack)) #days[i]
                        output = mwcs(cur, ref, f.mwcs_low, f.mwcs_high, goal_sampling_rate, -maxlag, f.mwcs_wlen, f.mwcs_step)
                        outfolder = os.path.join(
                            'MWCS', "%02i" % filterid, "%03i_DAYS" % mov_stack, components, ref_name1)

                        if outfolder not in outfolders:
                            if not os.path.isdir(outfolder):
                                os.makedirs(outfolder)
                            outfolders.append(outfolder)
                        np.savetxt(os.path.join(outfolder, "%s.txt" % str(data_idx[i])).replace(':', '_'), output)
                        del output
                    clean_scipy_cache()
                    del data, cur
        # THIS SHOULD BE IN THE API
        massive_update_job(db, jobs, "D")
        if not params.hpc:
            for job in jobs:
                update_job(db, job.day, job.pair, 'DTT', 'T')

    logger.info('*** Finished: Compute MWCS ***')
#################################################################################################
    
from obspy.signal.regression import linear_regression


import logbook


def wavg_wstd(data, errors):
    d = data
    errors[errors == 0] = 1e-6
    w = 1. / errors
    wavg = (d * w).sum() / w.sum()
    N = len(np.nonzero(w)[0])
    wstd = np.sqrt(np.sum(w * (d - wavg) ** 2) / ((N - 1) * np.sum(w) / N))
    return wavg, wstd

def get_data_first(df,sta1,sta2, interstations, params, first, logger):
    pair = "%s_%s" % (sta1, sta2)
    
    n1, s1, l1 = sta1.split(".")
    n2, s2, l2 = sta2.split(".")
    dpair = "%s_%s_%s_%s" % (n1, s1, n2, s2)
    dist = interstations[dpair] if dpair in interstations else 0.0
    if dist == 0. and params.dtt_lag == "dynamic":
        logger.debug('%s: Distance is Zero?!' % pair)
    tArray = df.index.values
    if params.dtt_lag == "static":
        lmlag = -params.dtt_minlag
        rmlag = params.dtt_minlag
    else:
        lmlag = -dist / params.dtt_v
        rmlag = dist / params.dtt_v
    lMlag = lmlag - params.dtt_width
    rMlag = rmlag + params.dtt_width
    
    if params.dtt_sides == "both":
        tindex = np.where(((tArray >= lMlag) & (tArray <= lmlag)) | ((tArray >= rmlag) & (tArray <= rMlag)))[0]
    elif params.dtt_sides == "left":
        tindex = np.where((tArray >= lMlag) & (tArray <= lmlag))[0]
    else:
        tindex = np.where((tArray >= rmlag) & (tArray <= rMlag))[0]
    
    tmp = np.setdiff1d(np.arange(len(tArray)),tindex)
    df.iloc[tmp, df.columns.get_indexer(['err', ])] = 1.0
    df.iloc[tmp, df.columns.get_indexer(['coh', ])] = 0.0
    

    tArray = df.index.values
    dtArray = df['dt']
    errArray = df['err']
    cohArray = df['coh']
    pairArray = [pair, ]
    first = False

    del df
    return tArray, dtArray, errArray, cohArray, pairArray, first

def get_data_second(df,sta1,sta2, interstations, params,logger, dtArray, errArray, cohArray, pairArray):
    pair = "%s_%s" % (sta1, sta2)
    
    n1, s1, l1 = sta1.split(".")
    n2, s2, l2 = sta2.split(".")
    dpair = "%s_%s_%s_%s" % (n1, s1, n2, s2)
    dist = interstations[dpair] if dpair in interstations else 0.0
    if dist == 0. and params.dtt_lag == "dynamic":
        logger.debug('%s: Distance is Zero?!' % pair)
    tArray = df.index.values
    if params.dtt_lag == "static":
        lmlag = -params.dtt_minlag
        rmlag = params.dtt_minlag
    else:
        lmlag = -dist / params.dtt_v
        rmlag = dist / params.dtt_v
    lMlag = lmlag - params.dtt_width
    rMlag = rmlag + params.dtt_width
    
    if params.dtt_sides == "both":
        tindex = np.where(((tArray >= lMlag) & (tArray <= lmlag)) | ((tArray >= rmlag) & (tArray <= rMlag)))[0]
    elif params.dtt_sides == "left":
        tindex = np.where((tArray >= lMlag) & (tArray <= lmlag))[0]
    else:
        tindex = np.where((tArray >= rmlag) & (tArray <= rMlag))[0]
    
    tmp = np.setdiff1d(np.arange(len(tArray)),tindex)
    df.iloc[tmp, df.columns.get_indexer(['err', ])] = 1.0
    df.iloc[tmp, df.columns.get_indexer(['coh', ])] = 0.0
    
    
    dtArray = np.vstack((dtArray, df['dt']))
    errArray = np.vstack((errArray, df['err']))
    cohArray = np.vstack((cohArray, df['coh']))
    pairArray.append(pair)
    del df
    return tArray, dtArray, errArray, cohArray, pairArray

def compute_dtt(tArray, dtArray, errArray, cohArray, pairArray, params, current, logger):
    Dates = []
    Pairs = []
    M = []
    EM = []
    A = []
    EA = []
    M0 = []
    EM0 = []
    if len(pairArray) != 1:
        # first stack all pairs to a ALL mean pair, using
        # indexes of selected values:
        new_dtArray = np.zeros(len(tArray))
        new_errArray = np.zeros(len(tArray)) + 9999
        new_cohArray = np.zeros(len(tArray))
        for i in range(len(tArray)):
            #~ if i in tindex:
            if 1:
                cohindex = np.where(
                    cohArray[:, i] >= params.dtt_mincoh)[0]
                errindex = np.where(
                    errArray[:, i] <= params.dtt_maxerr)[0]
                dtindex = np.where(
                    np.abs(dtArray[:, i]) <= params.dtt_maxdt)[0]

                index = np.intersect1d(cohindex, errindex)
                index = np.intersect1d(index, dtindex)

                wavg, wstd = wavg_wstd(
                    dtArray[:, i][index],
                    errArray[:, i][index])
                new_dtArray[i] = wavg
                new_errArray[i] = wstd
                new_cohArray[i] = 1.0

        dtArray = np.vstack((dtArray, new_dtArray))
        errArray = np.vstack((errArray, new_errArray))
        cohArray = np.vstack((cohArray, new_cohArray))
        pairArray.append("ALL")
        del new_cohArray, new_dtArray, new_errArray,\
            cohindex, errindex, dtindex, wavg, wstd
        
        # then stack selected pais to GROUPS:
        groups = {}
        npairs = len(pairArray)-1
        for group in groups.keys():
            new_dtArray = np.zeros(len(tArray))
            new_errArray = np.zeros(len(tArray)) + 9999
            new_cohArray = np.zeros(len(tArray))
            pairindex = []
            for j, pair in enumerate(pairArray[:npairs]):
                net1, sta1, net2, sta2 = pair.split('_')
                if sta1 in groups[group] and \
                                sta2 in groups[group]:
                    pairindex.append(j)
            pairindex = np.array(pairindex)

            for i in range(len(tArray)):
                #~ if i in tindex:
                if 1:
                    cohindex = np.where(
                        cohArray[:, i] >= params.dtt_mincoh)[0]
                    errindex = np.where(
                        errArray[:, i] <= params.dtt_maxerr)[0]
                    dtindex = np.where(
                        np.abs(dtArray[:, i]) <= params.dtt_maxdt)[0]

                    index = np.intersect1d(cohindex,
                                           errindex)
                    index = np.intersect1d(index, dtindex)
                    index = np.intersect1d(index, pairindex)
                    

                    wavg, wstd = wavg_wstd(
                        dtArray[:, i][index],
                        errArray[:, i][index])
                    new_dtArray[i] = wavg
                    new_errArray[i] = wstd
                    new_cohArray[i] = 1.0

            dtArray = np.vstack((dtArray, new_dtArray))
            errArray = np.vstack((errArray, new_errArray))
            cohArray = np.vstack((cohArray, new_cohArray))
            pairArray.append(group)
            del new_cohArray, new_dtArray, new_errArray,\
                cohindex, errindex, dtindex, wavg, wstd
            # END OF GROUP HANDLING

    # then process all pairs + the ALL
    if len(dtArray.shape) == 1:  # if there is only one pair:
        dtArray = dtArray.values.reshape((1, dtArray.shape[0]))
        cohArray = cohArray.values.reshape((1, cohArray.shape[0]))
        errArray = errArray.values.reshape((1, errArray.shape[0]))

    used = np.zeros(dtArray.shape)

    for i, pair in enumerate(pairArray):
        cohindex = np.where(cohArray[i] >= params.dtt_mincoh)[0]
        errindex = np.where(errArray[i] <= params.dtt_maxerr)[0]
        dtindex = np.where(np.abs(dtArray[i]) <= params.dtt_maxdt)[0]

        #~ index = np.intersect1d(tindex, cohindex)
        index = np.intersect1d(cohindex, errindex)
        index = np.intersect1d(index, dtindex)

        used[i][index] = 1.0

        w = 1.0 / errArray[i][index]
        w[~np.isfinite(w)] = 1.0
        VecXfilt = tArray[index]
        VecYfilt = dtArray[i][index]
        if len(VecYfilt) >= 2:
            m, a, em, ea = linear_regression(
                VecXfilt, VecYfilt, w,
                intercept_origin=False)
            m0, em0 = linear_regression(
                VecXfilt, VecYfilt, w,
                intercept_origin=True)
            M.append(m)
            EM.append(em)
            A.append(a)
            EA.append(ea)

            M0.append(m0)
            EM0.append(em0)

            Dates.append(current)
            Pairs.append(pair)

            del m, a, em, ea, m0, em0

        del VecXfilt, VecYfilt, w
        del index, cohindex, errindex, dtindex

    logger.debug(
        "%s: exporting: %i pairs" % (current,
                                      len(pairArray)))
    df = pd.DataFrame(
        {'Pairs': Pairs, 'M': M, 'EM': EM, 'A': A, 'EA': EA,
         'M0': M0, 'EM0': EM0},
        index=Dates)
    del M, EM, A, EA, M0, EM0, Pairs, Dates, used
    del tArray, dtArray, errArray, cohArray, pairArray
    return df


def save_dtt(df, current, filterid, components, mov_stack, logger):
    # TODO Needs to be changed to save via the API
    output = os.path.join(
        'DTT', "%02i" % filterid, "%03i_DAYS" % mov_stack,
        components)
    if not os.path.isdir(output):
        os.makedirs(output)
    fn = os.path.join(output, '%s.txt' % current)
    if os.path.isfile(fn):
        existing = pd.read_csv(fn, index_col="Pairs", parse_dates=True)
        for id, row in df.iterrows():
            if row.Pairs in existing.index.values:
                existing.drop(row.Pairs, inplace=True)
                logger.debug("Pair: %s is already in the output file, overwriting" % row.Pairs)
        existing["Pairs"] = existing.index.values
        existing.set_index("Date", inplace=True)
        df = pd.concat([df, existing])
    df.to_csv(fn, index_label='Date')
    del df, output
    
def get_zoom_mwcs(session, station1, station2, filterid, components, date,
                mov_stack=1):
    """
    TODO
    """
    file = os.path.join('MWCS', "%02i" % filterid, "%03i_DAYS" % mov_stack,
                        components, "%s_%s" % (station1, station2),
                        '%s.txt' % str(date).replace(':','_'))
    print(date)
    print(file)
    if os.path.isfile(file):
        df = pd.read_csv(
            file, delimiter=' ', header=None, index_col=0,
            names=['t', 'dt', 'err', 'coh'])
        print(df)
        return df
    else:
        return pd.DataFrame()
    
def zoom_dtt(interval=1, loglevel="INFO"):
    logger = logbook.Logger(__name__)
    # Reconfigure logger to show the pid number in log records
    logger = get_logger('msnoise.compute_dtt_child', loglevel,
                        with_pid=True)
    logger.info('*** Starting: Compute DT/T ***')
    db = connect()
    params = get_params(db)

    start, end, datelist = build_movstack_datelist(db)
    print(start, end)
    daylist = [date.strftime("%Y-%m-%d") for date in datelist]

    mov_stack = get_config(db, "mov_stack")
    if mov_stack.count(',') == 0:
        mov_stacks = [int(mov_stack), ]
    else:
        mov_stacks = [int(mi) for mi in mov_stack.split(',')]

    components_to_compute = get_components_to_compute(db)
    updated_dtt = updated_days_for_dates(
        db, start, end, '%', jobtype='DTT', returndays=True,
        interval=datetime.timedelta(days=interval))
    interstations = {}
    for sta1, sta2 in get_station_pairs(db):
        s1 = "%s_%s" % (sta1.net, sta1.sta)
        s2 = "%s_%s" % (sta2.net, sta2.sta)
        if s1 == s2:
            interstations["%s_%s" % (s1, s2)] = 0.0
        else:
            interstations["%s_%s"%(s1,s2)] = get_interstation_distance(sta1,
                                                                       sta2,
                                                                       sta1.coordinates)
        
    filters = get_filters(db, all=False)
    while is_next_job(db, jobtype='DTT'):
        jobs = get_next_job(db, jobtype='DTT')
        if jobs[0].day in daylist:
            current = jobs[0].day
            dday= datetime.datetime.strptime(str(jobs[0].day)[:10], '%Y-%m-%d')
            print('DDAY', dday)
            ### FIND BETTER WAY
            netsta1, netsta2 = jobs[0].pair.split(':')
            filterid = int(filters[0].ref)
            components = params.all_components[0]
            data_date = get_results_all(db, netsta1, netsta2, filterid, components, [dday]) 
            hours = data_date.index
            
            ###
            for current in hours:
                stations = []
                pairs = []
                refs = []
                for f in filters:
                    filterid = int(f.ref)
                    for components in params.all_components:
                        for mov_stack in mov_stacks:
                            logger.info('Loading mov=%i days for filter=%02i' %
                                          (mov_stack, filterid))
                            first = True
                            for job in jobs: #pair change
                                refs.append(job.ref)
                                pairs.append(job.pair)
                                netsta1, netsta2 = job.pair.split(':')
                                stations.append(netsta1)
                                stations.append(netsta2)
                                sta1 = netsta1+'.'
                                sta2 = netsta2+'.'
                                df = get_zoom_mwcs(db, netsta1, netsta2, filterid, components, current, mov_stack)
                                if not len(df):
                                    continue
                                if first:
                                    tArrayh, dtArrayh, errArrayh, cohArrayh, pairArrayh, first = get_data_first(
                                        df, sta1, sta2, interstations, params, first, logger)
                                else:
                                    tArrayh, dtArrayh, errArrayh, cohArrayh, pairArrayh = get_data_second(
                                        df, sta1, sta2, interstations, params,logger, dtArrayh, errArrayh, cohArrayh, pairArrayh)  
                                                              
                            if not first:
                                df = compute_dtt(tArrayh, dtArrayh, errArrayh, cohArrayh, pairArrayh, params, current, logger)
                                save_dtt(df, str(current).replace(':','_'), filterid, components, mov_stack, logger)
                        
        # THIS SHOULD BE IN THE API
        massive_update_job(db, jobs, "D")

    logger.info('*** Finished: Compute DT/T ***')


#if __name__ == "__zoom_dtt__":
#    zoom_dtt()


#######################################################################################
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt

from matplotlib.dates import DateFormatter


def wavg(group, dttname, errname):
    d = group[dttname]
    group[errname][group[errname] == 0] = 1e-6
    w = 1. / group[errname]
    wavg = (d * w).sum() / w.sum()
    return wavg


def wstd(group, dttname, errname):
    d = group[dttname]
    group[errname][group[errname] == 0] = 1e-6
    w = 1. / group[errname]
    wavg = (d * w).sum() / w.sum()
    N = len(np.nonzero(w)[0])
    wstd = np.sqrt(np.sum(w * (d - wavg) ** 2) / ((N - 1) * np.sum(w) / N))
    return wstd


def get_wavgwstd(data, dttname, errname):
    grouped = data.groupby(level=0)
    g = grouped.apply(wavg, dttname=dttname, errname=errname)
    h = grouped.apply(wstd, dttname=dttname, errname=errname)
    return g, h


def zoomerrdvv(mov_stack=None, dttname="M", components='ZZ', filterid=1,
         pairs=[], showALL=False, show=False, outfile=None):
    db = connect()

    start, end, datelist = build_movstack_datelist(db)
  
    start = datetime.datetime.combine(start, datetime.time())
    end = datetime.datetime.combine(end, datetime.time())

    if mov_stack != 0:
        mov_stacks = [mov_stack, ]
    else:
        mov_stack = get_config(db, "mov_stack")
        if mov_stack.count(',') == 0:
            mov_stacks = [int(mov_stack), ]
        else:
            mov_stacks = [int(mi) for mi in mov_stack.split(',')]

    if components.count(","):
        components = components.split(",")
    else:
        components = [components, ]

    low = high = 0.0
    for filterdb in get_filters(db, all=True):
        if filterid == filterdb.ref:
            low = float(filterdb.low)
            high = float(filterdb.high)
            break

    gs = gridspec.GridSpec(len(mov_stacks), 1)
    fig = plt.figure(figsize=(12, 9))
    plt.subplots_adjust(bottom=0.06, hspace=0.3)
    first_plot = True
    for i, mov_stack in enumerate(mov_stacks):
        current = start
        first = True
        alldf = []
        while current <= end:
            for comp in components:
                day = os.path.join('DTT', "%02i" % filterid, "%03i_DAYS" %
                                   mov_stack, comp, '%s.txt' % str(current).replace(':','_'))
                if os.path.isfile(day):
                    df = pd.read_csv(day, header=0, index_col=0,
                                     parse_dates=True)
                    alldf.append(df)
            current += datetime.timedelta(minutes=30)
        if len(alldf) == 0:
            print("No Data for %s m%i f%i" % (components, mov_stack, filterid))
            continue

        alldf = pd.concat(alldf)
        print(mov_stack, alldf.head())
        if 'alldf' in locals():
            errname = "E" + dttname
            alldf.to_csv("tt.csv")
            alldf[dttname] *= -100
            alldf[errname] *= -100

            ALL = alldf[alldf['Pairs'] == 'ALL'].copy()
            allbut = alldf[alldf['Pairs'] != 'ALL'].copy()

            if first_plot == 1:
                ax = plt.subplot(gs[i])
            else:
                plt.subplot(gs[i], sharex=ax)
            # x = {}
            # for group in groups.keys():
            #     pairindex = []
            #     for j, pair in enumerate(allbut['Pairs']):
            #         net1, sta1, net2, sta2 = pair.split('_')
            #         if sta1 in groups[group] and sta2 in groups[group]:
            #             pairindex.append(j)
            #     tmp = allbut.iloc[np.array(pairindex)]
            #     tmp = tmp.resample('D', how='mean')
            #     #~ plt.plot(tmp.index, tmp[dttname], label=group)
            #     x[group] = tmp
            #
            # tmp = x["CRATER"] - x["VOLCAN"]
            # plt.plot(tmp.index, tmp[dttname], label="Crater - Volcan")

            for pair in pairs:
                print(pair)
                pair1 = alldf[alldf['Pairs'] == pair].copy()
                print(pair1.head())
                plt.plot(pair1.index, pair1[dttname], label=pair)
                plt.fill_between(pair1.index, pair1[dttname]-pair1[errname],
                                 pair1[dttname]+pair1[errname], zorder=-1,
                                 alpha=0.5)
                pair1.to_csv('%s-m%i-f%i.csv'%(pair, mov_stack, filterid))

            if showALL:
                plt.plot(ALL.index, ALL[dttname], c='r',
                         label='ALL: $\delta v/v$ of the mean network')
            
            tmp2 = allbut[dttname].resample('30T').median()
            
            plt.plot(tmp2.index, tmp2, label="median")
            # tmp2.plot(label='mean',)
            etmp2 = allbut[errname].resample('30T').median()
            plt.fill_between(tmp2.index, tmp2-etmp2,
                                 tmp2+etmp2, zorder=-1,
                                 alpha=0.5, label='error')
            #tmp3 = allbut[dttname].resample('D').median()
            # tmp3.plot(label='median')
            #plt.plot(tmp3.index, tmp3, label="median")
            plt.ylabel('dv/v (%)')
            if first_plot == 1:
                plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=4,
                           ncol=2, borderaxespad=0.)
                left, right = tmp2.index[0], tmp2.index[-1]
                if mov_stack == 1:
                    plt.title('No Moving Window')
                else:
                    plt.title('%i Days Moving Window' % mov_stack)
                first_plot = False
            else:
                plt.xlim(left, right)
                plt.title('%i Days Moving Window' % mov_stack)

            plt.grid(True)
            plt.gca().xaxis.set_major_formatter(DateFormatter("%Y-%m-%d %H:%M"))
            fig.autofmt_xdate()
            title = '%s, Filter %d (%.2f - %.2f Hz)' % \
                    (",".join(components), filterid, low, high)
            plt.suptitle(title)
            del alldf
    if outfile:
        if outfile.startswith("?"):
            if len(mov_stacks) == 1:
                outfile = outfile.replace('?', '%s-f%i-m%i-M%s' % (components,
                                                                   filterid,
                                                                   mov_stack,
                                                                   dttname))
            else:
                outfile = outfile.replace('?', '%s-f%i-M%s' % (components,
                                                               filterid,
                                                               dttname))
        outfile = "dvv " + outfile
        print("output to:", outfile)
        plt.savefig(outfile)
    if show:
        plt.show()


#if __name__ == "__zoomerrdvv__":
#    zoomerrdvv()
        


#################################################################################################################
################   ##########      ##########   ####             ####                  ##########################
#################   ########   ##   ########   #####   #####################   ##################################
##################   ######   ####   ######   ######   #####################   ##################################
###################   ####   ######   ####   #######   #####################   ##################################
####################   ##   ########   ##   ########   #####################   ##################################
#####################      ##########      #########             ###########   ##################################
#################################################################################################################
        #########   WCT main function : compute_wct, rest from wxs_dvv of Alec Yates   #########
#################################################################################################################
#################################################################################################################

## Disable Warnings
warnings.filterwarnings('ignore')
## conv2 function
# Returns the two-dimensional convolution of matrices x and y
def conv2(x, y, mode='same'):
    return np.rot90(convolve2d(np.rot90(x, 2), np.rot90(y, 2), mode=mode), 2)

## nextpow2 function
# Returns the exponents p for the smallest powers of two that satisfy the relation  : 2**p >= abs(x)
def nextpow2(x):
    res = np.ceil(np.log2(x))
    return res.astype('int')

## Smoothing function
# Smooth the dataset
def smoothCFS(cfs, scales, dt, ns, nt):
    """
    Smoothing function
    """
    N = cfs.shape[1]
    npad = int(2 ** nextpow2(N))
    omega = np.arange(1, np.fix(npad / 2) + 1, 1).tolist()
    omega = np.array(omega) * ((2 * np.pi) / npad)
    omega_save = -omega[int(np.fix((npad - 1) / 2)) - 1:0:-1]
    omega_2 = np.concatenate((0., omega), axis=None)
    omega_2 = np.concatenate((omega_2, omega_save), axis=None)
    omega = np.concatenate((omega_2, -omega[0]), axis=None)
    # Normalize scales by DT because we are not including DT in the angular frequencies here.
    # The smoothing is done by multiplication in the Fourier domain.
    normscales = scales / dt

    for kk in range(0, cfs.shape[0]):
        
        F = np.exp(-nt * (normscales[kk] ** 2) * omega ** 2)
        smooth = np.fft.ifft(F * np.fft.fft(cfs[kk - 1], npad))
        cfs[kk - 1] = smooth[0:N]
    # Convolve the coefficients with a moving average smoothing filter across scales.
    H = 1 / ns * np.ones((ns, 1))

    cfs = conv2(cfs, H)
    return cfs




## xwt function
def xwt(trace_ref, trace_current, fs, ns=3, nt=0.25, vpo=12, freqmin=0.1, freqmax=8.0, nptsfreq=100):
    """
    Wavelet coherence transform (WCT).

    The WCT finds regions in time frequency space where the two time
    series co-vary, but do not necessarily have high power.
    
    Modified from https://github.com/Qhig/cross-wavelet-transform

    Parameters
    ----------
    trace_ref, trace_current : numpy.ndarray, list
        Input signals.
    fs : float
        Sampling frequency.
    ns : smoothing parameter. 
        Default value is 3
    nt : smoothing parameter. 
        Default value is 0.25
    vpo : float,
        Spacing parameter between discrete scales. Default value is 12.
        Higher values will result in better scale resolution, but
        slower calculation and plot.
        
    freqmin : float,
        Smallest frequency
        Default value is 0.1 Hz
    freqmax : float,
        Highest frequency
        Default value is 8.0 Hz
    nptsfreq : int,
        Number of frequency points between freqmin and freqmax.
        Default value is 100 points
    
    ----------        
    TODO.    
    normalize (boolean, optional) :
        If set to true, normalizes CWT by the standard deviation of
        the signals.
    Phase unwrapping

    Returns
    """
    # Choosing a Morlet wavelet with a central frequency w0 = 6
    mother = wavelet.Morlet(6.)
    # nx represent the number of element in the trace_current array
    nx = np.size(trace_current)
    x_reference = np.transpose(trace_ref)
    x_current = np.transpose(trace_current)
    # Sampling interval
    dt = 1 / fs
    # Spacing between discrete scales, the default value is 1/12
    dj = 1 / vpo 
    # Number of scales less one, -1 refers to the default value which is J = (log2(N * dt / so)) / dj.
    J = -1
    # Smallest scale of the wavelet, default value is 2*dt
    s0 = 2 * dt  # Smallest scale of the wavelet, default value is 2*dt

    # Creation of the frequency vector that we will use in the continuous wavelet transform 
    freqlim = np.linspace(freqmax, freqmin, num=nptsfreq, endpoint=True, retstep=False, dtype=None, axis=0)

    # Calculation of the two wavelet transform independently
    # scales are calculated using the wavelet Fourier wavelength
    # fft : Normalized fast Fourier transform of the input trace
    # fftfreqs : Fourier frequencies for the calculated FFT spectrum.
    ###############################################################################################################
    cwt_reference, scales, freqs, coi, fft, fftfreqs = wavelet.cwt(x_reference, dt, dj, s0, J, mother, freqs=freqlim)
    cwt_current, _, _, _, _, _ = wavelet.cwt(x_current, dt, dj, s0, J, mother, freqs=freqlim)
    ###############################################################################################################

    scales = np.array([[kk] for kk in scales])
    invscales = np.kron(np.ones((1, nx)), 1 / scales)
    
    cfs2 = smoothCFS(invscales * abs(cwt_current) ** 2, scales, dt, ns, nt)
    cfs1 = smoothCFS(invscales * abs(cwt_reference) ** 2, scales, dt, ns, nt)
    
    crossCFS = cwt_reference * np.conj(cwt_current)
    WXamp = abs(crossCFS)
    # cross-wavelet transform operation with smoothing
    crossCFS = smoothCFS(invscales * crossCFS, scales, dt, ns, nt)
    WXspec = crossCFS / (np.sqrt(cfs1) * np.sqrt(cfs2))
    WXangle = np.angle(WXspec)
    Wcoh = abs(crossCFS) ** 2 / (cfs1 * cfs2)
    pp = 2 * np.pi * freqs
    pp2 = np.array([[kk] for kk in pp])
    WXdt = WXangle / np.kron(np.ones((1, nx)), pp2)


    return WXamp, WXspec, WXangle, Wcoh, WXdt, freqs, coi



## dv/v measurement function
def get_dvv(freqs, tvec, WXamp, Wcoh, delta_t, lag_min=5, coda_cycles=20, mincoh=0.5, maxdt=0.2, min_nonzero=0.25, freqmin=0.1, freqmax=2.0):
    """
    Measure velocity variations (dv/v) from the Wavelet coherence transform (WCT).
    
    Parameters
    ----------
    freqs :
    
    tvec :
    
    WXamp :
    
    Wcoh :
    
    delta_t :
    
    lag_min :
    
    coda_cycles :

    mincoh :

    maxdt :

    min_nonzero : % of weights needed to be non-zero to perform regression, otherwise nan (from 0 to 1)
    
    freqmin :
    
    freqmax :
    
    RETURNS:
    ------------------
    dvv*100 : estimated dv/v in %
    err*100 : error of dv/v estimation in %
    wf : weighting function used for the linear regressions
    
    """
    inx = np.where((freqs>=freqmin) & (freqs<=freqmax)) #TODO don't hardcode frequency range
    dvv, err = np.zeros(inx[0].shape), np.zeros(inx[0].shape) # Create empty vectors vor dvv and err
    
    t=tvec
    
    #print ('t',t)
    #print('dvv', dvv)

    ## Better weight function
    weight_func = np.log(np.abs(WXamp))/np.log(np.abs(WXamp)).max()
    zero_idx = (np.where((Wcoh<mincoh) | (delta_t>maxdt))) #TODO get values from db
    wf = weight_func+abs(np.nanmin(weight_func))
    wf = wf/wf.max()
    wf[zero_idx] = 0


    #print ('weight_func',weight_func)
    #print ('zero_idx',zero_idx)
    #print ('wf',wf)

    ## Coda selection
    #tindex = np.where(((t >= -lag_max) & (t <= -lag_min)) | ((t >= lag_min) & (t <= lag_max)))[0] # Index of the coda
    # loop through freq for linear regression
    
    #print('tindex',tindex)

    for ii, ifreq in enumerate(inx[0]): # Loop through frequencies index
     #   print ('tvec:', len(tvec))

        #get coda for specific freq
        period = 1.0/freqs[ifreq]
        lag_max = lag_min + (period*coda_cycles) 

        tindex = np.where(((t >= -lag_max) & (t <= -lag_min)) | ((t >= lag_min) & (t <= lag_max)))[0] # Index of the coda

        if len(tvec)>2: # check time vector size

            if not np.any(delta_t[ifreq]): # check non-empty dt array
                continue
            delta_t[ifreq][tindex]=np.nan_to_num(delta_t[ifreq][tindex])
            w = wf[ifreq] # weighting function for the specific frequency
            w[~np.isfinite(w)] = 1.0
            #m, a, em, ea = linear_regression(time_axis[indx], delta_t[indx], w, intercept_origin=False) # if note forcing through origin
            #m, a, em, ea = linear_regression(tvec[tindex], delta_t[ifreq][tindex], w[tindex], intercept_origin=False)

            #get perc of non-zero weights
            nzc_perc = (np.count_nonzero(w[tindex]))/len(tindex)
            
            #print(nzc)
            
            if nzc_perc >= min_nonzero:

                m, em = linear_regression(tvec[tindex], delta_t[ifreq][tindex], w[tindex], intercept_origin=True) #Forcing through origin
                dvv[ii], err[ii] = -m, em
            else:
                dvv[ii], err[ii] = np.nan, np.nan
                #print (' ' )
            #print ('m',m)

        else:
            print('not enough points to estimate dv/v for wct')
            #dvv[ii], err[ii]=np.nan, np.nan         
    
    return dvv*100, err*100, wf


def get_avgcoh(freqs, tvec, wcoh, freqmin, freqmax, lag_min=5, coda_cycles=20):

    inx = np.where((freqs>=freqmin) & (freqs<=freqmax)) #TODO don't hardcode frequency range
    coh = np.zeros(inx[0].shape) # Create empty vector for coherence
    t = tvec

    ## Coda selection
    #tindex = np.where(((t >= -lag_max) & (t <= -lag_min)) | ((t >= lag_min) & (t <= lag_max)))[0] # Index of the coda

    for ii, ifreq in enumerate(inx[0]): # Loop through frequencies index
     #   print ('tvec:', len(tvec))
     
        period = 1.0/freqs[ifreq]
        lag_max = lag_min + (period*coda_cycles) 
        tindex = np.where(((t >= -lag_max) & (t <= -lag_min)) | ((t >= lag_min) & (t <= lag_max)))[0] # Index of the coda


        if len(tvec)>2: # check time vector size
            if not np.any(wcoh[ifreq]): # check non-empty dt array
                continue
            #wcoh[ifreq][tindex]=np.nan_to_num(delta_t[ifreq][tindex])
           
            c = np.nanmean(wcoh[ifreq][tindex])
            #print(c)
            coh[ii] = c

        else:
            print('not enough points to compute average coherence') #not sure why it would ever get here, but just in case.
            coh[ii] = np.nan

    return coh

def do_plot(time, WXamp, WXspec, WXangle, Wcoh, WXdt, freqs, coi, w, sta, date, comp):   
    save_dir = "WCT/Figure"
    if not os.path.isdir(save_dir):
        os.makedirs(save_dir)
    cmap = "plasma"

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2)

    dt = ax1.pcolormesh(time, freqs, WXdt, cmap="seismic_r", edgecolors='none', vmin=-0.2,vmax=0.2)
    plt.colorbar(dt, ax=ax1)
    ax1.plot(time, 1/coi, 'w--', linewidth=2)
    ax1.set_ylim(freqs[-1], freqs[0])
    ax1.set_xlim(-40,40)
    ax1.set_title('Smoothed Time difference', fontsize=13)
    ax1.set_xlabel('Time (s)')
    ax1.set_ylabel('Frequency (Hz)')

    wc = ax2.pcolormesh(time, freqs, Wcoh, cmap=cmap, edgecolors='none', vmin=0.2, vmax=1)
    plt.colorbar(wc, ax=ax2)
    ax2.plot(time, 1/coi, 'w--', linewidth=2)
    ax2.set_ylim(freqs[-1], freqs[0])
    ax2.set_xlim(-40,40)
    ax2.set_title('Wavelet Coherence', fontsize=13)
    ax2.set_xlabel('Time (s)')
    ax2.set_ylabel('Frequency (Hz)')

    la = ax3.pcolormesh(time, freqs, np.log(WXamp), cmap=cmap, edgecolors='none')
    #plt.clim([-40, 0])
    plt.colorbar(la, ax = ax3)
    ax3.plot(time, 1/coi, 'w--', linewidth=2)
    ax3.set_ylim(freqs[-1], freqs[0])
    ax3.set_xlim(-40,40)
    ax3.set_title('(Logarithmic) Amplitude', fontsize=13)
    ax3.set_xlabel('Time (s)')
    ax3.set_ylabel('Frequency (Hz)')

    weigh = ax4.pcolormesh(time, freqs, w, cmap=cmap, edgecolors='none')
    plt.colorbar(weigh, ax=ax4)
    ax4.plot(time, 1/coi, 'w--', linewidth=2)
    ax4.set_ylim(freqs[-1], freqs[0])
    ax4.set_xlim(-40,40)
    ax4.set_title('Weighting function', fontsize=13)
    ax4.set_xlabel('Time (s)')
    ax4.set_ylabel('Frequency (Hz)')



    plt.tight_layout()
    plt.savefig(os.path.join(save_dir,"{} {}_{}.png".format(sta, comp, date)),dpi=300)
#    plt.show()
    plt.close(fig)
    
    return


def compute_wct_fct(loglevel="INFO"):
    logger = logbook.Logger(__name__)
    # Reconfigure logger to show the pid number in log records
    logger = get_logger('msnoise.compute_wct_child', loglevel,
                        with_pid=True)
    logger.info('*** Starting: Compute WCT ***')
    
    db = connect()

    params = get_params(db)
    export_format = get_config(db, 'export_format')
    if export_format == "BOTH":
        extension = ".MSEED"
    else:
        extension = "."+export_format
    mov_stacks = params.mov_stack

    goal_sampling_rate = float(get_config(db, "cc_sampling_rate"))
    #maxlag = float(get_config(db, "maxlag"))
    ns = float(get_config(db, "wct_ns"))
    nt = float(get_config(db, "wct_nt"))
    vpo = float(get_config(db, "wct_vpo"))
    freqmin_xwt = float(get_config(db, "wct_freqmin"))
    freqmax_xwt = float(get_config(db, "wct_freqmax"))
    nptsfreq = float(get_config(db, "wct_nptsfreq"))
    freqmin_dtt = float(get_config(db, "dtt_freqmin"))
    freqmax_dtt = float(get_config(db, "dtt_freqmax"))
    lag_min = float(get_config(db, "dtt_minlag")) 
    coda_cycles = int(get_config(db, "dtt_codacycles"))
    min_nonzero = float(get_config(db, "dvv_min_nonzero"))
    mincoh = float(get_config(db, "dtt_mincoh"))
    maxdt = float(get_config(db, "dtt_maxdt"))

    params = get_params(db)
    
    logger.debug('Ready to compute')
    # Then we compute the jobs
    outfolders = []
    filters = get_filters(db, all=False)
    time.sleep(np.random.random() * 5)
    
    while is_dtt_next_job(db, flag='T', jobtype='WCT'):
        #TODO would it be possible to make the next 8 lines in the API ?
        jobs = get_dtt_next_job(db, flag='T', jobtype='WCT')
        #print('LEN JOB =', len(jobs))
        if not len(jobs):
            # edge case, should only occur when is_next returns true, but
            # get_next receives no jobs (heavily parallelised calls).
            time.sleep(np.random.random())
            continue
        pair = jobs[0].pair
        refs, days = zip(*[[job.ref, job.day] for job in jobs])

        #logger.info(
        #    "There are WCT jobs for some days to recompute for %s" % pair)
        for f in filters:
            filterid = int(f.ref)

            for components in params.all_components:
                ref_name = pair.replace(':', '_')
                station1, station2 = pair.split(":")
                #ref = get_ref(db, station1, station2, filterid, components, params)
                ref_name = pair.replace('.', '_').replace(':', '_')
                rf = os.path.join("STACKS", "%02i" %
                                  filterid, "REF", components,
                                  ref_name + extension)
                if not os.path.isfile(rf):
                    logging.debug("No REF file named %s, skipping." % rf)
                    continue
                ref = read(rf)[0].data
                if not len(ref):
                    logging.debug("No REF file found for %s.%i.%s, skipping." %
                                  (ref_name, filterid, components))
                    continue
                #ref = ref.data
                
                for mov_stack in mov_stacks:
                    new_datelist = []
                    dvv_list = []
                    err_list = []
                    coh_list = []
                    mov_stack = int(mov_stack)
                    for day in days:
                        df = os.path.join(
                            "STACKS", "%02i" % filterid, "%03i_DAYS" %
                            mov_stack, components, ref_name,
                            str(day) + extension)
                        if os.path.isfile(df):
                            try:
                                cur = read(df)[0].data
                                if np.all(np.isnan(cur)):
                                    print('NaN found in data, skipping')
                                    continue
                            except:
                                logging.debug("Error reading %s, skipping." %
                                              df)
                                continue
                            logger.debug(
                                'Processing WCT for: %s.%s.%02i - %s - %02i days' %
                                (ref_name, components, filterid, day, mov_stack))
                            """"
                            n, data = get_results(db, station1, station2, filterid,
                                                components, days, mov_stack,
                                                format="matrix", params=params)
                            
                            for i, cur in enumerate(data):"""
                            #print(len(cur))
                            if np.all(np.isnan(cur)):
                                print('NaN found in data, skipping')
                                continue

                            #logger.debug(
                            #    'Processing WCT for: %s.%s.%02i - %s - %02i days' %
                            #    (ref_name, components, filterid, day, mov_stack))
                            #output = mwcs(cur, ref, f.mwcs_low, f.mwcs_high, goal_sampling_rate, -maxlag, f.mwcs_wlen, f.mwcs_step)

                            ########### WCT ##########
                            if params.wct_norm:
                                ori_waveform = (ref/ref.max()) 
                                new_waveform = (cur/cur.max())
                            else:
                                ori_waveform = ref
                                new_waveform = cur
                            t = get_t_axis(db)
                        
                            #WXamp, WXspec, WXangle, Wcoh, WXdt, freqs, coi = xwt(ori_waveform, new_waveform, fs, ns, nt, vpo, freqmin_xwt, freqmax_xwt, nptsfreq)
                            WXamp, WXspec, WXangle, Wcoh, WXdt, freqs, coi = xwt(ori_waveform, new_waveform, goal_sampling_rate, int(ns), int(nt), int(vpo), freqmin_xwt, freqmax_xwt, int(nptsfreq))

                            dvv, err, wf = get_dvv(freqs, t, WXamp, Wcoh, WXdt, lag_min=int(lag_min), coda_cycles=coda_cycles, mincoh=mincoh, maxdt=maxdt, min_nonzero=min_nonzero, freqmin=freqmin_dtt, freqmax=freqmax_dtt)
                            coh = get_avgcoh(freqs, t, Wcoh, freqmin_dtt, freqmax_dtt, lag_min=int(lag_min), coda_cycles=coda_cycles)

                            dvv_list.append(dvv)
                            err_list.append(err)
                            coh_list.append(coh)
                            new_datelist.append(day)

                    if len(dvv_list)>1: # Check if the list has more than 1 measurement to save it
                        inx = np.where((freqs>=freqmin_dtt) & (freqs<=freqmax_dtt)) # Select a new frequency range
                        dvv_df = pd.DataFrame(columns=freqs[inx], index=new_datelist)
                        err_df = pd.DataFrame(columns=freqs[inx], index=new_datelist)
                        coh_df = pd.DataFrame(columns=freqs[inx], index=new_datelist)
                        for i, date2 in enumerate(new_datelist): # create the corresponding f_t DataFramei
                            #print(i,date2)
                            dvv_df.iloc[i]=dvv_list[i]
                            err_df.iloc[i]=err_list[i]
                            coh_df.iloc[i]=coh_list[i]


                        outfolder = os.path.join(
                            'WCT', "%02i" % filterid, "%03i_DAYS" % mov_stack, components, ref_name)
                        if outfolder not in outfolders:
                            if not os.path.isdir(outfolder):
                                os.makedirs(outfolder)
                            outfolders.append(outfolder)
                            
                        dfn = "{}_{}_{}-{}.csv".format(pair.replace(":","_"),components,str(dvv_df.index[0]),str(dvv_df.index[-1])) #labeling
                        efn = "{}_{}_{}-{}_error.csv".format(pair.replace(":","_"),components,str(err_df.index[0]),str(err_df.index[-1])) 
                        cfn = "{}_{}_{}-{}_coh.csv".format(pair.replace(":","_"),components,str(coh_df.index[0]),str(coh_df.index[-1])) 

                        pathd = os.path.join(outfolder,dfn)
                        pathe = os.path.join(outfolder,efn)
                        pathc = os.path.join(outfolder,cfn)
                        dvv_df.to_csv(pathd)    # Save dvv to .csv
                        err_df.to_csv(pathe)    #Save err to another csv
                        coh_df.to_csv(pathc)    #Save coh to another csv
                        print(pathd)
                        #np.savetxt(os.path.join(outfolder, "%s.txt" % str(days[i])), output)
                        del dvv_df, err_df, coh_df
                    del cur

        # THIS SHOULD BE IN THE API
        massive_update_job(db, jobs, "D")
        if not params.hpc:
            for job in jobs:
                update_job(db, job.day, job.pair, 'DTT', 'D')
                update_job(db, job.day, job.pair, 'DVV', 'T')

    logger.info('*** Finished: Compute WCT ***')




def read_and_resample(filename):
    dvv_df = pd.read_csv(filename, parse_dates=True, index_col=0)
    dvv_df.columns = dvv_df.columns.astype(float)
    dvv_df = dvv_df.resample("D").mean()
    return dvv_df

def zoom_read_and_resample(filename):
    dvv_df = pd.read_csv(filename, parse_dates=True, index_col=0)
    dvv_df.columns = dvv_df.columns.astype(float)
    #dvv_df = dvv_df.resample("D").mean()
    return dvv_df

def read_and_resample_files(folder_path, files, suffix='_ZZ_'):
    dvv_dfs = []
    for f in files:
        if suffix in f:
            filename = f
            dvv_df = read_and_resample(os.path.join(folder_path, filename))
            dvv_dfs.append(dvv_df)
    return pd.concat(dvv_dfs) if len(dvv_dfs) > 0 else None

def zoom_read_and_resample_files(folder_path, files, suffix='_ZZ_'):
    dvv_dfs = []
    for f in files:
        if suffix in f:
            filename = f
            dvv_df = zoom_read_and_resample(os.path.join(folder_path, filename))
            dvv_dfs.append(dvv_df)
    return pd.concat(dvv_dfs) if len(dvv_dfs) > 0 else None

def create_folder(folder_path):
    try:
        os.makedirs(folder_path)
        print(f"Folder '{folder_path}' created successfully.")
    except FileExistsError:
        pass

def plot_events(ax, event_list, start, end):
    for e in event_list:
        if start < e < end:
            ax.axvline(e, color="k", ls="--", lw=2.2)    

def save_figure(fig, folder_path, filename, plot_all_period=False, start=None, end=None):
    fig_path = os.path.join(folder_path, 'Figures' if plot_all_period else 'Figures/Zooms')
    create_folder(fig_path)

    if start is not None and end is not None:
        filename = f'{filename}_{str(start)[:10]}_{str(end)[:10]}'

    filepath = os.path.join(fig_path, f'{filename}.png')
    fig.savefig(filepath, dpi=300, bbox_inches='tight', transparent=True)

def freq_time_plot(data_type, df_means, stat1, networks, stat2, current_config):
    fig, ax2 = plt.subplots(figsize=(16, 10))
    start = pd.to_datetime(current_config['startdateplot'])
    end = pd.to_datetime(current_config['enddateplot'])
    pair = '/'.join(networks)+"."+stat1+"_"+'/'.join(networks)+"."+stat2
    
    if data_type == 'dvv':
        clim = 0.5
        span = 30
        df1 = df_means.ewm(span=span).mean().astype(float).values
        df1[df_means.isnull()] = np.nan
        cmap = mpl.cm.seismic
        cmap.set_bad('lightgray', 1.)
        vmin = -clim
        cbar_label = 'dv/v (in percent)'
        fig_title = current_config['freqmin']+'-'+current_config['freqmax']+' Hz, dv/v, Average, %s' % pair
        save_title = f"{stat1}_{stat2}_{current_config['freqmin']}_{current_config['freqmax']}_wctmap"
    else:
        clim = 1.0
        span = 1
        df1 = df_means.ewm(span=span).mean().astype(float).values
        cmap='RdYlGn'
        vmin = 0
        cbar_label = 'coherence'
        fig_title = current_config['freqmin']+'-'+current_config['freqmax']+' Hz, coherence, Average, %s'%pair
        save_title = f"{stat1}_{stat2}_{current_config['freqmin']}_{current_config['freqmax']}_coh"
        
    
    ax2.pcolormesh(np.asarray(df_means.ewm(span=span).mean().index),
                   np.asarray(df_means.ewm(span=span).mean().columns),
                   df1.T,
                   cmap=cmap, edgecolors='none', vmin=vmin, vmax=clim)
    
    if current_config['plot_event']:
        plot_events(ax2, current_config['event_list'], start, end)

    #ax2.text(pd.to_datetime("2008-06-05 00:00:00"), 0.8, 'Ölfus doublet', horizontalalignment='center', rotation='vertical')
    #ax2.axvline(pd.to_datetime("2015-07-11 16:50:00"),color="b", ls="--", lw=1)
    
    norm = mpl.colors.Normalize(vmin=vmin, vmax=clim)
    cbar2 = plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax2)
    cbar2.ax.tick_params(labelsize=16,width=2, length=5)
    cbar2.set_label(cbar_label, rotation=270, fontsize=18)
    
    # ax2.set_ylim(freqmin, freqmax)
    ax2.set_title(fig_title, fontsize=18)
    ax2.set_xlim(start, end)
    ax2.set_ylabel('Frequency (Hz)', fontsize=18)
    ax2.tick_params(axis='both',which='both', labelsize=16,width=2, length=5)
    plt.gca().xaxis.set_major_formatter(DateFormatter("%Y-%m-%d"))
    fig.autofmt_xdate()
    
    if current_config['save_fig']:
        save_figure(fig, current_config['folder_path'], save_title, current_config['plot_all_period'], start, end)

    plt.show()
    
def var_time_plot(data_type, df_means, stat1, stat2, current_config):
    if data_type == 'dvv':
        cbar_label = 'dv/v (%)'
        fig_title = stat1+'-'+stat2+' dv/v curve'
        save_title = f"{stat1}_{stat2}_{current_config['freqmin']}_{current_config['freqmax']}_wctcurve"
    else:
        cbar_label = 'Coherence'
        fig_title = stat1+'-'+stat2+' coherence curve'
        save_title =f"{stat1}_{stat2}_{current_config['freqmin']}_{current_config['freqmax']}_cohcurve"
        
    start = pd.to_datetime(current_config['startdateplot'])
    end = pd.to_datetime(current_config['enddateplot'])
    
    smooth=30
    freqs = np.asarray(df_means.columns)
    df_freqslice_dvv = pd.DataFrame()
    fig, ax = plt.subplots(figsize=(12,8))
    
    for freqrange in current_config['freqranges']:
        
        freqmin_dvv=freqrange[0]
        freqmax_dvv=freqrange[1]
    
        #note freqs go from high to low
        minidx = (np.abs(freqs-freqmax_dvv)).argmin()
        maxidx = (np.abs(freqs-freqmin_dvv)).argmin()
    
        df_freqslice = df_means.iloc[:, minidx:maxidx].copy()
        df_freqslice['avg'] = df_freqslice.mean(axis=1).rolling(smooth,min_periods=10).mean()
        y = df_freqslice['avg']
        #print(np.nanmean(y[3000:]) , np.nanmin(y[4000:]), '% difference')
        """
        if current_config['seasonal_rm']:
            import statsmodels.api as sm
            signal = df_freqslice['avg'].values
            signal = pd.Series(signal).fillna(method='ffill').fillna(method='bfill').values # fill nan values with the previous ones
            res = sm.tsa.seasonal_decompose(signal, period=365) # detect seasonal component
            signal_no_seasonal = signal - res.seasonal # remove seasonal component
            signal_no_seasonal[df_freqslice['avg'].isnull()] = np.nan # put nan values back
            ax.plot(df_freqslice.index, signal_no_seasonal, label=str(freqrange[0])+'-'+str(freqrange[1])+' Hz')
        else:"""
        ax.plot(df_freqslice.index, df_freqslice['avg'].values, label=str(freqrange[0])+'-'+str(freqrange[1])+' Hz')
            
        df_freqslice_dvv[str(freqrange[0])+'-'+str(freqrange[1])+' Hz'] = df_freqslice['avg'].values
        
    if current_config['plot_event']:
        plot_events(ax, current_config['event_list'], start, end)
            
    ax.set_xlim(start,end)
    if current_config['same_dvv_scale']:
        ax.set_ylim(current_config['dvv_min'],current_config['dvv_max'])
    ax.set_ylabel(cbar_label,fontsize=18)
    #ax.set_title(freqmin+'-'+freqmax+' Hz, dv/v, Average, %s'%pair, fontsize=18)
    ax.set_title(fig_title, fontsize=22)
    ax.tick_params(axis='both',which='both', labelsize=16,width=2, length=5)
    ax.legend(fontsize=18)
    
    ax.xaxis.set_minor_locator(mdates.MonthLocator())
    plt.gca().xaxis.set_major_formatter(DateFormatter("%Y-%m"))
    fig.autofmt_xdate()
    
    if current_config['save_fig']:    
        save_figure(fig, current_config['folder_path'], save_title, current_config['plot_all_period'], start, end)

    plt.show()
    
    return df_freqslice.index, df_freqslice_dvv
    
def dvv_coh_plot(index, dvv_freqslice, coh_freqslice, stat1, stat2, current_config):
    curve_index = index
    start = pd.to_datetime(current_config['startdateplot'])
    end = pd.to_datetime(current_config['enddateplot'])
            
    color = ['Blues', 'Reds','Greens', 'Purples','Greys']
    color2 = ['blue', 'red', 'green', 'purple']
    freq_names = []
    legend_handles = []
    
    fig, ax = plt.subplots(figsize=(16,10))
    
    for i, freqrange in enumerate(current_config['freqranges']):
        freq_name = f"{freqrange[0]}-{freqrange[1]} Hz"
        dvv_curve = dvv_freqslice[freq_name]
        coh_curve = coh_freqslice[freq_name]
        
        # Filter data between start and end
        # mask = (curve_index >= start) & (curve_index <= end)
        # curve_index_filtered = curve_index[mask]
        # dvv_curve_filtered = dvv_curve[mask]
        # coh_curve_filtered = coh_curve[mask]
            
        freq_names.append(freq_name)
        ax.scatter([0,1], [0,1], c=[0,1], cmap=color[-1])
        norm1 = plt.Normalize(vmin=0, vmax=1)
        # sc1 = ax.scatter(curve_index_filtered, dvv_curve_filtered, c=coh_curve_filtered, cmap=color[i], norm=norm1, label=freq_name)

        sc1 = ax.scatter(curve_index, dvv_curve, c=coh_curve, cmap=color[i], norm=norm1, label=freq_name)
        # Add custom legend handle with the specified color
        legend_handles.append(Line2D([0], [0], marker='o', color=color2[i], markerfacecolor=color2[i], markersize=10, label=freq_name))
        
    if current_config['plot_event']:
        plot_events(ax, current_config['event_list'], start, end)
                
    ax.set_xlim(start,end)
    if current_config['same_dvv_scale']:
        ax.set_ylim(current_config['dvv_min'],current_config['dvv_max'])
    ax.set_ylabel('dv/v (%)',fontsize=18)
    #ax.set_title(freqmin+'-'+freqmax+' Hz, dv/v, Average, %s'%pair, fontsize=18)
    ax.set_title(f"{stat1}-{stat2} dvv and their coherence", fontsize=22)
    
    legend1 = ax.legend(handles=legend_handles, fontsize=22, loc='upper left')
    ax.add_artist(legend1)
    
    cbar1 = plt.colorbar(ax.collections[0], ax=ax, pad=0.02)
    cbar1.set_label('Coherence value', fontsize=18)
    norm1 = Normalize(vmin=0, vmax=1)
    #sm1 = plt.cm.ScalarMappable(cmap=color[-1])#, norm=norm1)
    #sm1.set_array([])
    # new_ticks = [0, 1.0]  # Adjust these values as needed
    # cbar1.set_ticks(new_ticks)
    
    ax.tick_params(axis='both',which='both', labelsize=16,width=2, length=5)
    ax.xaxis.set_minor_locator(mdates.MonthLocator())
    plt.gca().xaxis.set_major_formatter(DateFormatter("%Y-%m"))
    fig.autofmt_xdate()
    
    # Adjust the layout to make room for the colorbar legends
    fig.subplots_adjust(right=0.85)

    if current_config['save_fig']:    
        save_figure(fig,  current_config['folder_path'], f"{stat1}_{stat2}_{current_config['freqmin']}_{current_config['freqmax']}_wctcohcurve",  current_config['plot_all_period'], start, end)

    plt.show()

def get_dvv_and_coh(data_type, stat1, stat2, folder_path, comps):
    all_files = os.listdir(folder_path)
    networks = [n.split('.')[0] for n in all_files if stat1 in n and stat2 in n and '0.csv' in n]
    if data_type == 'dvv':
        dvv_files = [f for f in all_files if stat1 in f and stat2 in f and '0.csv' in f]
    if data_type == 'coh':
        dvv_files = [f for f in all_files if stat1 in f and stat2 in f and '0_coh.csv' in f]

    dvv_df_zz = read_and_resample_files(folder_path, dvv_files, '_ZZ_')
    dvv_df_nz = read_and_resample_files(folder_path, dvv_files, '_NZ_')
    dvv_df_en = read_and_resample_files(folder_path, dvv_files, '_EN_')
    dvv_df_ez = read_and_resample_files(folder_path, dvv_files, '_EZ_')
    #print(len(comps))
    # Average of NZ,EZ,EN for all frequencies
    if len(comps) <= 2:
        dvv_concat = dvv_df_zz
    else:
        dvv_concat = pd.concat((dvv_df_nz, dvv_df_ez, dvv_df_en), join='inner')

    try:
        dvv_by_row_index = dvv_concat.groupby(dvv_concat.index)
        # Average of cross-components-correct
        dvv_means = dvv_by_row_index.mean().astype(float)
        return dvv_means, networks

    except Exception as e:
        print('{data_type}_concat empty, skipping ', stat1, stat2)
        print(dvv_files)
        print(e)
    
    

def plot_wct(config_wct, wanted_plot, mov_stack=None, components='ZZ', filterid=1, pairs=[], showALL=False, show=False, outfile=None):
    db = connect()
    
    start, end, datelist = build_movstack_datelist(db)
    
    start = datetime.datetime.combine(start, datetime.time())
    end = datetime.datetime.combine(end, datetime.time())
    
    if config_wct['startdateplot'] is not None:
        start = datetime.datetime.strptime(config_wct['startdateplot'], "%Y-%m-%d")
    if config_wct['enddateplot'] is not None:
        end = datetime.datetime.strptime(config_wct['enddateplot'], "%Y-%m-%d")

    if mov_stack != 0:
        mov_stacks = [mov_stack, ]
    else:
        mov_stack = get_config(db, "mov_stack")
        if mov_stack.count(',') == 0:
            mov_stacks = [int(mov_stack), ]
        else:
            mov_stacks = [int(mi) for mi in mov_stack.split(',')]
    
    if components.count(","):
        components = components.split(",")
    else:
        components = [components, ]

    for i, mov_stack in enumerate(mov_stacks):
               
        for pair in pairs:
            stat1 = pair.replace('.', '_').split("_")[1]
            stat2 = pair.replace('.', '_').split("_")[4]
            ref_name = pair.replace('_', '').replace('.', '_').replace(':', '_')[:-1]
            for comp in components:
                folder_path = os.path.join('WCT', "%02i" % filterid, "%03i_DAYS" % mov_stack, comp, ref_name)

                config_wct = {**config_wct, 'folder_path': folder_path}
                

                # DATA
                dvv_means, networks = get_dvv_and_coh('dvv', stat1, stat2, folder_path, comp)
                coh_means, n = get_dvv_and_coh('coh', stat1, stat2, folder_path, comp)    
                
                # Plot: Freq, time, average dv/v  #####################################################  
                if 'dvv' in wanted_plot:
                    freq_time_plot('dvv', dvv_means, stat1, networks, stat2, config_wct)
                if 'coh' in wanted_plot:
                    freq_time_plot('coh', coh_means, stat1, networks, stat2, config_wct)   
                
                # Plot: dv/v time series within frequency band  #########################################
                if 'dvv_curve' in wanted_plot or 'dvv_coh' in wanted_plot:
                    index, dvv_freqslice = var_time_plot('dvv', dvv_means, stat1, stat2, config_wct)
                if 'coherence_curve' in wanted_plot or 'dvv_coh' in wanted_plot:
                    index, coh_freqslice = var_time_plot('coh', coh_means, stat1, stat2, config_wct)
                if 'dvv_coh' in wanted_plot:
                    dvv_coh_plot(index, dvv_freqslice, coh_freqslice, stat1, stat2, config_wct)




def compute_zoom_wct_fct(loglevel="INFO"):
    logger = logbook.Logger(__name__)
    # Reconfigure logger to show the pid number in log records
    logger = get_logger('msnoise.compute_wct_child', loglevel,
                        with_pid=True)
    logger.info('*** Starting: Compute WCT ***')
    
    db = connect()

    params = get_params(db)
    export_format = get_config(db, 'export_format')
    if export_format == "BOTH":
        extension = ".MSEED"
    else:
        extension = "."+export_format
    mov_stacks = params.mov_stack

    goal_sampling_rate = float(get_config(db, "cc_sampling_rate"))
    #maxlag = float(get_config(db, "maxlag"))
    ns = float(get_config(db, "wct_ns"))
    nt = float(get_config(db, "wct_nt"))
    vpo = float(get_config(db, "wct_vpo"))
    freqmin_xwt = float(get_config(db, "wct_freqmin"))
    freqmax_xwt = float(get_config(db, "wct_freqmax"))
    nptsfreq = float(get_config(db, "wct_nptsfreq"))
    freqmin_dtt = float(get_config(db, "dtt_freqmin"))
    freqmax_dtt = float(get_config(db, "dtt_freqmax"))
    lag_min = float(get_config(db, "dtt_minlag")) 
    coda_cycles = int(get_config(db, "dtt_codacycles"))
    min_nonzero = float(get_config(db, "dvv_min_nonzero"))
    mincoh = float(get_config(db, "dtt_mincoh"))
    maxdt = float(get_config(db, "dtt_maxdt"))

    params = get_params(db)
    
    logger.debug('Ready to compute')
    # Then we compute the jobs
    outfolders = []
    filters = get_filters(db, all=False)
    time.sleep(np.random.random() * 5)
    
    while is_dtt_next_job(db, flag='T', jobtype='WCT'):
        #TODO would it be possible to make the next 8 lines in the API ?
        jobs = get_dtt_next_job(db, flag='T', jobtype='WCT')
        #print('LEN JOB =', len(jobs))
        if not len(jobs):
            # edge case, should only occur when is_next returns true, but
            # get_next receives no jobs (heavily parallelised calls).
            time.sleep(np.random.random())
            continue
        pair = jobs[0].pair
        refs, days = zip(*[[job.ref, job.day] for job in jobs])

        startdate =  datetime.datetime.strptime(params['startdate'],"%Y-%m-%d")
        enddate = datetime.datetime.strptime(params['enddate'], "%Y-%m-%d")
        date_range = pd.date_range(startdate, enddate, freq='d')

        #logger.info(
        #    "There are WCT jobs for some days to recompute for %s" % pair)
        for f in filters:
            filterid = int(f.ref)

            for components in params.all_components:
                ref_name = pair.replace(':', '_')
                station1, station2 = pair.split(":")
                #ref = get_ref(db, station1, station2, filterid, components, params)
                ref_name = pair.replace('.', '_').replace(':', '_')
                rf = os.path.join("STACKS", "%02i" %
                                  filterid, "REF", components,
                                  ref_name + extension)
                if not os.path.isfile(rf):
                    logging.debug("No REF file named %s, skipping." % rf)
                    continue
                ref = read(rf)[0].data
                if not len(ref):
                    logging.debug("No REF file found for %s.%i.%s, skipping." %
                                  (ref_name, filterid, components))
                    continue
                #ref = ref.data
                days2 = []
                for day in date_range:
                    day= datetime.datetime.strptime(str(day)[:10], '%Y-%m-%d')
                    days2.append(day)
                days = days2

                for mov_stack in mov_stacks:
                    new_datelist = []
                    dvv_list = []
                    err_list = []
                    coh_list = []
                    mov_stack = int(mov_stack)


                    data = get_results_all(db, station1, station2, filterid, components, days)
                    data_idx = data.index
                    #print(data_idx)
                    data=data.to_numpy()
                    for i, cur in enumerate(data):
                        if np.all(np.isnan(cur)):
                            print('NaN found in data, skipping')
                            continue
                        logger.debug(
                            'Processing MWCS for: %s.%s.%02i - %s - %02i days' %
                            (ref_name, components, filterid, data_idx[i], mov_stack))
                        
                        day = data_idx[i]
                        ########### WCT ##########
                        if params.wct_norm:
                            ori_waveform = (ref/ref.max()) 
                            new_waveform = (cur/cur.max())
                        else:
                            ori_waveform = ref
                            new_waveform = cur
                        t = get_t_axis(db)
                    
                        #WXamp, WXspec, WXangle, Wcoh, WXdt, freqs, coi = xwt(ori_waveform, new_waveform, fs, ns, nt, vpo, freqmin_xwt, freqmax_xwt, nptsfreq)
                        WXamp, WXspec, WXangle, Wcoh, WXdt, freqs, coi = xwt(ori_waveform, new_waveform, goal_sampling_rate, int(ns), int(nt), int(vpo), freqmin_xwt, freqmax_xwt, int(nptsfreq))

                        dvv, err, wf = get_dvv(freqs, t, WXamp, Wcoh, WXdt, lag_min=int(lag_min), coda_cycles=coda_cycles, mincoh=mincoh, maxdt=maxdt, min_nonzero=min_nonzero, freqmin=freqmin_dtt, freqmax=freqmax_dtt)
                        coh = get_avgcoh(freqs, t, Wcoh, freqmin_dtt, freqmax_dtt, lag_min=int(lag_min), coda_cycles=coda_cycles)

                        dvv_list.append(dvv)
                        err_list.append(err)
                        coh_list.append(coh)
                        new_datelist.append(day)

                    if len(dvv_list)>1: # Check if the list has more than 1 measurement to save it
                        inx = np.where((freqs>=freqmin_dtt) & (freqs<=freqmax_dtt)) # Select a new frequency range
                        dvv_df = pd.DataFrame(columns=freqs[inx], index=new_datelist)
                        err_df = pd.DataFrame(columns=freqs[inx], index=new_datelist)
                        coh_df = pd.DataFrame(columns=freqs[inx], index=new_datelist)
                        for i, date2 in enumerate(new_datelist): # create the corresponding f_t DataFramei
                            #print(i,date2)
                            dvv_df.iloc[i]=dvv_list[i]
                            err_df.iloc[i]=err_list[i]
                            coh_df.iloc[i]=coh_list[i]


                        outfolder = os.path.join(
                            'WCT', "%02i" % filterid, "%03i_DAYS" % mov_stack, components, ref_name)
                        if outfolder not in outfolders:
                            if not os.path.isdir(outfolder):
                                os.makedirs(outfolder)
                            outfolders.append(outfolder)
                            
                        dfn = "{}_{}_{}-{}.csv".format(pair.replace(":","_"),components,str(dvv_df.index[0]).replace(":","_"),str(dvv_df.index[-1]).replace(":","_")) #labeling
                        efn = "{}_{}_{}-{}_error.csv".format(pair.replace(":","_"),components,str(err_df.index[0]).replace(":","_"),str(err_df.index[-1]).replace(":","_")) 
                        cfn = "{}_{}_{}-{}_coh.csv".format(pair.replace(":","_"),components,str(coh_df.index[0]).replace(":","_"),str(coh_df.index[-1]).replace(":","_")) 

                        pathd = os.path.join(outfolder,dfn)
                        pathe = os.path.join(outfolder,efn)
                        pathc = os.path.join(outfolder,cfn)
                        dvv_df.to_csv(pathd)    # Save dvv to .csv
                        err_df.to_csv(pathe)    #Save err to another csv
                        coh_df.to_csv(pathc)    #Save coh to another csv
                        print(pathd)
                        #np.savetxt(os.path.join(outfolder, "%s.txt" % str(days[i])), output)
                        del dvv_df, err_df, coh_df
                    del cur

        # THIS SHOULD BE IN THE API
        massive_update_job(db, jobs, "D")
        if not params.hpc:
            for job in jobs:
                update_job(db, job.day, job.pair, 'DTT', 'D')
                update_job(db, job.day, job.pair, 'DVV', 'T')

    logger.info('*** Finished: Compute WCT ***')

def get_zoom_dvv_and_coh(data_type, stat1, stat2, folder_path, comps):
    all_files = os.listdir(folder_path)
    networks = [n.split('.')[0] for n in all_files if stat1 in n and stat2 in n and '0.csv' in n]
    if data_type == 'dvv':
        dvv_files = [f for f in all_files if stat1 in f and stat2 in f and '0.csv' in f]
    if data_type == 'coh':
        dvv_files = [f for f in all_files if stat1 in f and stat2 in f and '0_coh.csv' in f]

    dvv_df_zz = zoom_read_and_resample_files(folder_path, dvv_files, '_ZZ_')
    dvv_df_nz = zoom_read_and_resample_files(folder_path, dvv_files, '_NZ_')
    dvv_df_en = zoom_read_and_resample_files(folder_path, dvv_files, '_EN_')
    dvv_df_ez = zoom_read_and_resample_files(folder_path, dvv_files, '_EZ_')

    # Average of NZ,EZ,EN for all frequencies
    if len(comps) <= 2:
        dvv_concat = dvv_df_zz
    else:
        dvv_concat = pd.concat((dvv_df_nz, dvv_df_ez, dvv_df_en), join='inner')

    try:
        dvv_by_row_index = dvv_concat.groupby(dvv_concat.index)
        # Average of cross-components-correct
        dvv_means = dvv_by_row_index.mean().astype(float)
        return dvv_means, networks

    except Exception as e:
        print('{data_type}_concat empty, skipping ', stat1, stat2)
        print(dvv_files)
        print(e)

def plot_zoom_wct(config_wct, wanted_plot, mov_stack=None, components='ZZ', filterid=1, pairs=[], showALL=False, show=False, outfile=None):
    db = connect()
    
    start, end, datelist = build_movstack_datelist(db)
    
    start = datetime.datetime.combine(start, datetime.time())
    end = datetime.datetime.combine(end, datetime.time())
    
    if config_wct['startdateplot'] is not None:
        start = datetime.datetime.strptime(config_wct['startdateplot'], "%Y-%m-%d")
    if config_wct['enddateplot'] is not None:
        end = datetime.datetime.strptime(config_wct['enddateplot'], "%Y-%m-%d")

    if mov_stack != 0:
        mov_stacks = [mov_stack, ]
    else:
        mov_stack = get_config(db, "mov_stack")
        if mov_stack.count(',') == 0:
            mov_stacks = [int(mov_stack), ]
        else:
            mov_stacks = [int(mi) for mi in mov_stack.split(',')]
    
    if components.count(","):
        components = components.split(",")
    else:
        components = [components, ]

    for i, mov_stack in enumerate(mov_stacks):
               
        for pair in pairs:
            stat1 = pair.replace('.', '_').split("_")[1]
            stat2 = pair.replace('.', '_').split("_")[4]
            ref_name = pair.replace('_', '').replace('.', '_').replace(':', '_')[:-1]
            for comp in components:
                folder_path = os.path.join('WCT', "%02i" % filterid, "%03i_DAYS" % mov_stack, comp, ref_name)

                config_wct = {**config_wct, 'folder_path': folder_path}
                

                # DATA
                dvv_means, networks = get_zoom_dvv_and_coh('dvv', stat1, stat2, folder_path, comp)
                coh_means, n = get_zoom_dvv_and_coh('coh', stat1, stat2, folder_path, comp)    
                
                # Plot: Freq, time, average dv/v  #####################################################  
                if 'dvv' in wanted_plot:
                    freq_time_plot('dvv', dvv_means, stat1, networks, stat2, config_wct)
                if 'coh' in wanted_plot:
                    freq_time_plot('coh', coh_means, stat1, networks, stat2, config_wct)   
                # Plot: dv/v time series within frequency band  #########################################
                if 'dvv_curve' in wanted_plot or 'dvv_coh' in wanted_plot:
                    index, dvv_freqslice = var_time_plot('dvv', dvv_means, stat1, stat2, config_wct)
                if 'coherence_curve' in wanted_plot or 'dvv_coh' in wanted_plot:
                    index, coh_freqslice = var_time_plot('coh', coh_means, stat1, stat2, config_wct)
                if 'dvv_coh' in wanted_plot:
                    dvv_coh_plot(index, dvv_freqslice, coh_freqslice, stat1, stat2, config_wct)