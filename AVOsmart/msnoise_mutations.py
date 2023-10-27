# -*- coding: utf-8 -*-
"""
Created on Wed Sep 27 15:45:59 2023

@author: laure
Contain get_ref(), zoom_s05compute_mwcs.py, zoom_s06compute_dtt.py, zoomerrdvv.py
"""

from .api import *
#from .api import get_config, get_extension, get_job_types, get_logger, get_params, get_filters, get_results_all, get_results
from .move2obspy import mwcs


import logbook
import subprocess

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
    rf = os.path.join("STACKS", "%02i" %
                      filterid, "REF", components,
                      ref_name + extension)
    if not os.path.isfile(rf):
        logging.debug("No REF file named %s, skipping." % rf)
        return Trace()

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

##########################################################################################
    
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
                ref_name1 = pair.replace(':', '_')
                ref_name = pair.replace('.', '_').replace(':', '_') #
                station1, station2 = pair.split(":")

                ref = get_ref(db, station1.replace('.', '_'), station2.replace('.', '_'), 
                              filterid, components, params)

                if not len(ref):
                    # print("error ref")
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
                # print(days)

                output_folder = get_config(db, 'output_folder')
                path = os.path.join(output_folder, "%02i" % int(filterid),
                                    station1, station2, components)
                print(path, days[0].strftime('%Y-%m-%d.h5'))
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
                    # print(data_idx)
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