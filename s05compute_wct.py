"""
Compute dt/t using wavelet transform based on 
Mao. el al 2020 theory and matlab code; 
Higueret python code https://github.com/Qhig
and Alec Yates python codes https://github.com/asyates

This code is modified from s05compute_mwcs.py 
with A. Yates code for MSNoise workflow https://github.com/ROBelgium/MSNoise
by L. Brenot https://github.com/LaureBrenot/msnoise_mutations
 
    Parallel Processing 
"""
# TODO check if still // while no more output per day

from .api import *
from .wxs_dvv import *

import logbook

from sqlalchemy.orm.exc import NoResultFound
import os
import pydoc
import pandas as pd
from obspy.core.util import AttribDict
from .msnoise_table_def import Config

def get_defaults(default_file_name):
    df = pd.read_csv(os.path.join(os.path.dirname(os.path.realpath(__file__)), default_file_name),
                     delimiter=",", encoding="latin-1", index_col=0)
    df["type"] = [pydoc.locate(t) for t in df["type"].values]
    df = df.fillna("")

    default = AttribDict()
    for id, row in df.iterrows():
        if len(row.possible_values):
            df.loc[id, "definition"] += " " + row.possible_values.replace(row.default, "[%s]"%row.default)
        # elif len(row.default):
        #     df.loc[id].definition += " " + "[%s]" % row.default
        default[id] = AttribDict(row.to_dict())
    return default

def get_params(session):
    """Get config parameters from the database.

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`

    :returns: a Param class containing the parameters
    """
    # TODO: this could be populated automatically from defauts iff defaults
    # would mention types
    from obspy.core.util.attribdict import AttribDict
    #from .default import default
    s = session
    params = AttribDict()
    default = get_defaults('default_mod.csv')
    for name in default.keys():
        itemtype = default[name].type
        if itemtype is bool:
            params[name] = get_config(s, name, isbool=True)
        else:
            params[name] = itemtype(get_config(s, name))

    # TODO remove reference to goal_sampling_rate
    params.goal_sampling_rate = params.cc_sampling_rate
    params.min30 = params.corr_duration * params.goal_sampling_rate
    params.components_to_compute = get_components_to_compute(s)
    params.components_to_compute_single_station = get_components_to_compute_single_station(s)
    params.all_components = np.unique(params.components_to_compute_single_station + \
                            params.components_to_compute)

    if params.mov_stack.count(',') == 0:
        params.mov_stack = [int(params.mov_stack), ]
    else:
        params.mov_stack = [int(mi) for mi in params.mov_stack.split(',')]

    return params

def add_config_entry(session, name, value):
    """
    Add a new configuration entry to the database.

    :param session: SQLAlchemy session
    :param name: Name of the configuration field to add
    :param value: Value for the configuration field
    """
    try:
        # Check if the configuration entry already exists
        existing_entry = session.query(Config).filter_by(name=name).one()
        #print(f"Configuration field '{name}' already exists in the database.")
    except NoResultFound:
        # If the entry doesn't exist, add a new row to the Config table
        new_entry = Config(name=name, value=value)
        session.add(new_entry)
        session.commit()
        #print(f"New configuration field '{name}' added successfully.")

def wavelet_default_update(session, new_default):
    defaultn = get_defaults(new_default)
    #ori_default = get_defaults('default.csv')
    i = 0
    for name in defaultn.keys():
        i += 1
        if i > 64:
            #ri_default[name]=defaultn[name]
            if defaultn[name].default:
                itemtype = defaultn[name].type
                add_config_entry(session, name, itemtype(defaultn[name].default))
                #ori_default[name]=defaultn[name]
            else:
                add_config_entry(session, name,"")
                #ori_default[name]=defaultn[name]
    #df_updated = pd.DataFrame(list(ori_default.values()), index=ori_default.keys())

    #file_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "default.csv")
    #df_original = pd.read_csv(file_path, delimiter=",", encoding="latin-1", index_col=0)
    #df_original.update(df_updated)
    #df_original.to_csv(file_path, sep=',', encoding='latin-1')
    #if not os.path.exists(os.path.join(os.path.dirname(os.path.realpath(__file__)),'default_save.csv')):
    default_path = os.path.join(os.path.dirname(os.path.realpath(__file__)),'default.csv')
    new_default = os.path.join(os.path.dirname(os.path.realpath(__file__)),'default_save.csv')
    mod_path =  os.path.join(os.path.dirname(os.path.realpath(__file__)),'default_mod.csv')
    new_mode = os.path.join(os.path.dirname(os.path.realpath(__file__)),'default.csv')
    os.rename(default_path,new_default)
    os.rename(mod_path, new_mode)

def main(loglevel="INFO"):
    logger = logbook.Logger(__name__)
    # Reconfigure logger to show the pid number in log records
    logger = get_logger('msnoise.compute_wct_child', loglevel,
                        with_pid=True)
    logger.info('*** Starting: Compute WCT ***')
    
    db = connect()
    #if not os.path.exists(os.path.join(os.path.dirname(os.path.realpath(__file__)),'default_save.csv')): 
    #    wavelet_default_update(db, 'default_mod.csv')


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
        print('LEN JOB =', len(jobs))
        if not len(jobs):
            # edge case, should only occur when is_next returns true, but
            # get_next receives no jobs (heavily parallelised calls).
            time.sleep(np.random.random())
            continue
        pair = jobs[0].pair
        refs, days = zip(*[[job.ref, job.day] for job in jobs])

        logger.info(
            "There are WCT jobs for some days to recompute for %s" % pair)
        for f in filters:
            filterid = int(f.ref)
            for components in params.all_components:
                ref_name = pair.replace(':', '_')
                station1, station2 = pair.split(":")
                ref = get_ref(db, station1, station2, filterid, components,
                              params)
                if not len(ref):
                    logging.debug("No REF file found for %s.%i.%s, skipping." %
                                  (ref_name, filterid, components))
                    continue
                ref = ref.data

                for mov_stack in mov_stacks:
                    
                    new_datelist = []
                    dvv_list = []
                    err_list = []
                    coh_list = []
                    mov_stack = int(mov_stack)
                    n, data = get_results(db, station1, station2, filterid,
                                          components, days, mov_stack,
                                          format="matrix", params=params)

                    for i, cur in enumerate(data):
                        if np.all(np.isnan(cur)):
                            continue

                        logger.debug(
                            'Processing WCT for: %s.%s.%02i - %s - %02i days' %
                            (ref_name, components, filterid, days[i], mov_stack))
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
                        new_datelist.append(days[i])


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

                    pathd = os.path.join("WCT",dfn)
                    pathe = os.path.join("WCT",efn)
                    pathc = os.path.join("WCT",cfn)
                    dvv_df.to_csv(pathd)    # Save dvv to .csv
                    err_df.to_csv(pathe)    #Save err to another csv
                    coh_df.to_csv(pathc)    #Save coh to another csv
                    print(pathd)
                    print(dfn)
                    #np.savetxt(os.path.join(outfolder, "%s.txt" % str(days[i])), output)
                    del dvv_df, err_df, coh_df
                del data

        # THIS SHOULD BE IN THE API
        massive_update_job(db, jobs, "D")
        if not params.hpc:
            for job in jobs:
                update_job(db, job.day, job.pair, 'DTT', 'D')
                update_job(db, job.day, job.pair, 'DVV', 'T')

    logger.info('*** Finished: Compute WCT ***')
