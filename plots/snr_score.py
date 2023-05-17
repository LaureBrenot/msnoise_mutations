"""
    Computes the minimum cross-correlation (CC) stack that meets a specified 
    SNR threshold (Clarke et al. 2011). The SNR threshold is 10% of the reference SNR. 
    Save the score of the processing for each station pair. The score is the 
    average number of CC stacked to reach 10% of the SNR reference.

"""
from ..api import * 

from scipy import signal
from scipy.signal import hilbert
import pandas as pd
from obspy.signal.filter import bandpass
from random import choice
import dask

def smooth(x,window_len=11,window='hanning'):
        """smooth the data using a window with requested size.
            
        This method is based on the convolution of a scaled window with the signal.
        The signal is prepared by introducing reflected copies of the signal 
        (with the window size) in both ends so that transient parts are minimized
         in the begining and end part of the output signal.
         		 By Alec usefulFuncs.py           
        """
        if x.ndim != 1:
            raise(ValueError, "smooth only accepts 1 dimension arrays.")
        if x.size < window_len:
            raise(ValueError, "Input vector needs to be bigger than window size.")
        if window_len<3:
            return(x)

        if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
            raise(ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")

        #print window_len                                         
        s=np.r_[x[int(window_len)-1:0:-1],x,x[-2:-int(window_len)-1:-1]]
        #print(len(s))
        if window == 'flat': #moving average
            w=np.ones(window_len,'d')
        else:
            w=eval('np.'+window+'(window_len)')

        y=np.convolve(w/w.sum(),s,mode='valid')

        #return to original size of x
        yfix = y[(int(window_len/2)-1):-(int(window_len/2))]

        if len(yfix) > len(x):
            yfix = yfix[1:]

        return yfix

def compute_ccf_snr(ccfarray, fs, smooth_win = 10, norm=False):
    """By Alec SNRanalysis.py"""
    #normalize ccfs
    if norm==True:
        #get max for each individual ccf
        ccfs_max = np.max(np.abs(ccfarray.tolist()), axis=1)
        ccfarray = ccfarray / ccfs_max

    #calculate mean of these, equivalent to the N-day stack
    ccf_mean = np.mean(ccfarray, axis=0)
    
    #define taper and taper signal for performing hilbert transform
    taper = signal.tukey(len(ccf_mean), alpha=0.2)
    ccf_mean_tapered = ccf_mean * taper

    #compute amplitude envelope
    analytic_signal = hilbert(ccf_mean_tapered)
    ampenv = np.abs(analytic_signal)
    
    #smooth with hanning window
    ampenv_smoothed = smooth(ampenv, window_len = int(fs*smooth_win))
    if len(ampenv_smoothed) != len(ccf_mean):
        ampenv_smoothed = signal.resample(ampenv_smoothed, len(ccf_mean))

    #compute noise signal, following Clarke et al. (2011):
    ccfarray_squared = np.square(ccfarray.tolist())      
    avg_squared = np.mean(ccfarray_squared,axis=0) #take avg of squares
    squared_avg = np.square(ccf_mean) #square the average

    noise = np.sqrt((avg_squared-squared_avg)/(len(ccfarray)-1))
    
    #smooth noise signal (with same length as earlier smooth)
    noise_smoothed = smooth(noise, window_len = int(fs*smooth_win))
    if len(noise_smoothed) != len(noise):
        noise_smoothed = signal.resample(noise_smoothed, len(noise))

    #compute SNR
    snr = ampenv_smoothed / noise_smoothed
    
    return snr, ampenv_smoothed, noise_smoothed  

@dask.delayed
def iterate_random_min_cc_stack(data, ref_snr_value, cc_sampling_rate):
    """ 
    Computes the minimum cross-correlation (CC) stack that meets a specified 
    SNR threshold. The SNR threshold is 10% of the reference SNR. 
    It randomly select a subset of CCs from the input data and stacking them 
    together. The SNR of the stacked CCs is then computed, and the process is 
    repeated by adding more CCs to the stack until the SNR threshold is met. 
    Finally, the index of the minimum CC stack that meets 10 % of the SNR 
    reference is returned.
    Return : index of the minimum CC stack reaching SNR threshold, the SNR 
    value at that index, the ref SNR value, and the threshold SNR value.
    """
    # TODO : 1) take x random cc to compute SNR, 2) SNR mean per x, 3) then cumul SNR -> more parallelization
    
    ten_ref_snr_value = (10*ref_snr_value)/100 #threshold SNR value as 10% of the reference SNR

    cc_stack_snr =  {} #store the SNR values for each CC stack
    #draw_list = list(range(np.shape(data)[0]))
    random_indices = np.random.choice(range(np.shape(data)[0]), size=np.shape(data)[0], replace=False)
    draw = random_indices[0]
    cumul_cc = data[draw,:]


    #Draw initialization
    #draw = choice(draw_list)
    #cumul_cc = data[draw,:]
    #draw_list.remove(draw)

    #stack cc
    for r in range(1,len(random_indices)):
    #for r in range(1,np.shape(data)[0]):
        draw = random_indices[r]
        #draw = choice(draw_list)
        cumul_cc = np.column_stack((cumul_cc, data[draw,:]))  #stack the CC selected 

        # Compute the mean SNR of the stacked CCs
        cumul_snr_cc, ampenv_smoothed, noise_smoothed = compute_ccf_snr(cumul_cc.T, cc_sampling_rate,10)
        cumul_snr_value = np.mean(cumul_snr_cc)

        cc_stack_snr[r] = cumul_snr_value
        #draw_list.remove(draw) #remove the used index from the list

    #Compute the difference between the ref SNR and the SNR of stacked CCs
    cc_stack_snr_values = pd.DataFrame(cc_stack_snr, index = ['snr']).T
    diff_ref_cumul = ref_snr_value - cc_stack_snr_values
    # Select the rows where the difference is less than the threshold SNR
    select_diff_ref_cumul = diff_ref_cumul[(diff_ref_cumul<ten_ref_snr_value).any(axis=1)]
    
    # Find the index of the minimum CC stack among the 10%
    if len(select_diff_ref_cumul) != 0:
        min_cc_stack = min(select_diff_ref_cumul.index)
        #print('minimum stack with a good snr at', min_cc_stack)
        
    # If empty find the index of the minimum SNR value
    else:
        print('error compute min')
        min_cc_stack= min(diff_ref_cumul.index)
        print( min_cc_stack, 'while target is', ref_snr_value)

    return min_cc_stack, cc_stack_snr_values['snr'][min_cc_stack], ref_snr_value, ten_ref_snr_value



def write_score_to_csv(param, station, score, filename):
    """
    Update the corresponding param and station (row and column) with the score.
    Return : CSV file
    """
    try:
        # Read existing csv into a DataFrame
        df = pd.read_csv(filename, sep='|', index_col=0)
    except FileNotFoundError:
        # If file doesn't exist, create a new DataFrame
        df = pd.DataFrame()
        
    if df.empty:
        df = pd.DataFrame({'parameters': [param], station: [score]})
        
    else:
        if param in df['parameters'].values and station in df.columns:
            # idx = df.loc[df['parameters'] == param].index[0]
            # df.at[idx, station] = score
            df.loc[df['parameters'] == param, station] = score
            
        if param not in df['parameters'].values and station in df.columns:
            new_row = {'parameters': param, station: score}
            df = df.append(pd.DataFrame([new_row]), ignore_index=True)
            
        if param in df['parameters'].values and station not in df.columns:
            df.loc[df['parameters'] == param, station] = score
            
        if param not in df['parameters'].values and station not in df.columns:
            df[station] = np.nan
            new_row = {'parameters': param, station : score}
            df = df.append( pd.DataFrame([new_row]), ignore_index=True)
    
    df.to_csv(filename, sep='|')


def checkcell(filename, param, station, overwrite=False):
    try:
        # Read existing csv into a DataFrame
        df = pd.read_csv(filename, sep='|', index_col=0)
        if param in df['parameters'].values and station in df.columns:
            #print('param and column exist')
            if df.loc[df['parameters'] == param, station].notna().any() :   #not empty
                #print('cell not empty')
                if overwrite == True :
                    #print('overwrite is True')
                    return True
                elif overwrite == False:
                    #print('overwrite is False')
                    return False      
            else:
                return True
        else:
           return True

    except FileNotFoundError:
        print('File not found')
        return True

  
def main(filterid, components, mov_stack=1, overwrite=False, show=True,
         outfile=None, refilter=None, **kwargs):
    """
    Save the score of the processing for each station pair. The score is the 
    average number of CC stacked to reach 10% of the SNR reference.
    """
    db = connect()

    jobtype = 'CC'
    jobs = db.query(Job).filter(Job.jobtype == jobtype)
    
    job_pair = []
    for job in jobs:
        job_pair.append(job.pair) 
    job_set = set(job_pair)
    print(job_set)
    print(len(job_set))
    
    #Loop on each station pair
    for pair in job_set:
        #pair = job.pair
        sta1, sta2 = pair.split(':')
        
        maxlag = float(get_config(db, 'maxlag'))
        cc_sampling_rate = float(get_config(db, 'cc_sampling_rate'))
        start, end, datelist = build_movstack_datelist(db)
        if refilter:
            freqmin, freqmax = refilter.split(':')
            freqmin = float(freqmin)
            freqmax = float(freqmax)

        if sta2 < sta1:
            print("Stations STA1 STA2 should be sorted alphabetically")
            return
    
        sta1 = check_stations_uniqueness(db, sta1)
        sta2 = check_stations_uniqueness(db, sta2)
        
        params = get_params(db)

        ##########
        PARAMS = list(params.items())[:66]
        if refilter:
            PARAMS.append(['freq_band', freqmin+'_'+freqmax])
        PARAMS = ''.join(str(x) for x in PARAMS) 
        
        filename = 'scores_stack_snr.csv'
        STATION = sta1+'_'+sta2
        doit = False
        overwrite = False
        #print(STATION, checkcell(filename, PARAMS, STATION, overwrite))
        if checkcell(filename, PARAMS, STATION, overwrite):
        #if doit:
            #########       
            pair = "%s:%s" % (sta1, sta2)
        
            print("New Data for %s-%s-%i-%i" % (pair, components, filterid,
                                                mov_stack))
            
            #Get MSNoise CCs of the station pair
            nstack, results = get_results(db, sta1, sta2, filterid, components,
                                              datelist, mov_stack, format="matrix")
            #Get MSNoise parameters cf 3msnoise_config.py
            
            if refilter:
                for i, d in enumerate(results):
                    results[i] = bandpass(results[i], freqmin, freqmax, cc_sampling_rate,
                                       zerophase=True)
            
            #Starting stuffs
            data_wn = results
            dates = datelist
    
            #Check for and remove any rows with missing values : cannot perform "compute_ccf_snr" otherwise
            mask = np.isnan(data_wn).any(axis=1)
            data = data_wn[~mask]
            print('nan remove from',np.shape(data_wn), 'to', np.shape(data))
              
            #compute ref SNR
            ref_snr_cc, ampenv_smoothed, noise_smoothed = compute_ccf_snr(data, cc_sampling_rate, 10)
            ref_snr_value = np.mean(ref_snr_cc)
            print('REF SNR VALUE', ref_snr_value)
    
            
            #compute cumul SNR
            delayed_tasks = []    
            repetition = 31
            for c in range(1,repetition):
                #print(c, 'over 50')
                delayed_task = iterate_random_min_cc_stack(data, ref_snr_value, cc_sampling_rate)
                delayed_tasks.append(delayed_task)    
            multi_min_nb_stack = pd.DataFrame(dask.compute(*delayed_tasks))
    
            # Compute the mean of the minimum CC stack reaching the SNR threshold
            stack_processing = np.mean(multi_min_nb_stack[0])
            # TODO : solution based au max freq appearence, mean, median ?
            print('solution', stack_processing, 'over', len(dates), 'cc \n')
        
            #Save the score of the processing
            # TODO: avoid to remove "all_components" of the list of params and dealing with dtype breaking line in csv

            SCORE = (stack_processing*100)/np.shape(data_wn)[0]
            
            write_score_to_csv(PARAMS, STATION, SCORE, filename)

        else:
            continue
