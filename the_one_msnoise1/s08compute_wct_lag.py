"""
Wavelet Coherence Transform (WCT) Computation
This script performs the computation of the Wavelet Coherence Transform (WCT), a tool used to analyze the correlation between two time series in the time-frequency domain. The script supports parallel processing and interacts with a database to manage job statuses.

Filter Configuration Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* |wct_minlag|
* |wct_maxdt|
* |wct_mincoh|
* |wct_codacycles|
* |wct_ns|
* |wct_nt|
* |wct_vpo|
* |wct_nptsfreq|
* |wct_min_nonzero|
* |wct_norm|
* |hpc|

This process is job-based, so it is possible to run several instances in
parallel.


To run this step:

.. code-block:: sh

    $ msnoise cc dvv compute_wct

This step also supports parallel processing/threading:

.. code-block:: sh

    $ msnoise -t 4 cc dvv compute_wct

will start 4 instances of the code (after 1 second delay to avoid database
conflicts). This works both with SQLite and MySQL but be aware problems
could occur with SQLite.

"""

import os
import time
import numpy as np
import pandas as pd
import xarray as xr
import scipy.signal
import pycwt as wavelet
from scipy.signal import convolve2d
from obspy.signal.regression import linear_regression
from .api import *
import logbook
import scipy
import scipy.optimize
import scipy.fft as sf 


# ============================================================================
# ADAPTIVE CODA WINDOW DETECTION (integrated with MSNoise API)
# ============================================================================

class CodaDetector:
    """Simple coda window detector with minimum time constraint"""
    
    def __init__(self, dt=0.05, velocity=2.0, safety_factor=1.3, min_coda_time=8.0):
        self.dt = dt
        self.velocity = velocity
        self.safety_factor = safety_factor
        self.min_coda_time = min_coda_time
        
    def detect_coda_window(self, t, ccf, distance_km):
        """Detect coda window using percentage method with minimum time"""
        
        if distance_km is None:
            distance_km = self._estimate_distance(t, ccf)
            if distance_km is None:
                return self.min_coda_time, -self.min_coda_time, {}
        
        # Core calculation
        surface_wave_time = distance_km / self.velocity
        calculated_time = surface_wave_time * self.safety_factor
        coda_start = max(calculated_time, self.min_coda_time)
        
        # Simple diagnostics
        diagnostics = {
            'distance_km': distance_km,
            'surface_wave_time': surface_wave_time,
            'calculated_time': calculated_time,
            'coda_start': coda_start,
            'constraint': 'minimum' if coda_start == self.min_coda_time else 'physics'
        }
        
        return coda_start, -coda_start, diagnostics
    
    def _estimate_distance(self, t, ccf):
        """Estimate distance from first arrival"""
        try:
            from scipy import signal
            envelope = np.abs(signal.hilbert(ccf))
            center_idx = len(t) // 2
            pos_envelope = envelope[center_idx:]
            pos_time = t[center_idx:]
            
            threshold = 0.5 * np.max(pos_envelope)
            peaks, _ = signal.find_peaks(pos_envelope, height=threshold, distance=int(5/self.dt))
            
            if len(peaks) > 0:
                first_peak_time = pos_time[peaks[0]]
                estimated_distance = first_peak_time * self.velocity * 1.1
                if 5 <= estimated_distance <= 500:
                    return estimated_distance
            return None
        except:
            return None

def get_station_distance(db, station1, station2):
    """
    Calculate distance between two stations using MSNoise's get_interstation_distance API.
    
    Parameters:
    -----------
    db : database session
        MSNoise database session
    station1, station2 : str
        Station names in format 'NET.STA' or just 'STA'
        
    Returns:
    --------
    float or None
        Distance in kilometers, or None if stations not found or coordinates not available
    """
    def parse_station_name(station_name):
        """Parse station name to extract network and station codes"""
        parts = station_name.split('.')
        if len(parts) >= 2:
            return parts[0], parts[1]  # network, station
        else:
            return None, parts[0]  # no network, just station
    
    # Parse station names
    net1, sta1 = parse_station_name(station1)
    net2, sta2 = parse_station_name(station2)
    
    logger.debug(f"Looking for distance between {net1}.{sta1} and {net2}.{sta2}")
    
    # Get all stations from database
    stations = get_stations(db, all=True)
    
    # Find the station objects
    station_obj1 = station_obj2 = None
    
    for station in stations:
        # Check for station1
        if ((net1 is None or station.net == net1) and station.sta == sta1):
            if station.X is not None and station.Y is not None:
                station_obj1 = station
                logger.debug(f"Found station object for {station.net}.{station.sta}")
        
        # Check for station2
        if ((net2 is None or station.net == net2) and station.sta == sta2):
            if station.X is not None and station.Y is not None:
                station_obj2 = station
                logger.debug(f"Found station object for {station.net}.{station.sta}")
    
    if station_obj1 is None:
        logger.debug(f"Station object not found for {station1}")
        return None
    if station_obj2 is None:
        logger.debug(f"Station object not found for {station2}")
        return None
    
    # Use MSNoise's built-in function to calculate distance
    try:
        distance = get_interstation_distance(station_obj1, station_obj2, coordinates="DEG")
        logger.debug(f"Calculated distance between {station1} and {station2}: {distance:.2f} km")
        return distance
    except Exception as e:
        logger.error(f"Error calculating distance between {station1} and {station2}: {e}")
        return None

# ============================================================================
# MODIFIED WCT FUNCTIONS WITH MSNoise API INTEGRATION
# ============================================================================

def compute_wct_dvv_adaptive(freqs, tvec, WXamp, Wcoh, delta_t, distance_km=None, 
                           coda_detector=None, ref_waveform=None, 
                           coda_cycles=20, mincoh=0.5, maxdt=0.2, 
                           min_nonzero=0.25, freqmin=0.1, freqmax=2.0):
    """
    Compute the dv/v values with adaptive coda window detection.
    
    Parameters
    ----------
    freqs : numpy.ndarray
        Frequency values corresponding to the wavelet transform.
    tvec : numpy.ndarray
        Time vector.
    WXamp : numpy.ndarray
        Amplitude of the cross-wavelet transform.
    Wcoh : numpy.ndarray
        Wavelet coherence.
    delta_t : numpy.ndarray
        Time delays between signals.
    distance_km : float, optional
        Distance between stations in km. If None, will be estimated.
    coda_detector : CodaDetector, optional
        Coda detector instance. If None, uses default parameters.
    ref_waveform : numpy.ndarray, optional
        Reference waveform for distance estimation if distance_km is None.
    coda_cycles : int, optional
        Number of coda cycles to consider. Default is 20.
    mincoh : float, optional
        Minimum coherence value for weighting. Default is 0.5.
    maxdt : float, optional
        Maximum time delay for weighting. Default is 0.2.
    min_nonzero : float, optional
        Minimum percentage of non-zero weights required for valid estimation. Default is 0.25.
    freqmin : float, optional
        Minimum frequency for calculation. Default is 0.1 Hz.
    freqmax : float, optional
        Maximum frequency for calculation. Default is 2.0 Hz.
    
    Returns
    -------
    tuple
        dvv values (percentage), errors (percentage), weighting function, and coda diagnostics.
    """   
    import warnings
    from scipy.optimize import OptimizeWarning
    
    # Initialize coda detector if not provided
    if coda_detector is None:
        dt = tvec[1] - tvec[0] if len(tvec) > 1 else 0.05
        coda_detector = CodaDetector(dt=dt)
    
    # Determine adaptive lag_min using coda detector
    if distance_km is not None or ref_waveform is not None:
        lag_min, _, coda_diagnostics = coda_detector.detect_coda_window(
            tvec, ref_waveform if ref_waveform is not None else np.zeros_like(tvec), 
            distance_km
        )
        logger.info(f"Adaptive coda window: lag_min = {lag_min:.2f}s, constraint = {coda_diagnostics.get('constraint', 'unknown')}")
    else:
        lag_min = coda_detector.min_coda_time  # Fallback to minimum
        coda_diagnostics = {'constraint': 'fallback', 'coda_start': lag_min}
        logger.warning("No distance or reference waveform provided, using minimum coda time")
    
    inx = np.where((freqs >= freqmin) & (freqs <= freqmax))  # Filter frequencies within the specified range
    dvv, err = np.zeros(len(inx[0])), np.zeros(len(inx[0])) # Initialize dvv and err arrays

    # Weighting function based on WXamp
    weight_func = np.log(np.abs(WXamp)) / np.log(np.abs(WXamp)).max()
    zero_idx = np.where((Wcoh < mincoh) | (delta_t > maxdt))
    wf = (weight_func + abs(np.nanmin(weight_func))) / weight_func.max()
    wf[zero_idx] = 0
    
    # Track frequencies with estimation issues
    problematic_freqs = []

    # Loop through frequency indices for linear regression
    for ii, ifreq in enumerate(inx[0]):
        period = 1.0 / freqs[ifreq]
        lag_max = lag_min + (period * coda_cycles)  # Use adaptive lag_min

        # Coda selection
        tindex = np.where(((tvec >= -lag_max) & (tvec <= -lag_min)) | ((tvec >= lag_min) & (tvec <= lag_max)))[0]

        if len(tvec) > 2:
            if not np.any(delta_t[ifreq]):
                continue

            delta_t[ifreq][tindex] = np.nan_to_num(delta_t[ifreq][tindex])
            w = wf[ifreq] # Weighting function for the specific frequency
            non_zero_mask = w > 0
            nzc_perc = np.count_nonzero(non_zero_mask[tindex]) / len(tindex)
            
            if nzc_perc >= min_nonzero:
                
                with warnings.catch_warnings(record=True) as w_catcher:
                    warnings.simplefilter("always", OptimizeWarning)

                    m, em = linear_regression(tvec[tindex], delta_t[ifreq][tindex], w[tindex], intercept_origin=True)
                    
                    if any(issubclass(warning.category, OptimizeWarning) for warning in w_catcher):
                        problematic_freqs.append(freqs[ifreq])
                        
                dvv[ii], err[ii] = -m, em
            else:
                dvv[ii], err[ii] = np.nan, np.nan
        else:
            logger.debug('Not enough points to estimate dv/v for WCT')
            
    if problematic_freqs:
        logger.warning(f"Covariance estimation issues between {min(problematic_freqs):.2f}-{max(problematic_freqs):.2f} Hz: Modify min_nonzero (current: {min_nonzero}), mincoh (current: {mincoh}), maxdt (current: {maxdt}) and/or coda_cycles (current: {coda_cycles})")
    
    return dvv * 100, err * 100, wf, coda_diagnostics

def get_avgcoh_adaptive(freqs, tvec, wcoh, freqmin, freqmax, distance_km=None, 
                       coda_detector=None, ref_waveform=None, coda_cycles=20):
    """
    Calculate the average wavelet coherence with adaptive coda window.
    
    Parameters similar to compute_wct_dvv_adaptive but focused on coherence calculation.
    """
    # Initialize coda detector if not provided
    if coda_detector is None:
        dt = tvec[1] - tvec[0] if len(tvec) > 1 else 0.05
        coda_detector = CodaDetector(dt=dt)
    
    # Determine adaptive lag_min
    if distance_km is not None or ref_waveform is not None:
        lag_min, _, _ = coda_detector.detect_coda_window(
            tvec, ref_waveform if ref_waveform is not None else np.zeros_like(tvec), 
            distance_km
        )
    else:
        lag_min = coda_detector.min_coda_time
    
    inx = np.where((freqs >= freqmin) & (freqs <= freqmax)) 
    coh = np.zeros(inx[0].shape)

    for ii, ifreq in enumerate(inx[0]):     
        period = 1.0/freqs[ifreq]
        lag_max = lag_min + (period * coda_cycles)  # Use adaptive lag_min
        tindex = np.where(((tvec >= -lag_max) & (tvec <= -lag_min)) | ((tvec >= lag_min) & (tvec <= lag_max)))[0]

        if len(tvec) > 2:
            if not np.any(wcoh[ifreq]) or wcoh[ifreq][tindex].size == 0:
                coh[ii] = np.nan
                continue
            c = np.nanmean(np.abs(wcoh[ifreq][tindex]))
            coh[ii] = c
        else:
            logger.debug('Not enough points to compute average coherence')
            coh[ii] = np.nan

    return coh

def process_wct_job_adaptive(pair, day, params, taxis, filters, db_session=None): 
    """
    Process a single WCT job with adaptive coda window detection using MSNoise API.
    
    Parameters:
    ----------
    pair : str
        Station pair string (station1:station2)
    day : str
        Day string
    params : dict-like
        Configuration parameters
    taxis : numpy.ndarray
        Time axis
    filters : list
        List of filters to use
    db_session : database session, optional
        MSNoise database session for getting station coordinates
        
    Returns:
    -------
    bool
        True if job was successfully processed, False otherwise
    """
    # Parse pair and day
    station1, station2 = pair.split(":")
    
    # Get distance between stations using MSNoise API
    distance_km = None
    if db_session:
        try:
            distance_km = get_station_distance(db_session, station1, station2)
            if distance_km:
                logger.debug(f"Distance between {station1} and {station2}: {distance_km:.2f} km")
            else:
                logger.warning(f"No coordinates found for {station1}-{station2}, will use minimum coda time")
        except Exception as e:
            logger.warning(f"Error calculating distance for {station1}-{station2}: {e}")
    
    # Initialize coda detector with parameters from config
    dt = taxis[1] - taxis[0] if len(taxis) > 1 else 0.05
    
    # Get coda detection parameters from config (with defaults)
    velocity = 2.0 #get_config(db_session, 'coda_velocity', 2.0) if db_session else 2.0
    safety_factor = 1.2 #get_config(db_session, 'coda_safety_factor', 1.3) if db_session else 1.3
    min_coda_time = 6.0 #params.dtt_minlag  # Use existing dtt_minlag parameter
    
    coda_detector = CodaDetector(
        dt=dt, 
        velocity=velocity, 
        safety_factor=safety_factor, 
        min_coda_time=min_coda_time
    )
    
    # Determine components to compute
    if station1 == station2:
        components_to_compute = params.components_to_compute_single_station
    else:
        components_to_compute = params.components_to_compute
        
    logger.info(f"Processing Pair: {pair}, Day: {day} for {len(filters)} filters, {len(components_to_compute)} components and {len(params.mov_stack)} moving windows")
    
    # Extract parameters
    ns = params.wct_ns
    nt = params.wct_nt 
    vpo = params.wct_vpo 
    nptsfreq = params.wct_nptsfreq
    coda_cycles = params.dtt_codacycles 
    min_nonzero = params.dvv_min_nonzero
    wct_norm = params.wct_norm
    wavelet_type = params.wavelet_type
    
    mov_stacks = params.mov_stack
    goal_sampling_rate = params.cc_sampling_rate
    maxdt = params.dtt_maxdt
    mincoh = params.dtt_mincoh
        
    # Convert date string to datetime
    try:
        date = pd.to_datetime(day)
    except:
        logger.error(f"Failed to parse date: {day}")
        return False
    
    operations_succeeded = 0

    # Process each filter
    for f in filters:
        filterid = int(f['ref'])
        freqmin = f['low']
        freqmax = f['high']
        
        # Process each component
        for component in components_to_compute:
            try:
                # Get reference waveform using MSNoise API
                ref = xr_get_ref(station1, station2, component, filterid, taxis, ignore_network=True)
                ref = ref.CCF.values
                if wct_norm:
                    ori_waveform = (ref/ref.max()) 
                else:
                    ori_waveform = ref
            except FileNotFoundError as fullpath:
                logger.error(f"FILE DOES NOT EXIST: {fullpath}, skipping")
                continue
            except Exception as e:
                logger.error(f"Error getting reference waveform: {str(e)}")
                continue
            
            if not len(ref):
                logger.warning(f"Empty reference waveform for {station1}:{station2} {component} f{filterid}")
                continue
            
            # Process each moving stack
            for mov_stack in mov_stacks:
                try:
                    # Get the data for this day and moving stack using MSNoise API
                    data = xr_get_ccf(station1, station2, component, filterid, mov_stack, taxis)
                    
                    # Filter data to get only this specific day
                    day_data = data[data.index == date]
                    
                    if day_data.empty:
                        logger.warning(f"No data for {date} in {station1}:{station2} {component} f{filterid} m{mov_stack}")
                        continue
                    
                    logger.debug(f"Processing {station1}:{station2} {component} f{filterid} m{mov_stack} - {date}")
                    
                    # Get the waveform for this day
                    waveform = day_data.iloc[0].values
                    if wct_norm:
                        new_waveform = (waveform/waveform.max()) 
                    else:
                        new_waveform = waveform
                    
                    # Calculate wavelet coherence transform
                    WXamp, WXspec, WXangle, Wcoh, WXdt, freqs, coi = xwt(
                        ori_waveform, new_waveform, goal_sampling_rate, 
                        int(ns), float(nt), int(vpo), 
                        freqmin, freqmax, int(nptsfreq), wavelet_type, maxdt
                    )
                    
                    # Compute dv/v values with adaptive coda window
                    dvv, err, wf, coda_diagnostics = compute_wct_dvv_adaptive(
                        freqs, taxis, WXamp, Wcoh, WXdt, 
                        distance_km=distance_km,
                        coda_detector=coda_detector,
                        ref_waveform=ori_waveform,
                        coda_cycles=coda_cycles, 
                        mincoh=mincoh, maxdt=maxdt, min_nonzero=min_nonzero, 
                        freqmin=freqmin, freqmax=freqmax
                    )
                    
                    # Calculate average coherence with adaptive window
                    coh = get_avgcoh_adaptive(
                        freqs, taxis, Wcoh, freqmin, freqmax, 
                        distance_km=distance_km,
                        coda_detector=coda_detector,
                        ref_waveform=ori_waveform,
                        coda_cycles=coda_cycles
                    )
                    
                    # Get frequency indices
                    inx = np.where((freqs >= freqmin) & (freqs <= freqmax))
                    
                    # Save results for this day including coda diagnostics
                    save_day_wct_results_adaptive(
                        station1, station2, date, component, 
                        filterid, mov_stack, taxis, 
                        dvv, err, coh, freqs, inx, coda_diagnostics
                    )
                    
                    operations_succeeded += 1
                
                except Exception as e:
                    logger.error(f"Error processing {component} f{filterid} m{mov_stack}: {str(e)}")
                    continue
    
    # Determine overall success status
    if operations_succeeded == 0:
        logger.error(f"Job failed completely: {pair} - {day}.")
        return False
    else:
        logger.info(f"Job processing complete: {pair} - {day}. {operations_succeeded} operations succeeded.")
        return True

def save_day_wct_results_adaptive(station1, station2, date, component, filterid, mov_stack, 
                                 taxis, dvv, err, coh, freqs, inx, coda_diagnostics):
    """
    Save the WCT results including coda diagnostics.
    
    Same as original save_day_wct_results but includes coda window information.
    """
    # Format filterid with leading zeros (e.g., 1 -> "01")
    filter_str = f"{int(filterid):02d}"
    
    # Format mov_stack (could be string like '10d_1d' or int)
    if isinstance(mov_stack, (list, tuple)):
        mov_stack_str = '_'.join(str(m) for m in mov_stack)
    else:
        mov_stack_str = str(mov_stack)
    
    # Station pair string
    pair_str = f"{station1}_{station2}"
    
    # Create directory structure
    dir_path = os.path.join('WCT', filter_str, mov_stack_str, component, pair_str)
    os.makedirs(dir_path, exist_ok=True)
    
    # Filename is just the date: YYYY-MM-DD.npz
    date_str = date.strftime('%Y-%m-%d')
    filename = f"{date_str}.npz"
    filepath = os.path.join(dir_path, filename)
    
    # Create DataFrames
    freqs_subset = freqs[inx]
    
    # Save to numpy compressed format including coda diagnostics
    np.savez_compressed(
        filepath,
        date=date,
        freqs=freqs_subset,
        dvv=dvv,
        err=err,
        coh=coh,
        taxis=taxis,
        station1=station1,
        station2=station2,
        component=component,
        filterid=filterid,
        mov_stack=mov_stack,
        # Add coda diagnostics
        coda_distance_km=coda_diagnostics.get('distance_km'),
        coda_surface_wave_time=coda_diagnostics.get('surface_wave_time'),
        coda_start=coda_diagnostics.get('coda_start'),
        coda_constraint=coda_diagnostics.get('constraint')
    )
    
    logger.debug(f"Saved WCT results to {filepath} with coda constraint: {coda_diagnostics.get('constraint')}")
    return filepath

# ============================================================================
# KEEP ALL OTHER ORIGINAL FUNCTIONS (smoothCFS, conv2, get_wavelet_type, xwt, etc.)
# ============================================================================


def main(loglevel="INFO", batch_size=5):
    """
    Main function to process WCT jobs with adaptive coda window detection using MSNoise API.
    
    Parameters:
    -----------
    loglevel : str
        Logging level
    batch_size : int
        Number of jobs to process in each batch
    """
    global logger
    logger = get_logger('msnoise.compute_wct_child', loglevel, with_pid=True)
    logger.info('*** Starting: Compute WCT with Adaptive Coda Detection ***')
    
    # Add small random delay to prevent all processes starting at exactly the same time
    time.sleep(np.random.random() * 2)
    
    try:
        # Initial database connection to get configuration
        db = connect()
        
        # Ensure no lingering transactions
        if hasattr(db, 'in_transaction') and db.in_transaction():
            db.rollback()
            
        params = get_params(db)
        taxis = get_t_axis(db)
        filters = get_filters(db, all=False)
        
        # Convert filters to dictionaries before closing the session
        filter_dicts = []
        for f in filters:
            filter_dicts.append({
                'ref': f.ref,
                'low': f.low,
                'high': f.high,
                'mwcs_low': getattr(f, 'mwcs_low', None),
                'mwcs_high': getattr(f, 'mwcs_high', None),
                'rms_threshold': getattr(f, 'rms_threshold', None),
                'mwcs_wlen': getattr(f, 'mwcs_wlen', None),
                'mwcs_step': getattr(f, 'mwcs_step', None)
            })
        
        # Add new configuration parameters for adaptive coda detection if they don't exist
#        try:
#            # Try to get existing config values, add defaults if not present
#            if not get_config(db, 'coda_velocity', None):
#                update_config(db, 'coda_velocity', '2.0')
#                logger.info("Added coda_velocity parameter with default value 2.0 km/s")
#            
#            if not get_config(db, 'coda_safety_factor', None):
#                update_config(db, 'coda_safety_factor', '1.3')
#                logger.info("Added coda_safety_factor parameter with default value 1.3")
#                
#            db.commit()
#        except Exception as e:
#            logger.warning(f"Could not add/check coda parameters: {e}")
#            db.rollback()
        
        # Process jobs in batches
        jobs_processed = 0
        max_jobs_per_process = 1000
        consecutive_empty_batches = 0
        max_consecutive_empty = 3
        
        while jobs_processed < max_jobs_per_process:
            # Use a fresh database connection for each batch
            db_batch = connect()
            
            if hasattr(db_batch, 'in_transaction') and db_batch.in_transaction():
                db_batch.rollback()
            
            job_exists = is_dtt_next_job(db_batch, flag='T', jobtype='WCT')
            if not job_exists:
                logger.info("No more WCT jobs to process")
                db_batch.close()
                break
                
            jobs = claim_available_jobs(db_batch, batch_size)
            
            if not jobs:
                consecutive_empty_batches += 1
                if consecutive_empty_batches >= max_consecutive_empty:
                    logger.info(f"Exiting after {consecutive_empty_batches} empty batches")
                    db_batch.close()
                    break
                db_batch.close()
                time.sleep(2)
                continue
            else:
                consecutive_empty_batches = 0
                
            # Process each job in the batch
            for job in jobs:
                pair = job.pair
                day = job.day

                # Process this job with adaptive coda detection, passing db session
                success = process_wct_job_adaptive(pair, day, params, taxis, filter_dicts, db_batch)
                                
                # Update job status
                if success:
                    update_job(db_batch, day, pair, 'WCT', 'D')
                else:
                    update_job(db_batch, day, pair, 'WCT', 'E')
                
                jobs_processed += 1 
                
            db_batch.close()
            time.sleep(0.5)
            
        logger.info(f"Process completed {jobs_processed} jobs")
        
        # Close initial database connection
        db.close()
                        
    except Exception as e:
        logger.error(f"Error in main processing: {str(e)}")
        import traceback
        logger.error(traceback.format_exc())
    finally:
        try:
            if 'db' in locals() and db is not None:
                if hasattr(db, 'in_transaction') and db.in_transaction():
                    db.rollback()
                db.close()
        except:
            pass

    logger.info(f'*** Finished: Compute WCT with Adaptive Coda Detection ***')

          
            
            
            
            
            
####################################
####################################
def smoothCFS(cfs, scales, dt, ns, nt):
    """
    Smooth the continuous wavelet transform coefficients using a Fourier domain approach.
    Parameters
    ----------
    cfs : numpy.ndarray
        Continuous wavelet transform coefficients.
    scales : numpy.ndarray
        Scales used in the wavelet transform.
    dt : float
        Sampling interval.
    ns : int
        Smoothing parameter for the moving average filter.
    nt : float
        Smoothing parameter for the Gaussian filter.
    Returns
    -------
    numpy.ndarray
        Smoothed continuous wavelet transform coefficients.
    """
    N = cfs.shape[1]
    npad = sf.next_fast_len(N, real=True)
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
    
def conv2(x, y, mode='same'):
    """
    Perform 2D convolution of matrices x and y
    """
    return np.rot90(convolve2d(np.rot90(x, 2), np.rot90(y, 2), mode=mode), 2)

def get_wavelet_type(wavelet_type):
    """
    return a wavelet object based on the specified wavelet type and associated parameter
    """
    # Default parameters for each wavelet type
    default_params = {
        'Morlet': 6,
        'Paul': 4,
        'DOG': 2,
        'MexicanHat': 2  # MexicanHat inherits from DOG with m=2
    }

    wavelet_name = wavelet_type[0]

    # If a second argument is provided, use it; otherwise, use the default value
    if len(wavelet_type) == 2:
        param = float(wavelet_type[1])
    else:
        param = default_params[wavelet_name]

    # Get the corresponding wavelet object
    if wavelet_name == 'Morlet':
        return wavelet.Morlet(param)
    elif wavelet_name == 'Paul':
        return wavelet.Paul(param)
    elif wavelet_name == 'DOG':
        return wavelet.DOG(param)
    elif wavelet_name == 'MexicanHat':
        return wavelet.MexicanHat()  # Uses m=2, so no need for param
    else:
        raise logger.error(f"Unknown wavelet type: {wavelet_name}")

def xwt(trace_ref, trace_current, fs, ns=3, nt=0.25, vpo=12, freqmin=0.1, freqmax=8.0, nptsfreq=100, wavelet_type=('Morlet',6.), maxdt=0.1):
    """
    Wavelet coherence transform (WCT) on two time series..
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
    wavelet_type: list,
        Wavelet type and associated parameter.
        Default Morlet wavelet with a central frequency w0 = 6
    maxdt:float,
        Maximum allowed time delay in seconds. Time delays with absolute values
        exceeding this threshold will be clipped. 
        Default value is 0.1 second.
    
    Returns
        ----------
    WXamp : numpy.ndarray
        Amplitude of the cross-wavelet transform.
    WXspec : numpy.ndarray
        Complex cross-wavelet transform, representing both magnitude and phase information.
    WXangle : numpy.ndarray
        Phase angles of the cross-wavelet transform, indicating the phase relationship between the input signals.
    Wcoh : numpy.ndarray
        Wavelet coherence, representing the degree of correlation between the two signals in time-frequency space.
    WXdt : numpy.ndarray
        Time delay between the signals, estimated from the phase angles.
    freqs : numpy.ndarray
        Frequencies corresponding to the scales of the wavelet transform.
    coi : numpy.ndarray
        Cone of influence, representing the region of the wavelet spectrum where edge effects become significant.
    
    """
    mother = get_wavelet_type(wavelet_type) # mother wavelet class: Morlet, Paul, DOG, MexicanHat param 
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

    # Creation of the frequency vector
    freqlim = np.linspace(freqmax, freqmin, num=nptsfreq, endpoint=True, retstep=False, dtype=None, axis=0)

    try:
        # Wavelet transforms     # fft : Normalized fast Fourier transform of the input trace, fftfreqs : Fourier frequencies for the calculated FFT spectrum.
        cwt_reference, scales, freqs, coi, fft, fftfreqs = wavelet.cwt(x_reference, dt, dj, s0, J, mother, freqs=freqlim)
        cwt_current, _, _, _, _, _ = wavelet.cwt(x_current, dt, dj, s0, J, mother, freqs=freqlim)
    except Exception as e:
        logger.error(f"Error in wavelet transform: {str(e)}")
        raise

    scales = np.array([[kk] for kk in scales])
    invscales = np.kron(np.ones((1, nx)), 1 / scales)

    # Apply smoothing
    # Before the smoothCFS calls, ensure complex data types:
    power_ref = (invscales * abs(cwt_reference) ** 2).astype(complex)
    power_cur = (invscales * abs(cwt_current) ** 2).astype(complex)
    
    crossCFS = cwt_reference * np.conj(cwt_current)
    WXamp = abs(crossCFS)
    cross_spectrum = (invscales * crossCFS).astype(complex)
    
    # Then call smoothCFS with complex arrays:
    cfs1 = smoothCFS(power_ref, scales, dt, ns, nt)
    cfs2 = smoothCFS(power_cur, scales, dt, ns, nt)
    crossCFS = smoothCFS(cross_spectrum, scales, dt, ns, nt)
    
    # Handle zeros to prevent divide-by-zero
    mask1 = cfs1 > 0
    mask2 = cfs2 > 0
    valid_mask = mask1 & mask2  # Areas where both have power
            
    # Initialize with NaNs
    WXspec = np.full_like(crossCFS, np.nan, dtype=complex)
    Wcoh = np.full_like(crossCFS, np.nan)
    
    # Division operations
    WXspec[valid_mask] = crossCFS[valid_mask] / (np.sqrt(cfs1[valid_mask]) * np.sqrt(cfs2[valid_mask]))
    Wcoh[valid_mask] = abs(crossCFS[valid_mask]) ** 2 / (cfs1[valid_mask] * cfs2[valid_mask])
    WXangle = np.angle(WXspec)

    
    # Clip coherence values to [0,1]
    Wcoh = np.clip(Wcoh, 0.0, 1.0)

    # Calculate time delay
    pp = 2 * np.pi * freqs
    pp2 = np.array([[kk] for kk in pp])
    WXdt = WXangle / np.kron(np.ones((1, nx)), pp2)
    
    #if maxdt is not None:
    #WXdt = np.clip(WXdt, -maxdt, maxdt) 
    
    return WXamp, WXspec, WXangle, Wcoh, WXdt, freqs, coi
    
    
def compute_wct_dvv(freqs, tvec, WXamp, Wcoh, delta_t, lag_min=5, coda_cycles=20, mincoh=0.5, maxdt=0.2, 
            min_nonzero=0.25, freqmin=0.1, freqmax=2.0):
    """
    Compute the dv/v values and associated errors from the wavelet transform results.
    Parameters
    ----------
    freqs : numpy.ndarray
        Frequency values corresponding to the wavelet transform.
    tvec : numpy.ndarray
        Time vector.
    WXamp : numpy.ndarray
        Amplitude of the cross-wavelet transform.
    Wcoh : numpy.ndarray
        Wavelet coherence.
    delta_t : numpy.ndarray
        Time delays between signals.
    lag_min : int, optional
        Minimum lag in seconds. Default is 5.
    coda_cycles : int, optional
        Number of coda cycles to consider. Default is 20.
    mincoh : float, optional
        Minimum coherence value for weighting. Default is 0.5.
    maxdt : float, optional
        Maximum time delay for weighting. Default is 0.2.
    min_nonzero : float, optional
        Minimum percentage of non-zero weights required for valid estimation. Default is 0.25.
    freqmin : float, optional
        Minimum frequency for calculation. Default is 0.1 Hz.
    freqmax : float, optional
        Maximum frequency for calculation. Default is 2.0 Hz.
    Returns
    -------
    tuple
        dvv values (percentage), errors (percentage), and weighting function used.
    """   
    import warnings
    from scipy.optimize import OptimizeWarning
    
    inx = np.where((freqs >= freqmin) & (freqs <= freqmax))  # Filter frequencies within the specified range
    dvv, err = np.zeros(len(inx[0])), np.zeros(len(inx[0])) # Initialize dvv and err arrays

    # Weighting function based on WXamp
    weight_func = np.log(np.abs(WXamp)) / np.log(np.abs(WXamp)).max()
    zero_idx = np.where((Wcoh < mincoh) | (delta_t > maxdt))
    wf = (weight_func + abs(np.nanmin(weight_func))) / weight_func.max()
    wf[zero_idx] = 0
    
    # Track frequencies with estimation issues
    problematic_freqs = []

    # Loop through frequency indices for linear regression
    for ii, ifreq in enumerate(inx[0]):
        period = 1.0 / freqs[ifreq]
        lag_max = lag_min + (period * coda_cycles)

        # Coda selection
        tindex = np.where(((tvec >= -lag_max) & (tvec <= -lag_min)) | ((tvec >= lag_min) & (tvec <= lag_max)))[0]

        if len(tvec) > 2:
            if not np.any(delta_t[ifreq]):
                continue

            delta_t[ifreq][tindex] = np.nan_to_num(delta_t[ifreq][tindex])
            w = wf[ifreq] # Weighting function for the specific frequency
            non_zero_mask = w > 0
            nzc_perc = np.count_nonzero(non_zero_mask[tindex]) / len(tindex)
            
            if nzc_perc >= min_nonzero:
                
                with warnings.catch_warnings(record=True) as w_catcher:
                    warnings.simplefilter("always", OptimizeWarning)

                    m, em = linear_regression(tvec[tindex], delta_t[ifreq][tindex], w[tindex], intercept_origin=True)
                    
                    if any(issubclass(warning.category, OptimizeWarning) for warning in w_catcher):
                        problematic_freqs.append(freqs[ifreq])
                        
                dvv[ii], err[ii] = -m, em
            else:
                dvv[ii], err[ii] = np.nan, np.nan
        else:
            logger.debug('Not enough points to estimate dv/v for WCT')
    if problematic_freqs:
        logger.warning(f"Covariance estimation issues between {min(problematic_freqs):.2f}-{max(problematic_freqs):.2f} Hz: Modify min_nonzero (current: {min_nonzero}), mincoh (current: {mincoh}), maxdt (current: {maxdt}) and/or coda_cycles (current: {coda_cycles})")
    
    return dvv * 100, err * 100, wf
 
def get_avgcoh(freqs, tvec, wcoh, freqmin, freqmax, lag_min=5, coda_cycles=20):
    """
    Calculate the average wavelet coherence over a specified frequency range and time lags.

    :param freqs: A numpy array that represents frequency values.
    :type freqs: numpy.ndarray
    :param tvec: A time vector represented as a numpy array.
    :type tvec: numpy.ndarray
    :param wcoh: The wavelet coherence array, represented as a numpy array.
    :type wcoh: numpy.ndarray
    :param freqmin: The minimum frequency for coherence calculation, represented as a floating-point number.
    :type freqmin: float
    :param freqmax: The maximum frequency for coherence calculation, represented as a floating-point number.
    :type freqmax: float
    :param lag_min: The minimum lag in seconds for coherence calculation. This is optional and it defaults to 5.
    :type lag_min: int, optional
    :param coda_cycles: The number of coda cycles to consider. This is optional and it defaults to 20.
    :type coda_cycles: int, optional
    :returns: A numpy array of average coherence values computed over the specified frequency range and time lags.
    :rtype: numpy.ndarray
    """
    inx = np.where((freqs>=freqmin) & (freqs<=freqmax)) 
    coh = np.zeros(inx[0].shape) # Create empty vector for coherence

    for ii, ifreq in enumerate(inx[0]): # Loop through frequencies index     
        period = 1.0/freqs[ifreq]
        lag_max = lag_min + (period*coda_cycles) 
        tindex = np.where(((tvec >= -lag_max) & (tvec <= -lag_min)) | ((tvec >= lag_min) & (tvec <= lag_max)))[0] # Index of the coda

        if len(tvec)>2: # check time vector size
            if not np.any(wcoh[ifreq]) or wcoh[ifreq][tindex].size == 0:
                coh[ii] = np.nan
                continue
            c = np.nanmean(np.abs(wcoh[ifreq][tindex]))
            coh[ii] = c

        else:
            logger.debug('Not enough points to compute average coherence') #not sure why it would ever get here, but just in case.
            coh[ii] = np.nan

    return coh

def save_day_wct_results(station1, station2, date, component, filterid, mov_stack, taxis, dvv, err, coh, freqs, inx):
    """
    Save the WCT results for a single day to a file in a hierarchical folder structure:
    WCT/[filter_id]/[mov_stack]/[component]/[station_pair]/[date].npz
    
    Parameters:
    ----------
    station1, station2 : str
        Names of the stations
    date : datetime
        The date of the measurement
    component : str
        Component code
    filterid : int
        Filter ID
    mov_stack : int or str
        Moving stack value
    taxis : numpy.ndarray
        Time axis
    dvv, err, coh : numpy.ndarray
        The dv/v values, errors, and coherence values
    freqs : numpy.ndarray
        Frequency values
    inx : tuple
        Indices of valid frequencies
        
    Returns:
    -------
    str
        Path to the saved file
    """
    # Format filterid with leading zeros (e.g., 1 -> "01")
    filter_str = f"{int(filterid):02d}"
    
    # Format mov_stack (could be string like '10d_1d' or int)
    if isinstance(mov_stack, (list, tuple)):
        mov_stack_str = '_'.join(str(m) for m in mov_stack)
    else:
        mov_stack_str = str(mov_stack)
    
    # Station pair string
    pair_str = f"{station1}_{station2}"
    
    # Create directory structure
    # WCT/[filter_id]/[mov_stack]/[component]/[station_pair]/
    dir_path = os.path.join('WCT', filter_str, mov_stack_str, component, pair_str)
    os.makedirs(dir_path, exist_ok=True)
    
    # Filename is just the date: YYYY-MM-DD.npz
    date_str = date.strftime('%Y-%m-%d')
    filename = f"{date_str}.npz"
    filepath = os.path.join(dir_path, filename)
    
    # Create DataFrames
    freqs_subset = freqs[inx]
    
    # Save to numpy compressed format
    np.savez_compressed(
        filepath,
        date=date,
        freqs=freqs_subset,
        dvv=dvv,
        err=err,
        coh=coh,
        taxis=taxis,
        station1=station1,
        station2=station2,
        component=component,
        filterid=filterid,
        mov_stack=mov_stack
    )
    
    logger.debug(f"Saved WCT results to {filepath}")
    return filepath
    
def process_wct_job(pair, day, params, taxis, filters): 
    """
    Process a single WCT job without database access.
    
    Parameters:
    ----------
    job_data : dict
        Dictionary with job information (ref, pair, day, flag)
    params : dict-like
        Configuration parameters
    taxis : numpy.ndarray
        Time axis
    filters : list
        List of filters to use
        
    Returns:
    -------
    bool
        True if job was successfully processed, False otherwise
    """
    # Parse pair and day
    station1, station2 = pair.split(":")
    
    # Determine components to compute
    if station1 == station2:
        components_to_compute = params.components_to_compute_single_station
    else:
        components_to_compute = params.components_to_compute
        
    logger.info(f"Processing Pair: {pair}, Day: {day} for {len(filters)} filters, {len(components_to_compute)} components and {len(params.mov_stack)} moving windows")
    
    # Extract parameters
    ns = params.wct_ns
    nt = params.wct_nt 
    vpo = params.wct_vpo 
    nptsfreq = params.wct_nptsfreq
    coda_cycles = params.dtt_codacycles 
    min_nonzero = params.dvv_min_nonzero
    wct_norm = params.wct_norm
    wavelet_type = params.wavelet_type
    
    mov_stacks = params.mov_stack
    goal_sampling_rate = params.cc_sampling_rate
    lag_min = params.dtt_minlag
    maxdt = params.dtt_maxdt
    mincoh = params.dtt_mincoh
        
    # Convert date string to datetime
    try:
        date = pd.to_datetime(day)
    except:
        logger.error(f"Failed to parse date: {day}")
        return False
    
    operations_succeeded = 0  # Track if all operations for this job succeed

    # Process each filter
    for f in filters:
        filterid = int(f['ref'])
        freqmin = f['low']
        freqmax = f['high']
        
        # Process each component
        for component in components_to_compute:
            try:
                # Get reference waveform
                ref = xr_get_ref(station1, station2, component, filterid, taxis, ignore_network=True)
                ref = ref.CCF.values
                if wct_norm:
                    ori_waveform = (ref/ref.max()) 
                else:
                    ori_waveform = ref
            except FileNotFoundError as fullpath:
                logger.error(f"FILE DOES NOT EXIST: {fullpath}, skipping")
                continue
            except Exception as e:
                logger.error(f"Error getting reference waveform: {str(e)}")
                continue
            
            if not len(ref):
                logger.warning(f"Empty reference waveform for {station1}:{station2} {component} f{filterid}")
                continue
            
            # Process each moving stack
            for mov_stack in mov_stacks:
                try:
                    # Get the data for this day and moving stack
                    data = xr_get_ccf(station1, station2, component, filterid, mov_stack, taxis)
                    
                    # Filter data to get only this specific day
                    day_data = data[data.index == date]
                    
                    if day_data.empty:
                        logger.warning(f"No data for {date} in {station1}:{station2} {component} f{filterid} m{mov_stack}")
                        continue
                    
                    logger.debug(f"Processing {station1}:{station2} {component} f{filterid} m{mov_stack} - {date}")
                    
                    # Get the waveform for this day
                    waveform = day_data.iloc[0].values
                    if wct_norm:
                        new_waveform = (waveform/waveform.max()) 
                    else:
                        new_waveform = waveform
                    
                    # Calculate wavelet coherence transform
                    WXamp, WXspec, WXangle, Wcoh, WXdt, freqs, coi = xwt(
                        ori_waveform, new_waveform, goal_sampling_rate, 
                        int(ns), float(nt), int(vpo), 
                        freqmin, freqmax, int(nptsfreq), wavelet_type, maxdt
                    )
                    
                    # Compute dv/v values
                    dvv, err, wf = compute_wct_dvv(
                        freqs, taxis, WXamp, Wcoh, WXdt, 
                        lag_min=int(lag_min), coda_cycles=coda_cycles, 
                        mincoh=mincoh, maxdt=maxdt, min_nonzero=min_nonzero, 
                        freqmin=freqmin, freqmax=freqmax
                    )
                    
                    # Calculate average coherence
                    coh = get_avgcoh(
                        freqs, taxis, Wcoh, freqmin, freqmax, 
                        lag_min=int(lag_min), coda_cycles=coda_cycles
                    )
                    
                    # Get frequency indices
                    inx = np.where((freqs >= freqmin) & (freqs <= freqmax))
                    
                    # Save results for this day
                    save_day_wct_results(
                        station1, station2, date, component, 
                        filterid, mov_stack, taxis, 
                        dvv, err, coh, freqs, inx
                    )
                    
                    operations_succeeded += 1
                
                except Exception as e:
                    logger.error(f"Error processing {component} f{filterid} m{mov_stack}: {str(e)}")
                    continue
    
    # Determine overall success status
    if operations_succeeded == 0:
        logger.error(f"Job failed completely: {pair} - {day}.")
        return False
    else:
        logger.info(f"Job processing complete: {pair} - {day}. {operations_succeeded} operations succeeded.")
        return True

def claim_available_jobs(db, batch_size=5, max_attempts=3):
    """
    Claims a batch of available jobs
    
    Parameters:
    ----------
    db : Session
        Database session
    batch_size : int
        Maximum number of jobs to claim at once
    max_attempts : int
        Number of retries if no jobs found
        
    Returns:
    -------
    list
        List of claimed Job objects
    """
    from sqlalchemy import text
    import time
        
    for attempt in range(max_attempts):
        try:
            # Make sure we're not in a transaction already
            if db.in_transaction():
                db.rollback()
            
            # Find available jobs
            available_jobs = db.query(Job).filter(
                Job.flag == 'T',
                Job.jobtype == 'WCT',
                Job.day != 'REF'
            ).limit(batch_size).all()
            
            if not available_jobs:
                if attempt < max_attempts - 1:
                    time.sleep(1 * (attempt + 1))
                    continue
                else:
                    logger.info("No jobs available after {max_attempts} retries")
                    return []
            
            # Make a copy of the job data before updating
            jobs = []
            job_ids = []
            
            for job in available_jobs:
                # Properly create job copies with all required attributes
                job_copy = Job(
                    day=job.day,
                    pair=job.pair,
                    jobtype=job.jobtype,
                    flag='I'
                )
                job_copy.ref = job.ref
                jobs.append(job_copy)
                job_ids.append(job.ref)
            
            # Start a transaction for the update
            if db.in_transaction():
                db.rollback()
            db.begin()
            
            # Update the jobs to mark them as in progress
            update_query = text("""
                UPDATE jobs
                SET flag = 'I', lastmod = NOW()
                WHERE ref IN :job_ids
                AND flag = 'T'  -- Double-check it hasn't been claimed
            """)
            
            result = db.execute(update_query, {'job_ids': tuple(job_ids)})
            updated_count = result.rowcount
            
            # Commit the transaction
            db.commit()
            
            # Check if we actually updated anything (race condition check)
            if updated_count > 0:
                # Only return jobs that were actually updated
                return jobs[:updated_count]
            else:
                if attempt < max_attempts - 1:
                    time.sleep(1 * (attempt + 1))
                    continue
                else:
                    return []
                
        except Exception as e:
            if db.in_transaction():
                db.rollback()
            logger.error(f"Error claiming jobs: {str(e)}")
            time.sleep(1)
    
    return []
  
def main_f(loglevel="INFO", batch_size=5):
    """
    Main function to process WCT jobs using day job distribution approach
    """
    global logger
    logger = get_logger('msnoise.compute_wct_child', loglevel, with_pid=True)
    logger.info('*** Starting: Compute WCT ***')
    
    # Add small random delay to prevent all processes starting at exactly the same time
    time.sleep(np.random.random() * 2)
    
    try:
        # Initial database connection to get configuration
        db = connect()
        
        # Ensure no lingering transactions
        if hasattr(db, 'in_transaction') and db.in_transaction():
            db.rollback()
            
        params = get_params(db)
        taxis = get_t_axis(db)
        filters = get_filters(db, all=False)
        
        # Convert filters to dictionaries before closing the session
        filter_dicts = []
        for f in filters:
            filter_dicts.append({
                'ref': f.ref,
                'low': f.low,
                'high': f.high,
                'mwcs_low': getattr(f, 'mwcs_low', None),
                'mwcs_high': getattr(f, 'mwcs_high', None),
                'rms_threshold': getattr(f, 'rms_threshold', None),
                'mwcs_wlen': getattr(f, 'mwcs_wlen', None),
                'mwcs_step': getattr(f, 'mwcs_step', None)
            })
        db.close()
        
        # Process jobs in batches to avoid too many database connections
        jobs_processed = 0
        max_jobs_per_process = 1000  # Safety limit
        consecutive_empty_batches = 0
        max_consecutive_empty = 3  # Exit after this many empty batches
        
        while jobs_processed < max_jobs_per_process:
            # Connect to the database and get a batch of jobs
            db = connect()
            
            # Ensure no lingering transactions
            if hasattr(db, 'in_transaction') and db.in_transaction():
                db.rollback()
            
            # Check if there are any jobs left
            job_exists = is_dtt_next_job(db, flag='T', jobtype='WCT')
            if not job_exists:
                logger.info("No more WCT jobs to process")
                db.close()
                break
                
            # Try to claim a batch of jobs
            jobs = claim_available_jobs(db, batch_size)
            db.close()
            
            if not jobs:
                consecutive_empty_batches += 1
                if consecutive_empty_batches >= max_consecutive_empty:
                    logger.info(f"Exiting after {consecutive_empty_batches} empty batches")
                    break
                time.sleep(2)  # Wait before next attempt
                continue
            else:
                consecutive_empty_batches = 0  # Reset counter when we get jobs
                
            # Process each job in the batch
            for job in jobs:
                pair = job.pair
                day = job.day

                # Process this job
                success = process_wct_job(pair, day, params, taxis, filter_dicts)
                                
                # Update job status
                db = connect()
                if success:
                    update_job(db, day, pair, 'WCT', 'D')
                else:
                    update_job(db, day, pair, 'WCT', 'E')  # Mark as error
                db.close()
                
                jobs_processed += 1 
                
            # Small delay between batches
            time.sleep(0.5)
            
        logger.info(f"Process completed {jobs_processed} jobs")
                        
    except Exception as e:
        logger.error(f"Error in main processing: {str(e)}")
        import traceback
        logger.error(traceback.format_exc())
    finally:
        # Final cleanup
        try:
            if 'db' in locals() and db is not None:
                # Ensure we don't leave transactions open
                if hasattr(db, 'in_transaction') and db.in_transaction():
                    db.rollback()
                db.close()
        except:
            pass

    logger.info(f'*** Finished: Compute WCT ***')