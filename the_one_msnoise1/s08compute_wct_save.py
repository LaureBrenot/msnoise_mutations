"""
Wavelet Coherence Transform (WCT) Computation
This script performs the computation of the Wavelet Coherence Transform (WCT), a tool used to analyze the correlation between two time series in the time-frequency domain. The script supports parallel processing and interacts with a database to manage job statuses.

Filter Configuration Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* |dtt_minlag|
* |dtt_maxdt|
* |dtt_mincoh|
* |dtt_codacycles|
* |wct_ns|
* |wct_nt|
* |wct_vpo|
* |wct_nptsfreq|
* |dvv_min_nonzero|
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
import scipy.optimize
import scipy.signal
import pycwt as wavelet
from scipy.signal import convolve2d
from obspy.signal.regression import linear_regression
from .api import *
import logbook
import scipy
import scipy.fft as sf

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
            if not np.any(wcoh[ifreq]): # check non-empty dt array
                continue
            c = np.nanmean(wcoh[ifreq][tindex])
            coh[ii] = c

        else:
            logger.debug('not enough points to compute average coherence') #not sure why it would ever get here, but just in case.
            coh[ii] = np.nan

    return coh

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
        cfs[kk - 1] = np.real(smooth[0:N])
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

def compute_wct_dvv2(freqs, tvec, WXamp, Wcoh, delta_t, lag_min=5, coda_cycles=20, mincoh=0.5, maxdt=0.2, 
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
    inx = np.where((freqs >= freqmin) & (freqs <= freqmax))  # Filter frequencies within the specified range
    dvv, err = np.zeros(len(inx[0])), np.zeros(len(inx[0])) # Initialize dvv and err arrays

    # Weighting function based on WXamp
    weight_func = np.log(np.abs(WXamp)) / np.log(np.abs(WXamp)).max()
    zero_idx = np.where((Wcoh < mincoh) | (delta_t > maxdt))
    wf = (weight_func + abs(np.nanmin(weight_func))) / weight_func.max()
    wf[zero_idx] = 0

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
            w = wf[ifreq]  # Weighting function for the specific frequency
            w[~np.isfinite(w)] = 1.0

            # Ensure no zero weights to prevent division by zero
            w_for_regression = np.copy(w[tindex])
            w_for_regression[w_for_regression < 1e-10] = 1e-10
            #
            
            # Percentage of non-zero weights
            nzc_perc = np.count_nonzero(w[tindex]) / len(tindex)
            if nzc_perc >= min_nonzero:
                try:
                    m, em = linear_regression(tvec[tindex], delta_t[ifreq][tindex], w[tindex], intercept_origin=True)
                    dvv[ii], err[ii] = -m, em
                except RuntimeError:
                    # Handle the case when regression fails
                    logger.debug(f"Regression failed for frequency {freqs[ifreq]} Hz")
                    dvv[ii], err[ii] = np.nan, np.nan
            else:
                dvv[ii], err[ii] = np.nan, np.nan
        else:
            logger.debug('Not enough points to estimate dv/v for WCT')
    
    return dvv * 100, err * 100, wf

def compute_wct_dvv3(freqs, tvec, WXamp, Wcoh, delta_t, lag_min=5, coda_cycles=20, mincoh=0.5, maxdt=0.2, 
            min_nonzero=0.25, freqmin=0.1, freqmax=2.0):
    """
    Compute the dv/v values and associated errors from the wavelet transform results.
    """   
    inx = np.where((freqs >= freqmin) & (freqs <= freqmax))  # Filter frequencies within the specified range
    dvv, err = np.zeros(len(inx[0])), np.zeros(len(inx[0])) # Initialize dvv and err arrays

    # Prevent NaN and inf values in input arrays
    WXamp = np.nan_to_num(WXamp, nan=0.0, posinf=1.0, neginf=-1.0)
    Wcoh = np.nan_to_num(Wcoh, nan=0.0, posinf=1.0, neginf=0.0)
    delta_t = np.nan_to_num(delta_t, nan=0.0, posinf=maxdt, neginf=-maxdt)
    
    # Ensure all values of Wcoh are between 0 and 1
    Wcoh = np.clip(Wcoh, 0.0, 1.0)
    
    # Ensure delta_t values are within reasonable bounds
    delta_t = np.clip(delta_t, -maxdt, maxdt)

    # Weighting function based on WXamp
    # Prevent log(0) by adding a small value
    eps = 1e-10
    weight_func = np.log(np.abs(WXamp) + eps) / np.log(np.abs(WXamp).max() + eps)
    zero_idx = np.where((Wcoh < mincoh) | (np.abs(delta_t) > maxdt))
    wf = (weight_func + abs(np.nanmin(weight_func))) / (weight_func.max() + eps)
    wf[zero_idx] = 0

    # Loop through frequency indices for linear regression
    for ii, ifreq in enumerate(inx[0]):
        period = 1.0 / freqs[ifreq]
        lag_max = lag_min + (period * coda_cycles)

        # Coda selection
        tindex = np.where(((tvec >= -lag_max) & (tvec <= -lag_min)) | ((tvec >= lag_min) & (tvec <= lag_max)))[0]

        if len(tvec) > 2 and len(tindex) > 2:  # Need at least 3 points for regression
            if not np.any(delta_t[ifreq]):
                dvv[ii], err[ii] = np.nan, np.nan
                continue

            # Get weights and ensure they are valid
            w = np.copy(wf[ifreq])
            # Ensure weights are finite and positive
            w[~np.isfinite(w)] = 0.0
            w = np.clip(w, 0.0, 1.0)  # Ensure weights are between 0 and 1
            
            # Prepare data for regression
            x_data = tvec[tindex]
            y_data = delta_t[ifreq][tindex]
            weights = w[tindex]
            
            # Percentage of non-zero weights
            nzc_perc = np.count_nonzero(weights) / len(tindex)
            
            if nzc_perc >= min_nonzero and np.count_nonzero(weights) > 2:
                # Create modified weights that won't cause division by zero
                reg_weights = np.copy(weights)
                reg_weights[reg_weights < 1e-6] = 1e-6  # Minimum weight
                
                try:
                    # Try linear regression with safeguards
                    m, em = linear_regression(x_data, y_data, reg_weights, intercept_origin=True)
                    dvv[ii], err[ii] = -m, em
                except Exception as e:
                    logger.debug(f"Regression failed for frequency {freqs[ifreq]} Hz: {str(e)}")
                    dvv[ii], err[ii] = np.nan, np.nan
            else:
                dvv[ii], err[ii] = np.nan, np.nan
        else:
            logger.debug('Not enough points to estimate dv/v for WCT')
            dvv[ii], err[ii] = np.nan, np.nan
    
    return dvv * 100, err * 100, wf
 
def compute_wct_dvv(freqs, tvec, WXamp, Wcoh, delta_t, lag_min=5, coda_cycles=20, mincoh=0.5, maxdt=0.2, 
            min_nonzero=0.25, freqmin=0.1, freqmax=2.0):
    """
    Compute the dv/v values and associated errors from the wavelet transform results with robust regression.
    """   
    import numpy as np
    from numpy.linalg import lstsq
    
    inx = np.where((freqs >= freqmin) & (freqs <= freqmax))  # Filter frequencies within the specified range
    dvv, err = np.zeros(len(inx[0])), np.zeros(len(inx[0])) # Initialize dvv and err arrays

    # Prevent NaN and inf values in input arrays
    WXamp = np.nan_to_num(WXamp, nan=0.0, posinf=1.0, neginf=-1.0)
    Wcoh = np.nan_to_num(Wcoh, nan=0.0, posinf=1.0, neginf=0.0)
    delta_t = np.nan_to_num(delta_t, nan=0.0, posinf=maxdt, neginf=-maxdt)
    
    # Ensure all values of Wcoh are between 0 and 1
    Wcoh = np.clip(Wcoh, 0.0, 1.0)
    
    # Ensure delta_t values are within reasonable bounds
    delta_t = np.clip(delta_t, -maxdt, maxdt)

    # Weighting function based on WXamp
    # Prevent log(0) by adding a small value
    eps = 1e-10
    weight_func = np.log(np.abs(WXamp) + eps) / np.log(np.abs(WXamp).max() + eps)
    zero_idx = np.where((Wcoh < mincoh) | (np.abs(delta_t) > maxdt))
    wf = (weight_func + abs(np.nanmin(weight_func))) / (weight_func.max() + eps)
    wf[zero_idx] = 0

    # Loop through frequency indices for robust linear regression
    for ii, ifreq in enumerate(inx[0]):
        period = 1.0 / freqs[ifreq]
        lag_max = lag_min + (period * coda_cycles)

        # Coda selection
        tindex = np.where(((tvec >= -lag_max) & (tvec <= -lag_min)) | ((tvec >= lag_min) & (tvec <= lag_max)))[0]

        if len(tvec) > 2 and len(tindex) > 2:  # Need at least 3 points for regression
            if not np.any(delta_t[ifreq]):
                dvv[ii], err[ii] = np.nan, np.nan
                continue

            # Get weights and ensure they are valid
            w = np.copy(wf[ifreq][tindex])
            # Ensure weights are finite and positive
            w[~np.isfinite(w)] = 0.0
            w = np.clip(w, 0.0, 1.0)  # Ensure weights are between 0 and 1
            
            # Prepare data for regression
            x_data = tvec[tindex]
            y_data = delta_t[ifreq][tindex]
            
            # Percentage of non-zero weights
            valid_idx = w > 1e-6
            nzc_perc = np.sum(valid_idx) / len(tindex)
            
            if nzc_perc >= min_nonzero and np.sum(valid_idx) > 2:
                try:
                    # Use only valid data points
                    x_valid = x_data[valid_idx]
                    y_valid = y_data[valid_idx]
                    w_valid = w[valid_idx]
                    
                    # Apply weights to data
                    x_weighted = x_valid * np.sqrt(w_valid)
                    y_weighted = y_valid * np.sqrt(w_valid)
                    
                    # Simple linear regression through origin
                    # For y = mx (no intercept), using matrix algebra:
                    # m = (x'x)^-1 * x'y where x' is the transpose of x
                    # For weighted regression, we use the weighted x and y
                    
                    # Reshape for matrix operations
                    X = x_weighted.reshape(-1, 1)
                    Y = y_weighted.reshape(-1, 1)
                    
                    # Use numpy's least squares to solve for m
                    # This is more numerically stable than direct formula
                    m, residuals, rank, s = lstsq(X, Y, rcond=None)
                    
                    if rank == 0:  # Check if regression failed
                        dvv[ii], err[ii] = np.nan, np.nan
                        continue
                        
                    # m is now a 1x1 array, extract scalar
                    m_value = m[0, 0] if m.shape == (1, 1) else m[0]
                    
                    # Calculate standard error
                    n = len(x_valid)
                    if n > 2:
                        # Unweighted predictions for error calculation
                        y_pred = m_value * x_valid
                        
                        # Weighted residuals
                        weighted_residuals = (y_valid - y_pred) * np.sqrt(w_valid)
                        
                        # Sum of squared weighted residuals
                        ss_resid = np.sum(weighted_residuals**2)
                        
                        # Degrees of freedom
                        df = n - 1  # 1 parameter (slope) being estimated
                        
                        # Mean squared error
                        mse = ss_resid / df if df > 0 else np.nan
                        
                        # Standard error of the regression coefficient
                        # SE(m) = sqrt(MSE / sum(w_i * x_i^2))
                        se_m = np.sqrt(mse / np.sum(w_valid * x_valid**2)) if np.sum(w_valid * x_valid**2) > 0 else np.nan
                        
                        dvv[ii], err[ii] = -m_value, se_m
                    else:
                        dvv[ii], err[ii] = -m_value, np.nan
                
                except Exception as e:
                    logger.debug(f"Direct regression failed for frequency {freqs[ifreq]} Hz: {str(e)}")
                    
                    # Final fallback: simple non-weighted linear regression
                    try:
                        # Just use the formula directly
                        m_simple = np.sum(x_valid * y_valid) / np.sum(x_valid**2) if np.sum(x_valid**2) > 0 else np.nan
                        dvv[ii], err[ii] = -m_simple, np.nan
                    except:
                        dvv[ii], err[ii] = np.nan, np.nan
            else:
                dvv[ii], err[ii] = np.nan, np.nan
        else:
            logger.debug('Not enough points to estimate dv/v for WCT')
            dvv[ii], err[ii] = np.nan, np.nan
    
    return dvv * 100, err * 100, wf
      
def xwt(trace_ref, trace_current, fs, ns=3, nt=0.25, vpo=12, freqmin=0.1, freqmax=8.0, nptsfreq=100, wavelet_type=('Morlet',6.)):
    """
    Wavelet coherence transform (WCT) on two time series.
    """
    mother = get_wavelet_type(wavelet_type)
    nx = np.size(trace_current)
    x_reference = np.transpose(trace_ref)
    x_current = np.transpose(trace_current)
    dt = 1 / fs
    dj = 1 / vpo 
    J = -1
    s0 = 2 * dt

    # Creation of the frequency vector
    freqlim = np.linspace(freqmax, freqmin, num=nptsfreq, endpoint=True, retstep=False, dtype=None, axis=0)

    try:
        # Wavelet transforms
        cwt_reference, scales, freqs, coi, fft, fftfreqs = wavelet.cwt(x_reference, dt, dj, s0, J, mother, freqs=freqlim)
        cwt_current, _, _, _, _, _ = wavelet.cwt(x_current, dt, dj, s0, J, mother, freqs=freqlim)
    except Exception as e:
        logger.error(f"Error in wavelet transform: {str(e)}")
        raise

    scales = np.array([[kk] for kk in scales])
    invscales = np.kron(np.ones((1, nx)), 1 / scales)

    # Apply smoothing
    cfs2 = smoothCFS(invscales * abs(cwt_current) ** 2, scales, dt, ns, nt)
    cfs1 = smoothCFS(invscales * abs(cwt_reference) ** 2, scales, dt, ns, nt)
    
    # Handle zeros to prevent divide-by-zero
    eps = 1e-10  # Small constant to avoid division by zero
    cfs1_safe = np.maximum(cfs1, eps)
    cfs2_safe = np.maximum(cfs2, eps)
    
    crossCFS = cwt_reference * np.conj(cwt_current)
    WXamp = abs(crossCFS)
    
    # Smoothing
    crossCFS = smoothCFS(invscales * crossCFS, scales, dt, ns, nt)
    
    # Safe division operations
    WXspec = crossCFS / (np.sqrt(cfs1_safe) * np.sqrt(cfs2_safe))
    WXangle = np.angle(WXspec)
    Wcoh = abs(crossCFS) ** 2 / (cfs1_safe * cfs2_safe)
    
    # Clip coherence values to [0,1]
    Wcoh = np.clip(Wcoh, 0.0, 1.0)
    
    # Calculate time delay
    pp = 2 * np.pi * freqs
    pp2 = np.array([[kk] for kk in pp])
    WXdt = WXangle / np.kron(np.ones((1, nx)), pp2)
    
    # Limit extreme values in WXdt
    max_dt = 1.0  # Maximum allowed time delay
    WXdt = np.clip(WXdt, -max_dt, max_dt)

    return WXamp, WXspec, WXangle, Wcoh, WXdt, freqs, coi
    
      
def xwt2(trace_ref, trace_current, fs, ns=3, nt=0.25, vpo=12, freqmin=0.1, freqmax=8.0, nptsfreq=100, wavelet_type=('Morlet',6.)):
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
        
    cfs1_safe = np.copy(cfs1)
    cfs2_safe = np.copy(cfs2)
    
    # Set a small positive value instead of zeros
    min_positive = 1e-10
    cfs1_safe[cfs1_safe < min_positive] = min_positive
    cfs2_safe[cfs2_safe < min_positive] = min_positive
    
    crossCFS = cwt_reference * np.conj(cwt_current)
    WXamp = abs(crossCFS)
    # cross-wavelet transform operation with smoothing
    crossCFS = smoothCFS(invscales * crossCFS, scales, dt, ns, nt)
    
    # Use the safe versions for division
    WXspec = crossCFS / (np.sqrt(cfs1_safe) * np.sqrt(cfs2_safe))
    WXangle = np.angle(WXspec)
    Wcoh = abs(crossCFS) ** 2 / (cfs1_safe * cfs2_safe)
#    crossCFS = cwt_reference * np.conj(cwt_current)
#    WXamp = abs(crossCFS)
#    # cross-wavelet transform operation with smoothing
#    crossCFS = smoothCFS(invscales * crossCFS, scales, dt, ns, nt)
#    WXspec = crossCFS / (np.sqrt(cfs1) * np.sqrt(cfs2))
#    WXangle = np.angle(WXspec)
#    Wcoh = abs(crossCFS) ** 2 / (cfs1 * cfs2)
    pp = 2 * np.pi * freqs
    pp2 = np.array([[kk] for kk in pp])
    WXdt = WXangle / np.kron(np.ones((1, nx)), pp2)


    return WXamp, WXspec, WXangle, Wcoh, WXdt, freqs, coi

def main2(loglevel="INFO"):
    # Reconfigure logger to show the pid number in log records
    global logger
    logger = get_logger('msnoise.compute_wct_child', loglevel,
                        with_pid=True)
    logger.info('*** Starting: Compute WCT ***')

    db = connect()
    params = get_params(db)
    taxis = get_t_axis(db)

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

    logger.debug('Ready to compute')
    # Then we compute the jobs
    filters = get_filters(db, all=False)
    time.sleep(np.random.random() * 5)

    while is_dtt_next_job(db, flag='T', jobtype='WCT'):
        # TODO would it be possible to make the next 8 lines in the API ?
        jobs = get_dtt_next_job(db, flag='T', jobtype='WCT')

        if not len(jobs):
            # edge case, should only occur when is_next returns true, but
            # get_next receives no jobs (heavily parallelised calls).
            time.sleep(np.random.random())
            continue
        pair = jobs[0].pair
        print(pair)
        refs, days = zip(*[[job.ref, job.day] for job in jobs])

        logger.info(
            "There are WCT jobs for some days to recompute for %s" % pair)
        for f in filters:
            filterid = int(f.ref)
            freqmin = f.low
            freqmax = f.high

            station1, station2 = pair.split(":")
            if station1 == station2:
                components_to_compute = params.components_to_compute_single_station
            else:
                components_to_compute = params.components_to_compute

            for components in components_to_compute:
                try:
                    ref = xr_get_ref(station1, station2, components, filterid, taxis, ignore_network=True)
                    #print(ref)
                    ref = ref.CCF.values
                    if wct_norm:
                        ori_waveform = (ref/ref.max()) 
                    else:
                        ori_waveform = ref
                except FileNotFoundError as fullpath:
                    logger.error("FILE DOES NOT EXIST: %s, skipping" % fullpath)
                    continue
                if not len(ref):
                    continue
                #print(mov_stacks)
                for mov_stack in mov_stacks:
                    dvv_list = []
                    err_list = []
                    coh_list = []
                    data_dates=[]
                    try:
                        data = xr_get_ccf(station1, station2, components, filterid, mov_stack, taxis)
                    except FileNotFoundError as fullpath:
                        logger.error("FILE DOES NOT EXIST: %s, skipping" % fullpath)
                        continue
                    logger.debug("Processing %s:%s f%i m%s %s" % (station1, station2, filterid, mov_stack, components))

                    to_search = pd.to_datetime(days)
                    to_search = to_search.append(pd.DatetimeIndex([to_search[-1]+pd.Timedelta("1d"),]))
                    data = data[data.index.floor('d').isin(to_search)]
                    data = data.dropna()

                    cur = data#.CCF.values
                    if wct_norm:
                        new_waveform = (cur/cur.max()) 
                    else:
                        new_waveform = cur

                    for date, row in new_waveform.iterrows(): 
                        try:
                            WXamp, WXspec, WXangle, Wcoh, WXdt, freqs, coi = xwt(ori_waveform, row.values, goal_sampling_rate, int(ns), int(nt), int(vpo), freqmin, freqmax, int(nptsfreq), wavelet_type)
                            dvv, err, wf = compute_wct_dvv(freqs, taxis, WXamp, Wcoh, WXdt, lag_min=int(lag_min), coda_cycles=coda_cycles, mincoh=mincoh, maxdt=maxdt, min_nonzero=min_nonzero, freqmin=freqmin, freqmax=freqmax)
                            coh = get_avgcoh(freqs, taxis, Wcoh, freqmin, freqmax, lag_min=int(lag_min), coda_cycles=coda_cycles)
                            dvv_list.append(dvv)
                            err_list.append(err)
                            coh_list.append(coh)
                            data_dates.append(date)
                        except Exception as e:
                            logger.error(f"Error processing date {date}: {str(e)}")
                            # Add dummy values or skip this date
                            continue
                    if len(dvv_list) > 0:#1:
                        inx = np.where((freqs >= freqmin) & (freqs <= freqmax))

                        dvv_df = pd.DataFrame(dvv_list, columns=freqs[inx], index=data_dates)
                        err_df = pd.DataFrame(err_list, columns=freqs[inx], index=data_dates)
                        coh_df = pd.DataFrame(coh_list, columns=freqs[inx], index=data_dates)

                        # Saving
                        xr_save_wct(station1, station2, components, filterid, mov_stack, taxis, dvv_df, err_df, coh_df)

                        del dvv_df, err_df, coh_df
                    del cur

        massive_update_job(db, jobs, "D")

    logger.info('*** Finished: Compute WCT ***')
    
def get_fresh_db_connection():
    try:
        # Try to rollback any pending transactions
        try:
            db.rollback()
        except:
            pass
        
        # Dispose of the old connection
        try:
            db.dispose()
        except:
            pass
            
        # Create a new connection
        new_db = connect()
        return new_db
    except Exception as e:
        logger.error(f"Error getting fresh DB connection: {str(e)}")
        time.sleep(5)
        return connect()  # Last resort - try one more time


def get_wct_next_job2(session, flag='T', jobtype='DTT'):
    """
    Get the next DTT Jobs in the database for a randomly selected station pair,
    with flag=`flag` and jobtype=`jobtype`. This function atomically selects jobs
    and sets their flag to "I"n progress, preventing other processes from taking
    the same jobs.
    
    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object
    :type jobtype: str
    :param jobtype: Job type (e.g., 'WCT', 'STACK', 'DTT')
    :type flag: str
    :param flag: Status of the Job: "T"odo, "I"n Progress, "D"one.
    :rtype: list
    :returns: List of Job objects that were claimed
    """
    from sqlalchemy.sql.expression import func
    
    # Choose appropriate random function based on database type
    if read_db_inifile().tech == 2:  # MySQL
        rand = func.rand
    else:  # SQLite or PostgreSQL
        rand = func.random
        
    max_retries = 5
    retry_count = 0
    
    while retry_count < max_retries:
        try:
            # Start a transaction
            session.begin_nested()
            
            # First, find a random pair that has jobs to do
            random_pair_query = session.query(Job.pair).filter(Job.flag == flag).\
                filter(Job.jobtype == jobtype).filter(Job.day != 'REF').\
                order_by(rand()).limit(1)
                
            random_pair = random_pair_query.scalar()
            
            if random_pair is None:
                # No jobs found
                session.commit()
                return []
                
            # Get all jobs for this pair with FOR UPDATE lock
            jobs = session.query(Job).filter(Job.flag == flag).\
                filter(Job.jobtype == jobtype).filter(Job.day != 'REF').\
                filter(Job.pair == random_pair).\
                with_for_update().all()
                
            if not jobs:
                # Another process might have taken these jobs already
                session.rollback()
                retry_count += 1
                time.sleep(np.random.random())
                continue
                
            # Update all jobs to "I"n progress at once
            job_ids = [job.ref for job in jobs]
            session.query(Job).filter(Job.ref.in_(job_ids)).\
                update({Job.flag: "I"}, synchronize_session=False)
                
            # Commit the transaction
            session.commit()
            return jobs
            
        except Exception as e:
            session.rollback()
            if "deadlock" in str(e).lower():
                retry_count += 1
                delay = 1 * (2 ** retry_count) + np.random.random()
                print(f"Deadlock in get_dtt_next_job, retry {retry_count}/{max_retries} after {delay:.2f}s")
                time.sleep(delay)
            else:
                print(f"Error in get_dtt_next_job: {str(e)}")
                traceback.print_exc()
                return []
                
    # If we've exhausted all retries
    print(f"Failed to get next job after {max_retries} retries")
    return []


def get_wct_next_job3(session, flag='T', jobtype='WCT'):
    """
    Get the next WCT Jobs without using FOR UPDATE locks to prevent deadlocks.
    Follows the original pattern but with safer locking strategy.
    
    Parameters:
    ----------
    session : SQLAlchemy session
        Database session
    flag : str
        Job status flag (default 'T' for Todo)
    jobtype : str
        Type of job to retrieve (default 'WCT')
        
    Returns:
    -------
    list
        List of Job objects that were successfully claimed, or empty list if none
    """
    from sqlalchemy.sql.expression import func
    import traceback
    import time
    import numpy as np
    
    # Choose appropriate random function based on database type
    if read_db_inifile().tech == 2:  # MySQL
        rand = func.rand
    else:  # SQLite or PostgreSQL
        rand = func.random
    
    max_retries = 5
    for retry in range(max_retries):
        try:
            # Make sure we start fresh
            session.rollback()
            
            # Find a random pair (similar to original, but in a separate transaction)
            random_pair_query = session.query(Job.pair).filter(Job.flag == flag).\
                filter(Job.jobtype == jobtype).filter(Job.day != 'REF').\
                order_by(rand()).limit(1)
            
            random_pair = random_pair_query.scalar()
            
            if not random_pair:
                return []  # No jobs available
            
            # Get the jobs for this pair (WITHOUT using with_for_update())
            jobs_query = session.query(Job).filter(Job.flag == flag).\
                filter(Job.jobtype == jobtype).filter(Job.day != 'REF').\
                filter(Job.pair == random_pair)
            
            jobs = jobs_query.all()
            
            if not jobs:
                # No jobs for this pair (someone else might have taken them)
                time.sleep(0.1)
                continue
            
            # Prepare update data
            job_ids = [job.ref for job in jobs]
            mappings = [{'ref': job_id, 'flag': "I"} for job_id in job_ids]
            
            # Try to update with retries (similar to original)
            updated = False
            update_retries = 3
            
            for update_retry in range(update_retries):
                try:
                    # Use an update query that only affects rows still marked as 'T'
                    update_count = session.query(Job).\
                        filter(Job.ref.in_(job_ids)).\
                        filter(Job.flag == flag).\
                        update({Job.flag: "I"}, synchronize_session=False)
                    
                    if update_count == 0:
                        # All jobs were already taken
                        break
                    
                    session.commit()
                    updated = True
                    break
                except Exception as e:
                    session.rollback()
                    
                    if update_retry == update_retries - 1:
                        traceback.print_exc()
                    
                    # Determine wait time based on error type
                    if "deadlock" in str(e).lower() or "lock wait timeout" in str(e).lower():
                        wait_time = 0.5 + (update_retry * 0.5) + np.random.random()
                    else:
                        wait_time = 0.1 + np.random.random() * 0.2
                    
                    time.sleep(wait_time)
            
            if updated:
                # Fetch the jobs that were actually updated to "I"
                updated_jobs = session.query(Job).filter(Job.ref.in_(job_ids)).\
                    filter(Job.flag == "I").all()
                
                return updated_jobs
            else:
                # Couldn't update any jobs, try another pair
                continue
                
        except Exception as e:
            # If anything fails, retry with a different pair
            print(f"Error in get_wct_next_job (attempt {retry+1}/{max_retries}): {str(e)}")
            traceback.print_exc()
            
            wait_time = (1 + retry) + np.random.random()
            time.sleep(wait_time)
    
    # If we've exhausted all retries
    return []
    

def get_wct_next_job(session, flag='T', jobtype='WCT'):
    """
    Get the next WCT Jobs with improved transaction handling.
    
    Parameters:
    ----------
    session : SQLAlchemy session
        Database session
    flag : str
        Job status flag (default 'T' for Todo)
    jobtype : str
        Type of job to retrieve (default 'WCT')
        
    Returns:
    -------
    list
        List of Job objects that were successfully claimed, or empty list if none
    """
    from sqlalchemy.sql.expression import func
    import traceback
    import time
    import numpy as np
    
    # Choose appropriate random function based on database type
    if read_db_inifile().tech == 2:  # MySQL
        rand = func.rand
    else:  # SQLite or PostgreSQL
        rand = func.random
    
    # Always ensure we start with a clean session
    try:
        session.rollback()
    except:
        pass
    
    max_retries = 5
    for retry in range(max_retries):
        try:
            # First, find a random pair (in a separate, short transaction)
            random_pair = None
            try:
                random_pair_query = session.query(Job.pair).filter(Job.flag == flag).\
                    filter(Job.jobtype == jobtype).filter(Job.day != 'REF').\
                    order_by(rand()).limit(1)
                
                random_pair = random_pair_query.scalar()
                
                # Explicitly commit this query
                session.commit()
            except Exception as e:
                # If anything goes wrong, roll back
                try:
                    session.rollback()
                except:
                    pass
                raise
            
            if not random_pair:
                return []  # No jobs available
            
            # Get job IDs for this pair (in a separate, short transaction)
            job_ids = []
            try:
                job_id_query = session.query(Job.ref).filter(Job.flag == flag).\
                    filter(Job.jobtype == jobtype).filter(Job.day != 'REF').\
                    filter(Job.pair == random_pair)
                
                job_ids = [j[0] for j in job_id_query.all()]
                
                # Commit this query
                session.commit()
            except Exception as e:
                # If anything goes wrong, roll back
                try:
                    session.rollback()
                except:
                    pass
                raise
            
            if not job_ids:
                # No jobs for this pair (they were taken by another process)
                time.sleep(0.1)
                continue
            
            # Update job flags (in a separate transaction)
            update_success = False
            try:
                # Start a new transaction
                session.begin()
                
                # Update the flags (only for rows still with flag='T')
                update_count = session.query(Job).\
                    filter(Job.ref.in_(job_ids)).\
                    filter(Job.flag == flag).\
                    update({Job.flag: "I"}, synchronize_session=False)
                
                if update_count > 0:
                    # Commit the transaction if we updated some rows
                    session.commit()
                    update_success = True
                else:
                    # Roll back if no rows were updated (they were taken by another process)
                    session.rollback()
            except Exception as e:
                # If anything goes wrong, roll back
                try:
                    session.rollback()
                except:
                    pass
                raise
            
            if not update_success:
                # No jobs were updated, try another pair
                continue
            
            # Fetch the updated jobs (in a final, separate transaction)
            try:
                updated_jobs = session.query(Job).filter(Job.ref.in_(job_ids)).\
                    filter(Job.flag == "I").all()
                
                # Commit this query
                session.commit()
                
                if updated_jobs:
                    return updated_jobs
            except Exception as e:
                # If anything goes wrong, roll back
                try:
                    session.rollback()
                except:
                    pass
                raise
                
        except Exception as e:
            # Always ensure transaction is rolled back
            try:
                session.rollback()
            except:
                pass
            
            print(f"Error in get_wct_next_job (attempt {retry+1}/{max_retries}): {str(e)}")
            traceback.print_exc()
            
            wait_time = (2 ** retry) + np.random.random()
            time.sleep(wait_time)
    
    # Final cleanup to ensure no hanging transactions
    try:
        session.rollback()
    except:
        pass
    
    # If we've exhausted all retries
    return []

    
def main(loglevel="INFO"):
    # Reconfigure logger to show the pid number in log records
    global logger
    logger = get_logger('msnoise.compute_wct_child', loglevel, with_pid=True)
    logger.info('*** Starting: Compute WCT ***')

    # Function to safely get a fresh DB connection with transaction cleanup
    def get_fresh_db_connection():
        try:
            # If we have an existing DB object, try to clean it up
            try:
                if 'db' in globals():
                    try:
                        db.rollback()  # Ensure no open transactions
                    except:
                        pass
                    try:
                        db.close()  # Close the connection
                    except:
                        pass
                    try:
                        db.dispose()  # Dispose the engine
                    except:
                        pass
            except:
                pass
                
            # Create a new connection
            new_db = connect()
            
            # Test the connection with a simple query
            _ = get_params(new_db)
            
            return new_db
        except Exception as e:
            logger.error(f"Error getting fresh DB connection: {str(e)}")
            time.sleep(5)  # Wait before retrying
            return get_fresh_db_connection()  # Recursively retry

    # Initial connection
    db = get_fresh_db_connection()
    params = get_params(db)
    taxis = get_t_axis(db)
    
    # Get configuration parameters
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

    logger.debug('Ready to compute')
    
    # Get filters only once at the beginning
    try:
        filters = get_filters(db, all=False)
    except Exception as e:
        logger.error(f"Error getting filters: {str(e)}")
        db = get_fresh_db_connection()
        filters = get_filters(db, all=False)
        
    # Add small random delay to avoid all processes hitting DB at once
    time.sleep(np.random.random() * 5)

    # Track time to periodically refresh connection
    last_connection_refresh = time.time()
    connection_refresh_interval = 600  # Refresh connection every 10 minutes
    
    jobs_processed = 0
    
    while True:
        try:
            # Check if we need to refresh the connection
            if time.time() - last_connection_refresh > connection_refresh_interval:
                logger.info("Refreshing database connection")
                db = get_fresh_db_connection()
                last_connection_refresh = time.time()
            
            # Check if jobs exist
            if not is_dtt_next_job(db, flag='T', jobtype='WCT'):
                logger.info("No more jobs to process")
                break
            
            # Get jobs to process
            jobs = get_wct_next_job(db, flag='T', jobtype='WCT')
            
            if not jobs:
                # No jobs were claimed, wait briefly and try again
                time.sleep(0.5 + np.random.random())
                continue
            
            pair = jobs[0].pair
            logger.info(f"Processing {len(jobs)} jobs for pair: {pair}")
            jobs_processed += len(jobs)
              

                
            if not len(jobs):
                # edge case, should only occur when is_next returns true, but
                # get_next receives no jobs (heavily parallelised calls).
                time.sleep(np.random.random())
                continue
                
            pair = jobs[0].pair
            logger.info(f"Processing pair: {pair}")
            refs, days = zip(*[[job.ref, job.day] for job in jobs])

            logger.info(
                "There are WCT jobs for some days to recompute for %s" % pair)
                
            # Process each filter
            for f in filters:
                try:
                    filterid = int(f.ref)
                    freqmin = f.low
                    freqmax = f.high
                    
                    station1, station2 = pair.split(":")
                    if station1 == station2:
                        components_to_compute = params.components_to_compute_single_station
                    else:
                        components_to_compute = params.components_to_compute

                    # Process each component
                    for components in components_to_compute:
                        try:
                            # Get reference waveform
                            ref = xr_get_ref(station1, station2, components, filterid, taxis, ignore_network=True)
                            ref = ref.CCF.values
                            if wct_norm:
                                ori_waveform = (ref/ref.max()) 
                            else:
                                ori_waveform = ref
                        except FileNotFoundError as fullpath:
                            logger.error("FILE DOES NOT EXIST: %s, skipping" % fullpath)
                            continue
                        if not len(ref):
                            continue
                            
                        # Process each moving stack
                        for mov_stack in mov_stacks:
                            dvv_list = []
                            err_list = []
                            coh_list = []
                            data_dates = []
                            
                            # Get the data for this component and moving stack
                            try:
                                data = xr_get_ccf(station1, station2, components, filterid, mov_stack, taxis)
                            except FileNotFoundError as fullpath:
                                logger.error("FILE DOES NOT EXIST: %s, skipping" % fullpath)
                                continue
                                
                            logger.debug("Processing %s:%s f%i m%s %s" % (station1, station2, filterid, mov_stack, components))

                            # Filter data by dates
                            to_search = pd.to_datetime(days)
                            to_search = to_search.append(pd.DatetimeIndex([to_search[-1]+pd.Timedelta("1d"),]))
                            data = data[data.index.floor('d').isin(to_search)]
                            data = data.dropna()

                            cur = data
                            if wct_norm:
                                new_waveform = (cur/cur.max()) 
                            else:
                                new_waveform = cur

                            # Process each date
                            for date, row in new_waveform.iterrows():
                                try:
                                    # Calculate wavelet coherence transform
                                    WXamp, WXspec, WXangle, Wcoh, WXdt, freqs, coi = xwt(
                                        ori_waveform, row.values, goal_sampling_rate, 
                                        int(ns), float(nt), int(vpo), 
                                        freqmin, freqmax, int(nptsfreq), wavelet_type
                                    )
                                    
                                    # Compute dv/v values with improved error handling
                                    try:
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
                                        
                                        # Store results
                                        dvv_list.append(dvv)
                                        err_list.append(err)
                                        coh_list.append(coh)
                                        data_dates.append(date)
                                        
                                    except Exception as e:
                                        logger.error(f"Error in compute_wct_dvv for date {date}: {str(e)}")
                                        # Continue with next date
                                        continue
                                        
                                except Exception as e:
                                    logger.error(f"Error in xwt for date {date}: {str(e)}")
                                    # Continue with next date
                                    continue

                            # Save results if we have any
                            if len(dvv_list) > 0:
                                try:
                                    inx = np.where((freqs >= freqmin) & (freqs <= freqmax))

                                    dvv_df = pd.DataFrame(dvv_list, columns=freqs[inx], index=data_dates)
                                    err_df = pd.DataFrame(err_list, columns=freqs[inx], index=data_dates)
                                    coh_df = pd.DataFrame(coh_list, columns=freqs[inx], index=data_dates)

                                    # Saving
                                    xr_save_wct(station1, station2, components, filterid, mov_stack, taxis, dvv_df, err_df, coh_df)

                                    del dvv_df, err_df, coh_df
                                except Exception as e:
                                    logger.error(f"Error saving results: {str(e)}")
                                    
                            del cur
                
                except Exception as e:
                    logger.error(f"Error processing filter {f.ref}: {str(e)}")
                    # Try to reconnect if it's a database error
                    if "MySQL" in str(e) or "Connection" in str(e):
                        db = get_fresh_db_connection()
                    continue
                    
                        
            # Mark jobs as done
            try:
                massive_update_job(db, jobs, "D")
                logger.info(f"Successfully marked {len(jobs)} jobs as done")
            except Exception as e:
                logger.error(f"Error updating job status: {str(e)}")
                
                # Always roll back if there's an error
                try:
                    db.rollback()
                except:
                    pass
                    
                # Try with a fresh connection
                db = get_fresh_db_connection()
                try:
                    massive_update_job(db, jobs, "D")
                except Exception as e2:
                    logger.error(f"Failed to update job status after reconnection: {str(e2)}")
                    # Ensure we roll back
                    try:
                        db.rollback()
                    except:
                        pass
                
        except Exception as e:
            logger.error(f"Error in main job processing loop: {str(e)}")
            
            # Always roll back if there's an error
            try:
                db.rollback()
            except:
                pass
                
            # Handle different types of errors
            if "Can't reconnect until invalid transaction" in str(e):
                logger.warning("Transaction issue detected, getting fresh connection...")
                db = get_fresh_db_connection()
                time.sleep(2)
            elif "MySQL server has gone away" in str(e) or "Connection" in str(e):
                logger.warning("Database connection issue, reconnecting...")
                db = get_fresh_db_connection()
                time.sleep(2)
            elif "deadlock" in str(e).lower() or "lock wait timeout" in str(e).lower():
                logger.warning("Database contention detected, waiting...")
                
                # Try to clean up and then wait
                try:
                    db.rollback()
                except:
                    pass
                    
                time.sleep(5 + np.random.random() * 5)
            else:
                # Brief delay for other errors
                logger.error(f"Unexpected error: {str(e)}")
                time.sleep(1 + np.random.random())

    # Final cleanup
    try:
        db.rollback()  # Ensure no open transactions
    except:
        pass
    try:
        db.close()  # Close connection
    except:
        pass

    logger.info(f'*** Finished: Compute WCT - Processed {jobs_processed} jobs ***')
    
    
    
def main3(loglevel="INFO"):
    # Reconfigure logger to show the pid number in log records
    global logger
    logger = get_logger('msnoise.compute_wct_child', loglevel,
                        with_pid=True)
    logger.info('*** Starting: Compute WCT ***')

    # Function to establish a fresh database connection
    def get_fresh_db_connection():
        try:
            db = connect()
            # Test the connection
            _ = get_params(db)
            return db
        except Exception as e:
            logger.error(f"Error connecting to database: {str(e)}")
            time.sleep(5)  # Wait before retrying
            return get_fresh_db_connection()  # Recursively retry

    # Initial connection
    db = connect()
    params = get_params(db)
    taxis = get_t_axis(db)
    
    # Get configuration parameters
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

    logger.debug('Ready to compute')
    
    # Get filters only once at the beginning
    try:
        filters = get_filters(db, all=False)
    except Exception as e:
        logger.error(f"Error getting filters: {str(e)}")
        db = get_fresh_db_connection()
        filters = get_filters(db, all=False)
        
    time.sleep(np.random.random() * 5)

    # Track time to periodically refresh connection
    last_connection_refresh = time.time()
    connection_refresh_interval = 600  # Refresh connection every 10 minutes
    
    while True:
        # Check if we need to refresh the connection
        if time.time() - last_connection_refresh > connection_refresh_interval:
            logger.info("Refreshing database connection")
            try:
                db.dispose()  # Close the current connection
            except:
                pass
            db = get_fresh_db_connection()
            last_connection_refresh = time.time()
        
        # Try to get the next job with error handling
        try:
            next_job_exists = is_dtt_next_job(db, flag='T', jobtype='WCT')
            if not next_job_exists:
                logger.info("No more jobs to process")
                break
                
            jobs = get_wct_next_job(db, flag='T', jobtype='WCT')
              
#            if jobs:
#                # Print job structure and attributes
#                print(f"Job object structure: {vars(jobs[0])}")
#                print(f"Job attributes: {dir(jobs[0])}")
#                
#                # Print all fields for each job 
#                for job in jobs:
#                    print(f"Job ID: {job.ref}, Pair: {job.pair}, Day: {job.day}, Type: {job.jobtype}, Flag: {job.flag}")
#                    
#                # Print DB query info 
#                print(f"Total jobs: {len(jobs)}")
#                
#                # Print relationship with other entities (example)
#                station1, station2 = jobs[0].pair.split(":")
#                print(f"Related stations: {station1}, {station2}")
                
            if not len(jobs):
                # edge case, should only occur when is_next returns true, but
                # get_next receives no jobs (heavily parallelised calls).
                time.sleep(np.random.random())
                continue
                
            pair = jobs[0].pair
            logger.info(f"Processing pair: {pair}")
            refs, days = zip(*[[job.ref, job.day] for job in jobs])

            logger.info(
                "There are WCT jobs for some days to recompute for %s" % pair)
                
            # Process each filter
            for f in filters:
                try:
                    filterid = int(f.ref)
                    freqmin = f.low
                    freqmax = f.high
                    
                    station1, station2 = pair.split(":")
                    if station1 == station2:
                        components_to_compute = params.components_to_compute_single_station
                    else:
                        components_to_compute = params.components_to_compute

                    # Process each component
                    for components in components_to_compute:
                        try:
                            # Get reference waveform
                            ref = xr_get_ref(station1, station2, components, filterid, taxis, ignore_network=True)
                            ref = ref.CCF.values
                            if wct_norm:
                                ori_waveform = (ref/ref.max()) 
                            else:
                                ori_waveform = ref
                        except FileNotFoundError as fullpath:
                            logger.error("FILE DOES NOT EXIST: %s, skipping" % fullpath)
                            continue
                        if not len(ref):
                            continue
                            
                        # Process each moving stack
                        for mov_stack in mov_stacks:
                            dvv_list = []
                            err_list = []
                            coh_list = []
                            data_dates = []
                            
                            # Get the data for this component and moving stack
                            try:
                                data = xr_get_ccf(station1, station2, components, filterid, mov_stack, taxis)
                            except FileNotFoundError as fullpath:
                                logger.error("FILE DOES NOT EXIST: %s, skipping" % fullpath)
                                continue
                                
                            logger.debug("Processing %s:%s f%i m%s %s" % (station1, station2, filterid, mov_stack, components))

                            # Filter data by dates
                            to_search = pd.to_datetime(days)
                            to_search = to_search.append(pd.DatetimeIndex([to_search[-1]+pd.Timedelta("1d"),]))
                            data = data[data.index.floor('d').isin(to_search)]
                            data = data.dropna()

                            cur = data
                            if wct_norm:
                                new_waveform = (cur/cur.max()) 
                            else:
                                new_waveform = cur

                            # Process each date
                            for date, row in new_waveform.iterrows():
                                try:
                                    # Calculate wavelet coherence transform
                                    WXamp, WXspec, WXangle, Wcoh, WXdt, freqs, coi = xwt(
                                        ori_waveform, row.values, goal_sampling_rate, 
                                        int(ns), float(nt), int(vpo), 
                                        freqmin, freqmax, int(nptsfreq), wavelet_type
                                    )
                                    
                                    # Compute dv/v values with improved error handling
                                    try:
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
                                        
                                        # Store results
                                        dvv_list.append(dvv)
                                        err_list.append(err)
                                        coh_list.append(coh)
                                        data_dates.append(date)
                                        
                                    except Exception as e:
                                        logger.error(f"Error in compute_wct_dvv for date {date}: {str(e)}")
                                        # Continue with next date
                                        continue
                                        
                                except Exception as e:
                                    logger.error(f"Error in xwt for date {date}: {str(e)}")
                                    # Continue with next date
                                    continue

                            # Save results if we have any
                            if len(dvv_list) > 0:
                                try:
                                    inx = np.where((freqs >= freqmin) & (freqs <= freqmax))

                                    dvv_df = pd.DataFrame(dvv_list, columns=freqs[inx], index=data_dates)
                                    err_df = pd.DataFrame(err_list, columns=freqs[inx], index=data_dates)
                                    coh_df = pd.DataFrame(coh_list, columns=freqs[inx], index=data_dates)

                                    # Saving
                                    xr_save_wct(station1, station2, components, filterid, mov_stack, taxis, dvv_df, err_df, coh_df)

                                    del dvv_df, err_df, coh_df
                                except Exception as e:
                                    logger.error(f"Error saving results: {str(e)}")
                                    
                            del cur
                
                except Exception as e:
                    logger.error(f"Error processing filter {f.ref}: {str(e)}")
                    # Try to reconnect if it's a database error
                    if "MySQL" in str(e) or "Connection" in str(e):
                        db = get_fresh_db_connection()
                    continue
                    
                        
            # Mark jobs as done
            try:
                massive_update_job(db, jobs, "D")
                logger.info(f"Successfully marked {len(jobs)} jobs as done")
            except Exception as e:
                logger.error(f"Error updating job status: {str(e)}")
                # Try to reconnect and update again
                db = get_fresh_db_connection()
                try:
                    massive_update_job(db, jobs, "D")
                except Exception as e2:
                    logger.error(f"Failed to update job status after reconnection: {str(e2)}")
                
        except Exception as e:
            logger.error(f"Error in main job processing loop: {str(e)}")
            if "MySQL server has gone away" in str(e) or "Connection reset by peer" in str(e):
                logger.warning("Database connection lost, reconnecting...")
                db = get_fresh_db_connection()
                time.sleep(10)  # Give some time before continuing
            elif "Deadlock" in str(e):
                logger.warning("Database deadlock detected, retrying after delay...")
                time.sleep(20)  # Longer delay for deadlocks
            else:
                # For other errors, log and continue
                logger.error(f"Unexpected error: {str(e)}")
                time.sleep(5)

    logger.info('*** Finished: Compute WCT ***')
    
if __name__ == "__main__":
    main()
