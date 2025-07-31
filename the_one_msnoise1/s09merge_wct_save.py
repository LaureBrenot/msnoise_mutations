#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
s09merge_wct.py

This module merges daily WCT files for each station pair and parameter set
into a consolidated file with hierarchical output directories:
WCT_MERGED/[filter_id]/[mov_stack]/[component]/[station_pair].nc

This script is designed to be called by MSNoise through the command line
interface using the 'msnoise cc dvv merge_wct' command.

Usage:
msnoise cc dvv merge_wct
msnoise cc dvv merge_wct --wct-dir /path/to/wct --output-dir /path/to/output
msnoise -t 8 cc dvv merge_wct
"""

import os
import glob
import numpy as np
import pandas as pd
import xarray as xr
import datetime
from .api import (get_logger, xr_create_or_open, xr_save_and_close, 
                 xr_insert_or_update, validate_stack_data, xr_save_wct)

def parse_pair_name(pair_folder):
    return pair_folder.split('_', 1)


def get_filter_directories(wct_dir):
    if not os.path.exists(wct_dir):
        return []
    
    # Look for filter directories (numbered folders like "01", "02", etc.)
    filters = []
    for item in os.listdir(wct_dir):
        if os.path.isdir(os.path.join(wct_dir, item)) and item.isdigit():
            filters.append(item)
    
    return sorted(filters)


def get_movstack_directories(wct_dir, filter_dir):
    filter_path = os.path.join(wct_dir, filter_dir)
    movstacks = []
    for item in os.listdir(filter_path):
        if os.path.isdir(os.path.join(filter_path, item)):
            movstacks.append(item)
    
    return sorted(movstacks)


def get_component_directories(wct_dir, filter_dir, movstack_dir):
    path = os.path.join(wct_dir, filter_dir, movstack_dir)
    components = []
    for item in os.listdir(path):
        if os.path.isdir(os.path.join(path, item)):
            components.append(item)
    
    return sorted(components)


def get_station_pair_directories(wct_dir, filter_dir, movstack_dir, component_dir):
    path = os.path.join(wct_dir, filter_dir, movstack_dir, component_dir)
    pairs = []
    for item in os.listdir(path):
        if os.path.isdir(os.path.join(path, item)) and '_' in item:
            pairs.append(item)
    
    return sorted(pairs)


def merge_daily_files(wct_dir, output_dir, filter_dir, movstack_dir, component_dir, pair_dir, logger):
    """
    Merge daily files for a specific station pair and parameter set
    
    :type wct_dir: str
    :param wct_dir: Path to the WCT directory
    :type output_dir: str
    :param output_dir: Path to the output directory
    :type filter_dir: str
    :param filter_dir: Name of the filter directory
    :type movstack_dir: str
    :param movstack_dir: Name of the moving stack directory
    :type component_dir: str
    :param component_dir: Name of the component directory
    :type pair_dir: str
    :param pair_dir: Name of the station pair directory
    :type logger: logging.Logger
    :param logger: Logger object
    
    :rtype: str or None
    :returns: Path to the merged file if successful, None otherwise
    """
    full_path = os.path.join(wct_dir, filter_dir, movstack_dir, component_dir, pair_dir)
    station1, station2 = parse_pair_name(pair_dir)
    
    # Find all date files (*.npz)
    files = glob.glob(os.path.join(full_path, "*.npz"))
    
    if not files:
        # Silently skip if no files found - don't log as it clutters output
        logger.warning(f"No .npz files found in {full_path}")
        return None
    
    # Initialize lists to store data
    dates = []
    dvv_values = []
    err_values = []
    coh_values = []
    
    # Pre-define variables that will be used later
    filterid = int(filter_dir)
    component = component_dir
    
    # Load each file and extract data
    for file in sorted(files):
        try:
            data = np.load(file, allow_pickle=True)
            
            # Try to extract date from filename first (more reliable)
            filename = os.path.basename(file)
            date_str = filename.split('.')[0]  # Remove .npz extension
            
            try:
                date = pd.to_datetime(date_str)
            except ValueError:
                # If filename doesn't contain a valid date, try to get it from the data
                if 'date' in data:
                    date_val = data['date']
                    # Safely convert various date formats
                    if isinstance(date_val, np.ndarray):
                        # If it's an array of size 1, extract the item
                        if date_val.size == 1:
                            try:
                                date = pd.to_datetime(date_val.item())
                            except:
                                date = pd.to_datetime(str(date_val.item()))
                        else:
                            # If it's a larger array, use the first element
                            date = pd.to_datetime(date_val[0])
                    else:
                        # If it's already a scalar
                        date = pd.to_datetime(date_val)
                else:
                    # Silently skip problematic files
                    continue
            
            # Extract data values, ensuring they're in the correct format
            if 'dvv' in data:
                dvv = data['dvv']
                if isinstance(dvv, np.ndarray) and dvv.ndim > 0:
                    dvv_values.append(dvv)
                else:
                    continue
            else:
                continue
                
            if 'err' in data:
                err = data['err']
                if isinstance(err, np.ndarray) and err.ndim > 0:
                    err_values.append(err)
                else:
                    # Use zeros array of same shape as dvv if err is invalid
                    err_values.append(np.zeros_like(dvv))
            else:
                # Use zeros array of same shape as dvv if err is missing
                err_values.append(np.zeros_like(dvv))
                
            if 'coh' in data:
                coh = data['coh']
                if isinstance(coh, np.ndarray) and coh.ndim > 0:
                    coh_values.append(coh)
                else:
                    # Use ones array of same shape as dvv if coh is invalid
                    coh_values.append(np.ones_like(dvv))
            else:
                # Use ones array of same shape as dvv if coh is missing
                coh_values.append(np.ones_like(dvv))
            
            dates.append(date)
            
            # Get additional metadata from first file
            if len(dates) == 1:
                # Get taxis from first file
                if 'taxis' in data:
                    taxis = data['taxis']
                else:
                    taxis = np.linspace(-120, 120, 241)  # Default taxis
                    
                # Override metadata if available in the file
                if 'filterid' in data:
                    try:
                        filterid_val = data['filterid']
                        if isinstance(filterid_val, np.ndarray) and filterid_val.size == 1:
                            filterid = int(filterid_val.item())
                        else:
                            filterid = int(filterid_val)
                    except:
                        filterid = int(filter_dir)
                        
                if 'component' in data:
                    try:
                        component_val = data['component']
                        if isinstance(component_val, np.ndarray) and component_val.size == 1:
                            component = str(component_val.item())
                        else:
                            component = str(component_val)
                    except:
                        component = component_dir
        
        except Exception as e:
            # Skip silently on errors to avoid cluttering logs
            continue
    
    if not dates or len(dates) != len(dvv_values):
        return None
    
    # Get frequencies from the first file
    try:
        first_file = np.load(files[0], allow_pickle=True)
        if 'freqs' in first_file:
            freqs_data = first_file['freqs']
            if isinstance(freqs_data, np.ndarray) and freqs_data.size > 0:
                freqs = freqs_data
            else:
                freqs = np.linspace(0.1, 2.0, dvv_values[0].shape[0])
        else:
            freqs = np.linspace(0.1, 2.0, dvv_values[0].shape[0])
    except Exception:
        # Generate default frequencies based on data shape
        freqs = np.linspace(0.1, 2.0, dvv_values[0].shape[0])
    
    # Create DataFrames
    try:
        # Ensure all data arrays have the same shape
        common_shape = dvv_values[0].shape[0]
        valid_indices = []
        
        for i, (date, dvv, err, coh) in enumerate(zip(dates, dvv_values, err_values, coh_values)):
            if (dvv.shape[0] == common_shape and 
                err.shape[0] == common_shape and 
                coh.shape[0] == common_shape):
                valid_indices.append(i)
        
        if not valid_indices:
            return None
            
        # Filter to only use valid data points
        filtered_dates = [dates[i] for i in valid_indices]
        filtered_dvv = [dvv_values[i] for i in valid_indices]
        filtered_err = [err_values[i] for i in valid_indices]
        filtered_coh = [coh_values[i] for i in valid_indices]
        
        dvv_df = pd.DataFrame(filtered_dvv, index=filtered_dates, columns=freqs)
        err_df = pd.DataFrame(filtered_err, index=filtered_dates, columns=freqs)
        coh_df = pd.DataFrame(filtered_coh, index=filtered_dates, columns=freqs)
        
        # Sort by date
        dvv_df = dvv_df.sort_index()
        err_df = err_df.sort_index()
        coh_df = coh_df.sort_index()
    except Exception:
        return None
    
    # Create output directory matching the input hierarchy
    output_subdir = os.path.join(output_dir, filter_dir, movstack_dir, component_dir)
    os.makedirs(output_subdir, exist_ok=True)
    
    # Save to file - use the same station pair name for the output file
    filename = f"{pair_dir}.nc"
    output_path = os.path.join(output_subdir, filename)
    
    if '_' in movstack_dir:
        parts = movstack_dir.split('_')
        mov_stack = (parts[0], parts[1])
    else:
        mov_stack = (movstack_dir, movstack_dir)
        
    try:
        xr_save_wct(station1, station2, component, filterid, mov_stack, taxis, dvv_df, err_df, coh_df)
        return True  # or return the expected output path

    except Exception:
        return None


def process_combination(args):
    """
    Process a single combination of filter, movstack, component, and pair
    This function is designed to be used with multiprocessing

    :param args: Tuple containing (wct_dir, output_dir, filter_dir, movstack_dir, component_dir, pair_dir)
    :return: Tuple with (success, output_path)
    """
    # Unpack arguments
    wct_dir, output_dir, filter_dir, movstack_dir, component_dir, pair_dir = args
    
    # Create a logger that won't output anything (to avoid multiprocessing logging issues)
    class SilentLogger:
        def info(self, msg): pass
        def debug(self, msg): pass
        def warning(self, msg): pass
        def error(self, msg): pass
    
    logger = SilentLogger()
    
    # Process the combination
    result = merge_daily_files(wct_dir, output_dir, filter_dir, movstack_dir, component_dir, pair_dir, logger)
    
    return (result is not None, result)


def main(wct_dir='WCT', output_dir='WCT_MERGED', loglevel="INFO", num_processes=None, use_tqdm=True):
    """
    Main function to merge all daily WCT files
    
    :type wct_dir: str
    :param wct_dir: Directory containing WCT results [default: WCT]
    :type output_dir: str
    :param output_dir: Directory to save merged results [default: WCT_MERGED]
    :type loglevel: str
    :param loglevel: Logging level [default: INFO]
    :type num_processes: int
    :param num_processes: Number of parallel processes to use [default: None, which uses CPU count]
    :type use_tqdm: bool
    :param use_tqdm: Use tqdm progress bar (True) or logging-based progress (False) [default: True]
    
    :rtype: int
    :returns: 0 if successful, 1 if error
    """
    import multiprocessing as mp
    
    # Configure logger using MSNoise API
    logger = get_logger('msnoise.merge_wct', loglevel)
    logger.info(f"Starting merge process, reading from {wct_dir}, writing to {output_dir}")
    
    # Set the number of processes to use
    if num_processes is None:
        # Use 75% of available CPU cores, but at least 1
        num_processes = max(1, int(mp.cpu_count() * 0.75))
    logger.info(f"Using {num_processes} parallel processes")
    
    try:
        # Get list of filter directories
        filters = get_filter_directories(wct_dir)
        if not filters:
            logger.error(f"No filter directories found in {wct_dir}")
            return 1
        
        # Build the list of all work to do
        combinations = []
        for filter_dir in filters:
            movstacks = get_movstack_directories(wct_dir, filter_dir)
            for movstack_dir in movstacks:
                components = get_component_directories(wct_dir, filter_dir, movstack_dir)
                for component_dir in components:
                    pairs = get_station_pair_directories(wct_dir, filter_dir, movstack_dir, component_dir)
                    for pair_dir in pairs:
                        combinations.append((wct_dir, output_dir, filter_dir, movstack_dir, component_dir, pair_dir))
        
        total_combinations = len(combinations)
        logger.info(f"Found {len(filters)} filters, {total_combinations} total combinations to process")
        
        if total_combinations == 0:
            logger.warning("No data to merge")
            return 0
        
        # Start timing
        start_time = datetime.datetime.now()
        
        # Process combinations in parallel
        merged_count = 0
        failed_count = 0
        
        if num_processes > 1:
            # Use a Pool for parallel processing
            with mp.Pool(processes=num_processes) as pool:
                if use_tqdm:
                    try:
                        import tqdm
                        # Process combinations with a progress bar
                        results = list(tqdm.tqdm(
                            pool.imap(process_combination, combinations),
                            total=total_combinations,
                            desc="Merging files",
                            unit="combination"
                        ))
                    except ImportError:
                        # If tqdm is not available, fall back to logging
                        use_tqdm = False
                        logger.warning("tqdm not available, falling back to logging-based progress")
                
                if not use_tqdm:
                    # Use logging for progress
                    results = []
                    checkpoint_interval = max(1, total_combinations // 100)  # Log every 1%
                    
                    # Create a callback to log progress
                    def log_progress(result):
                        nonlocal merged_count, failed_count
                        if result[0]:  # If successful
                            merged_count += 1
                        else:
                            failed_count += 1
                        
                        current_count = merged_count + failed_count
                        if current_count % checkpoint_interval == 0 or current_count == total_combinations:
                            progress_pct = (current_count / total_combinations) * 100
                            elapsed = datetime.datetime.now() - start_time
                            
                            # Estimate time remaining
                            if current_count > 0:
                                avg_time_per_item = elapsed.total_seconds() / current_count
                                remaining_seconds = avg_time_per_item * (total_combinations - current_count)
                                remaining_time = datetime.timedelta(seconds=int(remaining_seconds))
                                eta = datetime.datetime.now() + remaining_time
                                logger.info(f"Progress: {current_count}/{total_combinations} ({progress_pct:.1f}%) - "
                                           f"ETA: {eta.strftime('%Y-%m-%d %H:%M:%S')} - "
                                           f"Elapsed: {str(elapsed).split('.')[0]} - "
                                           f"Success/Fail: {merged_count}/{failed_count}")
                    
                    # Start async processing with callback
                    async_results = []
                    for combination in combinations:
                        async_result = pool.apply_async(process_combination, (combination,), callback=log_progress)
                        async_results.append(async_result)
                    
                    # Wait for all processes to complete
                    for async_result in async_results:
                        async_result.wait()
                    
                    # Collect results
                    results = [async_result.get() for async_result in async_results]
                
                # Count if we didn't use the callback
                if use_tqdm:
                    merged_count = sum(1 for success, _ in results if success)
                    failed_count = total_combinations - merged_count
        else:
            # Process sequentially
            if use_tqdm:
                try:
                    import tqdm
                    # Process with a progress bar
                    results = []
                    for combination in tqdm.tqdm(combinations, desc="Merging files", unit="combination"):
                        result = process_combination(combination)
                        results.append(result)
                        if result[0]:  # If successful
                            merged_count += 1
                        else:
                            failed_count += 1
                except ImportError:
                    # If tqdm is not available, fall back to logging
                    use_tqdm = False
                    logger.warning("tqdm not available, falling back to logging-based progress")
            
            if not use_tqdm:
                # Use logging for progress in sequential mode
                results = []
                checkpoint_interval = max(1, total_combinations // 100)  # Log every 1%
                
                for i, combination in enumerate(combinations):
                    result = process_combination(combination)
                    results.append(result)
                    
                    if result[0]:  # If successful
                        merged_count += 1
                    else:
                        failed_count += 1
                    
                    if (i + 1) % checkpoint_interval == 0 or (i + 1) == total_combinations:
                        progress_pct = ((i + 1) / total_combinations) * 100
                        elapsed = datetime.datetime.now() - start_time
                        
                        # Estimate time remaining
                        if i > 0:
                            avg_time_per_item = elapsed.total_seconds() / (i + 1)
                            remaining_seconds = avg_time_per_item * (total_combinations - (i + 1))
                            remaining_time = datetime.timedelta(seconds=int(remaining_seconds))
                            eta = datetime.datetime.now() + remaining_time
                            logger.info(f"Progress: {i+1}/{total_combinations} ({progress_pct:.1f}%) - "
                                       f"ETA: {eta.strftime('%Y-%m-%d %H:%M:%S')} - "
                                       f"Elapsed: {str(elapsed).split('.')[0]} - "
                                       f"Success/Fail: {merged_count}/{failed_count}")
        
        # Calculate summary statistics
        elapsed = datetime.datetime.now() - start_time
        avg_time = elapsed.total_seconds() / total_combinations if total_combinations > 0 else 0
        
        logger.info(f"Merge process completed in {str(elapsed).split('.')[0]}")
        logger.info(f"Processed {total_combinations} combinations, successfully merged {merged_count} files, {failed_count} failed")
        logger.info(f"Average processing time: {avg_time:.2f} seconds per combination")
        
        return 0
    except Exception as e:
        logger.error(f"Error in merge process: {str(e)}")
        import traceback
        logger.debug(traceback.format_exc())
        return 1


def dvv_merge_wct(ctx, wct_dir, output_dir):
    """
    MSNoise command function to merge WCT files
    
    :type ctx: click.Context
    :param ctx: Click context containing MSNoise configuration
    :type wct_dir: str
    :param wct_dir: Directory containing WCT results
    :type output_dir: str
    :param output_dir: Directory to save merged results
    """
    loglevel = ctx.obj['MSNOISE_verbosity']
    threads = ctx.obj['MSNOISE_threads']
    
    # Determine if we're running in a batch environment
    # Check for common batch environment variables
    import os
    batch_environment = any(env in os.environ for env in 
                            ['SLURM_JOB_ID', 'PBS_JOBID', 'SGE_TASK_ID', 'LSB_JOBID'])
    
    # Run with specified number of processes, disabling tqdm if in batch mode
    main(wct_dir=wct_dir, output_dir=output_dir, loglevel=loglevel, 
         num_processes=threads, use_tqdm=not batch_environment)


if __name__ == "__main__":
    import sys
    sys.exit(main())