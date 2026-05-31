#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
merge_wct_daily_files.py

This script merges all daily WCT files for each station pair and parameter set
into a consolidated file with hierarchical output directories:
WCT_MERGED/[filter_id]/[mov_stack]/[component]/[station_pair].nc

Usage:
    python merge_wct_daily_files.py [--wct_dir WCT_DIR] [--output_dir OUTPUT_DIR] [--verbose]

Options:
    --wct_dir WCT_DIR      Directory containing WCT results [default: WCT]
    --output_dir OUTPUT_DIR Directory to save merged results [default: WCT_MERGED]
    --verbose              Enable verbose output
"""

import os
import sys
import glob
import argparse
import numpy as np
import pandas as pd
import xarray as xr
import logging
from pathlib import Path
from datetime import datetime

def setup_logger(verbose=False):
    """Set up logging configuration"""
    log_level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[logging.StreamHandler()]
    )
    return logging.getLogger('merge_wct')

def load_completed_merges(log_file):
    """Load the set of completed merges from log file"""
    completed = set()
    if os.path.exists(log_file):
        try:
            with open(log_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line and not line.startswith('#'):
                        completed.add(line)
        except Exception as e:
            print(f"Warning: Could not read log file {log_file}: {e}")
    return completed

def save_completed_merge(log_file, merge_key):
    """Save a completed merge to the log file"""
    try:
        with open(log_file, 'a') as f:
            f.write(f"{merge_key}\n")
    except Exception as e:
        print(f"Warning: Could not write to log file {log_file}: {e}")

def create_merge_key(filter_dir, movstack_dir, component_dir, pair_dir):
    """Create a unique key for a merge operation"""
    return f"{filter_dir}|{movstack_dir}|{component_dir}|{pair_dir}"

def initialize_log_file(log_file):
    """Initialize the log file with header if it doesn't exist"""
    if not os.path.exists(log_file):
        try:
            os.makedirs(os.path.dirname(log_file), exist_ok=True)
            with open(log_file, 'w') as f:
                f.write(f"# WCT Merge Progress Log - Created {datetime.now()}\n")
                f.write("# Format: filter|movstack|component|pair\n")
        except Exception as e:
            print(f"Warning: Could not create log file {log_file}: {e}")

def get_filter_directories(wct_dir):
    """Get list of filter directories in the WCT directory"""
    if not os.path.exists(wct_dir):
        raise FileNotFoundError(f"WCT directory not found: {wct_dir}")
    
    # Look for filter directories (numbered folders like "01", "02", etc.)
    filters = []
    for item in os.listdir(wct_dir):
        if os.path.isdir(os.path.join(wct_dir, item)) and item.isdigit():
            filters.append(item)
    
    return sorted(filters)

def get_movstack_directories(wct_dir, filter_dir):
    """Get list of moving stack directories for a filter"""
    filter_path = os.path.join(wct_dir, filter_dir)
    movstacks = []
    for item in os.listdir(filter_path):
        if os.path.isdir(os.path.join(filter_path, item)):
            print(item)
            movstacks.append(item)
    
    return sorted(movstacks)

def get_component_directories(wct_dir, filter_dir, movstack_dir):
    """Get list of component directories for a filter and moving stack"""
    path = os.path.join(wct_dir, filter_dir, movstack_dir)
    components = []
    for item in os.listdir(path):
        if os.path.isdir(os.path.join(path, item)):
            components.append(item)
    
    return sorted(components)

def get_station_pair_directories(wct_dir, filter_dir, movstack_dir, component_dir):
    """Get list of station pair directories"""
    path = os.path.join(wct_dir, filter_dir, movstack_dir, component_dir)
    pairs = []
    for item in os.listdir(path):
        if os.path.isdir(os.path.join(path, item)) and '_' in item:
            pairs.append(item)
    
    return sorted(pairs)

def parse_pair_name(pair_folder):
    """Parse station pair folder name into station1 and station2"""
    return pair_folder.split('_', 1)

def merge_daily_files(wct_dir, output_dir, filter_dir, movstack_dir, component_dir, pair_dir, logger):
    """Merge daily files for a specific station pair and parameter set"""
    full_path = os.path.join(wct_dir, filter_dir, movstack_dir, component_dir, pair_dir)
    station1, station2 = parse_pair_name(pair_dir)
    
    # Find all date files (*.npz)
    files = glob.glob(os.path.join(full_path, "*.npz"))
    
    if not files:
        logger.warning(f"No files found in {full_path}")
        return None
    
    logger.info(f"Found {len(files)} files for filter {filter_dir}, movstack {movstack_dir}, component {component_dir}, pair {pair_dir}")
    
    # Initialize lists to store data
    dates = []
    dvv_values = []
    err_values = []
    coh_values = []
    
    # Pre-define variables that will be used later
    filterid = int(filter_dir)
    component = component_dir

    curr_window_days_meta = None
    ref_window_days_meta  = None
    step_days_meta        = None
    ref_lag_days_meta     = None
            
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
                    logger.warning(f"Could not determine date for {file}, skipping")
                    continue
            
            # Extract data values, ensuring they're in the correct format
            if 'dvv' in data:
                dvv = data['dvv']
                if isinstance(dvv, np.ndarray) and dvv.ndim > 0:
                    dvv_values.append(dvv)
                else:
                    logger.warning(f"Invalid dvv data in {file}, skipping")
                    continue
            else:
                logger.warning(f"No dvv data in {file}, skipping")
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
                    logger.warning(f"No taxis found in {file}, will use default")
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
                        
                if 'curr_window_days' in data:
                    try:
                        v = data['curr_window_days']
                        curr_window_days_meta = int(v.item() if isinstance(v, np.ndarray) else v)
                    except:
                        curr_window_days_meta = None
                
                if 'ref_window_days' in data:
                    try:
                        v = data['ref_window_days']
                        ref_window_days_meta = int(v.item() if isinstance(v, np.ndarray) else v)
                    except:
                        ref_window_days_meta = None
                
                if 'step_days' in data:
                    try:
                        v = data['step_days']
                        step_days_meta = int(v.item() if isinstance(v, np.ndarray) else v)
                    except:
                        step_days_meta = None
                        
                if 'ref_lag_days' in data:
                    try:
                        v = data['ref_lag_days']
                        ref_lag_days_meta = int(v.item() if isinstance(v, np.ndarray) else v)
                    except:
                        ref_lag_days_meta = None
        
        except Exception as e:
            logger.error(f"Error loading file {file}: {str(e)}")
            continue
    
    if not dates or len(dates) != len(dvv_values):
        logger.warning(f"No valid data found or data mismatch for {pair_dir} - {component_dir} f{filter_dir} m{movstack_dir}")
        return None
    
    # Get frequencies from the first file
    try:
        first_file = np.load(files[0], allow_pickle=True)
        if 'freqs' in first_file:
            freqs_data = first_file['freqs']
            if isinstance(freqs_data, np.ndarray) and freqs_data.size > 0:
                freqs = freqs_data
            else:
                # Generate default frequencies if invalid
                logger.warning(f"Invalid freqs data in {files[0]}, using defaults")
                freqs = np.linspace(0.1, 2.0, dvv_values[0].shape[0])
        else:
            # Generate default frequencies if missing
            logger.warning(f"No freqs data in {files[0]}, using defaults")
            freqs = np.linspace(0.1, 2.0, dvv_values[0].shape[0])
    except Exception as e:
        logger.error(f"Error extracting frequencies: {str(e)}")
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
            else:
                logger.warning(f"Skipping data point {date} due to shape mismatch")
        
        if not valid_indices:
            logger.error(f"No valid data points with consistent shapes")
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
    except Exception as e:
        logger.error(f"Error creating DataFrames: {str(e)}")
        return None
    
    # Create xarray dataset
    coords = {
        'times': dvv_df.index,
        'freq': freqs
    }
    
    attrs = {
        'station1': station1,
        'station2': station2,
        'component': component,
        'filterid': filterid,
        'mov_stack': movstack_dir,
        'created': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
        'merged_from': f"{len(filtered_dates)} daily files",
        'curr_window_days': str(curr_window_days_meta) if curr_window_days_meta is not None else movstack_dir.split('_')[0].replace('d', ''),
        'ref_window_days':  str(ref_window_days_meta)  if ref_window_days_meta  is not None else movstack_dir.split('_')[0].replace('d', ''),
        'step_days':        str(step_days_meta)         if step_days_meta         is not None else movstack_dir.split('_')[1].replace('d', ''),
        'ref_lag_days': str(ref_lag_days_meta) if ref_lag_days_meta is not None else movstack_dir.split('_')[0].replace('d', ''),
    }
    
    dvv = xr.DataArray(
        dvv_df.values,
        dims=['times', 'freq'],
        coords=coords,
        name='dvv'
    )
    
    err = xr.DataArray(
        err_df.values,
        dims=['times', 'freq'],
        coords=coords,
        name='err'
    )
    
    coh = xr.DataArray(
        coh_df.values,
        dims=['times', 'freq'],
        coords=coords,
        name='coh'
    )
    
    ds = xr.Dataset(
        {'dvv': dvv, 'err': err, 'coh': coh},
        attrs=attrs
    )
    
    # Create output directory matching the input hierarchy
    output_subdir = os.path.join(output_dir, filter_dir, movstack_dir, component_dir)
    os.makedirs(output_subdir, exist_ok=True)
    
    # Save to file - use the same station pair name for the output file
    filename = f"{pair_dir}.nc"
    output_path = os.path.join(output_subdir, filename)
    
    try:
        ds.to_netcdf(output_path)
        logger.info(f"Saved merged file: {output_path}")
        
        # Also save as CSV for easy inspection (optional)
        csv_dir = os.path.join(output_dir, 'csv', filter_dir, movstack_dir, component_dir)
        os.makedirs(csv_dir, exist_ok=True)
        
        csv_base = pair_dir
        dvv_df.to_csv(os.path.join(csv_dir, f"{csv_base}_dvv.csv"))
        err_df.to_csv(os.path.join(csv_dir, f"{csv_base}_err.csv"))
        coh_df.to_csv(os.path.join(csv_dir, f"{csv_base}_coh.csv"))
        
        return output_path
    except Exception as e:
        logger.error(f"Error saving dataset: {str(e)}")
        return None

def main():
    """Main function to merge all daily WCT files"""
    parser = argparse.ArgumentParser(description='Merge daily WCT files into consolidated files')
    parser.add_argument('--wct_dir', default='WCT', help='Directory containing WCT results')
    parser.add_argument('--output_dir', default='WCT_MERGED_32', help='Directory to save merged results')
    parser.add_argument('--verbose', action='store_true', help='Enable verbose output')
    args = parser.parse_args()
    
    logger = setup_logger(args.verbose)
    logger.info(f"Starting merge process, reading from {args.wct_dir}, writing to {args.output_dir}")
    
    # Set up progress tracking
    log_file = os.path.join(args.output_dir, 'merge_progress.log')
    initialize_log_file(log_file)
    completed_merges = load_completed_merges(log_file)
    
    logger.info(f"Found {len(completed_merges)} previously completed merges")
    
    try:
        # Get list of filter directories
        filters = get_filter_directories(args.wct_dir)
        if not filters:
            logger.error(f"No filter directories found in {args.wct_dir}")
            return 1
        
        logger.info(f"Found {len(filters)} filter directories to process")
        
        # Process each parameter combination
        merged_count = 0
        skipped_count = 0
        
        for filter_dir in filters:
            logger.info(f"Processing filter: {filter_dir}")
            
            # Get moving stack directories
            movstacks = get_movstack_directories(args.wct_dir, filter_dir)
            for movstack_dir in movstacks:
                logger.info(f"Processing movstack: {movstack_dir}")
                
                # Get component directories
                components = get_component_directories(args.wct_dir, filter_dir, movstack_dir)
                for component_dir in components:
                    logger.info(f"Processing component: {component_dir}")
                    
                    # Get station pair directories
                    pairs = get_station_pair_directories(args.wct_dir, filter_dir, movstack_dir, component_dir)
                    for pair_dir in pairs:
                        # Check if this merge was already completed
                        merge_key = create_merge_key(filter_dir, movstack_dir, component_dir, pair_dir)
                        
                        if merge_key in completed_merges:
                            logger.info(f"Skipping already completed merge: {merge_key}")
                            skipped_count += 1
                            continue
                        
                        logger.info(f"Merging files for filter {filter_dir}, movstack {movstack_dir}, component {component_dir}, pair {pair_dir}")
                        
                        result = merge_daily_files(
                            args.wct_dir, args.output_dir, filter_dir, 
                            movstack_dir, component_dir, pair_dir, logger
                        )
                        
                        if result:
                            merged_count += 1
                            # Log the completed merge
                            save_completed_merge(log_file, merge_key)
                            logger.info(f"Completed merge {merged_count}: {merge_key}")
        
        logger.info(f"Merge process completed. {merged_count} new merged files created, {skipped_count} skipped (already completed)")
        logger.info(f"Results saved in {args.output_dir}")
        logger.info(f"Progress log: {log_file}")
        return 0
    except Exception as e:
        logger.error(f"Error in merge process: {str(e)}")
        return 1

if __name__ == "__main__":
    sys.exit(main())
