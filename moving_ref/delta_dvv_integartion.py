"""
DvvIntegrator for 10d_10d stack with ref(-19:-10)
--------------------------------------------------
Each dvv(t) = V(t) - V(t-1), where t is in 10-day units.
Reconstruction: v(t) = v(t-1) + dvv(t), anchored at v(0) = 0.
lag_days = 1 (one 10-day step between consecutive blocks).

NaN handling: if dvv(t) is NaN, carry v(t-1) forward silently.
The reconstructed v(t) therefore never contains NaN.
"""

import xarray as xr
import numpy as np
from pathlib import Path
from tqdm import tqdm


class DvvIntegrator:
    def __init__(self, input_base, output_base, start_date=None):
        """
        Parameters
        ----------
        curr_window_days : int
            Width of current stack in days (e.g. 10)
        ref_window_days : int
            Width of reference stack in days (e.g. 30)
        step_days : int
            Step between consecutive stacks in days (e.g. 10)
            Derived from mov_stack[1]
        """
        self.input_base  = Path(input_base)
        self.output_base = Path(output_base)
        self.start_date  = start_date
        
    def _reconstruct_dvv(self, dvv_array, times_array):
        n   = dvv_array.shape[0]
        v   = np.zeros_like(dvv_array)
    
        t_days = times_array.astype('datetime64[D]').astype(float)
        dt_days       = np.diff(t_days)
        gap_threshold = self.step_days * 3
        is_gap        = np.concatenate([[False], dt_days > gap_threshold])
    
        for t in range(1, n):
            delta = dvv_array[t]
    
            if is_gap[t]:
                v[t] = np.where(np.isnan(delta), v[t-1], v[t-1] + delta)
                continue
    
            # ref window in days: [t_days[t] - ref_lag_days - ref_window_days,
            #                       t_days[t] - ref_lag_days]
            ref_end_day   = t_days[t] - self.ref_lag_days
            ref_start_day = ref_end_day - self.ref_window_days
    
            # find indices of past reconstructed values whose timestamp
            # falls within the ref window
            in_ref_window = np.where(
                (t_days[:t] >= ref_start_day) &
                (t_days[:t] <= ref_end_day)
            )[0]
    
            if len(in_ref_window) == 0:
                # not enough history yet, carry forward
                v[t] = np.where(np.isnan(delta), v[t-1], v[t-1] + delta)
                continue
    
            v_ref = np.nanmean(v[in_ref_window], axis=0)
            v[t]  = np.where(np.isnan(delta), v_ref, v_ref + delta)
    
        return v    
        
    def _reconstruct_err(self, err_array, times_array):
        n   = err_array.shape[0]
        e   = np.zeros_like(err_array)
    
        t_days        = times_array.astype('datetime64[D]').astype(float)
        dt_days       = np.diff(t_days)
        gap_threshold = self.step_days * 3
        is_gap        = np.concatenate([[False], dt_days > gap_threshold])
    
        for t in range(1, n):
            err_t = err_array[t]
    
            if is_gap[t]:
                e[t] = np.where(
                    np.isnan(err_t),
                    e[t-1],
                    np.sqrt(e[t-1]**2 + err_t**2)
                )
                continue
    
            ref_end_day   = t_days[t] - self.ref_lag_days
            ref_start_day = ref_end_day - self.ref_window_days
    
            in_ref_window = np.where(
                (t_days[:t] >= ref_start_day) &
                (t_days[:t] <= ref_end_day)
            )[0]
    
            if len(in_ref_window) == 0:
                e[t] = np.where(np.isnan(err_t), e[t-1], np.sqrt(e[t-1]**2 + err_t**2))
                continue
    
            e_ref = np.sqrt(np.nanmean(e[in_ref_window]**2, axis=0))
            e[t]  = np.where(
                np.isnan(err_t),
                e_ref,
                np.sqrt(e_ref**2 + err_t**2)
            )
    
        return e
        
#    def _reconstruct_dvv(self, dvv_array, times_array):
#        n   = dvv_array.shape[0]
#        v   = np.zeros_like(dvv_array)
#    
#        dt_days = np.diff(times_array.astype('datetime64[D]').astype(float))
#        gap_threshold = self.step_days * 3
#        is_gap = np.concatenate([[False], dt_days > gap_threshold])
#    
#        for t in range(1, n):
#            delta = dvv_array[t]
#            lag   = self.ref_lag
#    
#            if is_gap[t]:
#                # find last valid v before the gap
#                v_before_gap = v[t - 1]          # always valid (zeros or reconstructed)
#                v[t] = np.where(np.isnan(delta), v_before_gap, v_before_gap + delta)
#                continue
#    
#            if t < lag:
#                continue
#    
#            v_ref = v[t - lag]
#            v[t]  = np.where(np.isnan(delta), v_ref, v_ref + delta)
#    
#        return v

#    def _reconstruct_err(self, err_array, times_array):
#        n   = err_array.shape[0]
#        e   = np.zeros_like(err_array)
#    
#        dt_days = np.diff(times_array.astype('datetime64[D]').astype(float))
#        gap_threshold = self.step_days * 3
#        is_gap = np.concatenate([[False], dt_days > gap_threshold])
#    
#        for t in range(1, n):
#            err_t = err_array[t]
#            lag   = self.ref_lag
#    
#            if is_gap[t]:
#                e_before_gap = e[t - 1]
#                e[t] = np.where(
#                    np.isnan(err_t),
#                    e_before_gap,
#                    np.sqrt(e_before_gap**2 + err_t**2)
#                )
#                continue
#    
#            if t < lag:
#                continue
#    
#            e_ref = e[t - lag]
#            e[t]  = np.where(
#                np.isnan(err_t),
#                e_ref,
#                np.sqrt(e_ref**2 + err_t**2)
#            )
#    
#        return e
    
    # ------------------------------------------------------------------
    # Dataset-level integration
    # ------------------------------------------------------------------

    def integrate_dataset(self, ds):
        """
        Integrate a single xarray Dataset.

        Parameters
        ----------
        ds : xarray.Dataset
            Must contain 'dvv', 'err', 'coh' with dimension 'times'.

        Returns
        -------
        xarray.Dataset
            Reconstructed dataset, same structure, no NaNs in dvv/err.
        """
        # Sort and slice
        ds_out = ds.sortby('times')

        if self.start_date is not None:
            ds_out = ds_out.sel(times=slice(self.start_date, None))

        if ds_out.sizes['times'] == 0:
            print(f"  Warning: no data after start_date={self.start_date}")
            return ds_out

        # Extract raw numpy arrays
        times_array = ds_out['times'].values 
        dvv_vals = ds_out['dvv'].values  # shape: (n_times, n_freqs)
        err_vals = ds_out['err'].values

        # Handle 1D edge case (single frequency)
        squeeze = dvv_vals.ndim == 1
        if squeeze:
            dvv_vals = dvv_vals[:, np.newaxis]
            err_vals = err_vals[:, np.newaxis]

        # Reconstruct — NaN-safe, output has no NaNs
        dvv_reconstructed = self._reconstruct_dvv(dvv_vals, times_array)   # <-- pass it
        err_reconstructed = self._reconstruct_err(err_vals, times_array)   # <-- pass it


        # Verify no NaNs leaked through
        assert not np.any(np.isnan(dvv_reconstructed)), \
            "BUG: NaN found in reconstructed dvv — should not happen"
        assert not np.any(np.isnan(err_reconstructed)), \
            "BUG: NaN found in reconstructed err — should not happen"

        if squeeze:
            dvv_reconstructed = dvv_reconstructed[:, 0]
            err_reconstructed = err_reconstructed[:, 0]

        # Build output dataset using assign to avoid in-place xarray issues
        ds_out = ds_out.assign(
            dvv=xr.DataArray(
                dvv_reconstructed,
                coords=ds_out['dvv'].coords,
                dims=ds_out['dvv'].dims
            ),
            err=xr.DataArray(
                err_reconstructed,
                coords=ds_out['err'].coords,
                dims=ds_out['err'].dims
            )
        )
        # coh is unchanged
        ds_out.attrs.update({
            'integrated':        'true',
            'integration_method': (
                'lag_inversion: v(t)=mean(v[t-ref_lag-ref_size:t-ref_lag])+dvv(t), '
                'NaN=carry-forward'
            ),
            'curr_size':         str(self.curr_size),
            'ref_size':          str(self.ref_size),
            'ref_lag':           str(self.ref_lag),
            'integration_start': str(self.start_date) if self.start_date else 'full',
        })

        return ds_out

    # ------------------------------------------------------------------
    # File-level processing
    # ------------------------------------------------------------------

    def get_output_path(self, input_path):
        relative_path = input_path.relative_to(self.input_base)
        return self.output_base / relative_path

    def process_file(self, input_path, overwrite=False):
        output_path = self.get_output_path(input_path)
    
        if output_path.exists() and not overwrite:
            print(f"  Skipping {input_path.name} (output exists)")
            return True
    
        try:
            ds = xr.open_dataset(input_path)
            
            curr_window_days = int(ds.attrs['curr_window_days'])
            ref_window_days  = int(ds.attrs['ref_window_days'])
            step_days        = int(ds.attrs['step_days'])
            ref_lag_days     = int(ds.attrs['ref_lag_days'])
            
            self.curr_window_days = curr_window_days
            self.ref_window_days  = ref_window_days
            self.ref_lag_days     = ref_lag_days
            self.step_days = step_days                        
            self.curr_size = curr_window_days / step_days
            self.ref_size  = ref_window_days  / step_days
            self.ref_lag   = ref_lag_days / step_days
            
            ds_integrated = self.integrate_dataset(ds)
            output_path.parent.mkdir(parents=True, exist_ok=True)
            ds_integrated.to_netcdf(output_path)
            ds.close()
            ds_integrated.close()
            return True
    
        except KeyError as e:
            print(f"  Missing stack config attribute in {input_path}: {e}")
            return False
        except Exception as e:
            print(f"  Error processing {input_path}: {e}")
            return False
            
    def process_all(self, overwrite=False, verbose=True):
        nc_files = list(self.input_base.rglob('*.nc'))

        if not nc_files:
            print(f"No NetCDF files found in {self.input_base}")
            return {'total': 0, 'success': 0, 'failed': 0}

        print(f"Found    : {len(nc_files)} NetCDF files")
        print(f"Input    : {self.input_base}")
        print(f"Output   : {self.output_base}")
        print(f"start    : {self.start_date}")
        print(f"Note     : curr_size/ref_size/ref_lag derived per file from ds.attrs")
        #print(f"start    : {self.start_date}")

        success_count = 0
        failed_count  = 0
        iterator = tqdm(nc_files) if verbose else nc_files

        for nc_file in iterator:
            if verbose:
                iterator.set_description(nc_file.name)
            if self.process_file(nc_file, overwrite=overwrite):
                success_count += 1
            else:
                failed_count += 1

        summary = {
            'total':   len(nc_files),
            'success': success_count,
            'failed':  failed_count,
        }
        print("\n" + "=" * 50)
        print(f"  Total:   {summary['total']}") 
        print(f"  Success: {summary['success']}")
        print(f"  Failed:  {summary['failed']}")
        print("=" * 50)
        return summary


# ----------------------------------------------------------------------
# Entry point
# ----------------------------------------------------------------------

def main():
    integrator = DvvIntegrator(
        input_base  = 'WCT_MERGED_32',
        output_base = 'WCT_INTEGRATED_32',
        start_date  = '2018-01-01'
    )
    summary = integrator.process_all(overwrite=True)
    return 0 if summary['failed'] == 0 else 1

if __name__ == "__main__":
    import sys
    sys.exit(main())
