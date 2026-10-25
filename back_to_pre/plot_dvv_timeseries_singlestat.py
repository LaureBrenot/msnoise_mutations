import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.dates import DateFormatter, YearLocator, MonthLocator
import numpy as np
from matplotlib.lines import Line2D
import seaborn as sns
import os
from pathlib import Path
from msnoise.api import connect, get_t_axis, xr_get_ccf, xr_get_ref
import datetime

class VolcanoConfig:
    def __init__(self):
        self.volcano = 'Askja'
        self.component = 'EN'
        self.window = '30d_1d'
        self.filter_id = '02'
        self.network = ["8K"]  # Changed to list like in second script
        
        # Station pairs - we'll focus on these four for autocorrelation
        self.stations = ["ASK", "JONS", "STAM", "KATT"]
        self.plot_all_pairs = True
        
        # Frequency ranges
        self.freqranges = [[0.5,1.0], [1.0,2.0], [2.0,4.0], [0.5,2.0]]
        
        # Time period
        self.startdate = '2018-01-01'
        self.enddate = '2024-09-01'
        
        # Ref to print
        self.titleref = "Dynamic ref -30-1"
        
        # Plot configuration
        self.plot_time_series = False 
        self.coherence_threshold = 0
        self.use_fixed_dvv_scale = False
        self.dvv_scale = {'min': -0.2, 'max': 0.2}
        self.plot_events = True
        self.event_list = [pd.to_datetime("2021-08-01 00:00:00")]
        self.event_names = [""] #01/08 inflation
        self.save_fig = True
        self.seasonal_rm = False
        self.show_plots = False

class DataLoader:
    def __init__(self, config):
        self.config = config
        self.base_path = Path("WCT_MERGEDc")
    
    def get_file_path(self, stat1, stat2, network):
        return (self.base_path / self.config.filter_id / self.config.window / 
                self.config.component / 
                f"{network}.{stat1}.--_{network}.{stat2}.--.nc")
    
    def read_data(self, stat1, stat2):
        """Modified to handle multiple networks like in the second script"""
        dvv_list, err_list, coh_list = [], [], []
    
        for net in self.config.network:
            file_path = self.get_file_path(stat1, stat2, net)
            print(f"Checking: {file_path}")
            if file_path.exists():
                ds = xr.open_dataset(file_path)
                dvv_list.append(ds['dvv'].to_pandas())
                err_list.append(ds['err'].to_pandas())
                coh_list.append(ds['coh'].to_pandas())
        
        if not dvv_list:
            return None, None, None
    
        # Concatenate and average across networks
        dvv_all = pd.concat(dvv_list).sort_index().groupby(level=0).mean()
        err_all = pd.concat(err_list).sort_index().groupby(level=0).mean()
        coh_all = pd.concat(coh_list).sort_index().groupby(level=0).mean()
        
        return dvv_all, err_all, coh_all

class Plotter:
    def __init__(self, config):
        self.config = config
        self.colorsp = ['Blues', 'Reds','Greens', 'Purples','Greys']
        self.colorsl = ['blue', 'red', 'green', 'purple']
        sns.set(style='ticks')
    
    def plot_multiple_pairs(self, pairs_data, station_pairs):
        """Plot four subplots with dv/V visualizations for autocorrelations"""
        fig, axs = plt.subplots(2, 2, figsize=(18, 12), sharex=True)  # Changed to 2x2
        axs = axs.flatten()
        
        # Use the same color scale for all plots, centered at zero
        if self.config.use_fixed_dvv_scale:
            vmin, vmax = self.config.dvv_scale['min'], self.config.dvv_scale['max']
        else:
            # Calculate global min/max from all datasets
            all_data = []
            for pair_data in pairs_data:
                if pair_data is not None:
                    dvv_df, _, coh_df = pair_data
                    coh_mask = coh_df.values > self.config.coherence_threshold
                    valid_data = dvv_df.copy()
                    valid_data.values[~coh_mask] = np.nan
                    all_data.append(valid_data.values)
            
            if all_data:
                all_data_combined = np.concatenate([data.flatten() for data in all_data])
                # Calculate percentiles but ensure the colormap is centered at zero
                abs_max = max(abs(np.nanpercentile(all_data_combined, 2)), 
                              abs(np.nanpercentile(all_data_combined, 98)))
                vmin, vmax = -abs_max, abs_max
            else:
                vmin, vmax = -0.2, 0.2  # Default if no data

        for i, (pair_data, ax) in enumerate(zip(pairs_data, axs)):
            if pair_data is None:
                ax.text(0.5, 0.5, "No data available", 
                        ha='center', va='center', transform=ax.transAxes, fontsize=16)
                continue
                
            dvv_df, err_df, coh_df = pair_data
            stat1, stat2 = station_pairs[i]
            
            if self.config.plot_time_series:
                # Plot time series instead of heatmap
                legend_elements = []  
                
                for j, freqrange in enumerate(self.config.freqranges):
                    freq_name = f"{freqrange[0]}-{freqrange[1]} Hz"
                    dvv_slice = self.get_freq_slice(dvv_df, freqrange)
                    coh_slice = self.get_freq_slice(coh_df, freqrange)
                    
                    mask = coh_slice['avg'] >= self.config.coherence_threshold
                    ax.scatter(dvv_slice.index[mask], 
                              dvv_slice['avg'][mask],
                              c=coh_slice['avg'][mask], 
                              cmap=self.colorsp[j], 
                              label=freq_name,
                              s=20)
                    
                    legend_elements.append(Line2D([0], [0], marker='o',
                                                color=self.colorsl[j],
                                                label=freq_name,
                                                markerfacecolor=self.colorsl[j],
                                                markersize=8))
                
                ax.legend(handles=legend_elements, fontsize=12)
                ax.set_ylabel('dv/v (%)', fontsize=24)
                
                if self.config.use_fixed_dvv_scale:
                    ax.set_ylim(self.config.dvv_scale['min'], self.config.dvv_scale['max'])
            else:

                # Create the freq-time plot
                coh_mask = coh_df.values > self.config.coherence_threshold
                dvv_masked = dvv_df.copy()
                #dvv_masked.values[~coh_mask] = np.nan
                
                # Using a diverging colormap centered at zero
                cmap = mpl.cm.RdBu_r
                
                full_index = pd.date_range(start=dvv_df.index.min(), end=dvv_df.index.max(), freq='D')
                data_reindexed = dvv_masked.reindex(full_index)
                time = pd.to_datetime(full_index)
                data_array = np.ma.masked_invalid(data_reindexed.values)
                
                freq = dvv_df.columns.astype(float)
                mesh = ax.pcolormesh(time, freq, data_array.T, cmap=cmap, vmin=vmin, vmax=vmax)
            
            ax.set_title(f"{stat1}-{stat2}", fontsize=28)
            
            # Add events if configured
            if self.config.plot_events:
                for event, name in zip(self.config.event_list, self.config.event_names):
                    ax.axvline(event, color='k', ls='--', lw=3)
            
            # Format date axis
            ax.xaxis.set_major_locator(YearLocator())
            ax.xaxis.set_minor_locator(MonthLocator())
            ax.xaxis.set_major_formatter(DateFormatter('%Y'))
            
            ax.tick_params(axis='x', which='major', length=10, width=2, labelsize=20, rotation=45,  labelleft=True, labelbottom=True)
            ax.tick_params(axis='x', which='minor', length=6, width=1)
            ax.tick_params(axis='y', which='major', length=8, width=1.5, labelsize=20)

            # Add x-label to bottom plots
            if i >= 2:  # Changed from 4 to 2 for 2x2 grid
                ax.set_xlabel('Date', fontsize=24)
            if i == 0 or i == 2:  # Left column
                if self.config.plot_time_series:
                    ax.set_ylabel('dv/v (%)', fontsize=24)
                else:
                    ax.set_ylabel('Frequency (Hz)', fontsize=24)
                    
        if not self.config.plot_time_series:
            # Add colorbar to the right of the subplots
            cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
            cbar = fig.colorbar(mesh, cax=cbar_ax)
            cbar.set_label('dv/v (%)', rotation=90, labelpad=20, fontsize=24)
            cbar.ax.tick_params(labelsize=20)

        
        plt.tight_layout(rect=[0, 0, 0.9, 1])
        
        if self.config.save_fig:
            if not os.path.exists('Figures'):
                os.makedirs('Figures')
            if not self.config.plot_time_series:
                fig.savefig(f'Figures/autocorr_dvv_{self.config.component}_{self.config.window}_{max(freq)}.png',  
                            dpi=300, bbox_inches='tight')
            else:
                fig.savefig(f'Figures/autocorr_dvv_timeserie_{self.config.component}_{self.config.window}.png',  
                            dpi=300, bbox_inches='tight')                            
        
        if self.config.show_plots:
            plt.show()
            
        return fig

    def get_freq_slice(self, data, freqrange):
        freqs = np.asarray(data.columns)
        minidx = (np.abs(freqs - freqrange[1])).argmin()
        maxidx = (np.abs(freqs - freqrange[0])).argmin()
        
        df_slice = data.iloc[:, minidx:maxidx].copy()
        df_slice['avg'] = df_slice.mean(axis=1).rolling(30, min_periods=10).mean()
        
        return df_slice
    
def main():
    config = VolcanoConfig()
    loader = DataLoader(config)
    print(f"Input path: {loader.base_path}")
    plotter = Plotter(config)
    
    # Generate autocorrelation pairs (same station with itself)
    pairs = [(station, station) for station in config.stations]
    
    # Store data for each pair
    pairs_data = []
    
    for stat1, stat2 in pairs:
        print(f"Processing {stat1}-{stat2}...")
        dvv_df, err_df, coh_df = loader.read_data(stat1, stat2)
        
        if dvv_df is None:
            print(f"No data found for {stat1}-{stat2}")
            pairs_data.append(None)
        else:
            db = connect()
            start = datetime.datetime.strptime(config.startdate, '%Y-%m-%d').date()
            end = datetime.datetime.strptime(config.enddate, '%Y-%m-%d').date()
            taxis = get_t_axis(db)
            mov = config.window.split('_')
            
            ccf_data_list = []
            ref_ccf_list = []
            
            for net in config.network:
                station1 = f"{net}.{stat1}.--"
                station2 = f"{net}.{stat2}.--"
                station1_ref = f"8K.{stat1}.--"
                station2_ref = f"8K.{stat2}.--"
                try:
                    ccf = xr_get_ccf(station1, station2, config.component, 
                                     int(config.filter_id), (mov[0],mov[1]), taxis)
                    ref = xr_get_ref(station1_ref, station2_ref, config.component, 
                                     int(config.filter_id), taxis)
                    
                    if ccf is not None:
                        ccf_data_list.append(ccf)
                    if ref is not None:
                        ref_ccf_list.append(ref)
                        
                except Exception as e:
                    print(f"Error with {net} network: {e}")
                    continue
            
            if ccf_data_list:
                ccf_data = pd.concat(ccf_data_list).sort_index().groupby(level=0).mean()
                ref_ccf = sum(ref_ccf_list) / len(ref_ccf_list)
                
                # Calculate correlations
                corr_df = pd.DataFrame(index=dvv_df.index, columns=dvv_df.columns)
                for time in ccf_data.index:
                    if time in dvv_df.index:
                        current_ccf = ccf_data.loc[time].values
                        if len(current_ccf.shape) > 1:
                            current_ccf = current_ccf.flatten()  
                        ref_data = ref_ccf.values.flatten()  
                        
                        if len(current_ccf) == len(ref_data):
                            corr = np.corrcoef(current_ccf, ref_data)[0,1]
                            corr_df.loc[time].values[:] = corr
            else:
                corr_df = coh_df  # Fallback to coherence if correlation calculation fails
            
            pairs_data.append((dvv_df, err_df, corr_df))
    
    # Plot all pairs in one figure
    plotter.plot_multiple_pairs(pairs_data, pairs)
    print(f"Output path: {os.path.abspath('Figures')}")

if __name__ == "__main__":
    main()