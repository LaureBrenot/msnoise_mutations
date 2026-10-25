import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.dates import DateFormatter, YearLocator, MonthLocator
import numpy as np
import seaborn as sns
import os
from pathlib import Path
from itertools import combinations
import geopandas as gpd

class VolcanoConfig:
    def __init__(self):
        self.volcano = 'Askja'
        self.component = 'ZZ'
        self.window = '10d_1d'
        self.filter_id = '02'
        self.network = '8K'
        
        # Station pairs - we'll focus on these four
        self.stations = ["ASK", "JONS", "KATT", "STAM"]
        self.plot_all_pairs = True
        
        # Frequency ranges
        self.freqranges = [[0.5,1.0], [1.0,2.0], [2.0,4.0], [0.5,2.0]]
        
        # Time period
        self.startdate = '2018-01-01'
        self.enddate = '2024-09-01'
        
        # Ref to print
        self.titleref = "Dynamic ref -30-1"
        
        # Plot configuration
        self.coherence_threshold = 0
        self.use_fixed_dvv_scale = False
        self.dvv_scale = {'min': -0.2, 'max': 0.2}
        self.plot_events = True
        self.event_list = [pd.to_datetime("2021-08-08 00:00:00")]
        self.event_names = [""] #01/08 inflation
        self.save_fig = True
        self.seasonal_rm = False
        self.show_plots = False
        
        # Map configuration
        self.shapefile_paths = {
            'lake': "/globalscratch/ulb/gtime/lbrenot/copy_last/msnoise__copy_static/velocities/oskjuvatn.shp",
            'caldera': "/globalscratch/ulb/gtime/lbrenot/copy_last/msnoise__copy_static/velocities/caldera.shp",
            'fissure': "/globalscratch/ulb/gtime/lbrenot/copy_last/msnoise__copy_static/velocities/fissure_swarm.shp"
        }

class DataLoader:
    def __init__(self, config):
        self.config = config
        self.base_path = Path("WCT_MERGEDc")
    
    def get_file_path(self, stat1, stat2):
        return (self.base_path / self.config.filter_id / self.config.window / 
                self.config.component / 
                f"{self.config.network}.{stat1}.--_{self.config.network}.{stat2}.--.nc")
    
    def read_data(self, stat1, stat2):
        file_path = self.get_file_path(stat1, stat2)
        if not file_path.exists():
            # Try the reverse order
            file_path = self.get_file_path(stat2, stat1)
            if not file_path.exists():
                return None, None, None
        
        ds = xr.open_dataset(file_path)
        return ds['dvv'].to_pandas(), ds['err'].to_pandas(), ds['coh'].to_pandas()

class Plotter:
    def __init__(self, config):
        self.config = config
        self.colorsp = ['Blues', 'Reds','Greens', 'Purples','Greys']
        self.colorsl = ['blue', 'red', 'green', 'purple']
        sns.set(style='ticks')
    
    def get_coordinates(self, station):
        """Get coordinates for each station"""
        coordinates = {
            'ASK': (65.051941,-16.64806),
            'JONS': (65.077469,-16.8057),
            'KATT': (64.999008,-16.965389),
            'STAM': (64.99691,-16.809589),
            'SVAD': (65.117462,-16.57498),
            'BRUN': (65.204613,-16.865971),
            'FLUR':(64.843536,-17.02693),
            'KOLL':(65.290237,-16.567261),
            'LIND':(64.852783,-16.452299),
            'MOFO':(64.984398,-16.65119),
            'RIFR':(64.915329,-16.371269),
            'TOHR':(64.91658,-16.784731),
            'VEGG':(65.38205,-16.37467),
            'DDAL':(65.077393,-16.933411),
            'DALR':(65.077332,-16.936701)
        }
        return coordinates.get(station, (0, 0))
    
    def create_station_map_inset(self, fig, position='upper right'):
        """Create a station map inset on the main figure"""
        # Create inset axes in the top right corner
        # Position: [x0, y0, width, height] in figure coordinates
        if position == 'upper right':
            inset_ax = fig.add_axes([0.80, 0.85, 0.2, 0.15])  # Adjust these values as needed
        else:
            inset_ax = fig.add_axes([0.80, 0.85, 0.2, 0.15])  # Default to upper right
        
        try:
            # Load shapefiles
            gdf_lake = gpd.read_file(self.config.shapefile_paths['lake'])
            gdf_caldera = gpd.read_file(self.config.shapefile_paths['caldera'])
            # Uncomment if you want the fissure swarm
            # gdf_fissure = gpd.read_file(self.config.shapefile_paths['fissure'])
            
            # Get station coordinates and calculate bounds
            station_coords = [self.get_coordinates(station) for station in self.config.stations]
            lats = [coord[0] for coord in station_coords]
            lons = [coord[1] for coord in station_coords]
            
            # Set plot limits based on station extent with padding
            padding = 0.05
            min_lat, max_lat = min(lats) - padding, max(lats) + padding
            min_lon, max_lon = min(lons) - padding, max(lons) + padding
            
            # Plot the shapefiles
            gdf_lake.plot(ax=inset_ax, color='lightblue')
            gdf_caldera.plot(ax=inset_ax, linewidth=2, color='k')
            # gdf_fissure.plot(ax=inset_ax, linewidth=2, color='lightgrey')
            
            # Plot stations
            latitude_threshold = 65.03
            for station in self.config.stations:
                lat, lon = self.get_coordinates(station)
                inset_ax.scatter(lon, lat, marker='^', c='r', s=100, zorder=10)
                
                # Adjust text position based on latitude threshold
                if lat > latitude_threshold:
                    text_x_offset = -0.08
                    text_y_offset = +0.015
                    va = 'top'
                else:
                    text_x_offset = -0.01
                    text_y_offset = -0.02
                    va = 'bottom'
                
                inset_ax.text(lon + text_x_offset, lat + text_y_offset, station, 
                            fontsize=18, ha='left', va=va)
            
            # Set the plot limits
            inset_ax.set_xlim(min_lon, max_lon)
            inset_ax.set_ylim(min_lat+0.02, max_lat-0.02)
            
#            # Remove axes
            inset_ax.get_xaxis().set_visible(False)
            inset_ax.get_yaxis().set_visible(False)
#            inset_ax.axis('off')

            # Style the inset
#            inset_ax.tick_params(axis='both', which='major', labelsize=16)
#            inset_ax.set_xlabel('Longitude', fontsize=18)
#            inset_ax.set_ylabel('Latitude', fontsize=18)
            
        except Exception as e:
            print(f"Warning: Could not create station map inset: {e}")
            # If shapefiles fail to load, create a simple station plot
            station_coords = [self.get_coordinates(station) for station in self.config.stations]
            lats = [coord[0] for coord in station_coords]
            lons = [coord[1] for coord in station_coords]
            
            for i, station in enumerate(self.config.stations):
                lat, lon = self.get_coordinates(station)
                inset_ax.scatter(lon, lat, marker='^', c='r', s=50, zorder=10)
                inset_ax.text(lon + 0.01, lat + 0.01, station, fontsize=8)
            
            padding = 0.05
            min_lat, max_lat = min(lats) - padding, max(lats) + padding
            min_lon, max_lon = min(lons) - padding, max(lons) + padding
            inset_ax.set_xlim(min_lon, max_lon)
            inset_ax.set_ylim(min_lat, max_lat)
            inset_ax.tick_params(axis='both', which='major', labelsize=6)
            inset_ax.set_xlabel('Longitude', fontsize=8)
            inset_ax.set_ylabel('Latitude', fontsize=8)
        
        return inset_ax

    def plot_multiple_pairs(self, pairs_data, station_pairs):
        """Plot six subplots with dv/V visualizations for the specified station pairs"""
        fig, axs = plt.subplots(3, 2, figsize=(20, 24), sharex=True)
        axs = axs.flatten()
        
        # Convert config dates to datetime objects
        start_date = pd.to_datetime(self.config.startdate)
        end_date = pd.to_datetime(self.config.enddate)
        
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
            
            # Filter data to the specified date range
            dvv_filtered = dvv_df.loc[start_date:end_date]
            coh_filtered = coh_df.loc[start_date:end_date]
            
            # Create the freq-time plot
            coh_mask = coh_filtered.values > self.config.coherence_threshold
            dvv_masked = dvv_filtered.copy()
            #dvv_masked.values[~coh_mask] = np.nan
            
            # Using a diverging colormap centered at zero
            cmap = mpl.cm.RdBu_r
            
            full_index = pd.date_range(start=start_date, end=end_date, freq='D')
            data_reindexed = dvv_masked.reindex(full_index)
            time = pd.to_datetime(full_index)
            data_array = np.ma.masked_invalid(data_reindexed.values)
            
            freq = dvv_filtered.columns.astype(float)
            mesh = ax.pcolormesh(time, freq, data_array.T, cmap=cmap, vmin=vmin, vmax=vmax)
            
            ax.set_title(f"{stat1}-{stat2}", fontsize=28)
            
            # Set x-axis limits to the specified date range
            ax.set_xlim(start_date, end_date)
            
            # Add events if configured
            if self.config.plot_events:
                for event, name in zip(self.config.event_list, self.config.event_names):
                    if start_date <= event <= end_date:  # Only plot events within the date range
                        ax.axvline(event, color='k', ls='--', lw=3)
                        #ax.text(event, ax.get_ylim()[1], name, rotation=90, va='bottom', fontsize=12)
            
            # Format date axis
            ax.xaxis.set_major_locator(YearLocator())
            ax.xaxis.set_minor_locator(MonthLocator())
            ax.xaxis.set_major_formatter(DateFormatter('%Y'))
            
            ax.tick_params(axis='x', which='major', length=10, width=2 , labelsize=25, labelbottom=True)
            ax.tick_params(axis='x', which='minor', length=6, width=1, )
            ax.tick_params(axis='y', which='major', length=8, width=1.5 , labelsize=25)
    
            # Only add x-label to bottom plots
            if i >= 4:
                ax.set_xlabel('Date', fontsize=30)
            if i == 0 or i == 2 or i == 4:
                ax.set_ylabel('Frequency (Hz)', fontsize=30)
        
        # Add colorbar to the right of the subplots
        cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
        cbar = fig.colorbar(mesh, cax=cbar_ax)
        cbar.set_label('dv/v (%)', rotation=90, labelpad=20, fontsize=30)
        cbar.ax.tick_params(labelsize=25)
        
        # Add the station map inset
        self.create_station_map_inset(fig, position='upper right')
        
        plt.tight_layout(rect=[0, 0, 0.9, 1])
        
        if self.config.save_fig:
            if not os.path.exists('Figures'):
                os.makedirs('Figures')
            fig.savefig(f'Figures/multi_station_dvv_{self.config.component}_{self.config.window}_{max(freq)}.png',  
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
    
    # Generate pairs between ASK, JONS, STAM, and KATT
    pairs = list(combinations(config.stations, 2))
    
    # Store data for each pair
    pairs_data = []
    
    for stat1, stat2 in pairs:
        print(f"Processing {stat1}-{stat2}...")
        dvv_df, err_df, coh_df = loader.read_data(stat1, stat2)
        
        if dvv_df is None:
            print(f"No data found for {stat1}-{stat2}")
            pairs_data.append(None)
        else:
            # Only collect required data for the multi-plot
            pairs_data.append((dvv_df, err_df, coh_df))
    
    # Plot all pairs in one figure
    plotter.plot_multiple_pairs(pairs_data, pairs)
    print(f"Output path: {os.path.abspath('Figures')}")

if __name__ == "__main__":
    main()
