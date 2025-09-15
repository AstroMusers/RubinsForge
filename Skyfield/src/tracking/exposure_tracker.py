import numpy as np
from astropy.table import Table
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
import time

class ExposureTracker:
    def __init__(self, exposure_id, start_time, end_time, observer_location):
        """
        Initialize an ExposureTracker object

        Parameters:
        -----------
        exposure_id : int
            Unique identifier for this exposure
        start_time : Skyfield time
            Start time of exposure
        end_time : Skyfield time
            End time of exposure
        observer_location : EarthLocation
            Observer location for coordinate transformations
        """
        self.exposure_id = exposure_id
        self.start_time = start_time
        self.end_time = end_time
        self.observer_location = observer_location
        
        # Exposure-level statistics
        self.total_events = 0
        self.trail_candidates = 0
        self.contaminations = 0
        self.no_topo_skip = 0
        self.no_brightness_skip = 0
        self.no_targets_skip = 0
        
        # Per-timestep data storage
        self.timestep_data = []  # List of timestep dictionaries
        
        # Exposure-level aggregated data
        self.all_satellite_positions = []
        self.all_contaminated_targets = []
        self.all_contaminated_names = []
        self.all_intensities = []
        self.all_ab_magnitudes = []
        self.all_separations = []
        self.all_satellite_names = []
        self.all_times = []
        self.all_angular_velocities = []
        
    def add_timestep_data(self, ti, satellite_data, target_data, contamination_data):
        """
        Add data for a specific timestep
        
        Parameters:
        -----------
        ti : Skyfield time
            Current time
        satellite_data : dict
            Dictionary containing satellite position and brightness data
        target_data : dict
            Dictionary containing target position data
        contamination_data : dict
            Dictionary containing contamination results
        """
        timestep_entry = {
            'time': ti,
            'satellites': satellite_data,
            'targets': target_data,
            'contamination': contamination_data
        }
        
        self.timestep_data.append(timestep_entry)
        
        # Update exposure-level aggregated data
        if satellite_data:
            for sat_name, sat_info in satellite_data.items():
                self.all_satellite_positions.append({
                    'name': sat_name,
                    'time': ti,
                    'alt': sat_info['alt'],
                    'az': sat_info['az'],
                    'height': sat_info['height'],
                    'intensity': sat_info.get('intensity', np.nan),
                    'ab_magnitude': sat_info.get('ab_magnitude', np.nan),
                    'angular_velocity': sat_info.get('angular_velocity', np.nan)
                })
                
                if not np.isnan(sat_info.get('intensity', np.nan)):
                    self.all_intensities.append(sat_info['intensity'])
                    self.all_ab_magnitudes.append(sat_info['ab_magnitude'])
                    self.all_satellite_names.append(sat_name)
                    self.all_times.append(ti)
                if not np.isnan(sat_info.get('angular_velocity', np.nan)):
                    self.all_angular_velocities.append(sat_info['angular_velocity'])
        
        # Add contaminated targets
        if contamination_data and contamination_data.get('contaminated_targets'):
            for target in contamination_data['contaminated_targets']:
                self.all_contaminated_targets.append({
                    'time': ti,
                    'position': target,
                    'name': contamination_data.get('contaminated_names', []),
                    'intensity': contamination_data.get('intensities', []),
                    'ab_magnitude': contamination_data.get('ab_magnitudes', []),
                    'separation': contamination_data.get('separations', [])
                })
        
        # Add separations
        if contamination_data and contamination_data.get('separations'):
            self.all_separations.extend(contamination_data['separations'])

    def increment_counters(self, total_events=0, trail_candidates=0, contaminations=0, no_topo_skip=0, no_brightness_skip=0, no_targets_skip=0):
        """Update exposure counters"""
        self.total_events += total_events
        self.trail_candidates += trail_candidates
        self.contaminations += contaminations
        self.no_topo_skip += no_topo_skip
        self.no_brightness_skip += no_brightness_skip
        self.no_targets_skip += no_targets_skip

    def get_statistics(self):
        """Calculate and return exposure statistics"""
        stats = {
            'exposure_id': self.exposure_id,
            'start_time': self.start_time,
            'end_time': self.end_time,
            'duration_seconds': (self.end_time.utc_datetime() - self.start_time.utc_datetime()).total_seconds(),
            'total_events': self.total_events,
            'trail_candidates': self.trail_candidates,
            'contaminations': self.contaminations,
            'no_topo_skip': self.no_topo_skip,
            'no_brightness_skip': self.no_brightness_skip,
            'no_targets_skip': self.no_targets_skip,
            'total_satellite_observations': len(self.all_satellite_positions),
            'total_contamination_events': len(self.all_contaminated_targets),
            'unique_contaminated_targets': len(set(self.all_contaminated_names))
        }
        
        # Calculate intensity/magnitude statistics
        if self.all_intensities:
            stats.update({
                'peak_intensity': np.max(self.all_intensities),
                'average_intensity': np.mean(self.all_intensities),
                'peak_ab_magnitude': np.min(self.all_ab_magnitudes),  # Lower magnitude = brighter
                'average_ab_magnitude': np.mean(self.all_ab_magnitudes)
            })
        
        # Calculate separation statistics
        if self.all_separations:
            stats.update({
                'closest_approach_arcsec': np.min(self.all_separations),
                'average_separation_arcsec': np.mean(self.all_separations)
            })
        
        return stats
    
    def create_fits_extensions(self):
        """
        Create FITS table extensions for this exposure
        
        Returns:
        --------
        list : List of FITS HDU objects
        """
        extensions = []
        
        # Main satellite data table
        if self.all_satellite_positions:
            satellite_times = [pos['time'].to_astropy() for pos in self.all_satellite_positions]
            satellite_names = [pos['name'] for pos in self.all_satellite_positions]
            satellite_alts = [pos['alt'] for pos in self.all_satellite_positions]
            satellite_azs = [pos['az'] for pos in self.all_satellite_positions]
            satellite_heights = [pos['height'] for pos in self.all_satellite_positions]
            satellite_intensities = [pos['intensity'] for pos in self.all_satellite_positions]
            satellite_magnitudes = [pos['ab_magnitude'] for pos in self.all_satellite_positions]
            satellite_angular_velocities = [pos.get('angular_velocity', np.nan) for pos in self.all_satellite_positions]
            
            main_table = Table([
                satellite_times, satellite_names, satellite_azs, satellite_alts,
                satellite_heights, satellite_intensities, satellite_magnitudes, satellite_angular_velocities
            ], names=[
                'Time', 'Satellite_Name', 'Azimuth', 'Altitude',
                'Height', 'Intensity', 'AB_Magnitude', 'Angular_Velocity'
            ])
            
            main_hdu = fits.BinTableHDU(main_table)
            main_hdu.header['EXTNAME'] = 'SATELLITE_DATA'
            main_hdu.header['EXPID'] = self.exposure_id
            main_hdu.header['EXPSTART'] = self.start_time.utc_iso()
            main_hdu.header['EXPEND'] = self.end_time.utc_iso()
            main_hdu.header['EVENTS'] = self.total_events
            main_hdu.header['TRAILCAN'] = self.trail_candidates
            main_hdu.header['CONTAM'] = self.contaminations
            main_hdu.header['NSAT'] = len(set(satellite_names))
            
            # Add statistics to header
            stats = self.get_statistics()
            if 'peak_intensity' in stats:
                main_hdu.header['PEAKI'] = stats['peak_intensity']
                main_hdu.header['AVGI'] = stats['average_intensity']
                main_hdu.header['PEAKM'] = stats['peak_ab_magnitude']
                main_hdu.header['AVGM'] = stats['average_ab_magnitude']
            
            extensions.append(main_hdu)
        
        # Contaminated targets table
        if self.all_contaminated_targets:
            target_times = [ct['time'].to_astropy() for ct in self.all_contaminated_targets]
            target_names = [ct['name'] for ct in self.all_contaminated_targets]
            target_coords = [ct['position'] for ct in self.all_contaminated_targets]
            target_alts = [coord.alt.degree for coord in target_coords]
            target_azs = [coord.az.degree for coord in target_coords]
            
            contamination_table = Table([target_names,
                target_times, target_alts, target_azs
            ], names=['Target_Name', 'Time', 'Target_Alt', 'Target_Az'])
            
            contamination_hdu = fits.BinTableHDU(contamination_table)
            contamination_hdu.header['EXTNAME'] = 'CONTAMINATED_TARGETS'
            contamination_hdu.header['EXPID'] = self.exposure_id
            # contamination_hdu.header['EXPSTART'] = self.start_time.utc_iso()
            contamination_hdu.header['NCONTAM'] = len(self.all_contaminated_targets)
            contamination_hdu.header['NTARG'] = len(set(self.all_contaminated_names))

            
            extensions.append(contamination_hdu)
        
        # Separations table
        if self.all_separations:
            separations_table = Table([self.all_separations], names=['Separation_arcsec'])
            separations_hdu = fits.BinTableHDU(separations_table)
            separations_hdu.header['EXTNAME'] = 'SEPARATIONS'
            separations_hdu.header['EXPID'] = self.exposure_id
            separations_hdu.header['NSEP'] = len(self.all_separations)
            separations_hdu.header['MINSEP'] = np.min(self.all_separations)
            separations_hdu.header['AVGSEP'] = np.mean(self.all_separations)
            # separations_hdu.header['EXPSTART'] = self.start_time.utc_iso()
            
            extensions.append(separations_hdu)
        
        return extensions
    
    def print_summary(self):
        """Print a summary of this exposure"""
        stats = self.get_statistics()
        print(f"\nExposure {self.exposure_id} Summary:")
        print(f"  Duration: {stats['duration_seconds']:.0f}s")
        print(f"  Events: {stats['total_events']}")
        print(f"  Trail candidates: {stats['trail_candidates']}")
        print(f"  Contaminations: {stats['contaminations']}")
        print(f"  Satellite observations: {stats['total_satellite_observations']}")
        
        if 'peak_intensity' in stats:
            print(f"  Peak intensity: {stats['peak_intensity']:.2e} W/m²/sr")
            print(f"  Peak magnitude: {stats['peak_ab_magnitude']:.2f} AB mag")
            print(f"  Average magnitude: {stats['average_ab_magnitude']:.2f} AB mag")
        
        if 'closest_approach_arcsec' in stats:
            print(f"  Closest approach: {stats['closest_approach_arcsec']:.2f} arcsec")

 