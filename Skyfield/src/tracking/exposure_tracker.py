import numpy as np
from astropy.table import Table
from astropy.io import fits
from astropy.coordinates import SkyCoord, AltAz, EarthLocation
from astropy.visualization.wcsaxes import SphericalCircle
from skyfield.positionlib import ICRF, position_of_radec
from astropy import units as u
from skyfield.api import wgs84, load, Star, Angle, S, W

eph = load('de421.bsp')
earth = eph['earth']
sun = eph['sun']

class ExposureTracker:
    def __init__(self, exposure_id, start_time, end_time, visit_info=None):
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
        visit_info : dict, optional
            Additional visit information (e.g., field center, seeing, sky brightness)
        """
        self.exposure_id = exposure_id
        self.start_time = start_time
        self.end_time = end_time
        self.rubin_location_geocentric_astropy = EarthLocation.of_site('rubin')
        self.rubin_location_geocentric =  wgs84.latlon(-30.244633,  -70.749417)
        self.rubin_location_astrometric = earth +  wgs84.latlon(30.244633*S,  70.749417*W, elevation_m = 2647)
        self.visit_info = visit_info

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

    def get_visit_region_skycoord_object(self):
        """Return the visit region information if available"""
        if self.visit_info:
            center_ra = self.visit_info.get('field_ra', None)
            center_dec = self.visit_info.get('field_dec', None)
            if center_ra is not None and center_dec is not None:
                skycoord_object =  SkyCoord(ra=center_ra*u.deg, dec=center_dec*u.deg, frame='icrs').transform_to(AltAz(obstime=self.start_time.to_astropy(), location=self.rubin_location_geocentric_astropy))
                return skycoord_object
        return None
    
    def get_visit_region_apparent(self):
        """Return the visit region center as an apparent position if available"""
        if self.visit_info:

            apparent = self.rubin_location_astrometric.at(self.start_time).observe(self.visit_info.get('visit_center_as_star')).apparent()
        return apparent

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
        
        # Visit information table (add this first)
        if self.visit_info:
            # Create a single-row table with visit information
            visit_data = {}
            
            # Standard visit fields
            standard_fields = {
                'visit_id': 'VISIT_ID',
                'field_ra': 'FIELD_RA',
                'field_dec': 'FIELD_DEC',
                'filter': 'FILTER',
                'exposure_time': 'EXPOSURE_TIME',
                'num_exposures': 'NUM_EXPOSURES',
                'seeing': 'SEEING',
                'sky_brightness': 'SKY_BRIGHTNESS',

            }
            
            # Extract available visit info fields
            for key, column_name in standard_fields.items():
                if key in self.visit_info:
                    value = self.visit_info[key]
                    
                    # Handle different data types
                    if isinstance(value, str):
                        visit_data[column_name] = [value]
                    elif isinstance(value, (int, float)):
                        visit_data[column_name] = [float(value)]
                    else:
                        visit_data[column_name] = [str(value)]
            
            # Add exposure-specific fields
            visit_data['EXPOSURE_ID'] = [self.exposure_id]
            visit_data['EXP_START'] = [self.start_time.utc_iso()]
            visit_data['EXP_END'] = [self.end_time.utc_iso()]
            visit_data['EXP_DURATION'] = [(self.end_time.utc_datetime() - self.start_time.utc_datetime()).total_seconds()]
            
            # Calculate field center in alt/az if possible
            if 'field_ra' in self.visit_info and 'field_dec' in self.visit_info:
                try:
                    field_center = self.get_visit_region_skycoord_object()
                    if field_center:
                        visit_data['FIELD_ALT'] = [float(field_center.alt.degree)]
                        visit_data['FIELD_AZ'] = [float(field_center.az.degree)]
                except Exception as e:
                    print(f"Warning: Could not calculate field alt/az: {e}")
            
            # Create the table
            if visit_data:
                visit_table = Table(visit_data)
                
                # Set appropriate data types
                for col_name in visit_table.colnames:
                    if col_name in ['VISIT_ID', 'EXPOSURE_ID', 'NUM_EXPOSURES']:
                        visit_table[col_name] = visit_table[col_name].astype('int64')
                    elif col_name in ['EXP_START', 'EXP_END', 'FILTER']:
                        visit_table[col_name] = visit_table[col_name].astype('U25')
                    else:
                        visit_table[col_name] = visit_table[col_name].astype('float64')
                
                visit_hdu = fits.BinTableHDU(visit_table)
                visit_hdu.header['EXTNAME'] = 'VISIT_INFO'
                visit_hdu.header['EXPID'] = self.exposure_id
                
                # Add visit info to header as well
                if 'visit_id' in self.visit_info:
                    visit_hdu.header['VISITID'] = self.visit_info['visit_id']
                if 'field_ra' in self.visit_info:
                    visit_hdu.header['FIELDRA'] = self.visit_info['field_ra']
                if 'field_dec' in self.visit_info:
                    visit_hdu.header['FIELDDEC'] = self.visit_info['field_dec']
                if 'filter' in self.visit_info:
                    visit_hdu.header['FILTER'] = self.visit_info['filter']
                
                extensions.append(visit_hdu)
        
        # Main satellite data table
        if self.all_satellite_positions:
            if self.trail_candidates:
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
            
            # Add visit info to satellite data header as well
            if self.visit_info:
                if 'visit_id' in self.visit_info:
                    main_hdu.header['VISITID'] = self.visit_info['visit_id']
                if 'filter' in self.visit_info:
                    main_hdu.header['FILTER'] = self.visit_info['filter']
            
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
            contamination_hdu.header['NCONTAM'] = len(self.all_contaminated_targets)
            contamination_hdu.header['NTARG'] = len(set(self.all_contaminated_names))
            
            # Add visit info to contamination header
            if self.visit_info and 'visit_id' in self.visit_info:
                contamination_hdu.header['VISITID'] = self.visit_info['visit_id']
            
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
            
            # Add visit info to separations header
            if self.visit_info and 'visit_id' in self.visit_info:
                separations_hdu.header['VISITID'] = self.visit_info['visit_id']
            
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

 