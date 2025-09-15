import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord, EarthLocation
import time

class TargetTracker:
    def __init__(self, targets, names, rubin_location_astrometric, rubin_location_geocentric=None):
        """
        Initialize target tracker
        
        Parameters:
        -----------
        targets : list
            List of Skyfield Star objects
        rubin_location : EarthLocation
            Astropy EarthLocation for Rubin Observatory
        """
        self.targets = targets
        self.names = names
        self.rubin_location_geocentric = rubin_location_geocentric or EarthLocation.of_site('rubin')
        self.rubin_location_astrometric = rubin_location_astrometric 

        # Storage for target positions at each time
        self.time_cache = {}  # ti -> target positions
        self.exposure_times = []
        self.current_exposure_targets = None

    def calculate_targets_for_exposure(self, exposure_times, observer, zone, horizon_degrees=30.0):
        """
        Pre-calculate all target positions for an entire exposure
        
        Parameters:
        -----------
        exposure_times : list
            List of Skyfield time objects for the exposure
        observer : Skyfield observer position
        horizon_degrees : float
            Minimum altitude for targets to be considered "up"
        """
        print(f"Calculating targets for exposure from {exposure_times[0].astimezone(zone)} to {exposure_times[-1].astimezone(zone)}")
        calc_start = time.time()
        
        # Clear previous cache
        self.time_cache = {}
        self.exposure_times = exposure_times
        
        # Filter targets that are up during this exposure (check first time only)
        first_time = exposure_times[0]
        targets_up = []
        names_up = []

        for target, name in zip(self.targets, self.names):
            # Check if target is above horizon at start of exposure
            target_astrometric = observer.at(first_time).observe(target)
            target_apparent = target_astrometric.apparent()
            target_alt, target_az, target_distance = target_apparent.altaz()
            
            if target_alt.degrees > horizon_degrees:
                targets_up.append(target)
                names_up.append(name)

        print(f"Found {len(targets_up)} targets above {horizon_degrees}° horizon")
        
        # Calculate positions for each time in exposure
        for ti in exposure_times:
            rubin_at_time = observer.at(ti)
            
            # Calculate apparent positions for all up targets
            topocentric_positions = [rubin_at_time.observe(target) for target in targets_up]
            apparent_positions = [topo.apparent() for topo in topocentric_positions]
            
            # Get alt/az coordinates
            alt_az_data = [pos.altaz(temperature_C=12.0, pressure_mbar=750) for pos in apparent_positions]
            alt_t, az_t, height_t = zip(*alt_az_data) if alt_az_data else ([], [], [])
            
            # Create SkyCoord objects for easy manipulation
            skycoords = [
                SkyCoord(
                    alt=alt.degrees*u.degree, 
                    az=az.degrees*u.degree, 
                    frame='altaz', 
                    unit='deg', 
                    obstime=ti.to_astropy(), 
                    location=self.rubin_location_geocentric
                ) for alt, az in zip(alt_t, az_t)
            ] if alt_t else []
            
            # Store in cache
            self.time_cache[ti] = {
                'targets_up': targets_up,
                'names_up': names_up,
                'apparent_positions': apparent_positions,
                'skycoords': skycoords,
                'alt_az_tuples': list(zip(alt_t, az_t)) if alt_t else []
            }
        
        calc_time = time.time() - calc_start
        print(f"Target calculation completed in {calc_time:.2f} seconds for {len(exposure_times)} time steps")
        
        # Store current exposure info
        self.current_exposure_targets = targets_up
        
        return targets_up
    
    def get_targets_at_time(self, ti):
        """
        Get pre-calculated target positions for a specific time
        
        Parameters:
        -----------
        ti : Skyfield time object
        
        Returns:
        --------
        dict : Target position data for this time
        """
        if ti not in self.time_cache:
            raise ValueError(f"Target positions not calculated for time {ti}. Call calculate_targets_for_exposure first.")
        
        return self.time_cache[ti]
    
    def get_apparent_positions(self, ti):
        """Get apparent positions for time ti"""
        return self.time_cache[ti]['apparent_positions']
    
    def get_skycoords(self, ti):
        """Get SkyCoord objects for time ti"""
        return self.time_cache[ti]['skycoords']
    
    def get_targets_up(self, ti=None):
        """Get list of targets that are up (uses current exposure targets)"""
        if ti is not None and ti in self.time_cache:
            return self.time_cache[ti]['targets_up']
        return self.current_exposure_targets
    
    def get_names_up(self, ti=None):
        """Get list of target names that are up (uses current exposure targets)"""
        if ti is not None and ti in self.time_cache:
            return self.time_cache[ti]['names_up']
        return self.names
    def has_time_cached(self, ti):
        """Check if time is in cache"""
        return ti in self.time_cache
    
    def clear_cache(self):
        """Clear the time cache"""
        self.time_cache = {}
        self.exposure_times = []
        self.current_exposure_targets = None
    
    def get_cache_info(self):
        """Get information about current cache"""
        return {
            'cached_times': len(self.time_cache),
            'exposure_duration': len(self.exposure_times),
            'targets_up': len(self.current_exposure_targets) if self.current_exposure_targets else 0
        }
    
    def print_cache_status(self):
        """Print current cache status"""
        info = self.get_cache_info()
        print(f"TargetTracker cache: {info['cached_times']} times cached, {info['targets_up']} targets up")