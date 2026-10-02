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
        

    def calculate_targets_for_exposure(
        self,
        exposure_times,
        observer,
        zone,
    ):
        """
        Calculate target positions separately at every time step of an exposure.
        """

        print(
            f"Calculating targets for exposure from "
            f"{exposure_times[0].astimezone(zone)} to "
            f"{exposure_times[-1].astimezone(zone)}"
        )

        calc_start = time.time()

        self.time_cache = {}
        self.exposure_times = exposure_times

        # ---------------------------------------------------------
        # 1. Decide which targets should be followed in this exposure
        # ---------------------------------------------------------


        targets= list(self.targets)
        names = list(self.names)

        print(
            f"Found {len(targets)} targets. "
        )

        # ---------------------------------------------------------
        # 2. Actually calculate positions independently at each time
        # ---------------------------------------------------------

        for ti in exposure_times:

            rubin_at_t = observer.at(ti)

            apparent_positions = []
            skycoords = []
            alt_az_tuples = []

            for target in targets:

                topocentric = rubin_at_t.observe(target)
                apparent = topocentric.apparent()

                # IMPORTANT:
                # no refraction here, to stay consistent with your
                # satellite .altaz() calls
                alt, az, _ = apparent.altaz()

                apparent_positions.append(apparent)

                skycoords.append(
                    SkyCoord(
                        alt=alt.degrees * u.degree,
                        az=az.degrees * u.degree,
                        frame="altaz",
                        obstime=ti.to_astropy(),
                        location=self.rubin_location_geocentric,
                    )
                )

                alt_az_tuples.append(
                    (alt.degrees, az.degrees)
                )

            self.time_cache[ti] = {
                "targets": targets,
                "names": names,
                "apparent_positions": apparent_positions,
                "skycoords": skycoords,
                "alt_az_tuples": alt_az_tuples,
            }

        self.current_exposure_targets = targets

        print(
            f"Target calculation completed in "
            f"{time.time() - calc_start:.2f} seconds "
            f"for {len(exposure_times)} time steps"
        )

        return targets
    
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
    
    def get_targets(self, ti=None):
        """Get list of targets that are up (uses current exposure targets)"""
        if ti is not None and ti in self.time_cache:
            return self.time_cache[ti]['targets']
        return self.current_exposure_targets
    
    def get_names(self, ti=None):
        """Get list of target names that are up (uses current exposure targets)"""
        if ti is not None and ti in self.time_cache:
            return self.time_cache[ti]['names']
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
            'targets': len(self.current_exposure_targets) if self.current_exposure_targets else 0
        }
    
    def print_cache_status(self):
        """Print current cache status"""
        info = self.get_cache_info()
        print(f"TargetTracker cache: {info['cached_times']} times cached, {info['targets']} targets up")