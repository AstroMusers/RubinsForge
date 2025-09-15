import numpy as np
import time
from astropy import units as u
from ..utils.calculations import point_brightness_calculation
from ..satellites import starlink_sat
import lumos.calculator
import lumos.conversions

class SatelliteTracker:
    def __init__(self):
        # Current satellite position (accessible globally until redefined)
        self.current_alt = None
        self.current_az = None
        self.current_height = None
        self.current_topocentric = None
        self.current_time = None
        self.current_satellite_name = None
        
        # Previous satellite position (ti-1)
        self.previous_alt = None
        self.previous_az = None
        self.previous_height = None
        self.previous_topocentric = None
        self.previous_time = None
        self.previous_satellite_name = None
        
        # Position history for analysis
        self.position_history = []
    
    def _update_positions(self, satal, sataz, sat_height, topocentric, ti, satellite_name):
        """
        Internal method to update current and previous positions
        """
        # Move current to previous
        self.previous_alt = self.current_alt
        self.previous_az = self.current_az
        self.previous_height = self.current_height
        self.previous_topocentric = self.current_topocentric
        self.previous_time = self.current_time
        self.previous_satellite_name = self.current_satellite_name
        
        # Update current
        self.current_alt = satal.degrees
        self.current_az = sataz.degrees
        self.current_height = sat_height.to(u.m).value
        self.current_topocentric = topocentric
        self.current_time = ti
        self.current_satellite_name = satellite_name
        
    def trail_event(self, starlink, rubin_obs, ti, sat_t, sat_al, sat_az, sat_hght):
        """
        Modified trail_event that stores current and previous positions
        """
        start = time.time()
        max_duration = 60
        status = False
        
        while (((time.time() - start) % 60) < max_duration) & (status == False):
            difference = starlink - rubin_obs
            topocentric = difference.at(ti)
            satal, sataz, sat_height = topocentric.altaz()
            
            # Update positions (current becomes previous, new becomes current)
            self._update_positions(satal, sataz, sat_height, topocentric, ti, starlink.name)
            
            if (sat_height.to(u.km).value < 2000):
                sat_al.append(self.current_alt)
                sat_az.append(self.current_az)
                sat_hght.append(self.current_height)
                sat_t.append(ti)
                
                # Store in history
                self.position_history.append({
                    'alt': self.current_alt,
                    'az': self.current_az,
                    'height': self.current_height,
                    'time': ti,
                    'satellite': starlink.name,
                    'topocentric': topocentric
                })
                
                status = True
                return topocentric
            else:
                status = True
                return None
    
    def point_trail_event(self, starlink, rubin_obs, ti):
        """
        Modified point_trail_event that stores current and previous positions
        """
        start = time.time()
        max_duration = 60
        status = False
        
        while (((time.time() - start) % 60) < max_duration) & (status == False):
            difference = starlink - rubin_obs
            topocentric = difference.at(ti)
            satal, sataz, sat_height = topocentric.altaz()
            
            # Update positions
            self._update_positions(satal, sataz, sat_height, topocentric, ti, starlink.name)
            
            if (sat_height.to(u.km).value < 2000):
                status = True
                return topocentric
            else:
                status = True
                return None
    
    def brightness_calculation_current(self, rubinobs_astr, sun):
        """
        Calculate brightness using current satellite position
        """
        if self.current_alt is None:
            raise ValueError("No current satellite position available.")
        
        return point_brightness_calculation(
            rubinobs_astr, self.current_time, sun, 
            [self.current_alt], [self.current_az], self.current_height
        )
    
    def brightness_calculation_previous(self, rubinobs_astr, sun):
        """
        Calculate brightness using previous satellite position
        """
        if self.previous_alt is None:
            raise ValueError("No previous satellite position available.")
        
        return point_brightness_calculation(
            rubinobs_astr, self.previous_time, sun, 
            [self.previous_alt], [self.previous_az], self.previous_height
        )
    
    def calculate_angular_velocity(self):
        """
        Calculate angular velocity between previous and current positions
        Returns velocity in arcsec/second
        """
        if self.previous_alt is None or self.current_alt is None:
            return None
        
        from ..utils.calculations import manual_angular_separation
        # Check if both times exist (correct way for Skyfield Time objects)
        if self.previous_time is None or self.current_time is None:
            return None       
        
        # Calculate angular separation
        separation_arcsec = manual_angular_separation(
            self.previous_alt, self.previous_az,
            self.current_alt, self.current_az
        )
        
        # Calculate time difference using Skyfield Time methods
        try:
            time_diff = (self.current_time.utc_datetime() - self.previous_time.utc_datetime()).total_seconds()
            if time_diff > 0:
                return separation_arcsec / time_diff
            else:
                return None
        except Exception as e:
            print(f"Error calculating time difference: {e}")
            return None
    
    def get_current_position(self):
        """
        Get current satellite position
        """
        return {
            'alt': self.current_alt,
            'az': self.current_az,
            'height': self.current_height,
            'topocentric': self.current_topocentric,
            'time': self.current_time,
            'satellite_name': self.current_satellite_name
        }
    
    def get_previous_position(self):
        """
        Get previous satellite position (ti-1)
        """
        return {
            'alt': self.previous_alt,
            'az': self.previous_az,
            'height': self.previous_height,
            'topocentric': self.previous_topocentric,
            'time': self.previous_time,
            'satellite_name': self.previous_satellite_name
        }
    
    def get_both_positions(self):
        """
        Get both current and previous positions
        """
        return {
            'current': self.get_current_position(),
            'previous': self.get_previous_position()
        }
    
    def has_previous_position(self):
        """Check if previous position is available"""
        return self.previous_alt is not None
    
    def has_current_position(self):
        """Check if current position is available"""
        return self.current_alt is not None
    
    def has_both_positions(self):
        """Check if both current and previous positions are available"""
        return self.has_current_position() and self.has_previous_position()
    
    def get_position_history(self):
        """Get full position history"""
        return self.position_history
    
    def clear_history(self):
        """Clear position history"""
        self.position_history = []
    
    def reset_current_position(self):
        """Reset current position to None"""
        self.current_alt = None
        self.current_az = None
        self.current_height = None
        self.current_topocentric = None
        self.current_time = None
        self.current_satellite_name = None
    
    def reset_previous_position(self):
        """Reset previous position to None"""
        self.previous_alt = None
        self.previous_az = None
        self.previous_height = None
        self.previous_topocentric = None
        self.previous_time = None
        self.previous_satellite_name = None
    
    def reset_all_positions(self):
        """Reset both current and previous positions"""
        self.reset_current_position()
        self.reset_previous_position()
    
    def get_altitude_azimuth(self):
        """Quick access to current alt/az"""
        if self.current_alt is None:
            return None, None
        return self.current_alt, self.current_az
    
    def get_previous_altitude_azimuth(self):
        """Quick access to previous alt/az"""
        if self.previous_alt is None:
            return None, None
        return self.previous_alt, self.previous_az
    
    def get_height_km(self):
        """Get current height in kilometers"""
        if self.current_height is None:
            return None
        return self.current_height / 1000.0
    
    def get_previous_height_km(self):
        """Get previous height in kilometers"""
        if self.previous_height is None:
            return None
        return self.previous_height / 1000.0
    
    def print_current_status(self):
        """Print current satellite status"""
        if self.current_alt is None:
            print("No current satellite position available")
        else:
            print(f"Current satellite: {self.current_satellite_name}")
            print(f"  Current: Alt={self.current_alt:.2f}°, Az={self.current_az:.2f}°, Height={self.get_height_km():.1f}km")
            
            if self.previous_alt is not None:
                print(f"  Previous: Alt={self.previous_alt:.2f}°, Az={self.previous_az:.2f}°, Height={self.get_previous_height_km():.1f}km")
                
                velocity = self.calculate_angular_velocity()
                if velocity:
                    print(f"  Angular velocity: {velocity:.2f} arcsec/s")
            else:
                print("  No previous position available")
    
    @staticmethod
    def find_satellites_up(satellites, observer, start_time, end_time, altitude_degrees=30.0, check_sunlit=True, eph=None, return_periods=False):
        """
        Find satellites that are above horizon and optionally sunlit during time period
        
        Parameters:
        -----------
        satellites : list
            List of Skyfield satellite objects
        observer : Skyfield observer position
            Observer position (e.g., rubin_obs)
        start_time : Skyfield time
            Start time for checking
        end_time : Skyfield time  
            End time for checking
        altitude_degrees : float
            Minimum altitude in degrees (default 30.0)
        check_sunlit : bool
            Whether to check if satellite is sunlit (default True)
        eph : Skyfield ephemeris
            Ephemeris for sunlit check (required if check_sunlit=True)
        return_periods : bool
            If True, return detailed period information for each satellite
            
        Returns:
        --------
        If return_periods=False:
            list : Satellites that meet the criteria
        If return_periods=True:
            tuple : (satellites_up, satellite_periods)
                satellites_up: list of satellites
                satellite_periods: dict with satellite -> period info
        """
        satellites_up = []
        satellite_periods = {}
        
        for satellite in satellites:
            try:
                # Get rise/set events during the time period
                times, events = satellite.find_events(
                    observer, start_time, end_time, 
                    altitude_degrees=altitude_degrees
                )
                
                has_events = len(times) > 0
                
                # Check if satellite is sunlit at start time (if requested)
                is_sunlit_start = True
                if check_sunlit and eph is not None:
                    is_sunlit_start = satellite.at(start_time).is_sunlit(eph)
                
                # Check altitude at start time
                difference = satellite - observer
                topocentric = difference.at(start_time)
                alt, az, distance = topocentric.altaz()
                above_horizon_start = alt.degrees > altitude_degrees
                
                # If satellite meets basic criteria, calculate detailed periods
                if has_events and is_sunlit_start and above_horizon_start:
                    satellites_up.append(satellite)
                    
                    if return_periods:
                        # Calculate visible and sunlit periods
                        visible_periods = []
                        sunlit_periods = []
                        
                        # Process events to find visible periods
                        if len(times) > 0:
                            current_period_start = None
                            
                            for i, (time, event) in enumerate(zip(times, events)):
                                if event == 0:  # Rise event
                                    current_period_start = time
                                elif event == 2:  # Set event
                                    if current_period_start is not None:
                                        visible_periods.append((current_period_start, time))
                                    current_period_start = None
                            
                            # If satellite is still up at end of period
                            if current_period_start is not None:
                                visible_periods.append((current_period_start, end_time))
                        
                        # Check if satellite is above horizon at start of observation
                        if above_horizon_start:
                            if not visible_periods or visible_periods[0][0] > start_time:
                                visible_periods.insert(0, (start_time, visible_periods[0][0] if visible_periods else end_time))
                        
                        # Find sunlit periods within visible periods
                        if check_sunlit and eph is not None:
                            for vis_start, vis_end in visible_periods:
                                # Sample points to check sunlit status
                                if satellite.at(vis_start).is_sunlit(eph) & satellite.at(vis_end).is_sunlit(eph):
                                    sunlit_periods.append((vis_start, vis_end))
                                    continue
                                duration = (vis_end.utc_datetime() - vis_start.utc_datetime()).total_seconds()
                                num_samples = max(10, int(duration / 60))  # Sample every minute minimum
                                
                                sample_times = [vis_start + (vis_end - vis_start) * i / num_samples 
                                                for i in range(num_samples + 1)]
                                
                                sunlit_start = None
                                for sample_time in sample_times:
                                    is_sunlit = satellite.at(sample_time).is_sunlit(eph)
                                    
                                    if is_sunlit and sunlit_start is None:
                                        sunlit_start = sample_time
                                    elif not is_sunlit and sunlit_start is not None:
                                        sunlit_periods.append((sunlit_start, sample_time))
                                        sunlit_start = None
                                
                                # If still sunlit at end
                                if sunlit_start is not None:
                                    sunlit_periods.append((sunlit_start, vis_end))
                        else:
                            sunlit_periods = visible_periods
                        
                        satellite_periods[satellite] = {
                            'visible_periods': visible_periods,
                            'sunlit_periods': sunlit_periods,
                            'combined_periods': sunlit_periods if check_sunlit else visible_periods,
                            'total_visible_duration': sum([
                                (end.utc_datetime() - start.utc_datetime()).total_seconds() 
                                for start, end in visible_periods
                            ]),
                            'total_sunlit_duration': sum([
                                (end.utc_datetime() - start.utc_datetime()).total_seconds() 
                                for start, end in sunlit_periods
                            ]) if sunlit_periods else 0
                        }
                        
            except Exception as e:
                print(f"Error checking satellite {satellite.name}: {e}")
                continue
        
        if return_periods:
            return satellites_up, satellite_periods
        else:
            return satellites_up
    
    @staticmethod
    def find_satellites_up_simple(satellites, observer, check_time, altitude_degrees=30.0, check_sunlit=True, eph=None):
        """
        Simplified version - check satellites at a single time point
        
        Parameters:
        -----------
        satellites : list
            List of Skyfield satellite objects
        observer : Skyfield observer position
        check_time : Skyfield time
            Time to check satellite positions
        altitude_degrees : float
            Minimum altitude in degrees
        check_sunlit : bool
            Whether to check if satellite is sunlit
        eph : Skyfield ephemeris
            Ephemeris for sunlit check
            
        Returns:
        --------
        list : Satellites above horizon (and sunlit if requested)
        """
        satellites_up = []
        
        for satellite in satellites:
            try:
                # Check altitude
                difference = satellite - observer
                topocentric = difference.at(check_time)
                alt, az, distance = topocentric.altaz()
                
                if alt.degrees > altitude_degrees:
                    # Check if sunlit (if requested)
                    if check_sunlit and eph is not None:
                        if satellite.at(check_time).is_sunlit(eph):
                            satellites_up.append(satellite)
                    else:
                        satellites_up.append(satellite)
                        
            except Exception as e:
                print(f"Error checking satellite {satellite.name}: {e}")
                continue
        
        return satellites_up
    
    @staticmethod
    def get_satellite_count_info(satellites, observer, check_time, altitude_degrees=30.0, eph=None):
        """
        Get detailed count information about satellites
        
        Returns:
        --------
        dict : Satellite counts and statistics
        """
        total_satellites = len(satellites)
        above_horizon = 0
        sunlit = 0
        above_horizon_and_sunlit = 0
        
        for satellite in satellites:
            try:
                # Check altitude
                difference = satellite - observer
                topocentric = difference.at(check_time)
                alt, az, distance = topocentric.altaz()
                
                is_above_horizon = alt.degrees > altitude_degrees
                is_sunlit = satellite.at(check_time).is_sunlit(eph) if eph else False
                
                if is_above_horizon:
                    above_horizon += 1
                    
                if is_sunlit:
                    sunlit += 1
                    
                if is_above_horizon and is_sunlit:
                    above_horizon_and_sunlit += 1
                    
            except Exception as e:
                continue
        
        return {
            'total': total_satellites,
            'above_horizon': above_horizon,
            'sunlit': sunlit,
            'above_horizon_and_sunlit': above_horizon_and_sunlit,
            'time': check_time
        }
    
    @staticmethod
    def is_time_in_periods(check_time, periods):
        """
        Check if a time falls within any of the given periods
        
        Parameters:
        -----------
        check_time : Skyfield time
            Time to check
        periods : list of tuples
            List of (start_time, end_time) tuples
            
        Returns:
        --------
        bool : True if time falls within any period
        """
        check_datetime = check_time.utc_datetime()
        
        for start_time, end_time in periods:
            start_datetime = start_time.utc_datetime()
            end_datetime = end_time.utc_datetime()
            
            if start_datetime <= check_datetime <= end_datetime:
                return True
        
        return False