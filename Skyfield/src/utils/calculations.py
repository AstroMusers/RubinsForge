import os
import sys
import logging
import numpy as np
from astropy import units as u
from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord
import skyfield as sf
from skyfield.api import load, wgs84, EarthSatellite, N, W,S,E, Star, Angle
import datetime as dt
import time
from pytz import timezone
from skyfield import almanac
from skyfield.framelib import galactic_frame
from skyfield.data import hipparcos
import lumos.calculator
import random
from ..satellites import starlink_sat
import lumos.plot
import lumos.brdf.library



    
def brightness_calculation(rubinobs_astr, ti, sun, s, f, sat_al, sat_az, sat_hght):  
    start = time.time()
    max_duration = 60
    status = False
    while (((time.time() - start) % 60) < max_duration) & (status == False):
        astro = rubinobs_astr.at(ti).observe(sun)
        app = astro.apparent()
        sunalt, sunaz, sundistance = app.altaz()
        sun_al = sunalt.degrees
        sun_az = sunaz.degrees

        sat_altitudes, sat_azimuths = \
            np.meshgrid(
            sat_al[s:f], 
            sat_az[s:f])
        sat_heights = sat_hght[s:f] # in meters
        intensities = lumos.calculator.get_intensity_observer_frame(
                starlink_sat.SURFACES_LAB_BRDFS, sat_heights, sat_altitudes, sat_azimuths,
                sun_altitude = sun_al, sun_azimuth = sun_az,
                include_earthshine = False)

        # Convert intensity to AB Magnitude
        ab_magnitudes = lumos.conversions.intensity_to_ab_mag(intensities)
        #print(f'task brightness done in {((time.time() - start) % 60)} seconds')
        status = True
        return intensities, ab_magnitudes
    else:
        return None
    
def point_brightness_calculation(rubinobs_astr, ti, sun, al, az, height):  
    start = time.time()
    max_duration = 60
    status = False
    while (((time.time() - start) % 60) < max_duration) & (status == False):
        astro = rubinobs_astr.at(ti).observe(sun)
        app = astro.apparent()
        sunalt, sunaz, sundistance = app.altaz()
        sun_al = sunalt.degrees
        sun_az = sunaz.degrees
        # sun_al, sun_az = lumos.calculator.get_sun_alt_az(rubinobs_astr, ti.to_astropy())


        sat_altitudes, sat_azimuths = \
            np.meshgrid(
            al, 
            az)
        sat_heights = height # in meters
        intensity = lumos.calculator.get_intensity_observer_frame(
                starlink_sat.SURFACES_LAB_BRDFS, sat_heights, sat_altitudes, sat_azimuths,
                sun_altitude = sun_al, sun_azimuth = sun_az,
                include_earthshine = False)

        # Convert intensity to AB Magnitude
        ab_magnitude = lumos.conversions.intensity_to_ab_mag(intensity)
        #print(f'task brightness done in {((time.time() - start) % 60)} seconds')
        status = True
        return intensity[0][0], ab_magnitude[0][0]
    else:
        return None 


def trail_event(starlink, rubin_obs, ti, sat_t, sat_al, sat_az, sat_hght):
    start = time.time()
    max_duration = 60
    status = False
    while (((time.time() - start) % 60) < max_duration) & (status == False):
        difference = starlink - rubin_obs
        topocentric = difference.at(ti)
        # satra, satdec, distance = topocentric.radec() # ICRF ("J2000")
        satal, sataz, sat_height = topocentric.altaz()
        if (sat_height.to(u.km).value < 2000):
            # sat_ra.append(satra.hours)
            # sat_dec.append(satdec.degrees)
            # sat_dist.append(distance.to(u.m).value)
            # satal, sataz, sat_height = topocentric.altaz()
            sat_al.append(satal.degrees)
            sat_az.append(sataz.degrees)
            sat_hght.append(sat_height.to(u.m).value)

            # diff = ti.astimezone(zone) - midnight.astimezone(zone)
            # sat_t.append(diff.total_seconds() / 3600)
            sat_t.append(ti)
            status = True
            return topocentric
        else:
            status = True
            return None

def point_trail_event(starlink, rubin_obs, ti, return_type=None):
    start = time.time()
    max_duration = 60
    status = False
    while (((time.time() - start) % 60) < max_duration) & (status == False):
        difference = starlink - rubin_obs
        topocentric = difference.at(ti)
        # satra, satdec, distance = topocentric.radec() # ICRF ("J2000")
        satal, sataz, sat_height = topocentric.altaz()
        if (sat_height.to(u.km).value < 2000):
            # sat_ra.append(satra.hours)
            # sat_dec.append(satdec.degrees)
            # sat_dist.append(distance.to(u.m).value)
            # satal, sataz, sat_height = topocentric.altaz()
            # al = satal.degrees
            # az = sataz.degrees
            # height = sat_height.to(u.m).value

            # # diff = ti.astimezone(zone) - midnight.astimezone(zone)
            # # sat_t.append(diff.total_seconds() / 3600)
            # t=ti
            status = True
            return topocentric if return_type is None else (satal.degrees, sataz.degrees)
        else:
            status = True
            return None
        
def target_check(target_apr, target_altaz, topocentric, separation, treshold):
    start = time.time()
    max_duration = 3
    status = False
    while (((time.time() - start) % 60) < max_duration) & (status == False):

        difference_angle = topocentric.separation_from(target_apr)
        separation.append(difference_angle.arcseconds())

        if (difference_angle.arcseconds() <= treshold):
            print(f'Contaminated targets alt/az : {target_apr.alt.degrees} , {target_apr.az.degrees}')
            # L.append(f'Contaminated targets alt/az : {target_apr.alt.degrees} , {target_apr.az.degrees}')

            print(f'Contamination difference angle  : {difference_angle.arcseconds()}')
            # L.append(f'Contamination difference angle  : {difference_angle.arcseconds()}')
            print(f'target check done{((time.time() - start))} seconds')
            status = True
            return True
        else:
            status = True
            return False


def filter_targets(targets_altaz, targets_altaz_pre, topocentric, topocentric_pre):
    """
    Vectorized target filtering for better performance

    Parameters:
    -----------
    targets_altaz : list of SkyCoord
        Current target positions
    targets_altaz_pre : list of SkyCoord
        Previous target positions
    box_bounds : tuple
        (max_alt, min_alt, max_az, min_az)

    Returns:
    --------
    numpy.ndarray : Boolean mask for targets within box

    """
    pre_alt, pre_az, pre_height = topocentric_pre.altaz()
    alt, az, height = topocentric.altaz()
    pre_alt, pre_az = pre_alt.degrees, pre_az.degrees
    alt, az = alt.degrees, az.degrees
    #print(f'Box area: Altitude ({min(pre_alt, alt)} to {max(pre_alt, alt)}), Azimuth ({min(pre_az, az)} to {max(pre_az, az)})')
    max_alt, min_alt, max_az, min_az = (max(pre_alt, alt), min(pre_alt, alt), max(pre_az, az) , min(pre_az, az))

    # Extract coordinates as arrays
    current_alts = np.array([t.alt.degree for t in targets_altaz])
    current_azs = np.array([t.az.degree for t in targets_altaz])
    prev_alts = np.array([t.alt.degree for t in targets_altaz_pre])
    prev_azs = np.array([t.az.degree for t in targets_altaz_pre])

    # Check altitude bounds
    current_alt_mask = (current_alts >= min_alt) & (current_alts <= max_alt)
    prev_alt_mask = (prev_alts >= min_alt) & (prev_alts <= max_alt)

    # Check azimuth bounds (handle wraparound)
    if max_az >= min_az:
        # Normal case
        current_az_mask = (current_azs >= min_az) & (current_azs <= max_az)
        prev_az_mask = (prev_azs >= min_az) & (prev_azs <= max_az)
    else:
        # Wraparound case
        current_az_mask = (current_azs >= min_az) | (current_azs <= max_az)
        prev_az_mask = (prev_azs >= min_az) | (prev_azs <= max_az)

    # Combine conditions
    current_in_box = current_alt_mask & current_az_mask
    prev_in_box = prev_alt_mask & prev_az_mask

    # Target is included if either current or previous position is in box
    return current_in_box | prev_in_box

def check_visit_region(visit_skycoord, topocentric, topocentric_pre):
    """
    Check if a satellite's path intersects the visit region

    Parameters:
    -----------
    visit : dict
        Visit information containing skycoord object
    topocentric : Skyfield position
        Current satellite position
    topocentric_pre : Skyfield position
        Previous satellite position

    Returns:
    --------
    bool : True if satellite path intersects visit region, else False
    """
    try:
        pre_alt, pre_az, pre_height = topocentric_pre.altaz()
        alt, az, height = topocentric.altaz()
        pre_alt, pre_az = pre_alt.degrees, pre_az.degrees
        alt, az = alt.degrees, az.degrees

        radius = 1.75  # degrees

        # Check if either current or previous position is inside the visit region
        current_in_region = topocentric.separation_from(visit_skycoord).degrees <= radius
        previous_in_region = topocentric_pre.separation_from(visit_skycoord).degrees <= radius

        return current_in_region or previous_in_region

    except Exception as e:
        print(f"Error in check_visit_region: {e}")
        return False


def calculated_target_check(current_targets,previous_targets,topocentric,topocentric_pre, separation, treshold):

    max_duration = 10
    status = False
    contaminated = []
    contaminated_names = []
    i = 0

    targets_boxed_mask = filter_targets(current_targets['skycoords'], previous_targets['skycoords'], topocentric, topocentric_pre)
    #print(f'Number of targets in box: {np.sum(targets_boxed_mask)}')
    if np.sum(targets_boxed_mask) == 0:
        return False, contaminated, contaminated_names

    # Convert to numpy array for boolean indexing
    names = np.array(current_targets['names_up'], dtype=object)
    visible_targets_array = np.array(current_targets['apparent_positions'], dtype=object)
    pre_visible_targets_array = np.array(previous_targets['apparent_positions'], dtype=object)
    visible_targets_altaz_array = np.array(current_targets['skycoords'], dtype=object)
    pre_visible_targets_altaz_array = np.array(previous_targets['skycoords'], dtype=object)


    # Now you can use boolean indexing
    print(f'filtered visible targets : {len(visible_targets_array[targets_boxed_mask])} out of {len(visible_targets_array)}')

    # Use the filtered arrays for processing
    visible_targets_appr = visible_targets_array[targets_boxed_mask].tolist()
    pre_visible_targets_appr = pre_visible_targets_array[targets_boxed_mask].tolist()
    visible_targets_altaz = visible_targets_altaz_array[targets_boxed_mask].tolist()
    pre_visible_targets_altaz = pre_visible_targets_altaz_array[targets_boxed_mask].tolist()
    names_up = names[targets_boxed_mask].tolist()

    for name_up,target_apr, target_apr_pre, target_altaz, target_altaz_pre in zip(names_up,visible_targets_appr, pre_visible_targets_appr, visible_targets_altaz, pre_visible_targets_altaz):
        start = time.time()
        while (((time.time() - start) < max_duration) & (status == False)):


            difference_angle_sat_target = topocentric.separation_from(target_apr)

            difference_angle_presat_target = topocentric_pre.separation_from(target_apr)

            difference_angle_sat_pretarget = topocentric.separation_from(target_apr_pre)

            difference_angle_presat_pretarget = topocentric_pre.separation_from(target_apr_pre)

            difference_angle_sat_presat = topocentric.separation_from(topocentric_pre)

            difference_angle_target_pretarget = target_apr.separation_from(target_apr_pre)

            difference_angle_target_astropy = target_altaz.separation(target_altaz_pre)

            separation.append(np.mean([difference_angle_sat_target.arcseconds(), difference_angle_presat_pretarget.arcseconds()]))

            pre_cond = difference_angle_sat_pretarget.arcseconds() + difference_angle_presat_pretarget.arcseconds() <= difference_angle_sat_presat.arcseconds() + treshold
            post_cond = difference_angle_sat_target.arcseconds() + difference_angle_presat_target.arcseconds() <= difference_angle_sat_presat.arcseconds() + treshold
            print(f'Precond check : {difference_angle_sat_pretarget.arcseconds()} + {difference_angle_presat_pretarget.arcseconds()} <= {difference_angle_sat_presat.arcseconds()} + {treshold} is {pre_cond}')
            print(f'Postcond check : {difference_angle_sat_target.arcseconds()} + {difference_angle_presat_target.arcseconds()} <= {difference_angle_sat_presat.arcseconds()} + {treshold} is {post_cond}')

            if pre_cond & post_cond:
                print(f'Contaminated target name : {name_up} ra : {target_apr.radec()[0].degrees} , dec : {target_apr.radec()[1].degrees}')
                print(f'Contaminated target alt/az : {target_altaz.alt.degree} , {target_altaz.az.degree}, pre alt/az : {target_altaz_pre.alt.degree} , {target_altaz_pre.az.degree}')
                #L.append(f'Contaminated targets alt/az : {target_altaz.alt.degree} , {target_altaz.az.degree}, pre alt/az : {target_altaz_pre.alt.degree} , {target_altaz_pre.az.degree}')
                print(f'Contaminating satellite alt/az : {topocentric.altaz()[0].degrees} , {topocentric.altaz()[1].degrees}, pre alt/az : {topocentric_pre.altaz()[0].degrees} , {topocentric_pre.altaz()[1].degrees}')
                print(f'Precond check : {difference_angle_sat_pretarget.arcseconds()} + {difference_angle_presat_pretarget.arcseconds()} <= {difference_angle_sat_presat.arcseconds()} + {treshold} is {pre_cond}')
                print(f'Postcond check : {difference_angle_sat_target.arcseconds()} + {difference_angle_presat_target.arcseconds()} <= {difference_angle_sat_presat.arcseconds()} + {treshold} is {post_cond}')
                print(f'Target check done {((time.time() - start))} seconds')
                contaminated.append(target_altaz)
                contaminated_names.append(name_up)
                break

            else:

                break
            
        else:
            print(f'Target is taking too long, skipping: {time.time() - start}')

    return len(contaminated) > 0, contaminated, contaminated_names

def calculate_target_positions(fits_file, deltatime, deltaday,local_start, local_end, rubinobs_astr, pathToTargets):
    # Calculate target positions for the given time range
    targets = []
    names = []
    times = []
    ti = local_start
    while ti < local_end:
        next_day = ti + dt.timedelta(days=1)
        while ti < next_day:
            times.append(ti)
            ti += deltatime

        ti += deltaday

    n_times = len(times)
     
    with fits.open(pathToTargets) as hdul:
        data = hdul[1].data
        ra = data['Ra']
        dec = data['Dec']
        name = data['Name']
        for i in range(0,10):

            target_i = Star(Angle(degrees=ra[i]),Angle(degrees=dec[i]))
            target_i_name = name[i]

            targets.append(target_i)
            names.append(target_i_name)

    hdul.close()

    n_targets = len(targets)
    fields = ['alt', 'az', 'distance_m', 'visible']
    n_fields = len(fields)
    print(f'len times: {n_times}, len targets: {n_targets}')

    cube_a = np.zeros((n_times, n_targets, n_fields), dtype=float)
    cube_b = np.zeros((n_times, n_targets, 1), dtype=object)
    for t_idx, ti in enumerate(times):
        rubin = rubinobs_astr.at(ti)
        for s_idx, star in enumerate(targets):
            # print(f'Calculating target {s_idx+1}/{n_targets} for time {t_idx+1}/{n_times}')
            topocentric_t = rubin.observe(star)
            appr = topocentric_t.apparent()
            # print(appr)
            alt_t, az_t, height = appr.altaz(temperature_C=12.0, pressure_mbar=750)
            # print(f'Apparent {appr}: Altitude={alt_t.degrees}, Azimuth={az_t.degrees}, Height={height.to(u.m).value}')
            cube_a[t_idx, s_idx,:] = [alt_t.degrees, az_t.degrees, height.to(u.m).value, 0 if alt_t.degrees < 30 else 1]
            cube_b[t_idx, s_idx, :] = appr

    # Save the cube to a FITS file
    hdu_a = fits.PrimaryHDU(cube_a)
    hdu_b = fits.PrimaryHDU(cube_b)
    hdul = fits.HDUList([hdu_a, hdu_b])
    hdu_b.header["EXTNAME"] = "APPARENT_POSITIONS"
    hdu_b.header["FIELDS"] = "apparent"
    hdu_a.header["FIELDS"] = ",".join(fields)
    hdu_a.header["NTARGETS"] = n_targets
    hdu_a.header["TIMES"] = n_times
    hdu_a.header["TARGETS"] = ",".join(names)
    hdu_a.writeto(fits_file, overwrite=True)

    print(f'FITS file saved to {fits_file}')

    return fits_file

def manual_angular_separation(ra1, dec1, ra2, dec2):
    """
    Calculate the angular separation between two points on the celestial sphere.

    Parameters:
        ra1, dec1: Coordinates of the first point in degrees.
        ra2, dec2: Coordinates of the second point in degrees.

    Returns:
        Angular separation in arcseconds.
    """
    ra1_rad = np.radians(ra1)
    dec1_rad = np.radians(dec1)
    ra2_rad = np.radians(ra2)
    dec2_rad = np.radians(dec2)

    delta_ra = ra2_rad - ra1_rad
    delta_dec = dec2_rad - dec1_rad

    a = np.sin(delta_dec / 2)**2 + np.cos(dec1_rad) * np.cos(dec2_rad) * np.sin(delta_ra / 2)**2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))

    # Convert radians to arcseconds
    separation_arcsec = np.degrees(c) * 3600.0
    return separation_arcsec




# def calculate_target_motion(target, rubinobs_astr, ti, dt_seconds=1.0):
#     """
#     Calculate the angular motion of a target as seen from the observer
    
#     Parameters:
#         target: Skyfield Star object
#         rubinobs_astr: Observer location
#         ti: Skyfield time object (current time)
#         dt_seconds: Time step in seconds for motion calculation
    
#     Returns:
#         dict: Motion data including angular velocities and distances
#     """
    
#     # Calculate position at current time
#     pos1 = rubinobs_astr.at(ti).observe(target).apparent()
#     ra1, dec1, dist1 = pos1.radec()
#     alt1, az1, height1 = pos1.altaz()
    
#     # Calculate position at slightly later time
#     ts = ti.ts
#     ti_next = ts.tt_jd(ti.tt + dt_seconds / 86400.0)  # Add dt_seconds in days
#     pos2 = rubinobs_astr.at(ti_next).observe(target).apparent()
#     ra2, dec2, dist2 = pos2.radec()
#     alt2, az2, height2 = pos2.altaz()
    
#     # Calculate angular differences (in degrees)
#     delta_ra = (ra2.hours - ra1.hours) * 15.0  # Convert hours to degrees
#     delta_dec = dec2.degrees - dec1.degrees
#     delta_alt = alt2.degrees - alt1.degrees
#     delta_az = az2.degrees - az1.degrees
    
#     # Handle RA wraparound (crossing 0h/24h boundary)
#     if delta_ra > 180:
#         delta_ra -= 360
#     elif delta_ra < -180:
#         delta_ra += 360
        
#     # Handle azimuth wraparound (crossing 0°/360° boundary)
#     if delta_az > 180:
#         delta_az -= 360
#     elif delta_az < -180:
#         delta_az += 360
    
#     # Calculate angular velocities (arcseconds per second)
#     ra_velocity = (delta_ra * 3600.0) / dt_seconds  # arcsec/sec
#     dec_velocity = (delta_dec * 3600.0) / dt_seconds  # arcsec/sec
#     alt_velocity = (delta_alt * 3600.0) / dt_seconds  # arcsec/sec
#     az_velocity = (delta_az * 3600.0) / dt_seconds  # arcsec/sec
    
#     # Calculate total angular velocity in RA-Dec and Alt-Az
#     radec_angular_velocity = np.sqrt(ra_velocity**2 + dec_velocity**2)  # arcsec/sec
#     altaz_angular_velocity = np.sqrt(alt_velocity**2 + az_velocity**2)  # arcsec/sec
    
#     # Calculate physical distances if distance is available
#     distance_km = dist1.km
    
#     # Convert angular velocity to linear velocity at target distance
#     # Angular velocity in radians/sec * distance = linear velocity
#     radec_linear_velocity = np.radians(radec_angular_velocity / 3600.0) * distance_km  # km/sec
#     altaz_linear_velocity = np.radians(altaz_angular_velocity / 3600.0) * distance_km  # km/sec
    
#     return {
#         'time_mjd': ti.tt,
#         'time_iso': ti.utc_iso(),
#         'dt_seconds': dt_seconds,
        
#         # Current position
#         'ra_hours': ra1.hours,
#         'dec_degrees': dec1.degrees,
#         'altitude_degrees': alt1.degrees,
#         'azimuth_degrees': az1.degrees,
#         'distance_km': distance_km,
        
#         # Angular velocities (arcsec/sec)
#         'ra_velocity_arcsec_per_sec': ra_velocity,
#         'dec_velocity_arcsec_per_sec': dec_velocity,
#         'alt_velocity_arcsec_per_sec': alt_velocity,
#         'az_velocity_arcsec_per_sec': az_velocity,
        
#         # Total angular velocities
#         'radec_angular_velocity_arcsec_per_sec': radec_angular_velocity,
#         'altaz_angular_velocity_arcsec_per_sec': altaz_angular_velocity,
        
#         # Linear velocities (km/sec)
#         'radec_linear_velocity_km_per_sec': radec_linear_velocity,
#         'altaz_linear_velocity_km_per_sec': altaz_linear_velocity,
        
#         # Motion in different units
#         'radec_motion_arcsec_per_minute': radec_angular_velocity * 60,
#         'radec_motion_arcsec_per_hour': radec_angular_velocity * 3600,
#         'altaz_motion_arcsec_per_minute': altaz_angular_velocity * 60,
#         'altaz_motion_arcsec_per_hour': altaz_angular_velocity * 3600,
#     }

# def calculate_all_targets_motion(targets, rubinobs_astr, ti, dt_seconds=1.0):
#     """
#     Calculate motion for all targets at a given time
    
#     Parameters:
#         targets: List of Skyfield Star objects
#         rubinobs_astr: Observer location
#         ti: Skyfield time object
#         dt_seconds: Time step for motion calculation
    
#     Returns:
#         list: Motion data for all targets
#     """
    
#     all_motion_data = []
    
#     for i, target in enumerate(targets):
#         try:
#             motion_data = calculate_target_motion(target, rubinobs_astr, ti, dt_seconds)
#             motion_data['target_id'] = i
#             all_motion_data.append(motion_data)
            
#             if i % 1000 == 0:  # Progress indicator
#                 print(f"Calculated motion for {i+1}/{len(targets)} targets")
                
#         except Exception as e:
#             print(f"Error calculating motion for target {i}: {e}")
#             continue
    
#     return all_motion_data

# def save_motion_to_fits(motion_data_list, filename):
#     """
#     Save target motion data to a FITS file
    
#     Parameters:
#         motion_data_list: List of motion dictionaries from calculate_all_targets_motion
#         filename: Output FITS filename
#     """
    
#     if not motion_data_list:
#         print("No motion data to save")
#         return
    
#     # Extract data for FITS table
#     target_ids = [d['target_id'] for d in motion_data_list]
#     time_mjds = [d['time_mjd'] for d in motion_data_list]
#     time_isos = [d['time_iso'] for d in motion_data_list]
    
#     ra_hours = [d['ra_hours'] for d in motion_data_list]
#     dec_degrees = [d['dec_degrees'] for d in motion_data_list]
#     alt_degrees = [d['altitude_degrees'] for d in motion_data_list]
#     az_degrees = [d['azimuth_degrees'] for d in motion_data_list]
#     distances = [d['distance_km'] for d in motion_data_list]
    
#     ra_velocities = [d['ra_velocity_arcsec_per_sec'] for d in motion_data_list]
#     dec_velocities = [d['dec_velocity_arcsec_per_sec'] for d in motion_data_list]
#     alt_velocities = [d['alt_velocity_arcsec_per_sec'] for d in motion_data_list]
#     az_velocities = [d['az_velocity_arcsec_per_sec'] for d in motion_data_list]
    
#     radec_angular_vel = [d['radec_angular_velocity_arcsec_per_sec'] for d in motion_data_list]
#     altaz_angular_vel = [d['altaz_angular_velocity_arcsec_per_sec'] for d in motion_data_list]
    
#     radec_linear_vel = [d['radec_linear_velocity_km_per_sec'] for d in motion_data_list]
#     altaz_linear_vel = [d['altaz_linear_velocity_km_per_sec'] for d in motion_data_list]
    
#     # Create table
#     motion_table = Table([
#         target_ids, time_mjds, time_isos,
#         ra_hours, dec_degrees, alt_degrees, az_degrees, distances,
#         ra_velocities, dec_velocities, alt_velocities, az_velocities,
#         radec_angular_vel, altaz_angular_vel, radec_linear_vel, altaz_linear_vel
#     ], names=[
#         'TARGET_ID', 'TIME_MJD', 'TIME_ISO',
#         'RA_HOURS', 'DEC_DEGREES', 'ALT_DEGREES', 'AZ_DEGREES', 'DISTANCE_KM',
#         'RA_VEL_ARCSEC_SEC', 'DEC_VEL_ARCSEC_SEC', 'ALT_VEL_ARCSEC_SEC', 'AZ_VEL_ARCSEC_SEC',
#         'RADEC_ANGULAR_VEL_ARCSEC_SEC', 'ALTAZ_ANGULAR_VEL_ARCSEC_SEC',
#         'RADEC_LINEAR_VEL_KM_SEC', 'ALTAZ_LINEAR_VEL_KM_SEC'
#     ])
    
#     # Add metadata
#     motion_table.meta['DT_SECONDS'] = motion_data_list[0]['dt_seconds']
#     motion_table.meta['N_TARGETS'] = len(motion_data_list)
#     motion_table.meta['TIME_MJD'] = motion_data_list[0]['time_mjd']
#     motion_table.meta['PURPOSE'] = 'Target motion analysis'
    
#     # Save to FITS
#     hdu = fits.BinTableHDU(motion_table)
#     hdu.header['EXTNAME'] = 'TARGET_MOTION'
#     hdu.header['DT_SEC'] = motion_data_list[0]['dt_seconds']
#     hdu.header['UNITS_ANG'] = 'arcsec/sec'
#     hdu.header['UNITS_LIN'] = 'km/sec'
    
#     primary_hdu = fits.PrimaryHDU()
#     primary_hdu.header['PURPOSE'] = 'Target motion database'
#     primary_hdu.header['CREATED'] = dt.datetime.now().isoformat()
    
#     hdul = fits.HDUList([primary_hdu, hdu])
#     hdul.writeto(filename, overwrite=True)
    
#     print(f"Motion data saved to {filename}")
#     return filename

# Add compatibility functions to work with the tracker

def get_global_tracker():
    """
    Get or create global tracker instance
    """
    global _global_tracker
    if '_global_tracker' not in globals():
        from ..tracking import SatelliteTracker
        _global_tracker = SatelliteTracker()
    return _global_tracker

# Wrapper functions for backward compatibility
def trail_event_with_tracker(starlink, rubin_obs, ti, sat_t, sat_al, sat_az, sat_hght):
    """Wrapper that uses global tracker"""
    tracker = get_global_tracker()
    return tracker.trail_event(starlink, rubin_obs, ti, sat_t, sat_al, sat_az, sat_hght)

def point_trail_event_with_tracker(starlink, rubin_obs, ti):
    """Wrapper that uses global tracker"""
    tracker = get_global_tracker()
    return tracker.point_trail_event(starlink, rubin_obs, ti)

def get_current_satellite_position():
    """Get current satellite position from global tracker"""
    tracker = get_global_tracker()
    return tracker.get_current_position()
