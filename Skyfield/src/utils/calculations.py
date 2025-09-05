
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


def target_check(L,target_apr, target_altaz, topocentric, separation, treshold):
    start = time.time()
    max_duration = 3
    status = False
    while (((time.time() - start) % 60) < max_duration) & (status == False):

        difference_angle = topocentric.separation_from(target_apr)
        separation.append(difference_angle.arcseconds())

        if (difference_angle.arcseconds() <= treshold):
            print(f'Contaminated targets alt/az : {target_apr.alt.degrees} , {target_apr.az.degrees}')
            L.append(f'Contaminated targets alt/az : {target_apr.alt.degrees} , {target_apr.az.degrees}')

            print(f'Contamination difference angle  : {difference_angle.arcseconds()}')
            L.append(f'Contamination difference angle  : {difference_angle.arcseconds()}')
            print(f'target check done{((time.time() - start))} seconds')
            status = True
            return True
        else:
            status = True
            return False

def calculated_target_check(L,visible_targets_appr,topocentric, separation, treshold):
    start = time.time()
    max_duration = 3
    status = False

    for target_apr in visible_targets_appr:

        while (((time.time() - start) % 60) < max_duration) & (status == False):

            difference_angle = topocentric.separation_from(target_apr)
            separation.append(difference_angle.arcseconds())

            if (difference_angle.arcseconds() <= treshold):
                print(f'Contaminated targets alt/az : {target_apr.alt.degrees} , {target_apr.az.degrees}')
                L.append(f'Contaminated targets alt/az : {target_apr.alt.degrees} , {target_apr.az.degrees}')

                print(f'Contamination difference angle  : {difference_angle.arcseconds()}')
                L.append(f'Contamination difference angle  : {difference_angle.arcseconds()}')
                print(f'target check done{((time.time() - start))} seconds')
                status = True
                return True
            else:
                status = True
                return False

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
