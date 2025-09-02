
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
        sat_heights = sat_hght[s:f]
        intensities = lumos.calculator.get_intensity_observer_frame(
                starlink_sat.SURFACES_INFER_BRDFS, sat_heights, sat_altitudes, sat_azimuths,
                sun_altitude = sun_al, sun_azimuth = sun_az,
                include_earthshine = False)
            
        # Convert intensity to AB Magnitude
        ab_magnitudes = lumos.conversions.intensity_to_ab_mag(intensities)
        #print(f'task brightness done in {((time.time() - start) % 60)} seconds')
        status = True
        return intensities, ab_magnitudes
    else:
        return None
    
def trail_event(starlink, rubin_obs, ti, midnight, zone, sat_ra, sat_dec, sat_dist, sat_t, sat_al, sat_az, sat_hght):
    start = time.time()
    max_duration = 60
    status = False
    while (((time.time() - start) % 60) < max_duration) & (status == False):
        difference = starlink - rubin_obs
        topocentric = difference.at(ti)
        satra, satdec, distance = topocentric.radec() # ICRF ("J2000")
        if (distance.to(u.km).value < 2000):
            sat_ra.append(satra.hours)
            sat_dec.append(satdec.degrees)
            sat_dist.append(distance.to(u.m).value)
            satal, sataz, sat_height = topocentric.altaz()
            sat_al.append(satal.degrees)
            sat_az.append(sataz.degrees)
            sat_hght.append(sat_height.to(u.m).value)

            diff = ti.astimezone(zone) - midnight.astimezone(zone)
            sat_t.append(diff.total_seconds() / 3600)
            status = True           
            return topocentric         
        else:
            status = True
            return None
        
def target_check(target, rubin_obs, ti, topocentric, separation, treshold):
    start = time.time()
    max_duration = 60
    status = False
    while (((time.time() - start) % 60) < max_duration) & (status == False):
        difference = target - rubin_obs
        topocentric_t = difference.at(ti)
        difference_angle = topocentric.separation_from(topocentric_t)
        separation.append(difference_angle.arcseconds())

        if (difference_angle.arcseconds() <= treshold):
            ra_t, dec_t, distance_t = topocentric_t.radec()
            print(f'Contaminated targets ra/dec : {ra_t} , {dec_t}')
            L.append(f'Contaminated targets ra/dec : {ra_t} , {dec_t}')

            print(f'Contamination difference angle  : {difference_angle.arcseconds()}')
            L.append(f'Contamination difference angle  : {difference_angle.arcseconds()}')
            #print(f'target check done{((time.time() - start) % 60)} seconds')
            status = True
            return True
        else:
            status = True
            return False



print(os.getcwd())