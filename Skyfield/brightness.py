import lumos.calculator
import numpy as np
import lumos.plot
import lumos.brdf.library
import numpy as np
import starlink_sat
from astropy import units as u
from astropy.coordinates import SkyCoord
import skyfield as sf
from skyfield.api import load, wgs84, EarthSatellite, N, W, Star
import datetime as dt
from pytz import timezone
from skyfield import almanac
from skyfield.api import N,S,E,W, wgs84

eph = load('de421.bsp')
earth = eph['earth']
sun = eph['sun']
def calculate_brightness(satellite, observatory,ti):
    # Altitude and azimuth in the sky of a
    # specific geographic location
    difference = satellite - observatory
    topocentric = difference.at(ti)
    rubin = earth + observatory
    astro = rubin.at(ti).observe(sun)
    app = astro.apparent()

    sunalt, sunaz, sundistance = app.altaz()
    print(sunalt.dstr())
    print(sunaz.dstr())
    print(sundistance)

    satal, sataz, sat_height = topocentric.altaz()
    sat_height =  sat_height.to(u.m).value
    print(' sat height ', sat_height)

    sat_altitudes, sat_azimuths = \
        np.meshgrid(
        np.linspace(0, 90, 90), 
        np.linspace(0, 360, 180))


    intensities = lumos.calculator.get_intensity_observer_frame(
            starlink_sat.SURFACES, sat_height, sat_altitudes, sat_azimuths,
            sun_altitude = sunalt.degrees, sun_azimuth = sunaz.degrees,
            include_earthshine = False
        )
        
    # Convert intensity to AB Magnitude
    ab_magnitudes = lumos.conversions.intensity_to_ab_mag(intensities)

    peak_ab_mag = ab_magnitudes.min()
    idx = np.argmin(ab_magnitudes)
    peak_alt = sat_altitudes.flatten()[idx]
    peak_az = sat_azimuths.flatten()[idx]

    print( f"Brightest: {peak_ab_mag:0.1f} AB Magnitude")
    print( f"Altitude: {peak_alt:0.0f}°")
    print( f"Azimuth: {peak_az:0.0f}°")
    return peak_ab_mag