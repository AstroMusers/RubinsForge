import os
import numpy as np
import pytz
import cv2
from astropy import units as u
from astropy.io import fits
from astropy.coordinates import SkyCoord
import skyfield as sf
from skyfield.api import load, wgs84, EarthSatellite, N, W,S,E, Star, Angle
import datetime as dt
import time
from pytz import timezone
from skyfield import almanac
from skyfield.framelib import galactic_frame
from skyfield.data import hipparcos
from brightness import calculate_brightness
import lumos.calculator
import random
import starlink_sat
import lumos.plot
import lumos.brdf.library
from astropy.wcs import WCS
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.ticker import FuncFormatter
import matplotlib.dates as mdates
from datetime import datetime, timedelta
from matplotlib import colors as mcolors
from astropy.visualization import quantity_support, time_support
quantity_support()  
matplotlib.rcParams['figure.figsize'] = (16,12)
matplotlib.rcParams['font.size'] = 17
def ra_formatter(x, pos):
        hours = int(x)
        minutes = int((x - hours) * 60)
        seconds = int(((x - hours) * 60 - minutes) * 60)
        return f'{hours:02d}:{minutes:02d}:{seconds:02d}'

c = dt.datetime.now()
if not os.path.exists("/data/a.saricaoglu/Lumos-Sat/Files/" +  str(c.strftime("%m.%d"))): 
    directoryf = "/data/a.saricaoglu/Lumos-Sat/Files/"  +  str(c.strftime("%m.%d"))+"/" + c.strftime('%H%M')
    os.makedirs(directoryf) 
if not os.path.exists("/data/a.saricaoglu/Lumos-Sat/Plots/" +  str(c.strftime("%m.%d")) + "/" + c.strftime('%H%M')): 
    directoryp = "/data/a.saricaoglu/Lumos-Sat/Plots/"  +  str(c.strftime("%m.%d"))+ "/" + c.strftime('%H%M')
    os.makedirs(directoryp) 
directoryf = "/data/a.saricaoglu/Lumos-Sat/Files/"  +  str(c.strftime("%m.%d"))
directoryp = "/data/a.saricaoglu/Lumos-Sat/Plots/"  +  str(c.strftime("%m.%d"))

ts = load.timescale()
zone = timezone('Chile/Continental')
eph = load('de421.bsp')
earth = eph['earth']
sun = eph['sun']
rubin_obs = wgs84.latlon(-30.244633,  -70.749417)

targets_dir = '/data/a.saricaoglu/lumos-sat/master.fits'


obs_start = '2025-01-24 00:00:00-03:00'    
 # Parse the datetime string
local_dt = datetime.fromisoformat(obs_start)
# Convert the local datetime to UTC
utc_dt = local_dt.astimezone(pytz.utc)
# Create a ts.utc object
obs_start = ts.utc(utc_dt.year, utc_dt.month, utc_dt.day, utc_dt.hour, utc_dt.minute, utc_dt.second)

deltatime = timedelta(days=60)
print(obs_start)   
t00 = obs_start
targets = []
with fits.open(targets_dir) as hdul:
    data = hdul[1].data
    ra_i = data['RAJ2000_deg']
    dec_i = data['DEJ2000_deg']
    # for i in range(0,len(data)):
    #     target_i = Star(Angle(degrees=ra_i[i]),Angle(degrees=dec_i[i]))
    #     targets.append(target_i)
    for i in range(0,len(data)):
        target_i = sf.positionlib.position_of_radec(ra_i[i], dec_i[i], t=obs_start, center=399)
        targets.append(wgs84.subpoint(target_i))


# Create a plot for for targets and Starlink positions for each day in Ra, Dec for targets above 30° altitude
for i in range(0,8):

    fig, ax = plt.subplots()
    above_30 = 0
    below_30 = 0
    for target in targets:
        # Calculate the target's topocentric position relative to the observing location
        difference_t = target - rubin_obs
        topocentric_t = difference_t.at(t00)     
        # Convert to Alt-Az
        altitude, azimuth, distance = topocentric_t.altaz()

        # Filter targets above 30° altitude
        if altitude.degrees > 30:
            # Convert RA and Dec for plotting
            above_30 = above_30 + 1
            ra_t, dec_t, distance_t = topocentric_t.radec()
            
            # Plot the target's RA and Dec
            ax.scatter(ra_t.hours, dec_t.degrees, marker="*", c='b', label=f'Target above 30: {target}')
        else:
            below_30 = below_30 + 1
            ra_t, dec_t, distance_t = topocentric_t.radec()
            ax.scatter(ra_t.hours, dec_t.degrees, marker="*", c='r', label=f'Target below 30: {target}')
    print(f'above 30: {above_30}')
    print(f'below 30: {below_30}')
    # sat_ra_above = [ra in sat_ra_all for ra, alt in zip(sat_ra_all, sat_alt_all) if alt > 30]
    # sat_dec_above = [dec in sat_dec_all for dec, alt in zip(sat_dec_all, sat_alt_all) if alt > 30]
    # sat_height_above = [height in sat_height_all for height, alt in zip(sat_height_all, sat_alt_all) if alt.degrees > 30]

    plt.title(f'target positions for {t00.astimezone(zone)}')
    plt.grid()
    ax.xaxis.set_major_formatter(FuncFormatter(ra_formatter))
    ax.set_xlabel('RA (hh:mm:ss)')

    # Dec is already in degrees, so just label it
    ax.set_ylabel('Dec (degrees)')    
    plt.tight_layout()
    plt.savefig(directoryp +  "/" + c.strftime('%H%M') + "/Target_positions_for_" +str(t00.astimezone(zone)) + str(i) + "_radec.png")
    plt.close()
    t00 = t00 + deltatime

t00 = obs_start
for i in range(0,8):
    fig, ax = plt.subplots()
    above_30 = 0
    below_30 = 0
    for target in targets:
        # Calculate the target's topocentric position relative to the observing location
        difference_t = target - rubin_obs
        topocentric_t = difference_t.at(t00)     
        # Convert to Alt-Az
        altitude, azimuth, distance = topocentric_t.altaz()

        # Filter targets above 30° altitude
        if altitude.degrees > 30:
            # Convert RA and Dec for plotting
            above_30 = above_30 + 1
            dec_t, ra_t, distance_t = topocentric_t.altaz()
            
            # Plot the target's RA and Dec
            ax.scatter(ra_t.degrees, dec_t.degrees, marker="*", c='b', label=f'Target above 30: {target}')
        else:
            below_30 = below_30 + 1
            dec_t, ra_t, distance_t = topocentric_t.altaz()
            ax.scatter(ra_t.degrees, dec_t.degrees, marker="*", c='r', label=f'Target below 30: {target}')
    print(f'above 30: {above_30}')
    print(f'below 30: {below_30}')
    # sat_ra_above = [ra in sat_ra_all for ra, alt in zip(sat_ra_all, sat_alt_all) if alt > 30]
    # sat_dec_above = [dec in sat_dec_all for dec, alt in zip(sat_dec_all, sat_alt_all) if alt > 30]
    # sat_height_above = [height in sat_height_all for height, alt in zip(sat_height_all, sat_alt_all) if alt.degrees > 30]

    plt.title(f'target positions for {t00.astimezone(zone)}')
    plt.grid()
    ax.xaxis.set_major_formatter(FuncFormatter(ra_formatter))
    ax.set_xlabel('Az (degrees)')

    # Dec is already in degrees, so just label it
    ax.set_ylabel('Alt (degrees)')    
    plt.tight_layout()
    plt.savefig(directoryp +  "/" + c.strftime('%H%M') + "/Target_positions_for_" +str(t00.astimezone(zone)) + str(i) + "_altaz.png")
    plt.close()
    t00 = t00 + deltatime

# Create a plot for for targets and Starlink positions for each day in Ra, Dec for targets above 30° altitude
t00 = obs_start
timestep = timedelta(minutes=0.1)
for i in range(0,8):
    separation = []
    fig, ax = plt.subplots()
    above_30 = 0
    below_30 = 0

    for target in targets:
        # Calculate the target's topocentric position relative to the observing location
        difference_t = target - rubin_obs
        topocentric_t = difference_t.at(t00)     

        # Convert to Alt-Az
        altitude, azimuth, distance = topocentric_t.altaz()

        # Convert to Alt-A

        # Filter targets above 30° altitude
        if (altitude.degrees > 30):
            # Convert RA and Dec for plotting
            above_30 = above_30 + 1

            ra_t, dec_t, distance_t = topocentric_t.radec()
            
            topocentric_pre = difference_t.at(t00 + timestep)
            ra_tpre, dec_tpre, distance_tpre = topocentric_pre.radec()

            separation.append(topocentric_t.separation_from(topocentric_pre).arcseconds())

            # Plot the target's RA and Dec
            ax.scatter(ra_t.hours, dec_t.degrees, marker=".", c='b', label=f'Now: {target}')
            ax.scatter(ra_tpre.hours, dec_tpre.degrees, marker=".", c='r', label=f'Previous: {target}')

    print(f'above 30: {above_30}')
    print(f'below 30: {below_30}')
    # sat_ra_above = [ra in sat_ra_all for ra, alt in zip(sat_ra_all, sat_alt_all) if alt > 30]
    # sat_dec_above = [dec in sat_dec_all for dec, alt in zip(sat_dec_all, sat_alt_all) if alt > 30]
    # sat_height_above = [height in sat_height_all for height, alt in zip(sat_height_all, sat_alt_all) if alt.degrees > 30]

    plt.title(f'target positions for {t00.astimezone(zone)}')
    plt.grid()
    ax.xaxis.set_major_formatter(FuncFormatter(ra_formatter))
    ax.set_xlabel('RA (hh:mm:ss)')

    # Dec is already in degrees, so just label it
    ax.set_ylabel('Dec (degrees)')    
    plt.tight_layout()
    plt.savefig(directoryp +  "/" + c.strftime('%H%M') + "/Target_positions_for_" +str(t00.astimezone(zone)) + '_' + str(timestep) +  "_iradec.png")
    plt.close()
    plt.figure(figsize=(14, 6))


    # Plot the bars for a1 and a2
    plt.hist(separation,bins=30, alpha=1, color='red', edgecolor='black')
    plt.axvline(x=15,  linestyle='--')

    # plt.hist(a2, bins=10, label='Average Magnitudes', color='salmon', edgecolor='black')

    # Label the axes
    plt.xlabel('Distance [arcseconds]')
    plt.ylabel('Number of Trail Points')
    #plt.yscale('log')

    plt.legend()

    plt.savefig(directoryp +  "/" + c.strftime('%H%M') + "/Separation_target_histogram_for_" + str(t00.astimezone(zone))+ '_' +  str(timestep) + ".png")

    timestep = timedelta(minutes=0.1*(i+1))

        
