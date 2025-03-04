# Parent : targets_2.py
# Creates a mock target distribution with a population of 10k sources above 30 degrees.
import os
import sys
import logging
import numpy as np
from astropy import units as u
from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord
import skyfield as sf
from skyfield.positionlib import Apparent
from skyfield.api import load, wgs84, EarthSatellite, N, W,S,E, Star, Angle
import datetime as dt
import time
from pytz import timezone
from skyfield import almanac
from skyfield.framelib import galactic_frame
from skyfield.data import hipparcos
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
start_time = time.time()

def generate_dataset():
    year = 2025
    day_month_i = [20, 1]
    day_month_f = [30, 1]
    zone = timezone('Chile/Continental')
    local = timezone('Etc/GMT-6')
    obs_start = ts.utc(year,day_month_i[1], day_month_i[0],3,0,0)
    obs_end = ts.utc(year,day_month_f[1], day_month_f[0],3,0,0)
    eph = load('de421.bsp')
    earth = eph['earth']
    sun = eph['sun']
    rubin_obs = wgs84.latlon(-30.244633,  -70.749417)
    rubinobs_astr = earth +  wgs84.latlon(30.244633*S,  70.749417*W, elevation_m = 2647)
    num_sources = 10000  # Number of random sources you want to generate
    ra_min, ra_max = 0, 360         # RA ranges from 0 to 360 degrees
    dec_min, dec_max = -90,90       # Dec ranges from 30 to 90 degrees

    ra_values = []
    dec_values = []
    while len(ra_values) < num_sources:
    # Generate random RA and Dec values
        ra_i = np.random.uniform(ra_min, ra_max, 1)[0] 
        dec_i = np.random.uniform(dec_min, dec_max, 1)[0] 
        target_i = sf.positionlib.position_of_radec(ra_i, dec_i, t=obs_start, center=399)
        target_i = wgs84.subpoint(target_i)
        difference_t = target_i - rubin_obs
        topocentric_t = difference_t.at(obs_start)
        alt_t, az_t, height = topocentric_t.altaz()
        
        if alt_t.degrees > 30:
            ra_t, dec_t, distance_t = topocentric_t.radec()
            ra_values.append(ra_i)
            dec_values.append(dec_i)


    print(len(ra_values))
    # Create a table with the data
    data_table = Table([ra_values, dec_values], names=('RA', 'Dec'))

    # Define the FITS file and save the data
    fits_filename = directoryf + "/radec_mock_targets_list_" + str(num_sources) + ".fits"
    data_table.write(fits_filename, format='fits', overwrite=True)   

c = dt.datetime.now()
if not os.path.exists("/data/a.saricaoglu/Lumos-Sat/Files/" +  str(c.strftime("%m.%d"))): 
    directoryf = "/data/a.saricaoglu/Lumos-Sat/Files/"  +  str(c.strftime("%m.%d"))
    os.makedirs(directoryf) 
if not os.path.exists("/data/a.saricaoglu/Lumos-Sat/Plots/" +  str(c.strftime("%m.%d")) + "/" + c.strftime('%H%M')): 
    directoryp = "/data/a.saricaoglu/Lumos-Sat/Plots/"  +  str(c.strftime("%m.%d"))+ "/" + c.strftime('%H%M')
    os.makedirs(directoryp) 
directoryf = "/data/a.saricaoglu/Lumos-Sat/Files/"  +  str(c.strftime("%m.%d"))
directoryp = "/data/a.saricaoglu/Lumos-Sat/Plots/"  +  str(c.strftime("%m.%d"))

# Load necessary data
ts = load.timescale()
year = 2025
day_month_i = [20, 1]
day_month_f = [30, 1]
zone = timezone('Chile/Continental')
local = timezone('Etc/GMT-6')
obs_start = ts.utc(year,day_month_i[1], day_month_i[0],3,0,0)
obs_end = ts.utc(year,day_month_f[1], day_month_f[0],3,0,0)
eph = load('de421.bsp')
earth = eph['earth']
sun = eph['sun']
rubin_obs = wgs84.latlon(-30.244633,  -70.749417)
rubinobs_astr = earth +  wgs84.latlon(30.244633*S,  70.749417*W, elevation_m = 2647)
num_sources = 10000
fits_filename = "/data/a.saricaoglu/Lumos-Sat/Files/01.28/radec_mock_targets_list_" + str(num_sources) + ".fits"

# Load the FITS file data (update the path as per your file location)
with fits.open(fits_filename) as hdul:
    data = hdul[1].data
    hdul.info()
    ra_values = data['RA']
    dec_values = data['DEC']
hdul.close()   
fig, ax = plt.subplots(figsize=(10, 8))

# Loop through each target and plot it
for i in range(0,num_sources):
    ra_i = ra_values[i]
    dec_i = dec_values[i]
    target_i = sf.positionlib.position_of_radec(ra_i, dec_i, t=obs_start, center=399)
    target_i = wgs84.subpoint(target_i)
    difference_t = target_i - rubin_obs
    topocentric_t = difference_t.at(obs_start)
    ra_t, dec_t, dist_t = topocentric_t.radec()
    # Plot the target's RA and Dec
    ax.scatter(ra_t.hours, dec_t.degrees, marker="*", c='b', label='Target')

# Add custom legend
#legend2 = plt.Line2D([0], [0], marker='*', color='w', markerfacecolor='b', markeredgecolor='k', markersize=10, label='Target')
#ax.legend(handles=[legend2], loc='lower center', bbox_to_anchor=(0.5, 1.1), ncol=2)

# Set axis labels
ax.set_xlabel('RA (degrees)')
ax.set_ylabel('Dec (degrees)')
ax.xaxis.set_major_formatter(FuncFormatter(ra_formatter))
ax.set_xlabel('RA (hh:mm:ss)')

# Dec is already in degrees, so just label it
ax.set_ylabel('Dec (degrees)')

# Plot title and adjustments
plt.title('Mock Strong Lens Candidate Distribution for ' + str(num_sources) + ' Targets (RA-Dec)')
plt.tight_layout()
plt.savefig(directoryp +  "/" + c.strftime('%H%M') + "/Random_target_positions_radec.png")


fig, ax = plt.subplots(figsize=(10, 8))

# Loop through each target and plot it
for i in range(0,num_sources):
    ra_i = ra_values[i]
    dec_i = dec_values[i]
    target_i = sf.positionlib.position_of_radec(ra_i, dec_i, t=obs_start, center=399)
    target_i = wgs84.subpoint(target_i)
    difference_t = target_i - rubin_obs
    topocentric_t = difference_t.at(obs_start)
    alt_t, az_t, dist_t = topocentric_t.altaz()
    # Plot the target's RA and Dec
    ax.scatter(az_t.degrees, alt_t.degrees, marker="*", c='b', label='Target')

# Add custom legend
#legend2 = plt.Line2D([0], [0], marker='*', color='w', markerfacecolor='b', markeredgecolor='k', markersize=10, label='Target')
#ax.legend(handles=[legend2], loc='lower center', bbox_to_anchor=(0.5, 1.1), ncol=2)


ax.set_xlabel('Az (degrees)')
ax.set_ylabel('Alt (degrees)')

ax.set_xlabel('Az (degrees)')

# Dec is already in degrees, so just label it
ax.set_ylabel('Alt (degrees)')

# Plot title and adjustments
plt.title('Mock Strong Lens Candidate Distribution for ' + str(num_sources) + ' Targets (Alt-Az)')
plt.tight_layout()
plt.savefig(directoryp +  "/" + c.strftime('%H%M') + "/Random_target_positions_altaz.png")
end_time = time.time()
total_seconds = end_time - start_time

# Convert seconds to hours, minutes, and seconds
hours = int(total_seconds // 3600)
minutes = int((total_seconds % 3600) // 60)
seconds = total_seconds % 60
current_time = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime(time.time()))
# Print the formatted runtime
print(f"Current time: {current_time} Runtime: {hours}h {minutes}m {seconds:.2f}s")

# # Parameters for the RA and Dec distribution
# num_sources = 10000  # Number of random sources you want to generate
# az_min, az_max = 0, 360         # RA ranges from 0 to 360 degrees
# alt_min, alt_max = 30,90       # Dec ranges from 30 to 90 degrees

# # Generate random RA and Dec values
# az_values = np.random.uniform(az_min, az_max, num_sources)  # Random RA values
# alt_values = np.random.uniform(alt_min, alt_max, num_sources)  # Random Dec values

# # Create a table with the data
# data_table = Table([az_values, alt_values], names=('Az', 'Alt'))

# # Define the FITS file and save the data
# fits_filename = directoryf + "/altaz_mock_targets_list_" + str(num_sources) + ".fits"
# data_table.write(fits_filename, format='fits', overwrite=True)

# print(f"FITS file {fits_filename} created with {num_sources} random Alt and Az values.")


# # Load the FITS file data (update the path as per your file location)
# with fits.open(fits_filename) as hdul:
#     data = hdul[1].data
#     hdul.info()
#     az_values = data['Az']
#     alt_values = data['Alt']
# hdul.close()   

# fig, ax = plt.subplots(figsize=(8, 8))
# # Loop through each target and plot it
# for i in range(0,num_sources):
#     az_i = az_values[i]
#     alt_i = alt_values[i]
    
#     # Observe the target to get accurate coordinates
#     t = obs_start
#     astrometric = rubin_obs.at(t).from_altaz(alt_degrees = alt_i, az_degrees = az_i)
#     alt, az , distance = astrometric.radec()
#     print(alt, az)
#     # Plot the target's RA and Dec
#     ax.scatter(alt.hours, az.degrees, marker="*", c='b', label='Target')

# # Add custom legend
# legend2 = plt.Line2D([0], [0], marker='*', color='w', markerfacecolor='b', markeredgecolor='k', markersize=10, label='Target')
# ax.legend(handles=[legend2], loc='lower center', bbox_to_anchor=(0.5, 1.1), ncol=2)

# # Set axis labels
# ax.set_xlabel('Az (degrees)')
# ax.set_ylabel('Alt (degrees)')

# ax.xaxis.set_major_formatter(FuncFormatter(ra_formatter))
# ax.set_xlabel('Az (degrees)')

# # Dec is already in degrees, so just label it
# ax.set_ylabel('Alt (degrees)')

# # Plot title and adjustments
# plt.title('Mock Strong Lens Candidate Distribution for ' + str(num_sources) + ' Targets (Alt-Az)')
# plt.tight_layout()
# plt.savefig(directoryp +  "/" + c.strftime('%H%M') + "/Random_target_positions_altaz.png")
# end_time = time.time()
# total_seconds = end_time - start_time

# # Convert seconds to hours, minutes, and seconds
# hours = int(total_seconds // 3600)
# minutes = int((total_seconds % 3600) // 60)
# seconds = total_seconds % 60
# current_time = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime(time.time()))
# # Print the formatted runtime
# print(f"Current time: {current_time} Runtime: {hours}h {minutes}m {seconds:.2f}s")