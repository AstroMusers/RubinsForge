import os
import numpy as np
import pytz
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
sat_ra = []
sat_dec = [] 
sat_dist = []
sat_time = []
sat_az = []
sat_alt = []
sat_height = []

read_dir = '/data/a.saricaoglu/Lumos-Sat/Files/01.27/1947'

with fits.open(read_dir + '/1947_satellite_data.fits') as hdul:
    hdul.info()
    days = hdul[0].header['DAYS']
    obs_start = hdul[0].header['OBSSTART']
    obs_end = hdul[0].header['OBSEND'] 
    Nsats = hdul[0].header['NSATS']
    targets_dir = hdul[0].header['TARGETS']
    print(days)
    print(obs_start)
    print(obs_end)
    print(Nsats)
    print(targets_dir)
    for i in range(1,days+1):
        print(i)
        data = hdul[i].data
        sat_ra.append(data['RA'])
        sat_dec.append(data['DEC'])
        sat_dist.append(data['Distance'])
        sat_time.append(data['Time'])
        sat_az.append(data['Azimuth'])
        sat_alt.append(data['Altitude'])
        sat_height.append(data['Height'])

sat_ra_all = [ra for sat in sat_ra for ra in sat]
sat_dec_all = [dec for sat in sat_dec for dec in sat]
sat_dist_all = [dist for sat in sat_dist for dist in sat]
sat_time_all = [time for sat in sat_time for time in sat]
sat_az_all = [az for sat in sat_az for az in sat]   
sat_alt_all = [alt for sat in sat_alt for alt in sat]
sat_height_all = [height for sat in sat_height for height in sat]

print(len(sat_ra))
print(len(sat_ra_all))


# Parse the datetime string
local_dt = datetime.fromisoformat(obs_end)
# Convert the local datetime to UTC
utc_dt = local_dt.astimezone(pytz.utc)
# Create a ts.utc object
obs_end = ts.utc(utc_dt.year, utc_dt.month, utc_dt.day, utc_dt.hour, utc_dt.minute, utc_dt.second)
print(obs_end)
# Parse the datetime string
local_dt = datetime.fromisoformat(obs_start)
# Convert the local datetime to UTC
utc_dt = local_dt.astimezone(pytz.utc)
# Create a ts.utc object
obs_start = ts.utc(utc_dt.year, utc_dt.month, utc_dt.day, utc_dt.hour, utc_dt.minute, utc_dt.second)

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

t00 = obs_start

# Create a plot for for targets and Starlink positions for each day in Ra, Dec
for i in range(0,days):
    fig, ax = plt.subplots()
    # for target in targets:
    #     # Calculate the target's topocentric position relative to the observing location
    #     difference_t = target - rubin_obs
    #     topocentric_t = difference_t.at(t00)     
    #     # Convert to Alt-Az
    #     ra_t, dec_t, distance_t = topocentric_t.radec()
    #     ax.scatter(ra_t.hours, dec_t.degrees, marker="*", c='b')

    scatter = ax.scatter(sat_ra[i], sat_dec[i], c = sat_dist[i], marker='.', cmap = 'viridis')  

    plt.title(f'Starlink positions for {t00.astimezone(zone)} in RA-Dec')
    plt.grid()

    # Format the x-axis to show RA in hh:mm:ss format
    #ax.xaxis.set_major_formatter(FuncFormatter(ra_formatter))
    ax.set_xlabel('RA (hh:mm:ss)')
    fig.colorbar(scatter).set_label("Distance [km]")

    # Dec is already in degrees, so just label it
    ax.set_ylabel('Dec (degrees)')    
    plt.tight_layout()
    plt.savefig(directoryp +  "/" + c.strftime('%H%M') + "/Satellite_positions_for_" +str(t00.astimezone(zone)) + str(i) + "_A.png")
    plt.close()
    t00 = t00 + timedelta(days=1)
t00 = obs_start

# Create a plot for for targets and Starlink positions for each day in Alt, Az
for i in range(0,days):
    fig, ax = plt.subplots()
    for target in targets:
        # Calculate the target's topocentric position relative to the observing location
        difference_t = target - rubin_obs
        topocentric_t = difference_t.at(t00)     
        # Convert to Alt-Az
        alt_t, az_t, distance_t = topocentric_t.altaz()
        ax.scatter(az_t.degrees, alt_t.degrees, marker="*", c='b')

    scatter = ax.scatter(sat_az[i], sat_alt[i], c = sat_height[i], marker='.', cmap = 'viridis')  

    plt.title(f'Starlink positions for {t00.astimezone(zone)} in Alt-Az')
    plt.grid()

    ax.set_xlabel('Azimuth (degrees)')
    fig.colorbar(scatter).set_label("Distance [km]")

    # Dec is already in degrees, so just label it
    ax.set_ylabel('Altitude (degrees)')    
    plt.tight_layout()
    plt.savefig(directoryp +  "/" + c.strftime('%H%M') + "/Satellite_positions_for_" +str(t00.astimezone(zone)) + str(i) + "_B.png")
    plt.close()
    t00 = t00 + timedelta(days=1)
t00 = obs_start

# Create a plot for for targets and Starlink positions for each day in Ra, Dec for targets above 30° altitude
for i in range(0,days):
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
    scatter = ax.scatter(sat_ra[i], sat_dec[i], marker='o', c = sat_dist[i], cmap = 'viridis')  

    plt.title(f'Starlink positions for {t00.astimezone(zone)} in RA-Dec')
    plt.grid()
    ax.xaxis.set_major_formatter(FuncFormatter(ra_formatter))
    ax.set_xlabel('RA (hh:mm:ss)')
    fig.colorbar(scatter).set_label("Distance [km]")

    # Dec is already in degrees, so just label it
    ax.set_ylabel('Dec (degrees)')    
    plt.tight_layout()
    plt.savefig(directoryp +  "/" + c.strftime('%H%M') + "/Satellite_positions_for_" +str(t00.astimezone(zone)) + str(i) + "_C.png")
    plt.close()
    t00 = t00 + timedelta(days=1)
t00 = obs_start

# Create a plot for for targets and Starlink positions for each day in Ra, Dec for only targets above 30° altitude
for i in range(0,days):
    fig, ax = plt.subplots()
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
    # sat_ra_above = [ra in sat_ra_all for ra, alt in zip(sat_ra_all, sat_alt_all) if alt > 30]
    # sat_dec_above = [dec in sat_dec_all for dec, alt in zip(sat_dec_all, sat_alt_all) if alt > 30]
    # sat_height_above = [height in sat_height_all for height, alt in zip(sat_height_all, sat_alt_all) if alt.degrees > 30]
    scatter = ax.scatter(sat_ra[i], sat_dec[i], marker='o', c = sat_dist[i], cmap = 'viridis')  

    plt.title(f'Starlink positions for {t00.astimezone(zone)} in RA-Dec')
    plt.grid()
    ax.xaxis.set_major_formatter(FuncFormatter(ra_formatter))
    ax.set_xlabel('RA (hh:mm:ss)')
    fig.colorbar(scatter).set_label("Distance [km]")

    # Dec is already in degrees, so just label it
    ax.set_ylabel('Dec (degrees)')    
    plt.tight_layout()
    plt.savefig(directoryp +  "/" + c.strftime('%H%M') + "/Satellite_positions_for_" +str(t00.astimezone(zone)) + str(i) + "_D.png")
    plt.close()
    t00 = t00 + timedelta(days=1)
t00 = obs_start

# Create a plot for for targets and Starlink positions for each day in Alt, Az for targets above 30° altitude
for i in range(0,days):
    fig, ax = plt.subplots()
    for target in targets:
        # Calculate the target's topocentric position relative to the observing location
        difference_t = target - rubin_obs
        topocentric_t = difference_t.at(t00)     
        # Convert to Alt-Az
        altitude, azimuth, distance = topocentric_t.altaz()

        # Filter targets above 30° altitude
        if altitude.degrees > 30:

            ax.scatter(azimuth.degrees, altitude.degrees, marker="*", c='b', label=f'Target above 30: {target}')
        else:
            ra_t, dec_t, distance_t = topocentric_t.radec()
            ax.scatter(azimuth.degrees, altitude.degrees, marker="*", c='r', label=f'Target below 30: {target}')
  
    # sat_alt_above = [alt in sat_ra_all for ra, alt in zip(sat_ra_all, sat_alt_all) if alt > 30]
    # sat_dec_above = [dec in sat_dec_all for dec, alt in zip(sat_dec_all, sat_alt_all) if alt > 30]
    # sat_height_above = [height in sat_height_all for height, alt in zip(sat_height_all, sat_alt_all) if alt > 30]
    scatter = ax.scatter(sat_az[i], sat_alt[i], marker='o', c = sat_height[i], cmap = 'viridis')  

    plt.title(f'Starlink positions for {t00.astimezone(zone)} in Alt-Az')
    plt.grid()

    ax.set_xlabel('Azimuth (degrees)')
    fig.colorbar(scatter).set_label("Distance [km]")

    # Dec is already in degrees, so just label it
    ax.set_ylabel('Altitude (degrees)')    
    plt.tight_layout()
    plt.savefig(directoryp +  "/" + c.strftime('%H%M') + "/Satellite_positions_for_" +str(t00.astimezone(zone)) +str(i) + "_E.png")
    plt.close()
    t00 = t00 + timedelta(days=1)
t00 = obs_start

# Create a plot for for targets and Starlink positions for each day in Ra, Dec for targets above 30° altitude
# for i in range(0,6):
#     fig, ax = plt.subplots()
#     above_30 = 0
#     below_30 = 0
#     for target in targets:
#         # Calculate the target's topocentric position relative to the observing location
#         difference_t = target - rubin_obs
#         topocentric_t = difference_t.at(t00)     
#         # Convert to Alt-Az
#         altitude, azimuth, distance = topocentric_t.altaz()

#         # Filter targets above 30° altitude
#         if altitude.degrees > 30:
#             # Convert RA and Dec for plotting
#             above_30 = above_30 + 1
#             ra_t, dec_t, distance_t = topocentric_t.radec()
            
#             # Plot the target's RA and Dec
#             ax.scatter(ra_t.hours, dec_t.degrees, marker="*", c='b', label=f'Target above 30: {target}')
#         else:
#             below_30 = below_30 + 1
#             ra_t, dec_t, distance_t = topocentric_t.radec()
#             ax.scatter(ra_t.hours, dec_t.degrees, marker="*", c='r', label=f'Target below 30: {target}')
#     print(f'above 30: {above_30}')
#     print(f'below 30: {below_30}')
#     # sat_ra_above = [ra in sat_ra_all for ra, alt in zip(sat_ra_all, sat_alt_all) if alt > 30]
#     # sat_dec_above = [dec in sat_dec_all for dec, alt in zip(sat_dec_all, sat_alt_all) if alt > 30]
#     # sat_height_above = [height in sat_height_all for height, alt in zip(sat_height_all, sat_alt_all) if alt.degrees > 30]

#     plt.title(f'Target positions for {t00.astimezone(zone)} in RA-Dec')
#     plt.grid()
#     ax.xaxis.set_major_formatter(FuncFormatter(ra_formatter))
#     ax.set_xlabel('RA (hh:mm:ss)')

#     # Dec is already in degrees, so just label it
#     ax.set_ylabel('Dec (degrees)')    
#     plt.tight_layout()
#     plt.savefig(directoryp +  "/" + c.strftime('%H%M') + "/Satellite_positions_for_" +str(t00.astimezone(zone)) + str(i) + "_F.png")
#     plt.close()
#     t00 = t00 + timedelta(hours=2)
t00 = obs_start
# Create an array of dates for each day in September (1 to 30)
dates = [i for i in range(1,days+1)]



with open(read_dir +"/totaltrail_logfile_for_" + str(utc_dt.day) + str(utc_dt.month) + "_" + str(utc_dt.year) + ".txt") as file:
    total_trail = [float(line.strip()) for line in file.readlines()]

with open(read_dir + "/totalvalid_logfile_for_" + str(utc_dt.day) + str(utc_dt.month) + "_" + str(utc_dt.year) + ".txt") as file:
    total_valid = [float(line.strip()) for line in file.readlines()]

with open(read_dir + "/totalevent_logfile_for_"+ str(utc_dt.day) + str(utc_dt.month) + "_" + str(utc_dt.year) + ".txt") as file:
    total_event = [float(line.strip()) for line in file.readlines()]

with open(read_dir +"/totalclosestapproach_logfile_for_"+ str(utc_dt.day) + str(utc_dt.month) + "_" + str(utc_dt.year) + ".txt") as file:
    total_closest_approach = [float(line.strip()) for line in file.readlines()]

with open(read_dir +"/totalcontamination_logfile_for_"+ str(utc_dt.day) + str(utc_dt.month) + "_" + str(utc_dt.year) + ".txt") as file:
    total_contamination= [float(line.strip()) for line in file.readlines()]   

with open(read_dir +"/totalpeakmag_logfile_for_"+ str(utc_dt.day) + str(utc_dt.month) + "_" + str(utc_dt.year) + ".txt") as file:
    total_peak_mag= [float(line.strip()) for line in file.readlines()]  

with open(read_dir +"/totalaveragemag_logfile_for_" + str(utc_dt.day) + str(utc_dt.month) + "_" + str(utc_dt.year) + ".txt") as file:
    total_ave_mag = [float(line.strip()) for line in file.readlines()]  

with open(read_dir +"/totalpeakintensity_logfile_for_"+ str(utc_dt.day) + str(utc_dt.month) + "_" + str(utc_dt.year) + ".txt") as file:
    total_peak_intensity= [float(line.strip()) for line in file.readlines()]  

with open(read_dir +"/totalaverageintensity_logfile_for_" + str(utc_dt.day) + str(utc_dt.month) + "_" + str(utc_dt.year) + ".txt") as file:
    total_ave_intensity = [float(line.strip()) for line in file.readlines()]  

# Define variable 'a1' and 'a2' with random values for each day (for illustration)
a1 = total_trail  # Random data for demonstration
a2 = total_valid # Another set of random data for comparison
a3 = total_contamination
a4 = [int(ev / 3) for ev in total_event]
a5 = total_closest_approach

# Set positions for each set of bars
dates = np.array(dates) # Shift dates slightly right
# Create the plot
plt.figure(figsize=(20, 6))

# Plot the bars for a1 and a2
plt.bar(dates, a4, alpha=1.0, label='Events', color='green', edgecolor='black')
# plt.bar(dates, a2, alpha=1.0, label='Validated Starlinks', color='salmon', edgecolor='black')
plt.bar(dates, a1,  alpha=1.0,label='Trail Candidates', color='skyblue', edgecolor='black')
plt.bar(dates, a3, alpha=1.0, label='Strong Lens Crossing', color='red', edgecolor='black')

# Label the axes
plt.xlabel('Days')
plt.ylabel('Numbers')
plt.yscale('log')
# Title
plt.title(f'Starlink Events between {obs_start.astimezone(zone)} - {obs_end.astimezone(zone)} for {Nsats} Satellites per day' )
# Add a legend
plt.legend()
# Show the plot
plt.savefig(directoryp +  "/" + c.strftime('%H%M') + "/Event_histogram_for_"  + str(obs_start.astimezone(zone)) +  ".png")

# if len(dates) == 1:
#     def func(pct, allvalues):
#         absolute = int(pct / 100.*np.sum(allvalues))
#         return "{:.1f}%\n({:d} g)".format(pct, absolute)

#     data = [a1[0], a2[0], a4[0]]

#     size = 1
#     # normalizing data to 2 pi
#     norm = data / np.sum(data)*2 * np.pi

#     # obtaining ordinates of bar edges
#     left = np.cumsum(np.append(0,
#                             norm.flatten()[:-1])).reshape(data.shape)

#     # Creating color scale
#     cmap = plt.get_cmap("viridis")
#     outer_colors = cmap(np.arange(6)*4)
#     inner_colors = cmap(np.array([1, 2, 5, 6, 9,
#                                 10, 12, 13, 15,
#                                 17, 18, 20]))

#     # Creating plot
#     fig, ax = plt.subplots(figsize=(10, 7),
#                         subplot_kw=dict(polar=True))

#     ax.bar(x=left[:, 0],
#         width=norm.sum(axis=1),
#         bottom=1-size,
#         height=size,
#         color=outer_colors,
#         edgecolor='w',
#         linewidth=1,
#         align="edge")

#     ax.bar(x=left.flatten(),
#         width=norm.flatten(),
#         bottom=1-2 * size,
#         height=size,
#         color=inner_colors,
#         edgecolor='w',
#         linewidth=1,
#         align="edge")

#     ax.set_axis_off()

#     plt.savefig(directoryp +  "/" + c.strftime('%H%M') + "/Piechart_for_"  + str(obs_start.astimezone(zone)) +  ".png")    


# Width of bars
# Set positions for each set of bar
# Define variable 'a1' and 'a2' with random values for each day (for illustration)


aa1 = total_peak_mag  # Random data for demonstration
aa2 = total_ave_mag # Another set of random data for comparison
# Width of bars
# Set positions for each set of bars
bb1 = total_peak_intensity
bb2 = total_ave_intensity

plt.figure(figsize=(14, 6))

aa = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
# Plot the bars for a1 and a2
plt.hist(aa1,bins=aa, label='Peak Magnitudes', alpha=0.7, color='blue', edgecolor='black')
plt.hist(aa2,bins=aa, label='Average Magnitudes', alpha=0.5, color='red', edgecolor='black')
# Label the axes
plt.xlabel('Magnitude [AB]')
plt.ylabel('Number of Trail Candidate Events')

plt.legend()

plt.savefig(directoryp +  "/" + c.strftime('%H%M') + "/Brightness_histogram_for_" + str(obs_start.astimezone(zone)) +  ".png")

plt.figure(figsize=(14, 6))

bb = np.arange(0,10e-11,1e-11)

plt.hist(bb1,bins=bb, label='Peak Intensities', alpha=0.7, color='blue', edgecolor='black')
plt.hist(bb2,bins=bb, label='Average Intensities', alpha=0.5, color='red', edgecolor='black')

plt.xlabel('Intensity [W/m^2]')
plt.ylabel('Number of Trail Candidate Events')

plt.legend()

plt.savefig(directoryp +  "/" + c.strftime('%H%M') + "/Intensity_histogram_for_" + str(obs_start.astimezone(zone)) +  ".png")


plt.figure(figsize=(14, 6))


# Plot the bars for a1 and a2
cc = np.arange(0,300,10)
plt.hist(a5,bins=cc, alpha=1, color='red', edgecolor='black')
plt.xlim(left=0,right=500)
plt.axvline(x=15,  linestyle='--')

# plt.hist(a2, bins=10, label='Average Magnitudes', color='salmon', edgecolor='black')

# Label the axes
plt.xlabel('Distance [arcseconds]')
plt.ylabel('Number of Trail Points')
plt.yscale('log')
# Title
#plt.title('Separation from the closest target at any trail point')

# # Format the x-axis to show dates in mm/dd/yy format
# plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%m/%d/%y'))
# plt.gca().xaxis.set_major_locator(mdates.DayLocator(interval=2))  # Show every 2 days

# Rotate date labels for better readability
# plt.gcf().autofmt_xdate()

# Add a legend
plt.legend()

# Show the plot
plt.savefig(directoryp +  "/" + c.strftime('%H%M') + "/Separation_histogram_for_" + str(obs_start.astimezone(zone)) +  ".png")



