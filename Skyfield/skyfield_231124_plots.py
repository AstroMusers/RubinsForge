# %%
# Parent version: skyfield_satellite_positions_B_021124.py version
# Runs with daily increments, samples every 0.1m of the satellites and checks masterfıts in the sky.
# Plots altitude, azimuth and distance. ONLY positions no target contamination check.
import os
import numpy as np
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

ts = load.timescale()
year = 2025
day_month_i = [20, 12]
day_month_f = [24, 12]
obs_start = ts.utc(year,day_month_i[1], day_month_i[0])
obs_end = ts.utc(year,day_month_f[1], day_month_f[0])
eph = load('de421.bsp')
earth = eph['earth']
sun = eph['sun']
rubin_obs = wgs84.latlon(-30.244633,  -70.749417)
# rubin_obs_astr = earth +  wgs84.latlon(30.244633*S,  70.749417*W, elevation_m = 2647)

# %%
c = dt.datetime.now()
if not os.path.exists("/data/a.saricaoglu/Lumos-Sat/Files/" +  str(c.strftime("%m.%d"))): 
    directoryf = "/data/a.saricaoglu/Lumos-Sat/Files/"  +  str(c.strftime("%m.%d"))+"/" + c.strftime('%H%M')
    os.makedirs(directoryf) 
if not os.path.exists("/data/a.saricaoglu/Lumos-Sat/Plots/" +  str(c.strftime("%m.%d")) + "/" + c.strftime('%H%M')): 
    directoryp = "/data/a.saricaoglu/Lumos-Sat/Plots/"  +  str(c.strftime("%m.%d"))+ "/" + c.strftime('%H%M')
    os.makedirs(directoryp) 
directoryf = "/data/a.saricaoglu/Lumos-Sat/Files/"  +  str(c.strftime("%m.%d"))
directoryp = "/data/a.saricaoglu/Lumos-Sat/Plots/"  +  str(c.strftime("%m.%d"))

fit_filename = "/data/a.saricaoglu/Lumos-Sat/Files/01.22/1926_satellite_data.fits"

targets = []
with fits.open('/data/a.saricaoglu/lumos-sat/master.fits') as hdul:
    hdul.info()
    data = hdul[1].data
    ra_i = data['RAJ2000_deg']
    dec_i = data['DEJ2000_deg']
    print(ra_i, dec_i)
    # for i in range(0,len(data)):
    #     target_i = Star(Angle(degrees=ra_i[i]),Angle(degrees=dec_i[i]))
    #     targets.append(target_i)
    fig, ax = plt.subplots(figsize=(8, 8))
    for i in range(0,len(data)):
        target_i = sf.positionlib.position_of_radec(ra_i[i], dec_i[i])
        targets.append(target_i)

start_time = time.time()

Nsats = 64
if Nsats != len(starlinks):
    random.shuffle(starlinks)
status_names = ['Rise','Culminate', 'Set']
t00 = obs_start
d = 0
eclipse = 0
rubinobs = wgs84.latlon(-30.244633,  -70.749417)
rubinobs_astr = earth +  wgs84.latlon(30.244633*S,  70.749417*W, elevation_m = 2647)

t00 = obs_start
while t00 < obs_end:
    d = d + 1
    fig, ax = plt.subplots()
    for target in targets:
        # astrometric = earth.at(ti).observe(target)
        ra_t, dec_t, distance_t = target.radec()
        # target_i = Star(Angle(degrees=ra_i), Angle(degrees=dec_i))
        # Plot the target's RA and Dec
        ax.scatter(ra_t.hours, dec_t.degrees, marker="*", c='b', label='Target')

    with fits.open(fit_filename) as hdul:
        hdul.info()

        for i in range(1,len(hdul)):
            data = hdul[i].data
            ra_i = data['RA']
            dec_i = data['Dec']
            dist_i = data['Distance']
            time_i = data['Time']


    
    scatter = ax.scatter(sat_ra, sat_dec, marker='o', c = sat_dist, cmap = 'viridis')  

    plt.title(f'Starlink positions for {t00.astimezone(zone)}')
    plt.grid()
    ax.xaxis.set_major_formatter(FuncFormatter(ra_formatter))
    ax.set_xlabel('RA (hh:mm:ss)')
    fig.colorbar(scatter).set_label("Distance [km]")

    # Dec is already in degrees, so just label it
    ax.set_ylabel('Dec (degrees)')    
    plt.tight_layout()
    plt.savefig(directoryp +  "/" + c.strftime('%H%M') + "/Satellite_positions_for_" +str(t00.astimezone(zone)) +  "_B.png")
    plt.close()

    
    end_time = time.time()
    total_seconds = end_time - start_time
    # Convert seconds to hours, minutes, and seconds
    hours = int(total_seconds // 3600)
    minutes = int((total_seconds % 3600) // 60)
    seconds = total_seconds % 60
    current_time = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime(time.time()))
    # Print the formatted runtime
    print(f"Current time: {current_time} Runtime: {hours}h {minutes}m {seconds:.2f}s , Day {d} completed.")
    L.append(f"\n Current time: {current_time} Runtime: {hours}h {minutes}m {seconds:.2f}s , Day {d} completed for {valid_starlink_counter} valid events")
    L.append('\n #####################################################################################')
    #print('Check point 5 : ',time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime(time.time())))
    t00 = t1

total_trail = np.array([total_trail])
total_valid = np.array([total_valid])
total_invalid = np.array([total_invalid])
total_contamination = np.array([total_contamination])
total_event = np.array([total_event])
total_closest_approach = np.array([total_closest_approach])
total_trail_ti = np.array([total_trail_ti])
L.append(f'\n Total trail candidate events: { np.sum(total_trail)}')
L.append(f'\n Validated Starlink observations :{  np.sum(total_valid)}')
L.append(f'\n Invalidated Starlink obervations:{  np.sum(total_invalid)}')
L.append(f'\n Target eclipse events: {  np.sum(eclipse)}')
L.append(f'\n peak mag. dim check:{ len(total_peak_mag)}')
L.append(f'\n ave mag. dim check:{ len(total_average_mag)}')
L.append(f'\n sun alt at ti len cehck:{ len(total_sun_alt_at_ti)}')
L.append(f'\n total ti len check:{ len(total_trail_ti)}')
L.append(f'\n sat alt at ti len cehck:{ len(total_sat_alt_at_ti)}')
L.append(f'\n Average magnitude of all trail candidate events{ np.mean(total_average_mag)}')




deltadays = int((t00 - obs_start))
# Create an array of dates for each day in September (1 to 30)
dates = [i for i in range(1,deltadays+1)]

# Define variable 'a1' and 'a2' with random values for each day (for illustration)
a1 = total_trail[0]  # Random data for demonstration
a2 = total_valid[0] # Another set of random data for comparison
a3 = total_contamination[0]
a4 = total_event[0]
a5 = total_closest_approach[0]
a6 = total_trail_ti[0]

os.makedirs(directoryf +  "/" + str(c.strftime('%H%M')))
np.savetxt(directoryf +  "/" + str(c.strftime('%H%M')) +"/totaltrail_logfile_for_" + str(day_month_i[0]) + str(day_month_i[1]) + "_" + str(year) + ".txt", a1)

# os.makedirs(directoryf +  "/" + str(c.strftime('%H%M')) +"/totalvalid_logfile_for_" + str(day_month_i[0]) + str(day_month_i[1]) + "_" + str(year) + ".txt")
np.savetxt(directoryf +  "/" + str(c.strftime('%H%M')) +"/totalvalid_logfile_for_" + str(day_month_i[0]) + str(day_month_i[1]) + "_" + str(year) + ".txt", a2)

# os.makedirs(directoryf +  "/" + str(c.strftime('%H%M')) +"/totalcontamination_logfile_for_" + str(day_month_i[0]) + str(day_month_i[1]) + "_" + str(year) +  ".txt")
np.savetxt(directoryf +  "/" + str(c.strftime('%H%M')) +"/totalcontamination_logfile_for_" + str(day_month_i[0]) + str(day_month_i[1]) + "_" + str(year) +  ".txt", a3)

# os.makedirs(directoryf +  "/" + str(c.strftime('%H%M')) +"/totalevent_logfile_for_" + str(day_month_i[0]) + str(day_month_i[1]) + "_" + str(year) +  ".txt")
np.savetxt(directoryf +  "/" + str(c.strftime('%H%M')) +"/totalevent_logfile_for_" + str(day_month_i[0]) + str(day_month_i[1]) + "_" + str(year) +  ".txt", a4)

# os.makedirs(directoryf +  "/" + str(c.strftime('%H%M')) +"/totalclosestapproach_logfile_for_" + str(day_month_i[0]) + str(day_month_i[1]) + "_" + str(year) + ".txt")
np.savetxt(directoryf +  "/" + str(c.strftime('%H%M')) +"/totalclosestapproach_logfile_for_" + str(day_month_i[0]) + str(day_month_i[1]) + "_" + str(year) + ".txt", a5)
#np.savetxt(directoryf + "/totaltrailti_logfile_for_" + str(day_month_i[0]) + str(day_month_i[1]) + "_" + str(year) + "_" + c.strftime('%H%M') + ".txt", a6)



# Set positions for each set of bars
dates = np.array(dates) # Shift dates slightly right
# Create the plot
plt.figure(figsize=(20, 6))

# Plot the bars for a1 and a2
plt.bar(dates, a4, alpha=1.0, label='Events', color='green', edgecolor='black')
plt.bar(dates, a2, alpha=1.0, label='Validated Starlinks', color='salmon', edgecolor='black')
plt.bar(dates, a1,  alpha=1.0,label='Trail Candidates', color='skyblue', edgecolor='black')
plt.bar(dates, a3, alpha=1.0, label='Strong Lens Crossing', color='red', edgecolor='black')

# Label the axes
plt.xlabel('Days')
plt.ylabel('Numbers')
plt.yscale('log')
# Title
plt.title(f'Starlink Events between {obs_start.astimezone(zone)} - {obs_end.astimezone(zone)} at beginning of the month for {Nsats} Satellites per day' )
# Add a legend
plt.legend()
# Show the plot
plt.savefig(directoryp +  "/" + c.strftime('%H%M') + "/Event_histogram_for_"  + str(obs_start.astimezone(zone)) +  ".png")



# Define variable 'a1' and 'a2' with random values for each day (for illustration)
aa1 = total_peak_mag  # Random data for demonstration
aa2 = total_average_mag # Another set of random data for comparison
# Width of bars
# Set positions for each set of bars

np.savetxt(directoryf +  "/" + str(c.strftime('%H%M')) +"/totalpeakmag_logfile_for_" + str(day_month_i[0]) + str(day_month_i[1]) + "_" + str(year) + "_" + c.strftime('%H%M') + ".txt", aa1)

np.savetxt(directoryf  +  "/" + str(c.strftime('%H%M')) +"/totalaveragemag_logfile_for_" + str(day_month_i[0]) + str(day_month_i[1]) + "_" + str(year) + "_" + c.strftime('%H%M') + ".txt", aa2)

                                                 

plt.figure(figsize=(14, 6))

bb = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
# Plot the bars for a1 and a2
plt.hist(aa1,bins=bb, label='Peak Magnitudes', alpha=0.7, color='blue', edgecolor='black')
plt.hist(aa2,bins=bb, label='Average Magnitudes', alpha=0.5, color='red', edgecolor='black')


# plt.hist(a2, bins=10, label='Average Magnitudes', color='salmon', edgecolor='black')

# Label the axes
plt.xlabel('Magnitude')
plt.ylabel('Number of Trail Candidate Events')

# Title
#plt.title('AB Magnitude Distrubutions')

# # Format the x-axis to show dates in mm/dd/yy format
# plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%m/%d/%y'))
# plt.gca().xaxis.set_major_locator(mdates.DayLocator(interval=2))  # Show every 2 days

# Rotate date labels for better readability
# plt.gcf().autofmt_xdate()

# Add a legend
plt.legend()

# Show the plot
plt.savefig(directoryp +  "/" + c.strftime('%H%M') + "/Brightness_histogram_for_" + str(obs_start.astimezone(zone)) +  ".png")

plt.figure(figsize=(14, 6))


# Plot the bars for a1 and a2
plt.hist(a5,bins=30, alpha=1, color='red', edgecolor='black')
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



end_time = time.time()
total_seconds = end_time - start_time

# Convert seconds to hours, minutes, and seconds
hours = int(total_seconds // 3600)
minutes = int((total_seconds % 3600) // 60)
seconds = total_seconds % 60

# Print the formatted runtime
print(f"Runtime: {hours}h {minutes}m {seconds:.2f}s for {d} days and {np.sum(total_valid)} valid satellite events and {np.sum(total_event)} trail candidates.")



end_time = time.time()
total_seconds = end_time - start_time

# Convert seconds to hours, minutes, and seconds
hours = int(total_seconds // 3600)
minutes = int((total_seconds % 3600) // 60)
seconds = total_seconds % 60
current_time = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime(time.time()))
starting_time = c.strftime("%m.%d")
# Print the formatted runtime
print(f"Start time: {starting_time} Current time: {current_time} Runtime: {hours}h {minutes}m {seconds:.2f}s , {d} days completed with {Nsats} satellites per day for {len(targets)} targets")
L.append(f"\n Start time: {starting_time} Current time: {current_time} Runtime: {hours}h {minutes}m {seconds:.2f}s , {d} days completed with {Nsats} satellites per dayday for {len(targets)} targets")
file.writelines(L)
file.close()



