# %%
# Parent version: skyfield_satellite_positions_B_021124.py version
# Runs with daily increments, samples every 0.1m of the satellites and checks masterfıts in the sky.
# Plots altitude, azimuth and distance. ONLY positions no target contamination check.
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

def trail_event(starlink, rubin_obs, ti, zone, sat_ra, sat_dec, sat_dist, sat_t, sat_al, sat_az, sat_hght):
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
        
  


def brightness_calculation(rubinobs_astr, sun, s, f, sat_al, sat_az, sat_hght):  
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

# Get the script name
script_name = os.path.basename(__file__)
# Configure logging
log_filename = f"/data/a.saricaoglu/Lumos-Sat/Files/{dt.datetime.now().strftime('%m.%d')}/{dt.datetime.now().strftime('%H%M')}/script.log"
os.makedirs(os.path.dirname(log_filename), exist_ok=True)
logging.basicConfig(filename=log_filename, level=logging.INFO, 
                    format='%(asctime)s - %(levelname)s - %(message)s')

# Redirect stdout and stderr to the log file
class StreamToLogger:
    def __init__(self, logger, log_level):
        self.logger = logger
        self.log_level = log_level
        self.linebuf = ''

    def write(self, buf):
        for line in buf.rstrip().splitlines():
            self.logger.log(self.log_level, line.rstrip())

    def flush(self):
        pass

sys.stdout = StreamToLogger(logging.getLogger('STDOUT'), logging.INFO)
sys.stderr = StreamToLogger(logging.getLogger('STDERR'), logging.ERROR)
# Log the start of the script with the script name
logging.info(f'Script {script_name} started')

ts = load.timescale()
year = 2025
day_month_i = [20, 1]
day_month_f = [24, 1]
zone = timezone('Chile/Continental')
local = timezone('Etc/GMT-6')
Nsats = 6412
treshold = 15 # tresold for contamination
deltatime = dt.timedelta(minutes=0.1)
obs_start = ts.utc(year,day_month_i[1], day_month_i[0],3,0,0)
obs_end = ts.utc(year,day_month_f[1], day_month_f[0],3,0,0)
starlinks = load.tle_file('/data/a.saricaoglu/Lumos-Sat/Files/Weekly_Starlink_Archive/starlinks_01.22.txt')
eph = load('de421.bsp')
earth = eph['earth']
sun = eph['sun']
rubin_obs = wgs84.latlon(-30.244633,  -70.749417)
rubinobs_astr = earth +  wgs84.latlon(30.244633*S,  70.749417*W, elevation_m = 2647)
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

fit_filename = directoryf + f"/{c.strftime('%H%M')}/satellite_data.fits"
hdu_pr = fits.PrimaryHDU()
hdu_pr.header['RUNDATE'] = c.strftime('%Y-%m-%d')
hdu_pr.header['OBSSTART'] = str(obs_start.astimezone(zone))
hdu_pr.header['OBSEND'] = str(obs_end.astimezone(zone))
hdu_pr.header['NSATS'] = Nsats
hdu_pr.header['DAYS'] = int(obs_end - obs_start)
hdu_pr.header['TARGETS'] = '/data/a.saricaoglu/lumos-sat/master.fits'
hdu_pr.header['TRESHOLD'] = treshold
hdu_pr.header['DTI'] = str(deltatime)
hdu_pr.writeto(fit_filename, overwrite=True)
with fits.open(fit_filename, mode='update') as hdu:
    print(hdu[0].header)

targets = []
with fits.open('/data/a.saricaoglu/lumos-sat/master.fits') as hdul:
    hdul.info()
    data = hdul[1].data
    ra_i = data['RAJ2000_deg']
    dec_i = data['DEJ2000_deg']
    # print(ra_i, dec_i)
    for i in range(0,len(data)):
        target_i = Star(Angle(degrees=ra_i[i]),Angle(degrees=dec_i[i]))
        targets.append(target_i)
    # for i in range(0,len(data)):
    #     target_i = sf.positionlib.position_of_radec(ra_i[i], dec_i[i], t=obs_start, center=399)
    #     targets.append(wgs84.subpoint(target_i))
hdul.close()
file = open(directoryf + "/satellite_pos_logfile_for_" + str(day_month_i[0]) + str(day_month_i[1]) + "_" + str(year) + "_" + c.strftime('%H%M') + ".txt", 'w')
L = []
description = 'DESCRIPTION: Plots altitude, azimuth and distance. Also checks contamination of n targets in every 0.1 minutes'
L.append(description)
start_time = time.time()


# if Nsats != len(starlinks):
#     random.shuffle(starlinks)
status_names = ['Rise','Culminate', 'Set']
t00 = obs_start
d = 0
eclipse = 0

total_trail = []
total_valid = []
total_event = []
total_invalid = []
total_contamination = []
total_peak_mag = []
total_separation = []
total_closest_approach = []
total_average_mag = []
total_sun_alt_at_ti = []
total_sat_heights = []
total_trail_ti = []
total_sat_alt_at_ti = []
total_peak_intensity = []
total_average_intensity = []
t00 = obs_start

hdu = fits.open(fit_filename, mode='update')
while t00 < obs_end:
    d = d + 1
    
    L.append(f'\n Day #{d}')
    fig, ax = plt.subplots()

    tar_plot= True
    sat_ra= []
    sat_dec= []
    sat_dist = []
    sat_t = []
    sat_al = []
    sat_az = []
    sat_hght = []
    sun_al = []
    sun_az = []
    sun_hght = []
    daily_separation = []
    filtered_targets = []
    trailcan_counter = 0
    valid_starlink_counter = 0
    event_counter = 0
    invalid_starlink_counter = 0
    contamination_counter = 0
    skips = 0
    t0 = t00
    t1 = t0 + dt.timedelta(days=1)
    midnight = t0 + dt.timedelta(hours=3)
    
    L.append(f'\n Observation period: {t0.astimezone(zone)} to {t1.astimezone(zone)} \n Time zone: {zone}')
    print(' midnight ', midnight.astimezone(zone))
    
    f = sf.almanac.dark_twilight_day(eph, rubin_obs)
    times, events = almanac.find_discrete(t0, t1, f)
    previous_e = f(t0).item()

    for t, e in zip(times, events):
        tstr = str(t.astimezone(zone))[:16]
        if previous_e < e:
            #print(tstr, ' ', almanac.TWILIGHTS[e], 'starts')
            if str(almanac.TWILIGHTS[e]) == 'Day':
                day_start = t
        else:
            #print(tstr, ' ', almanac.TWILIGHTS[previous_e], 'ends')
            if (almanac.TWILIGHTS[previous_e]) == 'Day':
                day_end = t
        previous_e = e


    print(f'\n Day starts at:  {day_start.astimezone(zone)}')
    print(f'\n Day ends at: {day_end.astimezone(zone)}')
    sat_num = 0
    for starlink in starlinks[:Nsats+1]:
        
        sat_num = sat_num + 1
        
        #try: 
        ev = 0
        rising_time = []
        setting_time = []
        peak_time = []
        timee, event = starlink.find_events(rubin_obs, t0, t1, altitude_degrees=30.0)
        for t, e in zip(timee, event):
            if e == 0:
                rising_time.append(t)
                event_counter = event_counter + 1
            if e == 1:
                peak_time.append(t)
            if e == 2:
                setting_time.append(t)
        
        #print('Check point 1 : ',time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime(time.time())))
        for rise_t, peak_t, set_t in zip(rising_time,peak_time, setting_time):
            sunlit = starlink.at(rise_t).is_sunlit(eph) | starlink.at(set_t).is_sunlit(eph) | starlink.at(peak_t).is_sunlit(eph)
            night = ((rise_t < day_start) | (rise_t > day_end)) & ((set_t < day_start) | (set_t > day_end)) 
            print(f'sunlit? {sunlit}, night? {night}')

            ti = rise_t
            if night & sunlit:
                L.append(f'\n Event time: {ti.astimezone(zone)}')
                trailcan_counter = trailcan_counter + 1
                print(f'Trail candidate event! Transit no: {trailcan_counter}, starlink no: {sat_num}')
                print(f'Rising time: {rise_t.astimezone(zone)}, setting time: {set_t.astimezone(zone)}')

                s = len(sat_ra)
                filt_targets = []
                skip = 0
                for target in targets:
                    topocentric_t = rubinobs_astr.at(ti).observe(target)
                    appr = topocentric_t.apparent()
                    alt_t, az_t, height = appr.altaz()
                        
                    if alt_t.degrees > 30:
                        ra, dec, dist = appr.radec()
                        target_i = sf.positionlib.position_of_radec(ra.hours, dec.degrees, t=ti, center=399)
                        filt_targets.append(wgs84.subpoint(target_i))

                filtered_targets.append(len(filt_targets))

                if tar_plot:

                    for target in filt_targets:
                        diff = target - rubin_obs
                        topo = diff.at(ti)
                        ra_t, dec_t, h_t = topo.radec()
                        # Plot the target's RA and Dec
                        ax.scatter(ra_t.hours, dec_t.degrees, marker="*", c='b', label=f'Target above 30')   
                        tar_plot = False                    

                while (ti < set_t) & ((ti < day_start) | (ti > day_end)) & starlink.at(ti).is_sunlit(eph) & (skip < 10):

                    separation = []
                    topocentric = trail_event(starlink, rubin_obs, ti, zone, sat_ra, sat_dec, sat_dist, sat_t, sat_al, sat_az, sat_hght)   
                    
                    if topocentric is None:
                        skip = skip + 1
                        print(f'Step at ti in event {trailcan_counter} for satellite {starlink.name} is taking too long, stopping. / higher than altitude limit of <2000m')
                        continue       

                    skips = skips + (skip // 10)    

                    for target in filt_targets:

                        contamination = target_check(target, rubin_obs, ti, topocentric, separation, treshold)

                        if contamination:
                            contamination_counter = contamination_counter + 1
                        elif contamination is None:
                            print(f'Event {trailcan_counter} for satellite {starlink.name} is taking too long, stopping.')

                    if len(separation) != 0 :      
                        separation = np.sort(separation)
                        total_closest_approach.append(separation[:100])
                        daily_separation.append(np.min(separation))
                    # print('Check point 2 : ',time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime(time.time())))
                    total_trail_ti.append(ti) 
                    ti = ti + deltatime
                f = len(sat_ra)

                    #print('Check point 3 : ',time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime(time.time())))
                print(f'indices {s}, {f}')    
                if f-s != 0:
                    intensities, ab_magnitudes = brightness_calculation(rubinobs_astr, sun, s, f, sat_al, sat_az, sat_hght)

                    if intensities is None:
                        print(f'Brightness calc for event {trailcan_counter} for satellite {starlink.name} is taking too long, stopping.')
                        trailcan_counter = trailcan_counter - 1
                        continue


                    sat_altitudes, sat_azimuths = \
                    np.meshgrid(
                        sat_al[s:f], 
                        sat_az[s:f])
                    sat_heights = sat_hght[s:f]
                    peak_intensity = intensities.max()
                    average_intensity = np.mean(intensities)
                    peak_ab_mag = ab_magnitudes.min()
                    average_ab_mag = np.mean(ab_magnitudes)
                    idx = np.argmin(ab_magnitudes)
                    idx2 = np.argmax(intensities)
                    peak_alt = sat_altitudes.flatten()[idx]
                    peak_az = sat_azimuths.flatten()[idx]


                    print( f"Brightest: {peak_ab_mag:0.1f} AB Magnitude, Average: {average_ab_mag:0.1f} AB Magnitude")
                    print(f'Peak intensity: {peak_intensity} W/m^2/sr, Average intensity: {average_intensity} W/m^2/sr')
                    print( f"Altitude: {peak_alt:0.0f}°, Azimuth: {peak_az:0.0f}°")
                    total_peak_mag.append(peak_ab_mag)
                    total_peak_intensity.append(peak_intensity)
                    total_average_intensity.append(average_intensity)   
                    total_average_mag.append(average_ab_mag)
                    #print('Check point 3 : ',time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime(time.time())))
            else:
                print(f'No trail candidate transit rising time {rise_t.astimezone(zone)}, setting time {set_t.astimezone(zone)}')

                            
        valid_starlink_counter = valid_starlink_counter + 1


    print(f'Total skips for day {d} : {skips}')
    print(f'Average target number in the sky for the day {d} : {np.mean(filtered_targets)}')
    total_trail.append(trailcan_counter)
    total_valid.append(valid_starlink_counter)
    total_invalid.append(invalid_starlink_counter)
    total_contamination.append(contamination_counter)
    total_event.append(event_counter)

    position_hdu = fits.BinTableHDU(Table(data=[sat_t, sat_ra, sat_dec, sat_dist, sat_az, sat_al, sat_hght], names=['Time', 'RA', 'Dec', 'Distance', 'Azimuth', 'Altitude', 'Height']))
    position_hdu.header['DAY'] = str(d)
    position_hdu.header['VALSAT'] = valid_starlink_counter
    position_hdu.header['INVSAT'] = invalid_starlink_counter
    position_hdu.header['EVENTS'] = event_counter
    position_hdu.header['TRAILCAN'] = trailcan_counter
    position_hdu.header['CONT'] = contamination_counter
    position_hdu.header['CLOSESTD'] = np.min(daily_separation)

    hdu.append(position_hdu)

    L.append(f'\n Day # : {d}')
    L.append(f'\n Validated Starlink number : {valid_starlink_counter}')
    L.append(f'\n Invalidated Starlink number : {invalid_starlink_counter}')
    L.append(f'\n Event number : {event_counter}')
    L.append(f'\n Trail candidate number : {trailcan_counter}')
    L.append(f'\n Contaminated target number: {contamination_counter}')
    if len(daily_separation) != 0:
        L.append(f'\n Closest target-satellite event of the day: {np.min(daily_separation)} arcseconds')
    print('Day #',d)
    print('Validated Starlink events :', valid_starlink_counter)
    print('Invalidated Starlink events :', invalid_starlink_counter)
    print('trail candidates', trailcan_counter)
    print('Contaminated targets : ', contamination_counter)
    if len(daily_separation) != 0:
        print(f'\n Closest target-satellite event of the day: {np.min(daily_separation)} arcseconds')
    #print('Check point 4 : ',time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime(time.time())))
    
    scatter = ax.scatter(sat_ra, sat_dec, marker='o', c = sat_dist, cmap = 'viridis')  

    plt.title(f'Starlink positions for {t00.astimezone(zone)}')
    plt.grid()
    ax.xaxis.set_major_formatter(FuncFormatter(ra_formatter))
    ax.set_xlabel('RA (hh:mm:ss)')
    #fig.colorbar(scatter).set_label("Distance [km]")

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
hdu.close()
total_trail = np.array([total_trail])
total_valid = np.array([total_valid])
total_invalid = np.array([total_invalid])
total_contamination = np.array([total_contamination])
total_event = np.array([total_event])
total_closest_approach = np.array([x for y in total_closest_approach for x in y ])
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

print(f'\n Total trail candidate events: { np.sum(total_trail)}')
print(f'\n Validated Starlink observations :{  np.sum(total_valid)}')
print(f'\n Invalidated Starlink obervations:{  np.sum(total_invalid)}')
print(f'\n Target eclipse events: {  np.sum(eclipse)}')
print(f'\n peak mag. dim check:{ len(total_peak_mag)}')
print(f'\n ave mag. dim check:{ len(total_average_mag)}')
print(f'\n sun alt at ti len cehck:{ len(total_sun_alt_at_ti)}')
print(f'\n total ti len check:{ len(total_trail_ti)}')
print(f'\n sat alt at ti len cehck:{ len(total_sat_alt_at_ti)}')
print(f'\n Average magnitude of all trail candidate events{ np.mean(total_average_mag)}')


deltadays = int((t00 - obs_start))
# Create an array of dates for each day in September (1 to 30)
dates = [i for i in range(1,deltadays+1)]

# Define variable 'a1' and 'a2' with random values for each day (for illustration)
a1 = total_trail[0]  # Random data for demonstration
a2 = total_valid[0] # Another set of random data for comparison
a3 = total_contamination[0]
a4 = total_event[0]
a5 = total_closest_approach
a6 = total_trail_ti[0]

if not os.path.exists(directoryf +  "/" + str(c.strftime('%H%M'))): 
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

# Define variable 'a1' and 'a2' with random values for each day (for illustration)
aa1 = total_peak_mag  # Random data for demonstration
aa2 = total_average_mag # Another set of random data for comparison
# Width of bars
# Set positions for each set of bars
bb1 = total_peak_intensity
bb2 = total_average_intensity

np.savetxt(directoryf +  "/" + str(c.strftime('%H%M')) +"/totalpeakmag_logfile_for_" + str(day_month_i[0]) + str(day_month_i[1]) + "_" + str(year) + ".txt", aa1)
np.savetxt(directoryf  +  "/" + str(c.strftime('%H%M')) +"/totalaveragemag_logfile_for_" + str(day_month_i[0]) + str(day_month_i[1]) + "_" + str(year) + ".txt", aa2)

np.savetxt(directoryf  +  "/" + str(c.strftime('%H%M')) +"/totalpeakintensity_logfile_for_" + str(day_month_i[0]) + str(day_month_i[1]) + "_" + str(year) + ".txt", bb1)                                                
np.savetxt(directoryf  +  "/" + str(c.strftime('%H%M')) +"/totalaverageintensity_logfile_for_" + str(day_month_i[0]) + str(day_month_i[1]) + "_" + str(year) +".txt", bb2)

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

plt.xlabel('Intensity [W/m^2/sr]')
plt.ylabel('Number of Trail Candidate Events')

plt.legend()

plt.savefig(directoryp +  "/" + c.strftime('%H%M') + "/Intensity_histogram_for_" + str(obs_start.astimezone(zone)) +  ".png")


plt.figure(figsize=(14, 6))


# Plot the bars for a1 and a2
cc = np.arange(0,500,10)
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
current_time = dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
starting_time = c.strftime("%Y-%m-%d %H:%M:%S")
# Print the formatted runtime
print(f"Start time: {starting_time} Current time: {current_time} Runtime: {hours}h {minutes}m {seconds:.2f}s , {d} days completed with {Nsats} satellites per day for {len(targets)} targets")
L.append(f"\n Start time: {starting_time} Current time: {current_time} Runtime: {hours}h {minutes}m {seconds:.2f}s , {d} days completed with {Nsats} satellites per dayday for {len(targets)} targets")
file.writelines(L)
file.close()
