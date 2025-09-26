# %%
# Parent version: skyfield_101224_positions_and_contamination_final_ver4.2
# Runs with daily increments, samples every 0.1m of the satellites and checks masterfıts in the sky.
# This version simulates whole sky per second (satellite and time incremation swapped)
import os
import sys
import logging
import numpy as np
from astropy import units as u
from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord, EarthLocation
import skyfield as sf
from skyfield.api import load, wgs84, EarthSatellite, N, W,S,E, Star, Angle
import datetime as dt
import time
from skyfield import almanac
from skyfield.framelib import galactic_frame
from skyfield.data import hipparcos
import lumos.calculator
import random
from src.satellites.starlink_sat import *
from src.utils.calculations import *
from src.utils.plots import *
from src.utils.utils import *
from src.tracking.satellite_tracker import SatelliteTracker
from src.tracking.target_tracker import TargetTracker
from src.tracking.exposure_tracker import ExposureTracker
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

now = datetime.now()
cwd = os.getcwd()
print(cwd)


pathToTargets = os.path.join(cwd, "in_LSST_footprint_coordinates.fits")
pathToCalculatedTargets = os.path.join(cwd, "Skyfield/calculated_targets")
pathToTLE = os.path.join(cwd, "Skyfield/tle_data")
pathToFiles = os.path.join(cwd, "Skyfield/files")
pathToLogs = os.path.join(cwd, "Skyfield/logs")
pathToPlots = os.path.join(cwd, "Skyfield/plots")
# Create output directory if it doesn't exist
os.makedirs(pathToFiles, exist_ok=True)
os.makedirs(pathToLogs, exist_ok=True)
os.makedirs(pathToPlots, exist_ok=True)

# Get the script name
script_name = os.path.basename(__file__)
# Configure logging
log_filename = f"{pathToLogs}/{script_name.strip('.py')}_{now.strftime('%m.%d')}_{now.strftime('%H%M')}.log"
print(log_filename)
os.makedirs(os.path.dirname(log_filename), exist_ok=True)
logging.basicConfig(filename=log_filename, level=logging.INFO, 
                    format='%(asctime)s - %(levelname)s - %(message)s')

# Redirect stdout and stderr to the log file
sys.stdout = StreamToLogger(logging.getLogger('STDOUT'), logging.INFO)
sys.stderr = StreamToLogger(logging.getLogger('STDERR'), logging.ERROR)
# Log the start of the script with the script name
logging.info(f'Script {script_name} started')

ts = load.timescale()
year = 2025
day_month_i = [1, 8]
day_month_f = [31, 8]
zone_stl = timezone('Etc/GMT-6')
zone = timezone('Chile/Continental')
utc = timezone('UTC')
obs_start = dt.datetime(year,day_month_i[1], day_month_i[0],0,0,0,0)
obs_end = dt.datetime(year,day_month_f[1], day_month_f[0],0,0,0,0)
local_start = ts.from_datetime(zone.localize(obs_start).astimezone(utc))
local_end = ts.from_datetime(zone.localize(obs_end).astimezone(utc))
# Nsats = 50 # number of starlinks to check
threshold = 15 # tresold for contamination
deltatime = dt.timedelta(seconds=1)
deltaday = dt.timedelta(days=4)
twilight_extension = dt.timedelta(hours=1.5)
postexposure_skip = dt.timedelta(minutes=10)
#starlinks = load.tle_file('https://celestrak.org/NORAD/elements/gp.php?GROUP=starlink&FORMAT=tle')
eph = load('de421.bsp')
earth = eph['earth']
sun = eph['sun']
rubin_obs = wgs84.latlon(-30.244633,  -70.749417)
rubinobs_astr = earth +  wgs84.latlon(30.244633*S,  70.749417*W, elevation_m = 2647)


targets, names = get_targets_as_stars(pathToTargets, region='In_Euclid_Lowdust_Region')

file = open(pathToLogs + "/satellite_pos_logfile_for_" + obs_start.strftime('%d%m_%H%M') + "_" +now.strftime('%d%m_%H%M') + ".txt", 'w')
L = []
description = f'DESCRIPTION: Plots altitude, azimuth and distance. Also checks contamination of n targets in every {deltatime.total_seconds()/60} minutes'
L.append(description)
L.append(f'\n Observation period: {local_start.astimezone(zone).strftime("%Y-%m-%d %H:%M")} to {local_end.astimezone(zone).strftime("%Y-%m-%d %H:%M")} \n Time zone: {zone}')
print(f'Observation period: {local_start.astimezone(zone).strftime("%Y-%m-%d %H:%M")} to {local_end.astimezone(zone).strftime("%Y-%m-%d %H:%M")} \n Time zone: {zone}')
start_time = time.time()



schedule = get_time_array(local_start, local_end, eph, rubin_obs, zone, twilight_extension, deltatime, deltaday, postexposure_skip)

trailcan_counter = 0
valid_starlink_counter = 0
event_counter = 0
contamination_counter = 0
exposure_counter = 0

persecond_plot= False
exposure_plot= False

target_tracker = TargetTracker(targets, names, rubinobs_astr) #if rubin_obs is not given, rubin_obs_geocentric is taken as EarthLocation.of_site('rubin')
satellite_tracker = SatelliteTracker()

for day in schedule:

    fit_filename = os.path.join(pathToFiles, f"skyfield_data_{day[0][0].astimezone(zone).strftime('%Y%m%d')}.fits")
    hdu_pr = fits.PrimaryHDU()
    hdu_pr.header['RUNDATE'] = now.strftime('%Y-%m-%d')
    hdu_pr.header['OBSSTART'] = str(local_start.astimezone(zone).day)
    hdu_pr.header['OBSEND'] = str(local_end.astimezone(zone).day)
    # hdu_pr.header['TARGETS'] = '/data/a.saricaoglu/lumos-sat/master.fits'
    hdu_pr.header['THRSHLD'] = threshold
    hdu_pr.header['DTIME'] = str(deltatime)
    hdu_pr.header['ZONE'] = str(zone)
    hdu_pr.header['DDAY'] = str(deltaday)
    hdu_pr.header['POSTEXP'] = str(postexposure_skip)
    hdu_pr.header['TWILIGHT'] = str(twilight_extension)
    hdu_pr.writeto(fit_filename, overwrite=True)
    hdu = fits.open(fit_filename, mode='update')
    hdr = hdu[0].header

   #L.append(f'\n Day #{day[0][0].astimezone(zone).day}')
    print(f'\n Event day starts at:  {day[0][0].astimezone(zone)}')
    print(f'\n Event day ends at: {day[-1][0].astimezone(zone)}')

    tle, starlinks = find_closest_TLE(day[0][0], pathToTLE, zone)
    Nsats = len(starlinks)


    for exp in day:
        e_start = time.time()
        exposure_counter += 1
        exposure = ExposureTracker(exposure_id=exposure_counter, start_time=exp[0], end_time=exp[-1], observer_location=rubin_obs)

        satellites_up, satellite_periods = SatelliteTracker.find_satellites_up(starlinks, rubin_obs, exposure.start_time, exposure.end_time, altitude_degrees=30.0, return_periods=True)
        print(f'Number of satellites up at throughout the exposure: {len(satellites_up)}')
        

        if len(satellites_up) == 0:
            print(f'No satellites are up during the exposure starting at {exposure.start_time.astimezone(zone)}, skipping to next exposure')
            continue

        
        satellite_trackers = {}
        for starlink in satellites_up:
            satellite_trackers[starlink] = SatelliteTracker()
           

        targets_up = target_tracker.calculate_targets_for_exposure(exp, rubinobs_astr, zone, horizon_degrees=30.0)
        print(f'Number of targets up at start of exposure: {len(targets_up)}')

        for starlink in satellites_up:
            

            satellite_tracker = satellite_trackers[starlink]
            satellite_period = satellite_periods[starlink]['combined_periods']
            satellite_visible_duration = satellite_periods[starlink]['total_visible_duration']
            period_start_pos = point_trail_event(starlink, rubin_obs, satellite_period[0][0])
            period_end_pos = point_trail_event(starlink, rubin_obs, satellite_period[0][1])
            target_box = filter_targets(target_tracker.get_targets_at_time(satellite_period[0][1])['skycoords'], target_tracker.get_targets_at_time(satellite_period[0][0])['skycoords'], period_end_pos, period_start_pos)
            exposure.increment_counters(total_events=int(satellite_visible_duration))

            if np.sum(target_box) == 0:
                print(f'No targets in the path of satellite {starlink.name} during its visible period, skipping to next satellite')
                exposure.increment_counters(no_targets_skip=1)
                continue
            ti_skip = 0

            for ti in exp:

                if not SatelliteTracker.is_time_in_periods(ti, satellite_period):
                    print(f'Satellite {starlink.name} is not above altitude limit at time {ti.astimezone(zone)}, skipping to next time step')
                    continue

                #print(f'\n Trail calculation for satellite {starlink.name} at time {ti.astimezone(zone)} for exposure starting at {exp[0].astimezone(zone)}')

                topocentric = satellite_tracker.point_trail_event(starlink, rubin_obs, ti)

                if topocentric is None:
                    exposure.increment_counters(no_topo_skip=1)
                    print(f'Satellite {starlink.name} could not be tracked at time {ti.astimezone(zone)}, skipping to next time step')
                    continue

                satellite_data = {}
                contamination_data = {}

                exposure.increment_counters(trail_candidates=1)

                if satellite_tracker.has_current_position():
                    current_position = satellite_tracker.get_current_position()
                    try:
                        intensity, ab_magnitude = satellite_tracker.brightness_calculation_current(rubinobs_astr, sun)
                        satellite_data[starlink.name] = {
                            'alt': current_position['alt'],
                            'az': current_position['az'],
                            'height': current_position['height'],
                            'time': current_position['time'],
                            'intensity': intensity,
                            'ab_magnitude': ab_magnitude
                        }
                    except Exception as e:
                        exposure.increment_counters(no_brightness_skip=1)
                        print(f'Brightness calculation failed for {starlink.name}: {e}')
                        satellite_data[starlink.name] = {
                            'alt': current_position['alt'],
                            'az': current_position['az'],
                            'height': current_position['height'],
                            'time': current_position['time'],
                            'intensity': np.nan,
                            'ab_magnitude': np.nan
                        }

                # Previous position is the same satellite at ti-1
                if satellite_tracker.has_both_positions():
                    # For sanity check, otherwise irrelevant. Delete later.

                    # # Calculate angular velocity for this specific satellite
                    velocity = satellite_tracker.calculate_angular_velocity()
                    satellite_data[starlink.name].update({'angular_velocity': velocity})
                    # print(f"Satellite {starlink.name} angular velocity: {velocity:.2f} arcsec/s")

                    # Use previous position for contamination analysis
                    prev_pos = satellite_tracker.get_previous_position()
                    curr_pos = satellite_tracker.get_current_position()

                    # Calculate angular separation between current and previous positions
                    # difference_angle_sat = curr_pos['topocentric'].separation_from(prev_pos['topocentric'])
                    # print(f"Satellite {starlink.name} moved {difference_angle_sat.arcseconds():.2f} arcseconds")
                    # if difference_angle_sat.arcseconds() < threshold:
                    #     print(f"Satellite {starlink.name} moved less than {threshold} arcseconds, skipping contamination check")
                    #     ti_skip += 1
                    #     continue
                    # Contamination check using consecutive positions
                    prev_ti_index = (exp.index(ti) - 1)
                    prev_ti = exp[prev_ti_index]
                    
                    if target_tracker.has_time_cached(prev_ti):
                        # Get target positions for current and previous times
                        current_targets = target_tracker.get_targets_at_time(ti)
                        previous_targets = target_tracker.get_targets_at_time(prev_ti)

                        # Contamination analysis with proper satellite tracking
                        separation = []

                        contamination, contaminated_targets, names = calculated_target_check(
                            L, 
                            current_targets,
                            previous_targets,
                            topocentric, prev_pos['topocentric'],
                            separation, threshold
                        )
                        if contamination:
                            exposure.increment_counters(contaminations=1)
                            contamination_data = {'contaminated_targets': contaminated_targets,
                                                  'contaminated_names': names,
                                                  'separations': separation
                                                  }
                            print(f'Contamination detected for satellite {starlink.name} at time {ti.astimezone(zone)} with {len(contaminated_targets)} targets')
                        else:
                            contamination_data = {'separations': separation}
                        #satellite_tracker.print_current_status()
                        
                # Update exposure with satellite and contamination data

                target_data = target_tracker.get_targets_at_time(ti) if target_tracker.has_time_cached(ti) else {}
                exposure.add_timestep_data(ti, satellite_data, target_data, contamination_data)
                
                # if persecond_plot:
                #     plt.figure(figsize=(10, 5))

                #     ax = plt.subplot(111, projection='aitoff')
                #     plt.grid(True)
                #     sat_time = [t.to_astropy() for t in satellite_data['time']]
                #     filtered_targets_altaz = SkyCoord(target_data['skycoords'])
                #     ax.scatter(-filtered_targets_altaz.az.wrap_at(180 * u.degree).radian, 
                #     filtered_targets_altaz.alt.wrap_at(180 * u.degree).radian, s=2, c='green', alpha=0.7)
                #     contaminated_tar = SkyCoord(contamination_data.get('contaminated_targets', []))
                #     if len(contaminated_tar) != 0:
                #         ax.scatter(-contaminated_tar.az.wrap_at(180 * u.degree).radian, 
                #         contaminated_tar.alt.wrap_at(180 * u.degree).radian, s=5, c='red', alpha=0.7, label='Contaminated targets')
                #         plt.legend(loc='upper right')
                #     sat_coord = SkyCoord(az=satellite_data['az'], alt=satellite_data['alt'], frame='altaz', unit='deg', obstime=sat_time, location=EarthLocation.of_site('rubin'))
                #     scatter = ax.scatter(-sat_coord.az.wrap_at(180 * u.degree).radian, sat_coord.alt.radian, marker='o', c = satellite_data['height'], cmap = 'plasma')
                #     plt.title(f'Starlink positions for {ti.astimezone(zone)}', y=1.08)
                #     ax.set_xlabel('Azimuth (degrees)')
                #     #fig.colorbar(scatter).set_label("Distance [km]")

                #     # Dec is already in degrees, so just label it
                #     ax.set_ylabel('Altitude (degrees)')    
                #     # ax.set_ylim(np.radians(30), np.radians(90))  # Set y-axis limits to show only altitudes above 30 degrees
                #     plt.tight_layout()
                #     plt.savefig(pathToPlots +   "/Satellite_positions_for_" + ti.utc_iso(delimiter='_') +  "_A.png")
                #     plt.close()
            
            satellite_tracker.reset_all_positions()
            satellite_tracker.clear_history()
        # After all time steps in the exposure are processed, compile and save exposure data
        target_tracker.print_cache_status()
        exposure.print_summary()

        exposure_extensions = exposure.create_fits_extensions()
        for ext in exposure_extensions:
            hdu.append(ext)

        print(f'Exposure HDUs appended for {exposure_counter}')
        print(f'Exposure processed in {e_start - time.time():.2f} seconds')
        # if exposure_plot and exposure.total_events > 0:
        #     plt.figure(figsize=(10, 5))

        #     ax = plt.subplot(111, projection='aitoff')
        #     plt.grid(True)
        #     sat_time = [pos['time'].to_astropy() for pos in exposure.all_satellite_positions]
        #     filtered_targets_altaz = SkyCoord(target_tracker.get_targets_at_time(exp[0])['skycoords'])
        #     ax.scatter(-filtered_targets_altaz.az.wrap_at(180 * u.degree).radian, 
        #     filtered_targets_altaz.alt.wrap_at(180 * u.degree).radian, s=2, c='green', alpha=0.7, label='Targets')
            

            
        #     sat_coord = SkyCoord(alt=[pos['alt'] for pos in exposure.all_satellite_positions],
        #                          az=[pos['az'] for pos in exposure.all_satellite_positions],
        #                          frame='altaz', unit='deg')
        #     scatter = ax.scatter(-sat_coord.az.wrap_at(180 * u.degree).radian, sat_coord.alt.radian, marker='o', c = [pos['height']/1000 for pos in exposure.all_satellite_positions], cmap = 'plasma', label='Satellites')
        #     if len(exposure.all_contaminated_targets) != 0:
        #         for contaminated_tar in [ct['position'] for ct in exposure.all_contaminated_targets]:
        #             ax.scatter(-contaminated_tar.az.wrap_at(180 * u.degree).radian, 
        #             contaminated_tar.alt.wrap_at(180 * u.degree).radian, marker='x', c='red', alpha=1, label='Contaminated targets')
        #     plt.title(f'Starlink positions for {exp[0].astimezone(zone)} to {exp[-1].astimezone(zone)}', y=1.08)
        #     ax.set_xlabel('Azimuth (degrees)')
        #     plt.colorbar(scatter).set_label("Distance [km]")

        #     # Dec is already in degrees, so just label it
        #     ax.set_ylabel('Altitude (degrees)')    
        #     # ax.set_ylim(np.radians(30), np.radians(90))  # Set y-axis limits to show only altitudes above 30 degrees
        #     plt.tight_layout()
        #     plt.legend(loc='upper right')
        #     plt.savefig(pathToPlots +   "/Satellite_exposure_for_" + exp[0].astimezone(zone).strftime('%Y-%m-%d_%H-%M-%S') +  ".png")
        #     plt.close()


    



    hdu.flush()
    hdu.close()
    end_time = time.time()
    total_seconds = end_time - start_time
    # Convert seconds to hours, minutes, and seconds
    hours = int(total_seconds // 3600)
    minutes = int((total_seconds % 3600) // 60)
    seconds = total_seconds % 60
    current_time = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime(time.time()))
    # Print the formatted runtime
    print(f"Current time by day end: {current_time} Runtime: {hours}h {minutes}m {seconds:.2f}s ")
   #L.append(f"\n Current time: {current_time} Runtime: {hours}h {minutes}m {seconds:.2f}s , Day {day[0][0].astimezone(zone).day} completed ")
   #L.append('\n #####################################################################################')
    #print('Check point 5 : ',time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime(time.time())))



