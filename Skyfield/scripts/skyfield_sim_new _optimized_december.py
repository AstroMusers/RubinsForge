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
from astropy.visualization import quantity_support, time_support, wcsaxes
quantity_support()  

now = datetime.now()
cwd = os.getcwd()

pathToTargets = os.path.join(cwd, "in_LSST_footprint_coordinates.fits")
pathToCalculatedTargets = os.path.join(cwd, "Skyfield/calculated_targets")
pathToTLE = "/data/a.saricaoglu/repo/satchecker/tle_snapshots_reconstructed"
pathToFiles = os.path.join(cwd, "Skyfield/files/cont_improvement_satchecker_final")
pathToLogs = os.path.join(cwd, "Skyfield/logs")
pathToPlots = os.path.join(cwd, "Skyfield/plots")
exposures_file_path = '/data/a.saricaoglu/repo/RubinsForge/OpSim/scripts/twilight_visits_december2025.fits'

# Create output directory if it doesn't exist
os.makedirs(pathToFiles, exist_ok=True)
os.makedirs(pathToLogs, exist_ok=True)
os.makedirs(pathToPlots, exist_ok=True)

# Get the script name
script_name = os.path.basename(__file__)
log_filename = (
    f"{pathToLogs}/"
    f"{os.path.splitext(script_name)[0]}_"
    f"{now.strftime('%m.%d_%H%M')}.log"
)

logging.basicConfig(
    filename=log_filename,
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)

logger = logging.getLogger(__name__)

logger.info("cwd: %s", cwd)
logger.info("log file: %s", log_filename)

# Redirect stdout and stderr to the log file
sys.stdout = StreamToLogger(logging.getLogger('STDOUT'), logging.INFO)
sys.stderr = StreamToLogger(logging.getLogger('STDERR'), logging.ERROR)
# Log the start of the script with the script name
logging.info(f'Script {script_name} started')

logger.info('tle snapshot path: %s', pathToTLE)

ts = load.timescale()
year = 2025
day_month_i = [1, 12]
day_month_f = [31, 12]
zone_stl = timezone('Etc/GMT-6')
zone = timezone('Chile/Continental')
utc = timezone('UTC')
obs_start = dt.datetime(year,day_month_i[1], day_month_i[0],0,0,0,0)
obs_end = dt.datetime(year,day_month_f[1], day_month_f[0],0,0,0,0)
local_start = ts.from_datetime(zone.localize(obs_start).astimezone(utc))
local_end = ts.from_datetime(zone.localize(obs_end).astimezone(utc))
# Nsats = 50 # number of starlinks to check
threshold = 67.5 # tresold for contamination
deltatime = dt.timedelta(seconds=1)
deltaday = dt.timedelta(days=4)
twilight_extension = dt.timedelta(hours=1.1)
postexposure_skip = dt.timedelta(minutes=10)
#starlinks = load.tle_file('https://celestrak.org/NORAD/elements/gp.php?GROUP=starlink&FORMAT=tle')
eph = load('de421.bsp')
earth = eph['earth']
sun = eph['sun']
rubin_obs = wgs84.latlon(-30.244633,  -70.749417, elevation_m = 2647)
rubinobs_astr = earth +  wgs84.latlon(30.244633*S,  70.749417*W, elevation_m = 2647)



exposure_star_dict = exposure_targets_to_stars(exposures_file_path)


# file = open(pathToLogs + "/satellite_pos_logfile_for_" + obs_start.strftime('%d%m_%H%M') + "_" +now.strftime('%d%m_%H%M') + ".txt", 'w')
# L = []
# description = f'DESCRIPTION: Plots altitude, azimuth and distance. Also checks contamination of n targets in every {deltatime.total_seconds()/60} minutes'
# L.append(description)
# L.append(f'\n Observation period: {local_start.astimezone(zone).strftime("%Y-%m-%d %H:%M")} to {local_end.astimezone(zone).strftime("%Y-%m-%d %H:%M")} \n Time zone: {zone}')
logger.info('Observation period: %s to %s | Time zone: %s',
            local_start.astimezone(zone).strftime("%Y-%m-%d %H:%M"),
            local_end.astimezone(zone).strftime("%Y-%m-%d %H:%M"),
            zone)
start_time = time.time()



# schedule = get_time_array(local_start, local_end, eph, rubin_obs, zone, twilight_extension, deltatime, deltaday, postexposure_skip)
exposure_data = get_exposure_windows(
    exposures_file_path,
    step_seconds=deltatime.total_seconds(),
)
schedule = list(exposure_data.keys())
logger.info('Number of days in the schedule: %s', len(schedule))
logger.info('Number of exposures in the schedule: %s', sum(len(exposure_data[day]) for day in schedule))
# print(f'Index of the day december 16th : {[(i, day) for i, day in enumerate(schedule)]}')
# quit()

trailcan_counter = 0
valid_starlink_counter = 0
event_counter = 0
contamination_counter = 0
exposure_counter = 0

persecond_plot= False
exposure_plot= False

satellite_tracker = SatelliteTracker()
# # Check if ids and positions match for two kind of target populations above
# matches = 0
# for et in target_tracker_exp.get_names_up(local_start):
#     if any(et == target for target in target_tracker.get_names_up(local_start)):
#         matches += 1
# print('Number of matching targets at start of observation between two target populations:', matches)
snapshot_index = build_snapshot_index(pathToTLE)
logger.info(
    "TLE snapshot range: MJD %.8f -- %.8f",
    snapshot_index[0][0],
    snapshot_index[-1][0],
)
tle_snapshot_cache = {}

for day in list(exposure_data.keys()):  
    
    # print(exposure_data.keys())
    # quit()
    # Stop script if the day is after April 21st
    ex = list(exposure_data[day].values())
    if ex[0]['window'][0].astimezone(zone).day <= 22:
        print(f'Quitting by day {day} as it is before december 5th')
        # quit()
        continue
    fit_filename = os.path.join(pathToFiles, f"skyfield_data_visit_{day}.fits")
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
    logger.info('Event day %s', day)



    for visit in list(exposure_data[day].values()):
        e_start = time.time()
        exposure_counter += 1
        exp = visit['window']

        snapshot_path, snapshot_mjd, snapshot_offset_hours, starlinks, tle_df = (
            load_starlink_snapshot(
                visit_time=exp[0],
                snapshot_index=snapshot_index,
                ts=ts,
                max_age_hours=4,
                cache=tle_snapshot_cache,
            )
        )
        tle_metadata = {
            sat: {
                "tle_quality": row.tle_quality,
                "tle_epoch_mjd": float(row.tle_epoch_mjd),
            }
            for sat, row in zip(
                starlinks,
                tle_df.itertuples(index=False),
            )
        }
        Nsats = len(starlinks)

        logger.info(
            "Visit %s using TLE snapshot %s | "
            "snapshot MJD %.8f | offset %.2f h | %d Starlinks",
            visit["visit_id"],
            snapshot_path.name,
            snapshot_mjd,
            snapshot_offset_hours,
            Nsats,
        )

        # if visit['visit_id'] != 78425:
        #     print(f'Skipping visit {visit["visit_id"]} for day {day} as it is not the target visit')
        #     continue
        # Skip dates before aril 21st
        # if exp[0].astimezone(zone).day >= 21:
        #     print(f'Skipping visit {visit["visit_id"]} for day {day} as it is before April 21st')
        #     continue
        logger.info('Visit ID: %s, Exposure from %s to %s, total %s time steps of %s seconds',
                visit["visit_id"], exp[0].astimezone(zone), exp[-1].astimezone(zone), len(exp), deltatime.total_seconds())
        logger.info('Targets for visit %s: %s', visit["visit_id"], exposure_star_dict[visit['visit_id']])

        exposure_stars = exposure_star_dict[visit['visit_id']]
        targets, names = exposure_stars['stars'], exposure_stars['names']
        target_tracker = TargetTracker(targets, names, rubinobs_astr)
        assert len(targets) == len(names)
        assert len(targets) == len(np.unique(visit["target_ids"]))

        visit["tle_snapshot"] = snapshot_path.name
        visit["tle_snapshot_mjd"] = snapshot_mjd
        visit["tle_snapshot_offset_hours"] = snapshot_offset_hours

        exposure = ExposureTracker(exposure_id=visit['visit_id'], start_time=exp[0], end_time=exp[-1], visit_info=visit)

        satellites_up, satellite_periods = (
            SatelliteTracker.find_satellites_up(
                starlinks,
                rubin_obs,
                exposure.start_time,
                exposure.end_time,
                altitude_degrees=15.0,
                check_sunlit=False,
                eph=eph,
                return_periods=True,
            )
        )
        logger.info('Number of satellites up throughout the exposure: %s', len(satellites_up))
        

        if len(satellites_up) == 0:
            logger.info('No satellites are up during exposure starting at %s, skipping', exposure.start_time.astimezone(zone))

            for ext in exposure.create_fits_extensions():
                hdu.append(ext)

            continue

        
        satellite_trackers = {}
        for starlink in satellites_up:
            satellite_trackers[starlink] = SatelliteTracker()

        visit_center_as_star = (
            visit["visit_center_as_star"]
        )

        starlinks_path_over_visit = []
        for starlink in satellites_up:
            
            satellite_period = satellite_periods[starlink]['combined_periods']
            satellite_visible_duration = satellite_periods[starlink]['total_visible_duration']
            exposure.increment_counters(total_events=int(satellite_visible_duration))

            for period in satellite_period:
                
                if not satellite_period_crosses_visit(
                        starlink,
                        rubin_obs,
                        rubinobs_astr,
                        period[0],
                        period[1],
                        exp,
                        visit_center_as_star,
                    ):
                    exposure.increment_counters(no_targets_skip=1)
                    logger.debug('Visit region does not cover path of %s for visible period index %s, skipping satellite',
                                 starlink.name, satellite_period.index(period))
                elif starlink in starlinks_path_over_visit:
                    logger.info('Satellite %s path already checked for visit crossing, proceeding', starlink.name)
                else:
                    starlinks_path_over_visit.append(starlink)
                    logger.info('Satellite %s path crosses visit region, proceeding', starlink.name)


        if len(starlinks_path_over_visit) == 0:
            logger.info('No satellites pass over visit region during exposure starting at %s, skipping', exposure.start_time.astimezone(zone))
            for ext in exposure.create_fits_extensions():
                hdu.append(ext)

            continue

        all_target_times = list(exp)

        for starlink in starlinks_path_over_visit:
            for period_start, period_end in satellite_periods[
                starlink
            ]["combined_periods"]:
                all_target_times.extend([
                    period_start,
                    period_end,
                ])

        all_target_times = sorted(
            set(all_target_times),
            key=lambda t: t.tt,
        )

        targets_for_visit = target_tracker.calculate_targets_for_exposure(
            all_target_times,
            rubinobs_astr,
            zone,
        )

        logger.info('Number of targets in the visit: %s names: %s', len(visit["target_ids"]), visit["target_ids"])
            
        for starlink in starlinks_path_over_visit:

            satellite_periods_for_starlink = satellite_periods[
                starlink
            ]["combined_periods"]

            satellite_tracker = satellite_trackers[starlink]
    
            for period_start, period_end in satellite_periods_for_starlink:

                satellite_tracker.reset_all_positions()
                # Build sampling times for this visibility period
                sample_times = [period_start]

                # Add all exposure-grid times strictly inside the period
                for ti in exp:
                    if period_start < ti < period_end:
                        sample_times.append(ti)

                # Always include exact end boundary
                sample_times.append(period_end)

                # Optional but good: remove duplicates and keep sorted
                sample_times = sorted(
                    set(sample_times),
                    key=lambda t: t.tt
                )

                for ti in sample_times:

                    topocentric = satellite_tracker.point_trail_event(
                        starlink,
                        rubin_obs,
                        ti
                    )

                    #print(f'\n Trail calculation for satellite {starlink.name} at time {ti.astimezone(zone)} for exposure starting at {exp[0].astimezone(zone)}')

                    if topocentric is None:
                        exposure.increment_counters(no_topo_skip=1)
                        logger.info('Satellite %s could not be tracked at time %s, skipping to next time step', starlink.name, ti.astimezone(zone))
                        continue

                    satellite_data = {}
                    contamination_data = {}

                    exposure.increment_counters(trail_candidates=1)

                    if satellite_tracker.has_current_position():
                        current_position = satellite_tracker.get_current_position()

                        satellite_data[starlink.name] = {
                            'alt': current_position['alt'],
                            'az': current_position['az'],
                            'height': current_position['height'],
                            'time': current_position['time'],
                            'intensity': np.nan,
                            'ab_magnitude': np.nan,
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

                        # Calculate angular separation between current and previous positions
                        # difference_angle_sat = curr_pos['topocentric'].separation_from(prev_pos['topocentric'])
                        # print(f"Satellite {starlink.name} moved {difference_angle_sat.arcseconds():.2f} arcseconds")
                        # if difference_angle_sat.arcseconds() < threshold:
                        #     print(f"Satellite {starlink.name} moved less than {threshold} arcseconds, skipping contamination check")
                        #     ti_skip += 1
                        #     continue
                        
                        # Get target positions for current and previous times
                        current_targets = target_tracker.get_targets_at_time(ti)
                        previous_time = prev_pos["time"]
                        previous_targets = target_tracker.get_targets_at_time(previous_time)

                        sunlit = []
                        event_times = []
                        contaminated_intensities = []
                        contaminated_ab_magnitudes = []

                        contamination, contaminated_targets, names, contaminated_t_closest, contaminated_separations, separations = calculated_target_check(
                            current_targets,
                            previous_targets,
                            topocentric,
                            prev_pos["topocentric"],
                            threshold,
                        )

                        if contamination:

                            for idx, fraction in enumerate(contaminated_t_closest):

                                event_tt = (
                                previous_time.tt
                                + fraction * (ti.tt - previous_time.tt))

                                event_time = ts.tt_jd(event_tt)
                                is_sunlit = starlink.at(event_time).is_sunlit(eph)

                                sunlit.append(is_sunlit)
                                event_times.append(event_time)
                                logger.info('Contamination event for satellite %s at time %s: sunlit=%s', starlink.name, event_time.astimezone(zone), is_sunlit)
                                if is_sunlit:
                                    try:
                                        intensity, ab_magnitude = (
                                            satellite_tracker
                                            .brightness_calculation_at_event_time(
                                                starlink,
                                                rubin_obs,
                                                rubinobs_astr,
                                                sun,
                                                event_time,
                                            )
                                        )
                                        satellite_data[starlink.name].update({
                                            'intensity': intensity,
                                            'ab_magnitude': ab_magnitude
                                        })
                                        contaminated_intensities.append(intensity)
                                        contaminated_ab_magnitudes.append(ab_magnitude)
                                        logger.info('Brightness calculation for %s at event time %s: intensity=%s, ab_magnitude=%s', starlink.name, event_time.astimezone(zone), intensity, ab_magnitude)
                                    except Exception as e:
                                        exposure.increment_counters(no_brightness_skip=1)
                                        logger.info('Brightness calculation failed for %s: %s', starlink.name, e)
                                        satellite_data[starlink.name].update({
                                            'intensity': np.nan,
                                            'ab_magnitude': np.nan
                                        })
                                        contaminated_intensities.append(np.nan)
                                        contaminated_ab_magnitudes.append(np.nan)
                                else:
                                    satellite_data[starlink.name].update({
                                            'intensity': np.nan,
                                            'ab_magnitude': np.nan
                                        })
                                    contaminated_intensities.append(np.nan)
                                    contaminated_ab_magnitudes.append(np.nan)

                            event_mjd = float(
                                event_time.to_astropy().utc.mjd
                            )

                            event_tle_age_h = (
                                abs(
                                    event_mjd
                                    - tle_metadata[starlink]["tle_epoch_mjd"]
                                )
                                * 24.0
                            )

                            exposure.increment_counters(contaminations=1)
                            contamination_data = {'contaminated_targets': contaminated_targets,
                                                'contaminated_names': names,
                                                'separations': separations,
                                                'contaminated_t_closest': contaminated_t_closest,
                                                'event_times': event_times,
                                                'intensities': contaminated_intensities,
                                                'ab_magnitudes': contaminated_ab_magnitudes,
                                                'sunlit': sunlit,
                                                'contaminating_satellite': starlink.name,
                                                'contaminated_separations': contaminated_separations,
                                                'tle_age_event_h': event_tle_age_h,
                                                'tle_quality': tle_metadata[starlink]["tle_quality"],
                                                }
                            
                            logger.info('Contamination detected for satellite %s at time %s with %s targets while %s sunlit, tle_age: %s', starlink.name, ti.astimezone(zone), len(contaminated_targets), sum(sunlit), event_tle_age_h)
                        else:
                            logger.info('No contamination detected for satellite %s at time %s with separation minimum %s', starlink.name, ti.astimezone(zone), np.min(separations) if separations else 'N/A')
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
            
            
            satellite_tracker.clear_history()
        # After all time steps in the exposure are processed, compile and save exposure data
        target_tracker.print_cache_status()
        logger.info('Exposure summary:')
        exposure.print_summary()

        exposure_extensions = exposure.create_fits_extensions()
        for ext in exposure_extensions:
            hdu.append(ext)

        logger.info('Exposure HDUs appended for %s', exposure_counter)
        logger.info('Exposure processed in %.2f seconds', (-e_start + time.time()))
        if exposure_plot and exposure.total_events > 0:
            # More sophisticated version with proper spherical geometry

            def plot_fov_circle_polar(ax, center_alt, center_az, radius_deg):
                """Plot accurate field of view circle on polar plot"""
                center_alt_rad = np.radians(center_alt)
                center_az_rad = np.radians(center_az)
                radius_rad = np.radians(radius_deg)
                
                # Generate points around the circle using spherical geometry
                n_points = 100
                angles = np.linspace(0, 2*np.pi, n_points)
                
                circle_alt = []
                circle_az = []
                
                for angle in angles:
                    # Use spherical trigonometry to find points on circle
                    # This is approximate - for exact solution need full spherical geometry
                    delta_alt = radius_rad * np.cos(angle)
                    delta_az = radius_rad * np.sin(angle) / np.cos(center_alt_rad)
                    
                    alt = center_alt + np.degrees(delta_alt)
                    az = center_az + np.degrees(delta_az)
                    
                    # Handle wrapping
                    if az < 0:
                        az += 360
                    if az >= 360:
                        az -= 360
                        
                    if alt >= 0:  # Only plot visible part
                        circle_alt.append(alt)
                        circle_az.append(az)
                
                if circle_alt:
                    circle_theta = np.radians(circle_az)
                    circle_r = 90 - np.array(circle_alt)
                    ax.plot(circle_theta, circle_r, 'b-', linewidth=0.5, label='LSST FOV')

            # Use in your plotting code
            set_plot_style()
            if exposure_plot and exposure.total_events > 0:
                fig = plt.figure(figsize=(3.5,5))
                ax = plt.subplot(111, projection='polar')
                
                # Configure polar plot
                ax.set_theta_zero_location('N')
                ax.set_theta_direction(-1)       # Clockwise direction
                ax.set_rlabel_position(135)      # Move radial labels away from plotted data
                ax.set_ylim(0, 90)

                

                # Plot exposure targets
                filtered_targets_altaz_exp = SkyCoord(target_tracker.get_targets_at_time(exp[0])['skycoords'])
                target_theta_exp = filtered_targets_altaz_exp.az.radian
                target_r_exp = 90 - filtered_targets_altaz_exp.alt.degree
                ax.scatter(target_theta_exp, target_r_exp, s=2, c='orange', alpha=0.2, label='Exposure Targets')

                # Plot accurate FOV circle
                if hasattr(exposure, 'get_visit_region_skycoord_object'):
                    visit_center = exposure.get_visit_region_skycoord_object()
                    plot_fov_circle_polar(ax, visit_center.alt.degree, visit_center.az.degree, 1.75)
                    
                    # # Mark field center
                    # center_theta = visit_center.az.radian
                    # center_r = 90 - visit_center.alt.degree
                    # ax.scatter(center_theta, center_r, s=100, c='blue', marker='+', 
                    #         linewidth=3, label='Field Center')
                
                # Plot satellites with trails
                if exposure.all_satellite_positions:
                    sat_alt = [pos['alt'] for pos in exposure.all_satellite_positions]
                    sat_az = [pos['az'] for pos in exposure.all_satellite_positions]
                    sat_heights = [pos['height']/1000 for pos in exposure.all_satellite_positions]
                    
                    sat_theta = np.radians(sat_az)
                    sat_r = 90 - np.array(sat_alt)
                    ab_magnitudes = [pos['ab_magnitude'] for pos in exposure.all_satellite_positions]

                    
                    # Plot satellites colored by height
                    scatter = ax.scatter(sat_theta, sat_r, c=ab_magnitudes, s=50,
                                    cmap='plasma', alpha=0.8, edgecolors='black',
                                    linewidth=0.5, label='Satellites')
                    
                    # Colorbar
                    cbar = plt.colorbar(scatter, ax=ax, orientation='horizontal', pad=0.1, location='bottom')
                    cbar.set_label('AB Magnitude')
                    cbar.set_ticks(np.linspace(np.nanmin(ab_magnitudes), np.nanmax(ab_magnitudes), num=6))  # Set colorbar ticks
                    cbar.set_ticklabels([f'{tick:.1f}' for tick in cbar.get_ticks()])  # Format colorbar ticks
                    
                    ax.set_theta_zero_location('N')  # Set 0 degrees to North
                    ax.set_theta_direction(-1)       # Clockwise direction
                    ax.set_rlabel_position(135)      # Move radial labels away from plotted data
                
                # Plot contaminated targets
                if len(exposure.all_contaminated_targets) != 0:
                    for contaminated_tar in [ct['position'] for ct in exposure.all_contaminated_targets]:
                        cont_theta = contaminated_tar.az.radian
                        cont_r = 90 - contaminated_tar.alt.degree
                        ax.scatter(cont_theta, cont_r, c='red', s=50, marker='X', edgecolors='black', label='Contaminated Target')

                # Formatting
                # ax.set_title(f'Polar Sky View\n{exp[0].astimezone(zone).strftime("%Y-%m-%d %H:%M")} to ' + 
                #             f'{exp[-1].astimezone(zone).strftime("%H:%M")} Chile Time', 
                #             y=1.1, fontsize=14)
                

                
                # Format time on x-axis
                def format_time(x, pos=None):
                    total_seconds = (x / (2 * np.pi)) * 86400  # Convert radians to seconds
                    hours = int(total_seconds // 3600)
                    minutes = int((total_seconds % 3600) // 60)
                    return f'{hours:02d}:{minutes:02d}'
                plt.legend(bbox_to_anchor=(0.5, 1.0), ncols=2, markerscale=2)
                ax.xaxis.set_major_formatter(FuncFormatter(format_time))
                ax.set_yticks(ticks=np.radians([30, 60, 90]))  # Set y-ticks at specific altitudes
                ax.set_yticklabels(['30°', '60°', '90°'], fontdict={'horizontalalignment': 'right'})  # Label y-ticks
                ax.xaxis.set_tick_params(pad=5)
                
                plt.tight_layout()
                plt.savefig(pathToPlots + "/Satellite_visit_exposure_polar_advanced_" + 
                        exp[0].astimezone(zone).strftime('%Y-%m-%d_%H-%M-%S') + ".png", 
                        dpi=150, bbox_inches='tight')
                plt.close()


    



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
    logger.info('Current time by day end: %s Runtime: %sh %sm %.2fs', current_time, hours, minutes, seconds)
   #L.append(f"\n Current time: {current_time} Runtime: {hours}h {minutes}m {seconds:.2f}s , Day {day[0][0].astimezone(zone).day} completed ")
   #L.append('\n #####################################################################################')
    #print('Check point 5 : ',time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime(time.time())))



