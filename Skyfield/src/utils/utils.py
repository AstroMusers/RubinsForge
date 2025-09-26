import logging
import sys, os
import datetime as dt
from pytz import timezone
from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord, EarthLocation
import skyfield as sf
from skyfield.api import load, wgs84, EarthSatellite, N, W,S,E, Star, Angle
from skyfield import almanac

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

def find_closest_TLE(obs_date,tle_data_path, zone):
    """
    Find the closest TLE (Two-Line Element) set to a specific date.

    Args:
        obs_date (datetime): The start time of the observation.
        tle_data (list): A list of TLE data strings.

    Returns:
        str: The closest TLE string or None if not found.
    """
    closest_tle = None
    closest_time_diff = 30

    tle_data = os.listdir(tle_data_path)
    ts = load.timescale()
    for tle in tle_data:
        tle_time = parse_tle_time(tle)
        if tle_time:
            tle_time = ts.from_datetime(tle_time.astimezone(zone))
            time_diff = abs(tle_time - obs_date)
            
            if time_diff < closest_time_diff:
                closest_time_diff = time_diff
                closest_tle = tle
    print(f'Closest TLE time: {closest_tle}, Obs time: {obs_date.astimezone(zone).strftime("%Y-%m-%d %H:%M")}, Time diff: {closest_time_diff}')
    return closest_tle, load.tle_file(os.path.join(tle_data_path, closest_tle))

def find_closest_TLE_lines(obs_date,tle_data_path, zone):
    """
    Find the closest TLE (Two-Line Element) set to a specific date.

    Args:
        obs_date (datetime): The start time of the observation.
        tle_data (list): A list of TLE data strings.

    Returns:
        str: The closest TLE string or None if not found.
    """
    closest_tle = None
    closest_time_diff = 30

    tle_data = os.listdir(tle_data_path)
    ts = load.timescale()
    for tle in tle_data:
        tle_time = parse_tle_time(tle)
        if tle_time:
            tle_time = ts.from_datetime(tle_time.astimezone(zone))
            time_diff = abs(tle_time - obs_date)
            
            if time_diff < closest_time_diff:
                closest_time_diff = time_diff
                closest_tle = tle
    print(f'Closest TLE time: {closest_tle}, Obs time: {obs_date.astimezone(zone).strftime("%Y-%m-%d %H:%M")}, Time diff: {closest_time_diff}')
    """Read a TLE file and return (name, l1, l2) triplets."""
    tle_triplets = []
    path = os.path.join(tle_data_path, closest_tle)
    with open(path, 'r') as f:
        lines = f.read().strip().splitlines()
    
    # Every 3 lines = 1 satellite (name, line1, line2)
    for i in range(0, len(lines), 3):
        name = lines[i].strip()
        l1 = lines[i+1].strip()
        l2 = lines[i+2].strip()
        tle_triplets.append((name, l1, l2))

    return closest_tle, tle_triplets

def parse_tle_time(tle):
    """
    Parse the time from a TLE (Two-Line Element) string.

    Args:
        tle (str): The TLE string.

    Returns:
        datetime: The parsed time or None if parsing failed.
    """
    try:
        lines = tle.strip('.txt').split('_')
        # print(lines)


        # Extract the epoch time from the second line
        epoch_str = lines[1]
        dt.format = "%m.%d.%y"
        epoch_time = dt.datetime.strptime(epoch_str, dt.format)
        return epoch_time
    except Exception as e:
        print(f"Error parsing TLE time: {e}")
        return None

def get_time_array(local_start, local_end, eph, rubin_obs, zone, twilight_extension, deltatime, deltaday, postexposure_skip):
    day_array = []
    today = local_start
    f = sf.almanac.dark_twilight_day(eph, rubin_obs)
    while today < local_end:
        time_array = []
        print(f'starting time loop at {today.astimezone(zone)}')
        next_day = today + dt.timedelta(days=1)
        times, events = almanac.find_discrete(today, next_day, f)
        previous_e = f(local_start).item()

        for t, e in zip(times, events):
            tstr = str(t.astimezone(zone))[:16]
            if previous_e < e:
                print(f'{tstr} {almanac.TWILIGHTS[e]} starts')
                if str(almanac.TWILIGHTS[e]) == 'Day':
                    day_start = t # Twilight ends for morning
                if str(almanac.TWILIGHTS[e]) == 'Astronomical twilight':
                    night_end = t # Twilight starts for morning
            else:
                print(f'{tstr} {almanac.TWILIGHTS[previous_e]} ends')
                if (almanac.TWILIGHTS[previous_e]) == 'Day':
                    day_end = t # Twilight starts for evening
                if (almanac.TWILIGHTS[previous_e]) == 'Astronomical twilight':
                    night_start = t # Twilight ends for evening
            previous_e = e

        ti = night_end - twilight_extension
        print(f'starting time array loop at {ti.astimezone(zone)}')
        exposure_array = []
        while (ti < day_start) :
            exposure_array.append(ti)
            ti += deltatime
            if len(exposure_array) == 15:
                time_array.append(exposure_array)
                ti += postexposure_skip - dt.timedelta(seconds=15)
                exposure_array = []

        ti = day_end + deltatime

        exposure_array = []
        while (ti > day_end) & (ti < night_start + twilight_extension):
            exposure_array.append(ti)
            ti += deltatime
            if len(exposure_array) == 15:
                time_array.append(exposure_array)
                ti += (postexposure_skip - dt.timedelta(seconds=15))
                exposure_array = []


        print(f'ending time loop at {ti.astimezone(zone)}')
        today  += deltaday
        
        day_array.append(time_array)

    return day_array

def get_targets_as_stars(pathToTargets, region=None):
    """
    Get target stars from a FITS file.

    Args:
        pathToTargets (str): Path to the FITS file containing target coordinates.

    Returns:
        list: A list of Star objects representing the targets.
    """
    targets = []
    names = []
    with fits.open(pathToTargets) as hdul:
        hdul.info()
        
        if region is not None:
            data = hdul[f'{region}'].data
        else:
            data = hdul[f'In_LSST_footprint'].data

        ra_i = data['RA']
        dec_i = data['DEC']
        name_i = data['Target_no']

        for i in range(0, len(data)):

            target_i = Star(Angle(degrees=ra_i[i]), Angle(degrees=dec_i[i]))
            targets.append(target_i) 
            names.append(name_i[i])

    return targets, names

def get_exposure_windows(exposures_file):
    """
    Get exposure windows from a FITS file.

    Args:
        exposures_file (str): Path to the FITS file containing exposure windows.

    Returns:
        list: A list of exposure windows.
    """
    import numpy as np
    from astropy.time import Time, TimeDelta
    from pytz import timezone   
    from skyfield.api import load
    chile_tz = timezone('Chile/Continental')
    ts = load.timescale()

    with fits.open(exposures_file) as hdul:
        # print(hdul.info())
        data = hdul[1].data
        print(len(data['visit_id']))
        print(len(np.unique(data['visit_id'])))
        
        exposure_dict = {}
        for visit in np.unique(data['visit_id']):
            visit_mask = data['visit_id'] == visit
            visit_data = data[visit_mask]
            day = Time(visit_data['mjd'][0], format='mjd').to_datetime(timezone=chile_tz).date()

            if day not in exposure_dict:
                exposure_dict[day] = {}
            start = Time(visit_data['mjd'][0], format='mjd').to_datetime(timezone=chile_tz)
            start = ts.from_datetime(start)            
            window = [start + dt.timedelta(seconds=i) for i in range(0,int(visit_data['exposure_time'][0]) +1)]
            exposure_dict[day][visit] = {
                'visit_id': visit,
                'mjd': visit_data['mjd'][0],
                'chile_time': visit_data['chile_time'][0],
                'exposure_time': visit_data['exposure_time'][0],# total duration of the visit
                'num_exposures': visit_data['num_exposures'][0], # number of exposures that adds up to exposure time
                'filter': visit_data['filter'][0],
                'field_ra': visit_data['field_ra'][0],
                'field_dec': visit_data['field_dec'][0],
                'seeing': visit_data['seeing'][0],
                'sky_brightness': visit_data['sky_brightness'][0],
                'target_ids': visit_data['target_id'],
                'target_ras': visit_data['target_ra'],
                'target_decs': visit_data['target_dec'],
                'visit_center_as_star': Star(Angle(degrees=visit_data['field_ra'][0]), Angle(degrees=visit_data['field_dec'][0])),
                'window': window,
            }

            # print(f"\nVisit ID: {visit}")
        return exposure_dict
    
def exposure_targets_to_stars(exposures_file):
    import numpy as np
    from astropy.time import Time, TimeDelta
    from pytz import timezone   
    from skyfield.api import load
    chile_tz = timezone('Chile/Continental')
    ts = load.timescale()

    with fits.open(exposures_file) as hdul:
        # print(hdul.info())
        data = hdul[1].data
        print(len(data['target_id']))
        print(len(np.unique(data['target_id'])))
        
        target_star = []
        target_name = []
        for target in np.unique(data['target_id']):
            target_mask = data['target_id'] == target
            target_data = data[target_mask]

            target_star.append(Star(Angle(degrees=target_data['target_ra'][0]), Angle(degrees=target_data['target_dec'][0])))
            target_name.append(f"Target {target}")

        return target_star, target_name

# def get_daily_twilight_periods_mjd(start_date, end_date, location=None, eph=None):
#     """
#     Get twilight periods for each day in MJD format
    
#     Parameters:
#     -----------
#     start_date : datetime
#         Start date for the period
#     end_date : datetime  
#         End date for the period
#     location : skyfield location object
#         Observatory location (e.g., rubin_obs)
#     eph : skyfield ephemeris
#         Ephemeris object for sun calculations
        
#     Returns:
#     --------
#     list of dict : Each dict contains twilight times for one day in MJD
#         Keys: 'date', 'night_start_mjd', 'night_end_mjd', 'day_start_mjd', 'day_end_mjd'
#     """
#     import skyfield as sf
#     from skyfield import almanac
#     from skyfield.api import load, wgs84, S, W
#     from pytz import timezone
#     from astropy.time import Time

#     if location is None or eph is None:
#         location = wgs84.latlon(-30.244633,  -70.749417) # Default to Rubin Observatory
#         eph = load('de421.bsp')  # Default ephemeris
    
#     ts = load.timescale()
#     twilight_periods = []
    
#     # Convert to skyfield times
#     current_date = start_date
    
#     while current_date <= end_date:
#         # Define day boundaries (from midnight to midnight to capture full night)
#         day_start_dt = current_date.replace(hour=0, minute=0, second=0, microsecond=0)
#         day_end_dt = (current_date + dt.timedelta(days=1)).replace(hour=0, minute=0, second=0, microsecond=0)
        
#         # Convert to skyfield times (assuming UTC input)
#         if hasattr(current_date, 'tzinfo') and current_date.tzinfo is not None:
#             sf_start = ts.from_datetime(day_start_dt.astimezone(timezone('UTC')))
#             sf_end = ts.from_datetime(day_end_dt.astimezone(timezone('UTC')))
#         else:
#             # Assume UTC if no timezone info
#             sf_start = ts.from_datetime(day_start_dt.replace(tzinfo=timezone('UTC')))
#             sf_end = ts.from_datetime(day_end_dt.replace(tzinfo=timezone('UTC')))
        
#         # Get twilight function
#         f = sf.almanac.dark_twilight_day(eph, location)
#         times, events = almanac.find_discrete(sf_start, sf_end, f)
        
#         # Initialize twilight times for this day
#         day_twilights = {
#             'date': current_date.strftime('%Y-%m-%d'),
#             'night_start_mjd': None,  # Evening astronomical twilight end
#             'night_end_mjd': None,    # Morning astronomical twilight start  
#             'day_start_mjd': None,    # Morning day start (end of twilight)
#             'day_end_mjd': None       # Evening day end (start of twilight)
#         }
        
#         # Process twilight events
#         previous_e = f(sf_start).item()
        
#         for t, e in zip(times, events):
#             # Convert skyfield time to MJD
#             mjd_time = Time(t.ut1, format='jd').mjd
            
#             if previous_e < e:
#                 # Transition to higher state (e.g., night -> twilight -> day)
#                 if str(almanac.TWILIGHTS[e]) == 'Day':
#                     day_twilights['day_start_mjd'] = mjd_time
#                 elif str(almanac.TWILIGHTS[e]) == 'Astronomical twilight':
#                     day_twilights['night_end_mjd'] = mjd_time
                    
#             else:
#                 # Transition to lower state (e.g., day -> twilight -> night)
#                 if str(almanac.TWILIGHTS[previous_e]) == 'Day':
#                     day_twilights['day_end_mjd'] = mjd_time
#                 elif str(almanac.TWILIGHTS[previous_e]) == 'Astronomical twilight':
#                     day_twilights['night_start_mjd'] = mjd_time
                    
#             previous_e = e
        
#         twilight_periods.append(day_twilights)
#         current_date += dt.timedelta(days=1)
    
#     return twilight_periods



# def get_observing_windows_mjd(start_date, end_date, location, eph):
#     """
#     Get nightly observing windows (astronomical twilight periods) in MJD
    
#     Returns:
#     --------
#     list of tuples : (night_start_mjd, night_end_mjd) for each night

#     """
#     from skyfield.api import load, wgs84
#     if location is None or eph is None:
#         location = wgs84.latlon(-30.244633,  -70.749417) # Default to Rubin Observatory
#         eph = load('de421.bsp')  # Default ephemeris

#     twilight_data = get_daily_twilight_periods_mjd(start_date, end_date, location, eph)
    
#     observing_windows = []
#     for day_data in twilight_data:
#         if day_data['night_start_mjd'] and day_data['night_end_mjd']:
#             observing_windows.append((day_data['night_start_mjd'], day_data['night_end_mjd']))
    
#     return observing_windows

