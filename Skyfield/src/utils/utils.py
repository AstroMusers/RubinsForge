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
                    day_start = t
                if str(almanac.TWILIGHTS[e]) == 'Astronomical twilight':
                    night_end = t
            else:
                print(f'{tstr} {almanac.TWILIGHTS[previous_e]} ends')
                if (almanac.TWILIGHTS[previous_e]) == 'Day':
                    day_end = t
                if (almanac.TWILIGHTS[previous_e]) == 'Astronomical twilight':
                    night_start = t
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
        while (ti > day_end) & (ti < day_end + twilight_extension):
            exposure_array.append(ti)
            ti += deltatime
            if len(exposure_array) == 15:
                time_array.append(exposure_array)
                ti += (postexposure_skip - dt.timedelta(seconds=15))
                exposure_array = []

        ti = day_end + deltatime

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