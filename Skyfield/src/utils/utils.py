import logging
import sys, os
import datetime as dt
from pytz import timezone
from skyfield.api import load
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