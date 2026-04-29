import logging
import sys, os
import datetime as dt
from pytz import timezone
from astropy.io import fits
from astropy.table import Table
from astropy import units as u
from astropy.coordinates import SkyCoord, EarthLocation
import skyfield as sf
from skyfield.api import load, wgs84, EarthSatellite, N, W,S,E, Star, Angle
from skyfield import almanac
from astropy.table import Table, Column
import numpy as np
import pandas as pd
from astropy.time import Time, TimeDelta

pathToFiles = '/data/a.saricaoglu/repo/RubinsForge/Skyfield/files'

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

def read_targets_file(targets_file, verbose=True):
    """
    Read target star coordinates from a CSV file.

    Args:
        targets_file (str): Path to the CSV file containing target coordinates.
        region (str, optional): Region name to select specific targets. Defaults to None.
        verbose (bool, optional): Whether to print verbose output. Defaults to True.

    Returns:
        tuple: A tuple containing the targets dictionary and the coordinates SkyCoord object.
    """
    # Read the CSV file and examine its structure
    df = pd.read_csv(targets_file)

    if verbose:
        print("CSV file structure:")
        print(f"Shape: {df.shape}")
        print(f"Columns: {list(df.columns)}")
        print(f"\nFirst few rows:")
        print(df.head())
        print(f"\nData types:")
        print(df.dtypes)

    # Create targets dictionary with all requested fields
    targets = {}

    for index, row in df.iterrows():
        # Create SkyCoord object for each target
        coord = SkyCoord(ra=row['ra'] * u.degree, dec=row['dec'] * u.degree, frame='icrs')
        
        # Create target ID (you can modify this naming scheme)
        target_id = f"target_{index:04d}"
        
        # Build the target dictionary with all requested fields
        targets[target_id] = {
            'ra': row['ra'],  # degrees
            'dec': row['dec'],  # degrees
            'coord': coord,  # SkyCoord object for convenience
            
            # Initialize other fields - update these based on your CSV columns
            'name': row.get('name',None),  # Use CSV value if exists, else default name
            'flag': row.get('flag', None),
            'n_img': row.get('n_img', None),  # Use CSV value if exists, else None
            'image_conf': row.get('image_conf', None),
            'lens_type': row.get('lens_type', None),
            'source_type': row.get('source_type', None),
            'mag': row.get('mag', None),
            'Dmag': row.get('Dmag', None), 
            'value': row.get('value', None),
            'dvalue_min': row.get('dvalue_min', None),
            'dvalue_max': row.get('dvalue_max', None),
            
        }

    # Also create a SkyCoord array for convenience (as in original code)
    
    print(f"\nCreated targets dictionary with {len(targets)} targets")
    return targets

def get_targets_as_coordinates(pathToTargets, region=None):
    """
    Get target coordinates from a FITS file.

    Args:
        pathToTargets (str): Path to the FITS file containing target coordinates.


    Returns:        SkyCoord: A SkyCoord object containing the target coordinates.  
    """
    with fits.open(pathToTargets) as hdul:
        hdul.info()
        
        if region is not None:
            data = hdul[f'{region}'].data
        else:
            data = hdul[f'In_LSST_footprint'].data

        ra_i = data['RA']
        dec_i = data['DEC']

    coordinates = SkyCoord(ra=ra_i*u.degree, dec=dec_i*u.degree, frame='icrs')
    print(f"Loaded {len(coordinates)} target coordinates from {pathToTargets}")
    return coordinates

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

def save_twilight_windows_numpy(local_start, local_end, twilight_events, filename_base):
    """Save twilight windows to numpy format"""
    
    # Extract MJD times for efficient storage
    twilight_mjds = []
    twilight_types = []
    
    for event_group in twilight_events:
        for event in event_group:
            twilight_types.append(event[0])  # event type
            twilight_mjds.append(event[2])   # MJD time
    
    # Save as structured array
    twilight_array = np.array(list(zip(twilight_types, twilight_mjds)), 
                             dtype=[('type', 'U20'), ('mjd', 'f8')])
    
    np.savez_compressed(f"{filename_base}.npz", 
                       twilight_data=twilight_array,
                       metadata=np.array([local_start, local_end]))
    
    print(f"Twilight data saved to {filename_base}.npz")

def find_twilight_window_visit_overlaps(twilight_windows, visit_exposure_groups, extension_time=1.5, verbose=True):
    """Find overlaps between twilight windows and visit/exposure groups"""
    counter = 0
    for visit in visit_exposure_groups.values():

        visit['twilight_overlap'] = False
        visit_mjd = Time(visit['mjd'], format='mjd')
        visit_chile_time = visit['chile_time']
        twilight_extension = extension_time * u.hr

        # Check if visit falls within any twilight window
        for window in windows:
            start_mjd, end_mjd = window
            start_mjd = Time(start_mjd, format='mjd')
            end_mjd = Time(end_mjd, format='mjd')
            if start_mjd - twilight_extension <= visit_mjd <= end_mjd + twilight_extension:
                counter += 1
                visit['twilight_overlap'] = True
                if verbose:
                    print(f"Visit {visit['visit_id']} has a visit at Chile Time {visit_chile_time} " +
                    f"which falls within twilight window {Time(start_mjd, format='mjd').to_datetime(timezone=chile_tz)} to {Time(end_mjd, format='mjd').to_datetime(timezone=chile_tz)}")
                break  # No need to check other windows for this visit
    print(f"Total twilight-overlapping visits on all days: {counter}")

    return visit_exposure_groups

def load_twilight_windows(filename, verbose=True):
    """Load twilight windows from numpy format"""
    data = np.load(filename)
    twilight_data = data['twilight_data']

    print(f"Loaded {len(twilight_data)} twilight events")
    windows = []
    window = []
    for i, event in enumerate(twilight_data):
        event_type = event['type']
        event_mjd = event['mjd']
        if len(window) != 2:
            window.append(event_mjd)
        if len(window) == 2:
            windows.append(window)
            window = []
        if verbose:    
            print(f"  {i+1}: {event_type} at MJD {event_mjd:.6f}")

    return windows

def create_twilight_windows(obs_start, obs_end):

    ts = load.timescale()
    zone_stl = timezone('Etc/GMT-6')
    zone = timezone('Chile/Continental')
    utc = timezone('UTC')
    local_start = ts.from_datetime(zone.localize(obs_start).astimezone(utc))
    local_end = ts.from_datetime(zone.localize(obs_end).astimezone(utc))
    eph = load('de421.bsp')
    earth = eph['earth']
    sun = eph['sun']
    rubin_obs = wgs84.latlon(-30.244633,  -70.749417)

    today = local_start
    next_day = local_end
    f = sf.almanac.dark_twilight_day(eph, rubin_obs)
    times, events = almanac.find_discrete(today, next_day, f)
    previous_e = f(local_start).item()
    twilight_events = []
    for t, e in zip(times, events):
        tstr = str(t.astimezone(zone))[:16]
        mjd_time = Time(t.ut1, format='jd').mjd
        twilight=[]
        if previous_e < e:
            print(f'{tstr} mjd {mjd_time} {almanac.TWILIGHTS[e]} starts')
            if str(almanac.TWILIGHTS[e]) == 'Day':
                day_start = t # Twilight ends for morning
                twilight.append(('am twilight end', tstr,mjd_time))
            if str(almanac.TWILIGHTS[e]) == 'Astronomical twilight':
                night_end = t # Twilight starts for morning
                twilight.append(('am twilight start', tstr,mjd_time))
        else:
            print(f'{tstr} mjd {mjd_time} {almanac.TWILIGHTS[previous_e]} ends')
            if (almanac.TWILIGHTS[previous_e]) == 'Day':
                day_end = t # Twilight starts for evening
                twilight.append(('pm twilight start', tstr,mjd_time))
            if (almanac.TWILIGHTS[previous_e]) == 'Astronomical twilight':
                night_start = t # Twilight ends for evening
                twilight.append(('pm twilight end', tstr,mjd_time))
        if twilight:
            twilight_events.append(twilight)
        previous_e = e

    # Save the data
    numpy_filename = f"{pathToFiles}/twilight_windows/{local_start.astimezone(zone).strftime('%Y%m%dT%H%M')}_to_{local_end.astimezone(zone).strftime('%Y%m%dT%H%M')}"
    save_twilight_windows_numpy(local_start, local_end, twilight_events, numpy_filename)

    return f'{numpy_filename}.npz'

def create_exposure_windows(visits_in_twilight):
    # Prepare data for FITS table
    fits_rows = []
    local_start = Time(visits_in_twilight[0]['mjd'], format='mjd').to_datetime(timezone=chile_tz)
    local_end = Time(visits_in_twilight[-1]['mjd'], format='mjd').to_datetime(timezone=chile_tz)
    for visit in visits_in_twilight:
        for target in visit['targets']:
            fits_rows.append({
                'visit_id': visit['visit_id'],
                'mjd': visit['mjd'],
                'chile_time': visit['chile_time'].strftime('%Y-%m-%d %H:%M:%S'),
                'exposure_time': visit['exposure_time'],
                'num_exposures': visit['num_exposures'],
                'filter': visit['filter'],
                'field_ra': visit['field_ra'],
                'field_dec': visit['field_dec'],
                'target_id': target['target_id'],
                'target_ra': target['target_ra'],
                'target_dec': target['target_dec'],
                'five_sigma_depth': target['five_sigma_depth'],
                'seeing': target['seeing'],
                'sky_brightness': target['sky_brightness']
            })

    fits_table = Table(rows=fits_rows)
    cols = [
        fits.Column(name='visit_id', format='K', array=fits_table['visit_id']),
        fits.Column(name='mjd', format='D', array=fits_table['mjd']),
        fits.Column(name='chile_time', format='20A', array=fits_table['chile_time']),
        fits.Column(name='exposure_time', format='E', array=fits_table['exposure_time']),
        fits.Column(name='num_exposures', format='I',   array=fits_table['num_exposures']), 
        fits.Column(name='filter', format='1A', array=fits_table['filter']),
        fits.Column(name='field_ra', format='D', array=fits_table['field_ra']),
        fits.Column(name='field_dec', format='D', array=fits_table['field_dec']),
        fits.Column(name='target_id', format='20A', array=fits_table['target_id']),
        fits.Column(name='target_ra', format='D', array=fits_table['target_ra']),
        fits.Column(name='target_dec', format='D', array=fits_table['target_dec']),
        fits.Column(name='five_sigma_depth', format='E', array=fits_table['five_sigma_depth']),
        fits.Column(name='seeing', format='E', array=fits_table['seeing']),
        fits.Column(name='sky_brightness', format='E', array=fits_table['sky_brightness'])
    ]
    hdu = fits.BinTableHDU.from_columns(cols)
    # Write to FITS file
    fits_filename = f"{pathToFiles}/twilight_visits/{local_start.strftime('%Y%m%d')}_to_{local_end.strftime('%Y%m%d')}.fits"
    hdu.writeto(fits_filename, overwrite=True)
    print(f"Written {len(fits_rows)} rows to {fits_filename}")
    print(f"FITS columns: {[col.name for col in hdu.columns]}")

    return fits_filename

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

def get_targets_by_visit_exposure(bundle_list, targets_dict, start_mjd, end_mjd):
    """Group targets by visit ID and exposure parameters"""
    
    visit_groups = {}
    target_ids = list(targets_dict.keys())

    # Process each target's observations
    for i, target_observations in enumerate(bundle_list[0].metric_values):
                    # Filter observations for the date range
        range_mask = (target_observations['observationStartMJD'] >= start_mjd) & \
                        (target_observations['observationStartMJD'] <= end_mjd)
            
        range_observations = target_observations[range_mask]
        if i >= len(target_ids):
            break
            
        target_id = target_ids[i]

        for obs in range_observations:
            visit_id = obs['observationId']
            exposure_time = obs['visitExposureTime']
            num_exposures = obs['numExposures']
            mjd = obs['observationStartMJD']
            
            # Create unique key for visit/exposure combination
            visit_key = f"{visit_id}"
            obs_time = Time(obs['observationStartMJD'], format='mjd')
            utc_datetime = obs_time.to_datetime(timezone=dt.timezone.utc)
            chile_datetime = utc_datetime.astimezone(chile_tz)
            if visit_key not in visit_groups:
                visit_groups[visit_key] = {
                    'visit_id': visit_id,
                    'mjd': mjd,
                    'chile_time': chile_datetime,
                    'exposure_time': exposure_time,
                    'num_exposures': num_exposures,
                    'filter': obs['filter'],
                    'field_ra': obs['fieldRA'],
                    'field_dec': obs['fieldDec'],
                    'targets': []
                }
            
            # Add target to this visit/exposure group
            visit_groups[visit_key]['targets'].append({
                'target_id': target_id,
                'target_ra': targets_dict[target_id]['ra'],
                'target_dec': targets_dict[target_id]['dec'],
                'five_sigma_depth': obs['fiveSigmaDepth'],
                'seeing': obs['seeingFwhmEff'],
                'sky_brightness': obs['skyBrightness']
            })
    
    return visit_groups    


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
    
def fits_sanitize_table(tab):
    """
    Sanitize an Astropy Table for FITS output.
    """
    out = Table()
    for name in tab.colnames:
        col = tab[name]
        # Skip columns without dtype (e.g., astropy Time objects)
        if not hasattr(col, 'dtype'):
            out[name] = col
            continue
        # Scalar strings -> fixed-width ASCII
        if (col.dtype.kind in ('U', 'S')) and col.ndim == 1:
            maxlen = max((len(x) for x in col.astype(str)), default=1)
            out[name] = np.asarray(col.astype(str), dtype=f'S{maxlen}')
            # Object dtype likely holds lists/arrays
        elif col.dtype == object:
            sample = next((v for v in col if v is not None), None)
            if isinstance(sample, (list, tuple, np.ndarray)):
                # If elements are strings -> join
                if len(sample) and all(isinstance(x, str) for x in sample):
                    joined = [','.join(map(str, v)) if v is not None else '' for v in col]
                    maxlen = max((len(s) for s in joined), default=1)
                    out[name] = np.asarray(joined, dtype=f'S{maxlen}')
                else:
                    # Assume numeric arrays: leave as object of ndarrays so astropy makes PE()/PJ()
                    out[name] = Column([None if v is None else np.asarray(v) for v in col], name=name)
            else:
                out[name] = col
        else:
            out[name] = col
    return out

def find_satellites_up(satellites, observer, start_time, end_time, altitude_degrees=30.0, check_sunlit=True, eph=None, return_periods=False):
    """
    Find satellites that are above horizon and optionally sunlit during time period
    
    Parameters:
    -----------
    satellites : list
        List of Skyfield satellite objects
    observer : Skyfield observer position
        Observer position (e.g., rubin_obs)
    start_time : Skyfield time
        Start time for checking
    end_time : Skyfield time  
        End time for checking
    altitude_degrees : float
        Minimum altitude in degrees (default 30.0)
    check_sunlit : bool
        Whether to check if satellite is sunlit (default True)
    eph : Skyfield ephemeris
        Ephemeris for sunlit check (required if check_sunlit=True)
    return_periods : bool
        If True, return detailed period information for each satellite
        
    Returns:
    --------
    If return_periods=False:
        list : Satellites that meet the criteria
    If return_periods=True:
        tuple : (satellites_up, satellite_periods)
            satellites_up: list of satellites
            satellite_periods: dict with satellite -> period info
    """
    satellites_up = []
    satellite_periods = {}

    for satellite in satellites:
        try:
            times, events = satellite.find_events(
                observer,
                start_time,
                end_time,
                altitude_degrees=altitude_degrees
            )

            difference = satellite - observer
            topocentric = difference.at(start_time)
            alt, az, distance = topocentric.altaz()
            above_horizon_start = alt.degrees > altitude_degrees

            visible_periods = []
            current_period_start = start_time if above_horizon_start else None

            for event_time, event in zip(times, events):
                if event == 0:
                    current_period_start = event_time
                elif event == 2 and current_period_start is not None:
                    visible_periods.append((current_period_start, event_time))
                    current_period_start = None

            if current_period_start is not None:
                visible_periods.append((current_period_start, end_time))

            if not visible_periods:
                continue

            sunlit_periods = []
            if check_sunlit and eph is not None:
                for vis_start, vis_end in visible_periods:
                    duration = (vis_end.utc_datetime() - vis_start.utc_datetime()).total_seconds()
                    num_samples = max(1, int(duration))
                    sample_times = [
                        vis_start + (vis_end - vis_start) * i / num_samples
                        for i in range(num_samples + 1)
                    ]

                    sunlit_start = None
                    for sample_time in sample_times:
                        is_sunlit = satellite.at(sample_time).is_sunlit(eph)
                        if is_sunlit and sunlit_start is None:
                            sunlit_start = sample_time
                        elif (not is_sunlit) and sunlit_start is not None:
                            sunlit_periods.append((sunlit_start, sample_time))
                            sunlit_start = None

                    if sunlit_start is not None:
                        sunlit_periods.append((sunlit_start, vis_end))
            else:
                sunlit_periods = visible_periods

            combined_periods = sunlit_periods if (check_sunlit and eph is not None) else visible_periods
            if not combined_periods:
                continue

            satellites_up.append(satellite)

            if return_periods:
                satellite_periods[satellite] = {
                    'visible_periods': visible_periods,
                    'sunlit_periods': sunlit_periods,
                    'combined_periods': combined_periods,
                    'total_visible_duration': sum([
                        (end.utc_datetime() - start.utc_datetime()).total_seconds()
                        for start, end in visible_periods
                    ]),
                    'total_sunlit_duration': sum([
                        (end.utc_datetime() - start.utc_datetime()).total_seconds()
                        for start, end in sunlit_periods
                    ]) if sunlit_periods else 0
                }

        except Exception as e:
            print(f"Error checking satellite {satellite.name}: {e}")
            continue
    
    if return_periods:
        return satellites_up, satellite_periods
    else:
        return satellites_up

def find_satellites_up_simple(satellites, observer, check_time, altitude_degrees=30.0, check_sunlit=True, eph=None):
    """
    Simplified version - check satellites at a single time point
    
    Parameters:
    -----------
    satellites : list
        List of Skyfield satellite objects
    observer : Skyfield observer position
    check_time : Skyfield time
        Time to check satellite positions
    altitude_degrees : float
        Minimum altitude in degrees
    check_sunlit : bool
        Whether to check if satellite is sunlit
    eph : Skyfield ephemeris
        Ephemeris for sunlit check
        
    Returns:
    --------
    list : Satellites above horizon (and sunlit if requested)
    """
    satellites_up = []
    
    for satellite in satellites:
        try:
            # Check altitude
            difference = satellite - observer
            topocentric = difference.at(check_time)
            alt, az, distance = topocentric.altaz()
            
            if alt.degrees > altitude_degrees:
                # Check if sunlit (if requested)
                if check_sunlit and eph is not None:
                    if satellite.at(check_time).is_sunlit(eph):
                        satellites_up.append(satellite)
                else:
                    satellites_up.append(satellite)
                    
        except Exception as e:
            print(f"Error checking satellite {satellite.name}: {e}")
            continue
    
    return satellites_up

def get_satellite_count_info(satellites, observer, check_time, altitude_degrees=30.0, eph=None):
    """
    Get detailed count information about satellites
    
    Returns:
    --------
    dict : Satellite counts and statistics
    """
    total_satellites = len(satellites)
    above_horizon = 0
    sunlit = 0
    above_horizon_and_sunlit = 0
    
    for satellite in satellites:
        try:
            # Check altitude
            difference = satellite - observer
            topocentric = difference.at(check_time)
            alt, az, distance = topocentric.altaz()
            
            is_above_horizon = alt.degrees > altitude_degrees
            is_sunlit = satellite.at(check_time).is_sunlit(eph) if eph else False
            
            if is_above_horizon:
                above_horizon += 1
                
            if is_sunlit:
                sunlit += 1
                
            if is_above_horizon and is_sunlit:
                above_horizon_and_sunlit += 1
                
        except Exception as e:
            continue
    
    return {
        'total': total_satellites,
        'above_horizon': above_horizon,
        'sunlit': sunlit,
        'above_horizon_and_sunlit': above_horizon_and_sunlit,
        'time': check_time
    }

def is_time_in_periods(check_time, periods):
    """
    Check if a time falls within any of the given periods
    
    Parameters:
    -----------
    check_time : Skyfield time
        Time to check
    periods : list of tuples
        List of (start_time, end_time) tuples
        
    Returns:
    --------
    bool : True if time falls within any period
    """
    check_datetime = check_time.utc_datetime()
    
    for start_time, end_time in periods:
        start_datetime = start_time.utc_datetime()
        end_datetime = end_time.utc_datetime()
        
        if start_datetime <= check_datetime <= end_datetime:
            return True
    
    return False


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

