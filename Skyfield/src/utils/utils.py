import logging
import sys, os
import datetime as dt
from pathlib import Path
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
from logging_stream import StreamToLogger

_CHILE_TIMEZONE = timezone('Chile/Continental')

pathToFiles = str(Path(__file__).resolve().parents[2] / 'files')


def _display_path(path):
    resolved_path = Path(path).resolve()
    repository_path = Path(__file__).resolve().parents[3]
    try:
        return str(resolved_path.relative_to(repository_path))
    except ValueError:
        return resolved_path.name


def _find_closest_tle_path(obs_date, tle_data_path, zone):
    """Return the nearest dated TLE filename and path, strictly within 30 days."""
    closest_tle = None
    closest_time_diff = 30  # [day]

    tle_data = os.listdir(tle_data_path)
    ts = load.timescale()
    for tle in tle_data:
        if not os.path.isfile(os.path.join(tle_data_path, tle)):
            continue
        tle_time = parse_tle_time(tle)
        if tle_time:
            tle_time = ts.from_datetime(tle_time.astimezone(zone))
            time_diff = abs(tle_time - obs_date)  # [day]
            
            if time_diff < closest_time_diff:
                closest_time_diff = time_diff
                closest_tle = tle
    if closest_tle is None:
        raise FileNotFoundError(f"No dated TLE file less than 30 days from the observation in {tle_data_path}")
    print(f'Closest TLE time: {closest_tle}, Obs time: {obs_date.astimezone(zone).strftime("%Y-%m-%d %H:%M")}, Time diff: {closest_time_diff}')
    return closest_tle, os.path.join(tle_data_path, closest_tle)


def find_closest_TLE(obs_date,tle_data_path, zone):
    """Return the nearest dated filename and its Skyfield satellite objects."""
    closest_tle, path = _find_closest_tle_path(obs_date, tle_data_path, zone)
    print(f'Reading from {path}...')
    return closest_tle, load.tle_file(path)

def find_closest_TLE_lines(obs_date,tle_data_path, zone):
    """Return the nearest dated filename and its (name, line1, line2) records."""
    closest_tle, path = _find_closest_tle_path(obs_date, tle_data_path, zone)
    tle_triplets = []
    print(f'Reading from {path}...')
    with open(path, 'r') as f:
        lines = f.read().strip().splitlines()
    if len(lines) % 3:
        raise ValueError(f"Expected complete three-line TLE records in {path}; found {len(lines)} lines")
    
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

def _sample_exposure_window(start, stop, sample_step, start_cadence):
    """Return complete 15-sample exposures and the final sampling cursor."""
    if sample_step <= dt.timedelta(0) or start_cadence <= dt.timedelta(0):
        raise ValueError("Sample step and exposure-start cadence must be positive")
    sample_count = 15
    exposures = []
    current_exposure = []
    current_time = start
    while current_time < stop:
        current_exposure.append(current_time)
        current_time += sample_step
        if len(current_exposure) == sample_count:
            exposures.append(current_exposure)
            current_time += start_cadence - sample_count * sample_step
            current_exposure = []
    return exposures, current_time


def get_time_array(local_start, local_end, eph, rubin_obs, zone, twilight_extension, deltatime, deltaday, postexposure_skip):
    """Build complete 15-sample dawn/dusk blocks within [local_start, local_end).

    Exposure cadence is start-to-start. Each processed day must contain the
    required dawn, day-start, day-end, and dusk transitions or ValueError is raised.
    """
    if deltaday <= dt.timedelta(0):
        raise ValueError("Day step must be positive")
    day_array = []
    today = local_start
    f = sf.almanac.dark_twilight_day(eph, rubin_obs)
    while today < local_end:
        time_array = []
        print(f'starting time loop at {today.astimezone(zone)}')
        next_day = today + dt.timedelta(days=1)
        times, events = almanac.find_discrete(today, next_day, f)
        previous_e = f(today).item()
        night_end = day_start = day_end = night_start = None

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

        boundaries = {
            'astronomical dawn': night_end, 'day start': day_start,
            'day end': day_end, 'astronomical dusk': night_start,
        }
        missing = [name for name, value in boundaries.items() if value is None]
        if missing:
            date = today.astimezone(zone).strftime('%Y-%m-%d')
            raise ValueError(f"Missing required twilight transitions on {date}: {', '.join(missing)}")
        ti = max(night_end - twilight_extension, local_start)
        print(f'starting time array loop at {ti.astimezone(zone)}')
        morning, _ = _sample_exposure_window(ti, min(day_start, local_end), deltatime, postexposure_skip)
        evening, ti = _sample_exposure_window(
            max(day_end + deltatime, local_start), min(night_start + twilight_extension, local_end),
            deltatime, postexposure_skip,
        )
        time_array.extend(morning)
        time_array.extend(evening)


        print(f'ending time loop at {ti.astimezone(zone)}')
        today  += deltaday
        
        day_array.append(time_array)

    return day_array

def read_targets_file(targets_file, verbose=True):
    """
    Read target star coordinates from a CSV file.

    Args:
        targets_file (str): Path to the CSV file containing target coordinates.
        verbose (bool, optional): Whether to print verbose output. Defaults to True.

    Returns:
        dict: Targets keyed by row ID, including an ICRS SkyCoord for each row.
    """
    # Read the CSV file and examine its structure
    print(f'Reading from {targets_file}...')
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

    print(f"\nCreated targets dictionary with {len(targets)} targets")
    return targets

def _read_target_catalog_columns(path, region, columns):
    """Copy requested columns from the selected target-catalog FITS extension."""
    print(f'Reading from {path}...')
    with fits.open(path) as hdul:
        hdul.info()
        extension = 'In_LSST_footprint' if region is None else f'{region}'
        data = hdul[extension].data
        return tuple(data[column].copy() for column in columns)


def get_targets_as_coordinates(pathToTargets, region=None):
    """
    Get target coordinates from a FITS file.

    Args:
        pathToTargets (str): Path to the FITS file containing target coordinates.
        region (str, optional): FITS extension name. Defaults to In_LSST_footprint.

    Returns:
        SkyCoord: ICRS target coordinates in catalog row order.
    """
    ra_i, dec_i = _read_target_catalog_columns(pathToTargets, region, ('RA', 'DEC'))
    coordinates = SkyCoord(ra=ra_i*u.degree, dec=dec_i*u.degree, frame='icrs')
    print(f"Loaded {len(coordinates)} target coordinates from {pathToTargets}")
    return coordinates

def get_targets_as_stars(pathToTargets, region=None):
    """
    Get target stars from a FITS file.

    Args:
        pathToTargets (str): Path to the FITS file containing target coordinates.
        region (str, optional): FITS extension name. Defaults to In_LSST_footprint.

    Returns:
        tuple: Skyfield stars and Target_no values in matching catalog row order.
    """
    targets = []
    names = []
    ra_i, dec_i, name_i = _read_target_catalog_columns(pathToTargets, region, ('RA', 'DEC', 'Target_no'))
    for right_ascension, declination, name in zip(ra_i, dec_i, name_i):
        targets.append(Star(Angle(degrees=right_ascension), Angle(degrees=declination)))
        names.append(name)

    return targets, names

def save_twilight_windows_numpy(local_start, local_end, twilight_events, filename_base):
    """Save labeled UTC-MJD twilight events with the existing metadata layout."""
    
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
    
    filename = f"{filename_base}.npz"
    os.makedirs(os.path.dirname(filename) or '.', exist_ok=True)
    display_filename = _display_path(filename)
    print(f'Writing to {display_filename}...')
    np.savez_compressed(filename,
                       twilight_data=twilight_array,
                       metadata=np.array([local_start, local_end]))
    
    print(f"Twilight data saved to {display_filename}")

def find_twilight_window_visit_overlaps(twilight_windows, visit_exposure_groups, extension_time=1.5, verbose=True):
    """Find overlaps between twilight windows and visit/exposure groups"""
    counter = 0
    windows = [(Time(start, format='mjd'), Time(end, format='mjd')) for start, end in twilight_windows]
    twilight_extension = extension_time * u.hr
    for visit in visit_exposure_groups.values():

        visit['twilight_overlap'] = False
        visit_mjd = Time(visit['mjd'], format='mjd')
        visit_chile_time = visit['chile_time']
        # Check if visit falls within any twilight window
        for start_mjd, end_mjd in windows:
            if start_mjd - twilight_extension <= visit_mjd <= end_mjd + twilight_extension:
                counter += 1
                visit['twilight_overlap'] = True
                if verbose:
                    print(f"Visit {visit['visit_id']} has a visit at Chile Time {visit_chile_time} " +
                    f"which falls within twilight window {start_mjd.to_datetime(timezone=_CHILE_TIMEZONE)} to {end_mjd.to_datetime(timezone=_CHILE_TIMEZONE)}")
                break  # No need to check other windows for this visit
    print(f"Total twilight-overlapping visits on all days: {counter}")

    return visit_exposure_groups

def load_twilight_windows(filename, verbose=True):
    """Load complete same-period start/end pairs; unmatched boundaries stay unpaired."""
    print(f'Reading from {filename}...')
    with np.load(filename) as data:
        twilight_data = data['twilight_data']
    print(f"Loaded {len(twilight_data)} twilight events")
    windows = []
    pending = None
    previous_mjd = None
    labels = {
        'am twilight start': ('am', 'start'), 'am twilight end': ('am', 'end'),
        'pm twilight start': ('pm', 'start'), 'pm twilight end': ('pm', 'end'),
    }
    for i, event in enumerate(twilight_data):
        event_type = event['type']
        event_mjd = event['mjd']
        label = event_type.decode() if isinstance(event_type, bytes) else str(event_type)
        if label not in labels:
            raise ValueError(f"Unknown twilight event type: {label!r}")
        if not np.isfinite(event_mjd) or (previous_mjd is not None and event_mjd < previous_mjd):
            raise ValueError("Twilight event MJDs must be finite and chronological")
        previous_mjd = event_mjd
        period, boundary = labels[label]
        if boundary == 'start':
            pending = (period, event_mjd)
        else:
            if pending is not None and pending[0] == period:
                windows.append([pending[1], event_mjd])
            pending = None
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
        mjd_time = t.to_astropy().utc.mjd
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
    """Write visit-target rows to FITS, preserving visit order and complete IDs."""
    if len(visits_in_twilight) == 0:
        raise ValueError("At least one visit is required")
    # Prepare data for FITS table
    fits_rows = []
    local_start = Time(min(visit['mjd'] for visit in visits_in_twilight), format='mjd').to_datetime(timezone=_CHILE_TIMEZONE)
    local_end = Time(max(visit['mjd'] for visit in visits_in_twilight), format='mjd').to_datetime(timezone=_CHILE_TIMEZONE)
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

    identifier_width = max(20, max((len(str(row['target_id'])) for row in fits_rows), default=0))
    column_formats = [
        ('visit_id', 'K'), ('mjd', 'D'), ('chile_time', '20A'),
        ('exposure_time', 'E'), ('num_exposures', 'I'), ('filter', '1A'),
        ('field_ra', 'D'), ('field_dec', 'D'), ('target_id', f'{identifier_width}A'),
        ('target_ra', 'D'), ('target_dec', 'D'), ('five_sigma_depth', 'E'),
        ('seeing', 'E'), ('sky_brightness', 'E'),
    ]
    cols = [
        fits.Column(name=name, format=column_format, array=[row[name] for row in fits_rows])
        for name, column_format in column_formats
    ]
    hdu = fits.BinTableHDU.from_columns(cols)
    # Write to FITS file
    fits_filename = f"{pathToFiles}/twilight_visits/{local_start.strftime('%Y%m%d')}_to_{local_end.strftime('%Y%m%d')}.fits"
    os.makedirs(os.path.dirname(fits_filename), exist_ok=True)
    display_filename = _display_path(fits_filename)
    print(f'Writing to {display_filename}...')
    hdu.writeto(fits_filename, overwrite=True)
    print(f"Written {len(fits_rows)} rows to {display_filename}")
    print(f"FITS columns: {[col.name for col in hdu.columns]}")

    return fits_filename

def _read_exposure_catalog(exposures_file, identifier_column):
    """Copy the visit-target FITS table and report its row/group counts."""
    print(f'Reading from {exposures_file}...')
    with fits.open(exposures_file) as hdul:
        data = hdul[1].data
        print(len(data[identifier_column]))
        print(len(np.unique(data[identifier_column])))
        return data.copy()


def get_exposure_windows(exposures_file):
    """
    Get exposure windows from a FITS file.

    Args:
        exposures_file (str): Path to the FITS file containing exposure windows.

    Returns:
        dict: Chile-local dates mapped to visit IDs and their exposure records.
    """
    data = _read_exposure_catalog(exposures_file, 'visit_id')
    exposure_dict = {}
    if len(data) == 0:
        return exposure_dict
    ts = load.timescale()
    for visit in np.unique(data['visit_id']):
        visit_mask = data['visit_id'] == visit
        visit_data = data[visit_mask]
        day = Time(visit_data['mjd'][0], format='mjd').to_datetime(timezone=_CHILE_TIMEZONE).date()
        if day not in exposure_dict:
            exposure_dict[day] = {}
        start = Time(visit_data['mjd'][0], format='mjd').to_datetime(timezone=_CHILE_TIMEZONE)
        start = ts.from_datetime(start)
        window = [start + dt.timedelta(seconds=i) for i in range(0, int(visit_data['exposure_time'][0]) + 1)]
        exposure_dict[day][visit] = {
            'visit_id': visit,
            'mjd': visit_data['mjd'][0],
            'chile_time': visit_data['chile_time'][0],
            'exposure_time': visit_data['exposure_time'][0],
            'num_exposures': visit_data['num_exposures'][0],
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
    return exposure_dict

def get_targets_by_visit_exposure(bundle_list, targets_dict, start_mjd, end_mjd):
    """Group targets by visit ID, rejecting conflicting visit-level metadata."""
    
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
            visit_metadata = {
                'visit_id': visit_id, 'mjd': mjd, 'exposure_time': exposure_time,
                'num_exposures': num_exposures, 'filter': obs['filter'],
                'field_ra': obs['fieldRA'], 'field_dec': obs['fieldDec'],
            }
            if visit_key not in visit_groups:
                obs_time = Time(mjd, format='mjd')
                utc_datetime = obs_time.to_datetime(timezone=dt.timezone.utc)
                chile_datetime = utc_datetime.astimezone(_CHILE_TIMEZONE)
                visit_groups[visit_key] = {
                    **visit_metadata, 'chile_time': chile_datetime, 'targets': [],
                }
            else:
                for field, value in visit_metadata.items():
                    existing = visit_groups[visit_key][field]
                    both_missing = pd.isna(existing) and pd.isna(value)
                    if not both_missing and (pd.isna(existing) or pd.isna(value) or existing != value):
                        raise ValueError(f"Conflicting metadata for visit {visit_id!r}: {field}")
            
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
    """Return one Skyfield star and display name per unique sorted target ID."""
    data = _read_exposure_catalog(exposures_file, 'target_id')
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
                sample = next((value for value in col if value is not None and len(value)), sample)
                # If elements are strings -> join
                if len(sample) and all(isinstance(x, str) for x in sample):
                    joined = [','.join(map(str, v)) if v is not None else '' for v in col]
                    maxlen = max((len(s) for s in joined), default=1)
                    out[name] = np.asarray(joined, dtype=f'S{maxlen}')
                else:
                    arrays = [np.asarray(value) if value is not None else np.asarray([]) for value in col]
                    populated = [array for array in arrays if array.size]
                    if not populated or len({array.shape for array in arrays}) > 1:
                        dtype = np.result_type(*(array.dtype for array in populated)) if populated else np.dtype(float)
                        values = np.empty(len(arrays), dtype=object)
                        values[:] = [array.astype(dtype, copy=False) for array in arrays]
                    else:
                        values = arrays
                    out[name] = Column(values, name=name)
            else:
                out[name] = col
        else:
            out[name] = col
    return out

def _satellite_above_horizon(satellite, observer, check_time, altitude_degrees):
    """Return whether instantaneous altitude strictly exceeds the degree limit."""
    altitude, _, _ = (satellite - observer).at(check_time).altaz()
    return altitude.degrees > altitude_degrees


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
    if check_sunlit and eph is None:
        raise ValueError("An ephemeris is required when check_sunlit=True")
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

            above_horizon_start = _satellite_above_horizon(satellite, observer, start_time, altitude_degrees)

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
    if check_sunlit and eph is None:
        raise ValueError("An ephemeris is required when check_sunlit=True")
    satellites_up = []
    
    for satellite in satellites:
        try:
            if _satellite_above_horizon(satellite, observer, check_time, altitude_degrees):
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
            is_above_horizon = _satellite_above_horizon(satellite, observer, check_time, altitude_degrees)
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
