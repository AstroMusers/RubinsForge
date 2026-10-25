import ast
import datetime as dt
import io
import logging
import os
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace
from unittest.mock import MagicMock, Mock

from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Column, Table
from astropy import units
import numpy as np
import pandas as pd
from pytz import timezone
import pytest

from logging_stream import StreamToLogger


class FakeTime:
    """Datetime-backed Skyfield time interface for offline file-selection tests."""

    def __init__(self, value):
        self.value = value

    def __sub__(self, other):
        return (self.value - other.value).total_seconds() / 86400

    def __add__(self, days):
        delta = days if isinstance(days, dt.timedelta) else dt.timedelta(days=days)
        return FakeTime(self.value + delta)

    def utc_datetime(self):
        return self.value

    def astimezone(self, zone):
        return self.value.astimezone(zone)


def test_path_to_files_is_relative_to_module_checkout():
    module_path = Path(__file__).parents[1] / 'src' / 'utils' / 'utils.py'
    print(f'Reading from {module_path}...')
    tree = ast.parse(module_path.read_text())
    assignment = next(
        node for node in tree.body
        if isinstance(node, ast.Assign)
        and any(isinstance(target, ast.Name) and target.id == 'pathToFiles' for target in node.targets)
    )
    namespace = {'Path': Path, '__file__': str(module_path)}

    exec(compile(ast.Module(body=[assignment], type_ignores=[]), str(module_path), 'exec'), namespace)

    assert Path(namespace['pathToFiles']) == module_path.resolve().parents[2] / 'files'


def test_display_path_uses_repository_relative_names(utilities, tmp_path):
    repository_output = Path(__file__).parents[1] / 'files' / 'twilight_windows' / 'example.npz'

    assert utilities['_display_path'](repository_output) == 'Skyfield/files/twilight_windows/example.npz'
    assert utilities['_display_path'](tmp_path / 'external-output.fits') == 'external-output.fits'


@pytest.fixture
def utilities():
    path = Path(__file__).parents[1] / 'src' / 'utils' / 'utils.py'
    print(f'Reading from {path}...')
    tree = ast.parse(path.read_text())
    names = {'_display_path', '_find_closest_tle_path', 'find_closest_TLE', 'find_closest_TLE_lines', 'parse_tle_time', '_satellite_above_horizon', 'find_satellites_up', 'find_satellites_up_simple', 'get_satellite_count_info', '_sample_exposure_window', 'get_time_array', 'find_twilight_window_visit_overlaps', '_read_target_catalog_columns', 'get_targets_as_coordinates', 'get_targets_as_stars', 'read_targets_file', 'get_targets_by_visit_exposure', 'create_exposure_windows', '_read_exposure_catalog', 'get_exposure_windows', 'exposure_targets_to_stars', 'load_twilight_windows', 'save_twilight_windows_numpy'}
    definitions = [
        node for node in tree.body
        if (isinstance(node, ast.FunctionDef) and node.name in names)
        or (isinstance(node, ast.Assign) and any(isinstance(target, ast.Name) and target.id in {'_CHILE_TIMEZONE', 'pathToFiles'} for target in node.targets))
    ]
    namespace = {
        'dt': dt, 'os': os, 'np': np, 'Path': Path, '__file__': str(path), 'Time': Time, 'u': units, 'timezone': timezone,
        'fits': fits, 'Table': Table, 'Column': Column, 'SkyCoord': SkyCoord, 'pd': pd,
        'Angle': lambda **kwargs: SimpleNamespace(**kwargs),
        'Star': lambda right_ascension, declination: SimpleNamespace(ra=right_ascension, dec=declination),
        'StreamToLogger': StreamToLogger,
        'load': SimpleNamespace(
            timescale=Mock(return_value=SimpleNamespace(from_datetime=FakeTime)),
            tle_file=Mock(return_value=['synthetic satellite']),
        ),
    }
    sanitizer = next(node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == 'fits_sanitize_table')
    definitions.append(sanitizer)
    exec(compile(ast.Module(body=definitions, type_ignores=[]), str(path), 'exec'), namespace)
    return namespace


@pytest.mark.parametrize('reader', ['find_closest_TLE', 'find_closest_TLE_lines'])
def test_empty_tle_directory_raises_clear_error_without_loading_satellites(utilities, tmp_path, reader):
    observed = FakeTime(dt.datetime(2026, 1, 31).astimezone(dt.timezone.utc))
    with pytest.raises(FileNotFoundError, match='less than 30 days'):
        utilities[reader](observed, tmp_path, dt.timezone.utc)
    utilities['load'].tle_file.assert_not_called()


@pytest.mark.parametrize('reader', ['find_closest_TLE', 'find_closest_TLE_lines'])
def test_both_tle_readers_choose_the_same_nearest_file(utilities, tmp_path, reader, capsys):
    for filename in ('tle_01.20.26.txt', 'tle_01.30.26.txt'):
        path = tmp_path / filename
        print(f'Writing to {path}...')
        path.write_text('SATELLITE\n1 synthetic\n2 synthetic\n')
    observed = FakeTime(dt.datetime(2026, 1, 31).astimezone(dt.timezone.utc))

    filename, satellites = utilities[reader](observed, tmp_path, dt.timezone.utc)

    assert filename == 'tle_01.30.26.txt'
    assert f'Reading from {tmp_path / filename}...' in capsys.readouterr().out
    if reader == 'find_closest_TLE':
        assert satellites == ['synthetic satellite']
        utilities['load'].tle_file.assert_called_once_with(str(tmp_path / filename))
    else:
        assert satellites == [('SATELLITE', '1 synthetic', '2 synthetic')]
        utilities['load'].tle_file.assert_not_called()


def test_older_pipeline_tle_call_supplies_the_existing_timezone(utilities, tmp_path):
    path = tmp_path / 'tle_01.30.26.txt'
    print(f'Writing to {path}...')
    path.write_text('SATELLITE\n1 synthetic\n2 synthetic\n')
    caller = Path(__file__).parents[1] / 'scripts' / 'skyfield_sim.py'
    print(f'Reading from {caller}...')
    tree = ast.parse(caller.read_text())
    call = next(node for node in ast.walk(tree) if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id == 'find_closest_TLE')
    namespace = {
        **utilities, 't00': FakeTime(dt.datetime(2026, 1, 31).astimezone(dt.timezone.utc)),
        'pathToTLE': tmp_path, 'zone': dt.timezone.utc,
    }

    filename, satellites = eval(compile(ast.Expression(call), str(caller), 'eval'), namespace)

    assert filename == path.name
    assert satellites == ['synthetic satellite']


@pytest.mark.parametrize('reader', ['find_closest_TLE', 'find_closest_TLE_lines'])
@pytest.mark.parametrize('days_apart', [29, 30, 31])
def test_tle_selection_keeps_the_existing_strict_thirty_day_limit(utilities, tmp_path, reader, days_apart):
    filename_date = dt.datetime(2026, 1, 31) - dt.timedelta(days=days_apart)
    path = tmp_path / f'tle_{filename_date:%m.%d.%y}.txt'
    print(f'Writing to {path}...')
    path.write_text('SATELLITE\n1 synthetic\n2 synthetic\n')
    observed = FakeTime(dt.datetime(2026, 1, 31).astimezone(dt.timezone.utc))

    if days_apart < 30:
        assert utilities[reader](observed, tmp_path, dt.timezone.utc)[0] == path.name
    else:
        with pytest.raises(FileNotFoundError, match='less than 30 days'):
            utilities[reader](observed, tmp_path, dt.timezone.utc)
        utilities['load'].tle_file.assert_not_called()


@pytest.mark.parametrize('reader', ['find_closest_TLE', 'find_closest_TLE_lines'])
def test_tle_selection_ignores_dated_directories_and_invalid_filenames(utilities, tmp_path, reader):
    (tmp_path / 'tle_01.31.26.txt').mkdir()
    invalid = tmp_path / 'README'
    print(f'Writing to {invalid}...')
    invalid.write_text('Not an orbit file')
    valid = tmp_path / 'tle_01.30.26.txt'
    print(f'Writing to {valid}...')
    valid.write_text('SATELLITE\n1 synthetic\n2 synthetic\n')
    observed = FakeTime(dt.datetime(2026, 1, 31).astimezone(dt.timezone.utc))

    filename, _ = utilities[reader](observed, tmp_path, dt.timezone.utc)

    assert filename == valid.name


@pytest.mark.parametrize('content', ['SATELLITE\n', 'SATELLITE\n1 synthetic\n'])
def test_tle_line_reader_rejects_incomplete_records_clearly(utilities, tmp_path, content):
    path = tmp_path / 'tle_01.30.26.txt'
    print(f'Writing to {path}...')
    path.write_text(content)
    observed = FakeTime(dt.datetime(2026, 1, 31).astimezone(dt.timezone.utc))

    with pytest.raises(ValueError, match='complete three-line TLE records'):
        utilities['find_closest_TLE_lines'](observed, tmp_path, dt.timezone.utc)


@pytest.mark.parametrize('empty', [False, True])
def test_tle_line_reader_preserves_complete_and_empty_record_lists(utilities, tmp_path, empty):
    path = tmp_path / 'tle_01.30.26.txt'
    print(f'Writing to {path}...')
    path.write_text('' if empty else 'FIRST\n1 first\n2 first\nSECOND\n1 second\n2 second\n')
    observed = FakeTime(dt.datetime(2026, 1, 31).astimezone(dt.timezone.utc))

    filename, records = utilities['find_closest_TLE_lines'](observed, tmp_path, dt.timezone.utc)

    assert filename == path.name
    assert records == ([] if empty else [('FIRST', '1 first', '2 first'), ('SECOND', '1 second', '2 second')])


@pytest.fixture
def visible_satellite():
    """Synthetic above-horizon satellite with an explicit illumination result."""
    satellite = MagicMock()
    satellite.name = 'SYNTHETIC'
    satellite.find_events.return_value = ([], [])
    satellite.__sub__.return_value.at.return_value.altaz.return_value = (
        SimpleNamespace(degrees=40.0), object(), object(),
    )
    satellite.at.return_value.is_sunlit.return_value = True
    return satellite


@pytest.fixture
def visibility_times():
    start = FakeTime(dt.datetime(2026, 1, 1, tzinfo=dt.timezone.utc))
    end = start + 2 / 86400
    return start, end


@pytest.mark.parametrize('selector', ['find_satellites_up', 'find_satellites_up_simple'])
def test_requested_sunlight_selection_requires_ephemeris(utilities, visible_satellite, visibility_times, selector):
    start, end = visibility_times
    arguments = ([visible_satellite], object(), start, end) if selector == 'find_satellites_up' else ([visible_satellite], object(), start)

    with pytest.raises(ValueError, match='ephemeris is required'):
        utilities[selector](*arguments, check_sunlit=True)

    visible_satellite.find_events.assert_not_called()
    visible_satellite.at.assert_not_called()


@pytest.mark.parametrize('selector', ['find_satellites_up', 'find_satellites_up_simple'])
@pytest.mark.parametrize('sunlit', [False, True])
def test_satellite_selection_applies_the_requested_illumination(utilities, visible_satellite, visibility_times, selector, sunlit):
    start, end = visibility_times
    visible_satellite.at.return_value.is_sunlit.return_value = sunlit
    arguments = ([visible_satellite], object(), start, end) if selector == 'find_satellites_up' else ([visible_satellite], object(), start)
    ephemeris = object()

    selected = utilities[selector](*arguments, eph=ephemeris)

    assert selected == ([visible_satellite] if sunlit else [])
    assert visible_satellite.at.return_value.is_sunlit.call_count > 0
    assert all(call.args == (ephemeris,) for call in visible_satellite.at.return_value.is_sunlit.call_args_list)


@pytest.mark.parametrize('selector', ['find_satellites_up', 'find_satellites_up_simple'])
def test_explicit_altitude_only_selection_does_not_require_ephemeris(utilities, visible_satellite, visibility_times, selector):
    start, end = visibility_times
    arguments = ([visible_satellite], object(), start, end) if selector == 'find_satellites_up' else ([visible_satellite], object(), start)

    assert utilities[selector](*arguments, check_sunlit=False) == [visible_satellite]
    visible_satellite.at.assert_not_called()


@pytest.mark.parametrize('filename', ['skyfield_sim_new _optimized.py', 'skyfield_sim_new _optimized_july.py'])
def test_optimized_callers_forward_their_ephemeris(utilities, visible_satellite, visibility_times, filename):
    path = Path(__file__).parents[1] / 'scripts' / filename
    print(f'Reading from {path}...')
    tree = ast.parse(path.read_text())
    call = next(node for node in ast.walk(tree) if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute) and node.func.attr == 'find_satellites_up')
    start, end = visibility_times
    ephemeris = object()
    namespace = {
        'SatelliteTracker': SimpleNamespace(find_satellites_up=utilities['find_satellites_up']),
        'starlinks': [visible_satellite], 'rubin_obs': object(),
        'exposure': SimpleNamespace(start_time=start, end_time=end), 'eph': ephemeris,
    }

    selected, periods = eval(compile(ast.Expression(call), str(path), 'eval'), namespace)

    assert selected == [visible_satellite]
    assert periods[visible_satellite]['total_sunlit_duration'] == 2.0
    assert all(item.args == (ephemeris,) for item in visible_satellite.at.return_value.is_sunlit.call_args_list)


@pytest.mark.parametrize('selector', ['find_satellites_up', 'find_satellites_up_simple'])
@pytest.mark.parametrize('altitude', [29.9, 30.0, 30.1])
def test_both_selectors_keep_the_strict_altitude_boundary(utilities, visible_satellite, visibility_times, selector, altitude):
    visible_satellite.__sub__.return_value.at.return_value.altaz.return_value = (
        SimpleNamespace(degrees=altitude), object(), object(),
    )
    start, end = visibility_times
    arguments = ([visible_satellite], object(), start, end) if selector == 'find_satellites_up' else ([visible_satellite], object(), start)

    selected = utilities[selector](*arguments, check_sunlit=False)

    assert selected == ([visible_satellite] if altitude > 30.0 else [])


@pytest.mark.parametrize(('altitude', 'sunlit'), [(29.9, True), (30.0, True), (30.1, True), (40.0, False)])
def test_count_summary_uses_the_same_altitude_boundary(utilities, visible_satellite, visibility_times, altitude, sunlit):
    visible_satellite.__sub__.return_value.at.return_value.altaz.return_value = (
        SimpleNamespace(degrees=altitude), object(), object(),
    )
    visible_satellite.at.return_value.is_sunlit.return_value = sunlit
    instant = visibility_times[0]

    counts = utilities['get_satellite_count_info']([visible_satellite], object(), instant, eph=object())

    assert counts == {
        'total': 1, 'above_horizon': int(altitude > 30.0), 'sunlit': int(sunlit),
        'above_horizon_and_sunlit': int(altitude > 30.0 and sunlit), 'time': instant,
    }


@pytest.mark.parametrize('selector', ['find_satellites_up', 'find_satellites_up_simple'])
def test_single_satellite_errors_do_not_discard_other_valid_satellites(utilities, visible_satellite, visibility_times, selector):
    broken = MagicMock()
    broken.name = 'BROKEN'
    broken.find_events.side_effect = ValueError('invalid orbital data')
    broken.__sub__.side_effect = ValueError('invalid orbital data')
    start, end = visibility_times
    arguments = ([broken, visible_satellite], object(), start, end) if selector == 'find_satellites_up' else ([broken, visible_satellite], object(), start)

    assert utilities[selector](*arguments, check_sunlit=False) == [visible_satellite]


def test_interval_selector_preserves_sampled_sunlit_periods_and_durations(utilities, visible_satellite, visibility_times):
    start, end = visibility_times
    visible_satellite.at.side_effect = lambda instant: SimpleNamespace(is_sunlit=lambda eph: instant.value.second == 1)

    selected, periods = utilities['find_satellites_up'](
        [visible_satellite], object(), start, end, eph=object(), return_periods=True,
    )

    assert selected == [visible_satellite]
    record = periods[visible_satellite]
    assert record['visible_periods'] == [(start, end)]
    assert len(record['sunlit_periods']) == 1
    beginning, ending = record['sunlit_periods'][0]
    assert beginning.utc_datetime() == start.utc_datetime() + dt.timedelta(seconds=1)
    assert ending.utc_datetime() == end.utc_datetime()
    assert record['combined_periods'] == record['sunlit_periods']
    assert record['total_visible_duration'] == 2.0
    assert record['total_sunlit_duration'] == 1.0


def test_empty_catalog_preserves_altitude_only_return_shapes(utilities, visibility_times):
    start, end = visibility_times
    assert utilities['find_satellites_up']([], object(), start, end, check_sunlit=False) == []
    assert utilities['find_satellites_up']([], object(), start, end, check_sunlit=False, return_periods=True) == ([], {})
    assert utilities['find_satellites_up_simple']([], object(), start, check_sunlit=False) == []


@pytest.fixture
def twilight_scheduler(utilities):
    """Use explicit synthetic dawn/dusk transitions, not ephemeris predictions."""
    states = Mock(return_value=np.int64(0))
    def transitions(start, end, function):
        return ([
            start + dt.timedelta(hours=6), start + dt.timedelta(hours=6, minutes=1),
            start + dt.timedelta(hours=18), start + dt.timedelta(hours=18, minutes=1),
        ], [1, 4, 1, 0])
    events = Mock(side_effect=transitions)
    utilities.update(
        sf=SimpleNamespace(almanac=SimpleNamespace(dark_twilight_day=lambda *args: states)),
        almanac=SimpleNamespace(
            find_discrete=events,
            TWILIGHTS={0: 'Dark', 1: 'Astronomical twilight', 4: 'Day'},
        ),
    )
    return SimpleNamespace(schedule=utilities['get_time_array'], states=states, events=events)


@pytest.mark.parametrize('step_days', [1, 4])
def test_twilight_scheduler_queries_each_processed_day_and_preserves_one_second_schedule(twilight_scheduler, step_days):
    start = dt.datetime(2026, 1, 1, tzinfo=dt.timezone.utc)
    day_step = dt.timedelta(days=step_days)

    schedule = twilight_scheduler.schedule(
        start, start + 2 * day_step, object(), object(), dt.timezone.utc,
        dt.timedelta(), dt.timedelta(seconds=1), day_step, dt.timedelta(seconds=30),
    )

    assert [call.args[0] for call in twilight_scheduler.states.call_args_list] == [start, start + day_step]
    assert [call.args[0] for call in twilight_scheduler.events.call_args_list] == [start, start + day_step]
    assert len(schedule) == 2
    for day_index, exposures in enumerate(schedule):
        day = start + day_index * day_step
        assert [exposure[0] for exposure in exposures] == [
            day + dt.timedelta(hours=6), day + dt.timedelta(hours=6, seconds=30),
            day + dt.timedelta(hours=18, seconds=1), day + dt.timedelta(hours=18, seconds=31),
        ]
        for exposure in exposures:
            assert len(exposure) == 15
            assert all(exposure[index + 1] - exposure[index] == dt.timedelta(seconds=1) for index in range(14))


def test_twilight_scheduler_keeps_exposure_start_cadence_for_two_second_samples(twilight_scheduler):
    start = dt.datetime(2026, 1, 1, tzinfo=dt.timezone.utc)

    schedule = twilight_scheduler.schedule(
        start, start + dt.timedelta(days=1), object(), object(), dt.timezone.utc,
        dt.timedelta(), dt.timedelta(seconds=2), dt.timedelta(days=1), dt.timedelta(seconds=30),
    )

    assert [exposure[0] for exposure in schedule[0]] == [
        start + dt.timedelta(hours=6), start + dt.timedelta(hours=6, seconds=30),
        start + dt.timedelta(hours=18, seconds=2),
    ]
    for exposure in schedule[0]:
        assert len(exposure) == 15
        assert all(exposure[index + 1] - exposure[index] == dt.timedelta(seconds=2) for index in range(14))


@pytest.mark.parametrize('sample_seconds', [0.5, 1.0, 2.0])
def test_shared_exposure_window_keeps_configured_steps_and_start_cadence(utilities, sample_seconds):
    start = dt.datetime(2026, 1, 1, tzinfo=dt.timezone.utc)
    stop = start + dt.timedelta(seconds=60)
    step = dt.timedelta(seconds=sample_seconds)

    exposures, cursor = utilities['_sample_exposure_window'](start, stop, step, dt.timedelta(seconds=30))

    assert [exposure[0] for exposure in exposures] == [start, start + dt.timedelta(seconds=30)]
    assert cursor >= stop
    for exposure in exposures:
        assert exposure == [exposure[0] + index * step for index in range(15)]


@pytest.mark.parametrize('seconds', [0.0, 10.0, 15.0])
def test_shared_exposure_window_discards_only_incomplete_blocks(utilities, seconds):
    start = dt.datetime(2026, 1, 1, tzinfo=dt.timezone.utc)

    exposures, cursor = utilities['_sample_exposure_window'](
        start, start + dt.timedelta(seconds=seconds), dt.timedelta(seconds=1), dt.timedelta(seconds=30),
    )

    assert len(exposures) == (1 if seconds == 15.0 else 0)
    assert cursor >= start + dt.timedelta(seconds=seconds)


@pytest.mark.parametrize(('step', 'cadence'), [(0, 30), (-1, 30), (1, 0), (1, -1)])
def test_shared_exposure_window_rejects_nonadvancing_steps(utilities, step, cadence):
    start = dt.datetime(2026, 1, 1, tzinfo=dt.timezone.utc)
    with pytest.raises(ValueError, match='must be positive'):
        utilities['_sample_exposure_window'](
            start, start + dt.timedelta(minutes=1), dt.timedelta(seconds=step), dt.timedelta(seconds=cadence),
        )


@pytest.mark.parametrize('day_step', [dt.timedelta(0), dt.timedelta(days=-1)])
def test_twilight_scheduler_rejects_nonadvancing_day_steps(twilight_scheduler, day_step):
    start = dt.datetime(2026, 1, 1, tzinfo=dt.timezone.utc)
    with pytest.raises(ValueError, match='Day step must be positive'):
        twilight_scheduler.schedule(
            start, start + dt.timedelta(days=1), object(), object(), dt.timezone.utc,
            dt.timedelta(), dt.timedelta(seconds=1), day_step, dt.timedelta(seconds=30),
        )
    twilight_scheduler.events.assert_not_called()


@pytest.mark.parametrize('failing_day', [0, 1])
@pytest.mark.parametrize('no_events', [False, True])
def test_twilight_scheduler_never_reuses_boundaries_when_a_day_lacks_transitions(twilight_scheduler, failing_day, no_events):
    start = dt.datetime(2026, 1, 1, tzinfo=dt.timezone.utc)
    failure_date = start + dt.timedelta(days=failing_day)
    complete_transitions = twilight_scheduler.events.side_effect
    def transitions(today, tomorrow, state):
        times, events = complete_transitions(today, tomorrow, state)
        if today == failure_date:
            return ([], []) if no_events else (times[1:], events[1:])
        return times, events
    twilight_scheduler.events.side_effect = transitions

    with pytest.raises(ValueError, match='Missing required twilight transitions') as error:
        twilight_scheduler.schedule(
            start, start + dt.timedelta(days=2), object(), object(), dt.timezone.utc,
            dt.timedelta(), dt.timedelta(seconds=1), dt.timedelta(days=1), dt.timedelta(seconds=30),
        )

    assert failure_date.strftime('%Y-%m-%d') in str(error.value)
    assert twilight_scheduler.events.call_count == failing_day + 1


@pytest.mark.parametrize('window_seconds', [14, 15, 40])
def test_twilight_scheduler_respects_observation_end_and_complete_exposures(twilight_scheduler, window_seconds):
    start = dt.datetime(2026, 1, 1, tzinfo=dt.timezone.utc)
    end = start + dt.timedelta(hours=6, seconds=window_seconds)

    schedule = twilight_scheduler.schedule(
        start, end, object(), object(), dt.timezone.utc, dt.timedelta(),
        dt.timedelta(seconds=1), dt.timedelta(days=1), dt.timedelta(seconds=30),
    )

    assert len(schedule) == 1
    assert len(schedule[0]) == (0 if window_seconds == 14 else 1)
    assert all(start <= instant < end for exposure in schedule[0] for instant in exposure)


def test_twilight_scheduler_clips_long_extensions_to_observation_bounds(twilight_scheduler):
    """A deliberately long synthetic extension isolates interval-bound behavior."""
    start = dt.datetime(2026, 1, 1, tzinfo=dt.timezone.utc)
    end = start + dt.timedelta(days=1)

    schedule = twilight_scheduler.schedule(
        start, end, object(), object(), dt.timezone.utc, dt.timedelta(hours=7),
        dt.timedelta(seconds=1), dt.timedelta(days=1), dt.timedelta(hours=1),
    )

    assert schedule[0]
    assert all(start <= instant < end for exposure in schedule[0] for instant in exposure)
    assert schedule[0][0][0] == start


def test_logging_stream_preserves_complete_print_lines(utilities):
    logger = Mock()
    stream = utilities['StreamToLogger'](logger, 20)

    print('satellite', 42, file=stream)

    logger.log.assert_called_once_with(20, 'satellite 42')
    assert stream.linebuf == ''


def test_logging_stream_buffers_partial_writes_until_newline_or_flush(utilities):
    logger = Mock()
    stream = utilities['StreamToLogger'](logger, 30)

    assert stream.write('first ') == 6
    assert stream.write('part') == 4
    logger.log.assert_not_called()
    assert stream.write('\nsecond') == 7
    logger.log.assert_called_once_with(30, 'first part')
    assert stream.linebuf == 'second'
    stream.flush()
    assert [item.args for item in logger.log.call_args_list] == [(30, 'first part'), (30, 'second')]
    assert stream.linebuf == ''
    stream.flush()
    assert logger.log.call_count == 2


def test_logging_stream_preserves_spacing_blank_lines_and_split_crlf(utilities):
    logger = Mock()
    stream = utilities['StreamToLogger'](logger, 40)

    stream.write('  indented  \r')
    stream.write('\n\nlast\n')

    assert [item.args for item in logger.log.call_args_list] == [
        (40, '  indented  '), (40, ''), (40, 'last'),
    ]


def test_logging_stream_accepts_empty_writes_and_rejects_nontext(utilities):
    logger = Mock()
    stream = utilities['StreamToLogger'](logger, 20)
    assert stream.write('') == 0
    stream.flush()
    logger.log.assert_not_called()
    with pytest.raises(TypeError):
        stream.write(b'bytes')
    assert stream.linebuf == ''


def test_existing_utility_reexports_the_lightweight_logging_stream():
    path = Path(__file__).parents[1] / 'src' / 'utils' / 'utils.py'
    print(f'Reading from {path}...')
    tree = ast.parse(path.read_text())
    imports = [node for node in tree.body if isinstance(node, ast.ImportFrom) and node.module == 'logging_stream']
    assert len(imports) == 1
    namespace = {}

    exec(compile(ast.Module(body=imports, type_ignores=[]), str(path), 'exec'), namespace)

    assert namespace['StreamToLogger'] is StreamToLogger


@pytest.mark.parametrize('path', sorted(
    Path(__file__).parents[1].glob('skyfield_101224_positions_and_contamination_final_ver*.py')
) + [
    Path(__file__).parents[1] / 'skyfield_101224_positions_and_contamination_targetpos.py',
    Path(__file__).parents[1] / 'skyfield_101224_positions_and_contamination_targetpos_funct.py',
], ids=lambda path: path.name)
def test_historical_stream_imports_share_the_repaired_class(path):
    print(f'Reading from {path}...')
    tree = ast.parse(path.read_text())
    assert not any(isinstance(node, ast.ClassDef) and node.name == 'StreamToLogger' for node in tree.body)
    imports = [node for node in tree.body if isinstance(node, ast.ImportFrom) and node.module == 'logging_stream']
    namespace = {}

    exec(compile(ast.Module(body=imports, type_ignores=[]), str(path), 'exec'), namespace)

    assert namespace['StreamToLogger'] is StreamToLogger
    logger = Mock()
    print('satellite', 42, file=namespace['StreamToLogger'](logger, 20))
    logger.log.assert_called_once_with(20, 'satellite 42')


def test_logging_stream_integrates_with_real_handlers_and_separate_stdout_stderr_buffers():
    output = io.StringIO()
    logger = logging.Logger('rubinsforge.stream.test', level=logging.DEBUG)
    handler = logging.StreamHandler(output)
    handler.setFormatter(logging.Formatter('%(levelname)s: %(message)s'))
    logger.addHandler(handler)
    stdout = StreamToLogger(logger, logging.INFO)
    stderr = StreamToLogger(logger, logging.ERROR)

    stdout.write('satellite ')
    stderr.write('error\n')
    print(42, file=stdout)

    assert output.getvalue() == 'ERROR: error\nINFO: satellite 42\n'


def test_logging_stream_import_does_not_load_satellite_or_numerical_packages():
    result = subprocess.run(
        [sys.executable, '-c',
         'import sys; before = set(sys.modules); from logging_stream import StreamToLogger; '
         'assert not {"src", "numpy", "skyfield", "lumos"}.intersection(set(sys.modules) - before); '
         'print(StreamToLogger.__name__)'],
        cwd=Path(__file__).parents[1], check=True, text=True, capture_output=True,
    )
    assert result.stdout == 'StreamToLogger\n'


@pytest.mark.parametrize('verbose', [False, True])
def test_twilight_visit_overlap_uses_supplied_windows_and_counts_each_visit_once(utilities, capsys, verbose):
    visits = {
        'before': {'mjd': 60000.0, 'visit_id': 1, 'chile_time': 'synthetic'},
        'inside': {'mjd': 60000.375, 'visit_id': 2, 'chile_time': 'synthetic'},
        'after': {'mjd': 60001.0, 'visit_id': 3, 'chile_time': 'synthetic'},
    }
    capsys.readouterr()

    returned = utilities['find_twilight_window_visit_overlaps'](
        [(60000.25, 60000.5), (60000.3, 60000.6)], visits, extension_time=0, verbose=verbose,
    )

    assert returned is visits
    assert [visit['twilight_overlap'] for visit in visits.values()] == [False, True, False]
    output = capsys.readouterr().out
    assert 'Total twilight-overlapping visits on all days: 1' in output
    assert ('Visit 2 has a visit' in output) == verbose


@pytest.mark.parametrize(
    ('mjd', 'expected'),
    [(60000.1875, True), (60000.1874, False), (60000.5625, True), (60000.5626, False)],
)
def test_twilight_visit_overlap_keeps_inclusive_hour_extension_bounds(utilities, mjd, expected):
    visits = {'visit': {'mjd': mjd, 'visit_id': 1, 'chile_time': 'synthetic'}}

    utilities['find_twilight_window_visit_overlaps']([(60000.25, 60000.5)], visits, extension_time=1.5, verbose=False)

    assert visits['visit']['twilight_overlap'] is expected


def test_twilight_visit_overlap_preserves_empty_inputs_and_resets_stale_flags(utilities):
    visits = {'visit': {'mjd': 60000.0, 'visit_id': 1, 'chile_time': 'synthetic', 'twilight_overlap': True}}
    assert utilities['find_twilight_window_visit_overlaps']([], visits, verbose=False) is visits
    assert not visits['visit']['twilight_overlap']
    assert utilities['find_twilight_window_visit_overlaps']([], {}, verbose=False) == {}


def test_twilight_visit_overlap_reuses_generator_windows_for_every_visit(utilities):
    windows = (window for window in [(60000.0, 60001.0)])
    visits = {
        key: {'mjd': 60000.5, 'visit_id': index, 'chile_time': 'synthetic'}
        for index, key in enumerate(['first', 'second'])
    }

    utilities['find_twilight_window_visit_overlaps'](windows, visits, extension_time=0, verbose=False)

    assert all(visit['twilight_overlap'] for visit in visits.values())


@pytest.fixture
def target_catalog(tmp_path):
    path = tmp_path / 'targets.fits'
    print(f'Writing to {path}...')
    default = Table({'RA': [10.0, 25.0], 'DEC': [-20.0, 30.0], 'Target_no': ['FIRST', 'SECOND']})
    selected = Table({'RA': [45.0], 'DEC': [15.0], 'Target_no': ['REGION']})
    fits.HDUList([
        fits.PrimaryHDU(), fits.BinTableHDU(default, name='In_LSST_footprint'),
        fits.BinTableHDU(selected, name='CUSTOM'),
    ]).writeto(path)
    return path


@pytest.mark.parametrize(('region', 'expected_ra', 'expected_dec', 'expected_names'), [
    (None, [10.0, 25.0], [-20.0, 30.0], ['FIRST', 'SECOND']),
    ('CUSTOM', [45.0], [15.0], ['REGION']),
])
def test_fits_catalog_adapters_share_region_selection_and_preserve_degree_order(utilities, target_catalog, capsys, region, expected_ra, expected_dec, expected_names):
    capsys.readouterr()

    coordinates = utilities['get_targets_as_coordinates'](target_catalog, region=region)
    stars, names = utilities['get_targets_as_stars'](target_catalog, region=region)

    np.testing.assert_array_equal(coordinates.ra.degree, expected_ra)
    np.testing.assert_array_equal(coordinates.dec.degree, expected_dec)
    assert [star.ra.degrees for star in stars] == expected_ra
    assert [star.dec.degrees for star in stars] == expected_dec
    assert names == expected_names
    assert capsys.readouterr().out.count(f'Reading from {target_catalog}...') == 2


def test_coordinate_only_catalog_does_not_require_star_identifiers(utilities, tmp_path):
    path = tmp_path / 'coordinates.fits'
    print(f'Writing to {path}...')
    fits.HDUList([
        fits.PrimaryHDU(), fits.BinTableHDU(Table({'RA': [10.0], 'DEC': [-20.0]}), name='In_LSST_footprint'),
    ]).writeto(path)

    coordinates = utilities['get_targets_as_coordinates'](path)

    np.testing.assert_array_equal(coordinates.ra.degree, [10.0])
    with pytest.raises(KeyError):
        utilities['get_targets_as_stars'](path)


def test_empty_fits_catalog_preserves_coordinate_and_star_return_shapes(utilities, tmp_path):
    path = tmp_path / 'empty.fits'
    print(f'Writing to {path}...')
    table = Table({'RA': np.array([], dtype=float), 'DEC': np.array([], dtype=float), 'Target_no': np.array([], dtype='U8')})
    fits.HDUList([fits.PrimaryHDU(), fits.BinTableHDU(table, name='In_LSST_footprint')]).writeto(path)

    coordinates = utilities['get_targets_as_coordinates'](path)
    stars, names = utilities['get_targets_as_stars'](path)

    assert coordinates.shape == (0,)
    assert stars == []
    assert names == []


def test_shared_catalog_columns_are_independent_copies_and_missing_regions_still_fail(utilities, target_catalog):
    right_ascension, declination = utilities['_read_target_catalog_columns'](target_catalog, None, ('RA', 'DEC'))
    assert right_ascension.flags.owndata
    assert declination.flags.owndata
    right_ascension[0] = 999.0

    reread = utilities['get_targets_as_coordinates'](target_catalog)

    assert reread.ra.degree[0] == 10.0
    with pytest.raises(KeyError):
        utilities['get_targets_as_coordinates'](target_catalog, region='MISSING')


@pytest.mark.parametrize('metadata', [False, True])
def test_csv_catalog_return_dictionary_coordinates_and_metadata_are_preserved(utilities, tmp_path, capsys, metadata):
    path = tmp_path / 'targets.csv'
    data = {'ra': [10.0, 25.0], 'dec': [-20.0, 30.0]}
    if metadata:
        data.update(name=['FIRST', 'SECOND'], mag=[20.0, 21.0])
    print(f'Writing to {path}...')
    pd.DataFrame(data).to_csv(path, index=False)
    capsys.readouterr()

    targets = utilities['read_targets_file'](path, verbose=False)

    assert list(targets) == ['target_0000', 'target_0001']
    for index, target in enumerate(targets.values()):
        assert target['ra'] == data['ra'][index]
        assert target['dec'] == data['dec'][index]
        assert target['coord'].ra.degree == data['ra'][index]
        assert target['coord'].dec.degree == data['dec'][index]
        assert target['name'] == (data['name'][index] if metadata else None)
        assert target['mag'] == (data['mag'][index] if metadata else None)
        assert target['flag'] is None
    assert f'Reading from {path}...' in capsys.readouterr().out


@pytest.fixture
def visit_record():
    return {
        'visit_id': 7, 'mjd': 60000.25,
        'chile_time': Time(60000.25, format='mjd').to_datetime(timezone=timezone('Chile/Continental')),
        'exposure_time': 15.0, 'num_exposures': 1, 'filter': 'r',
        'field_ra': 12.0, 'field_dec': -30.0,
        'targets': [{
            'target_id': 'TARGET', 'target_ra': 13.0, 'target_dec': -29.0,
            'five_sigma_depth': 24.0, 'seeing': 0.8, 'sky_brightness': 21.0,
        }],
    }


def test_visit_grouping_uses_the_declared_chile_timezone_and_keeps_target_metadata(utilities, target_observations):
    bundle = SimpleNamespace(metric_values=[target_observations, target_observations])
    targets = {'FIRST': {'ra': 13.0, 'dec': -29.0}, 'SECOND': {'ra': 14.0, 'dec': -28.0}}

    groups = utilities['get_targets_by_visit_exposure']([bundle], targets, 60000.25, 60000.25)

    assert list(groups) == ['7']
    visit = groups['7']
    assert visit['chile_time'] == Time(60000.25, format='mjd').to_datetime(timezone=timezone('Chile/Continental'))
    assert [target['target_id'] for target in visit['targets']] == ['FIRST', 'SECOND']
    assert [target['target_ra'] for target in visit['targets']] == [13.0, 14.0]
    assert visit['exposure_time'] == 15.0


def test_visit_fits_export_uses_the_shared_timezone_and_preserves_record_fields(utilities, tmp_path, visit_record):
    utilities['pathToFiles'] = str(tmp_path)
    (tmp_path / 'twilight_visits').mkdir()

    filename = utilities['create_exposure_windows']([visit_record])

    print(f'Reading from {filename}...')
    with fits.open(filename) as restored:
        table = restored[1].data
        assert len(table) == 1
        assert table['visit_id'][0] == 7
        assert table['mjd'][0] == 60000.25
        assert table['chile_time'][0] == visit_record['chile_time'].strftime('%Y-%m-%d %H:%M:%S')
        assert table['target_id'][0] == 'TARGET'
        assert table['target_ra'][0] == 13.0
        assert table['target_dec'][0] == -29.0
        assert table['five_sigma_depth'][0] == 24.0


def test_visit_fits_export_creates_directory_and_names_full_unsorted_date_range(utilities, tmp_path, visit_record):
    utilities['pathToFiles'] = str(tmp_path)
    earlier = {**visit_record, 'visit_id': 8, 'mjd': 60000.25}
    later = {**visit_record, 'visit_id': 9, 'mjd': 60002.25}

    filename = utilities['create_exposure_windows']([later, earlier])

    expected_start = Time(earlier['mjd'], format='mjd').to_datetime(timezone=timezone('Chile/Continental')).strftime('%Y%m%d')
    expected_end = Time(later['mjd'], format='mjd').to_datetime(timezone=timezone('Chile/Continental')).strftime('%Y%m%d')
    assert Path(filename).name == f'{expected_start}_to_{expected_end}.fits'
    print(f'Reading from {filename}...')
    with fits.open(filename) as restored:
        assert restored[1].data['visit_id'].tolist() == [9, 8]


def test_visit_fits_export_retains_complete_long_target_identifiers(utilities, tmp_path, visit_record):
    utilities['pathToFiles'] = str(tmp_path)
    (tmp_path / 'twilight_visits').mkdir()
    identifier = 'TARGET-' + 'X' * 40
    visit_record['targets'][0]['target_id'] = identifier

    filename = utilities['create_exposure_windows']([visit_record])

    print(f'Reading from {filename}...')
    with fits.open(filename) as restored:
        assert restored[1].data['target_id'][0] == identifier


def test_visit_fits_export_rejects_empty_visits_before_creating_outputs(utilities, tmp_path):
    utilities['pathToFiles'] = str(tmp_path)
    with pytest.raises(ValueError, match='At least one visit'):
        utilities['create_exposure_windows']([])
    assert list(tmp_path.iterdir()) == []


def test_visit_fits_export_keeps_a_declared_empty_schema_when_no_targets_match(utilities, tmp_path, visit_record):
    utilities['pathToFiles'] = str(tmp_path)
    visit_record['targets'] = []

    filename = utilities['create_exposure_windows']([visit_record])

    print(f'Reading from {filename}...')
    with fits.open(filename) as restored:
        assert len(restored[1].data) == 0
        assert restored[1].columns.names == [
            'visit_id', 'mjd', 'chile_time', 'exposure_time', 'num_exposures',
            'filter', 'field_ra', 'field_dec', 'target_id', 'target_ra',
            'target_dec', 'five_sigma_depth', 'seeing', 'sky_brightness',
        ]
        assert restored[1].columns['target_id'].format == '20A'


def test_exported_visits_round_trip_through_both_exposure_readers(utilities, tmp_path, visit_record, capsys):
    utilities['pathToFiles'] = str(tmp_path)
    visit_record['targets'].append({**visit_record['targets'][0], 'target_id': 'SECOND', 'target_ra': 14.0})
    filename = utilities['create_exposure_windows']([visit_record])
    utilities['load'].timescale.reset_mock()
    capsys.readouterr()

    stars, names = utilities['exposure_targets_to_stars'](filename)

    utilities['load'].timescale.assert_not_called()
    assert names == ['Target SECOND', 'Target TARGET']
    assert [star.ra.degrees for star in stars] == [14.0, 13.0]
    assert f'Reading from {filename}...' in capsys.readouterr().out
    visits = utilities['get_exposure_windows'](filename)
    day = visit_record['chile_time'].date()
    assert list(visits) == [day]
    record = visits[day][7]
    assert record['target_ids'].tolist() == ['TARGET', 'SECOND']
    assert record['target_ras'].tolist() == [13.0, 14.0]
    assert record['visit_center_as_star'].ra.degrees == 12.0
    assert record['visit_center_as_star'].dec.degrees == -30.0
    assert len(record['window']) == 16
    assert record['window'][0].utc_datetime() == visit_record['chile_time']
    assert (record['window'][-1].utc_datetime() - record['window'][0].utc_datetime()).total_seconds() == 15.0


def test_empty_exported_catalog_preserves_both_reader_return_shapes(utilities, tmp_path, visit_record):
    utilities['pathToFiles'] = str(tmp_path)
    visit_record['targets'] = []
    filename = utilities['create_exposure_windows']([visit_record])

    assert utilities['exposure_targets_to_stars'](filename) == ([], [])
    assert utilities['get_exposure_windows'](filename) == {}
    utilities['load'].timescale.assert_not_called()


def test_repeated_target_ids_across_visits_keep_first_file_row_and_date_groups(utilities, tmp_path, visit_record):
    """Synthetic differing coordinates isolate existing first-file-row behavior."""
    utilities['pathToFiles'] = str(tmp_path)
    later = {
        **visit_record, 'visit_id': 8, 'mjd': 60001.25,
        'chile_time': Time(60001.25, format='mjd').to_datetime(timezone=timezone('Chile/Continental')),
        'targets': [{**visit_record['targets'][0], 'target_ra': 14.0}],
    }
    filename = utilities['create_exposure_windows']([later, visit_record])

    stars, names = utilities['exposure_targets_to_stars'](filename)
    visits = utilities['get_exposure_windows'](filename)

    assert names == ['Target TARGET']
    assert [star.ra.degrees for star in stars] == [14.0]
    assert set(visits) == {visit_record['chile_time'].date(), later['chile_time'].date()}
    assert visits[visit_record['chile_time'].date()][7]['target_ras'].tolist() == [13.0]
    assert visits[later['chile_time'].date()][8]['target_ras'].tolist() == [14.0]


def test_shared_exposure_catalog_snapshot_does_not_change_original_fits(utilities, tmp_path, visit_record):
    utilities['pathToFiles'] = str(tmp_path)
    filename = utilities['create_exposure_windows']([visit_record])
    snapshot = utilities['_read_exposure_catalog'](filename, 'target_id')
    snapshot['target_ra'][0] = 99.0

    reread = utilities['_read_exposure_catalog'](filename, 'target_id')

    assert reread['target_ra'][0] == 13.0


@pytest.mark.parametrize('reader', ['get_exposure_windows', 'exposure_targets_to_stars'])
def test_exposure_readers_keep_missing_identifier_errors_explicit(utilities, tmp_path, reader):
    path = tmp_path / 'missing-identifiers.fits'
    print(f'Writing to {path}...')
    fits.HDUList([fits.PrimaryHDU(), fits.BinTableHDU(Table({'other': [1]}))]).writeto(path)

    with pytest.raises(KeyError):
        utilities[reader](path)
    utilities['load'].timescale.assert_not_called()


@pytest.mark.parametrize('scale', ['utc', 'tt'])
def test_twilight_writer_stores_utc_mjd_independently_of_synthetic_ut1_offset(scale):
    """A deliberately supplied DUT1 offset isolates time-scale conversion."""
    path = Path(__file__).parents[1] / 'src' / 'utils' / 'utils.py'
    print(f'Reading from {path}...')
    tree = ast.parse(path.read_text())
    function = next(node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == 'create_twilight_windows')
    assignment = next(node for node in ast.walk(function) if isinstance(node, ast.Assign) and any(isinstance(target, ast.Name) and target.id == 'mjd_time' for target in node.targets))
    utc = Time('2026-01-01T06:00:00', scale='utc')
    event = SimpleNamespace(
        ut1=utc.jd + 0.25 / 86400,
        to_astropy=lambda: getattr(utc, scale),
    )
    namespace = {'Time': Time, 't': event}

    exec(compile(ast.Module(body=[assignment], type_ignores=[]), str(path), 'exec'), namespace)

    assert namespace['mjd_time'] == utc.mjd


@pytest.mark.parametrize(
    ('events', 'expected'),
    [
        ([('am twilight start', 60000.2), ('am twilight end', 60000.25), ('pm twilight start', 60000.75), ('pm twilight end', 60000.8)], [[60000.2, 60000.25], [60000.75, 60000.8]]),
        ([('am twilight end', 60000.25), ('pm twilight start', 60000.75), ('pm twilight end', 60000.8)], [[60000.75, 60000.8]]),
        ([('am twilight start', 60000.2)], []),
        ([('am twilight start', 60000.2), ('pm twilight start', 60000.75), ('pm twilight end', 60000.8)], [[60000.75, 60000.8]]),
        ([], []),
    ],
)
def test_twilight_archive_loader_pairs_only_complete_matching_boundaries(utilities, tmp_path, events, expected):
    path = tmp_path / 'windows.npz'
    print(f'Writing to {path}...')
    data = np.array(events, dtype=[('type', 'U20'), ('mjd', 'f8')])
    np.savez_compressed(path, twilight_data=data)

    assert utilities['load_twilight_windows'](path, verbose=False) == expected


@pytest.mark.parametrize('events', [
    [('unknown', 60000.0)],
    [('am twilight start', np.nan)],
    [('am twilight start', 60000.5), ('am twilight end', 60000.4)],
])
def test_twilight_archive_loader_rejects_unknown_or_invalid_records(utilities, tmp_path, events):
    path = tmp_path / 'invalid.npz'
    print(f'Writing to {path}...')
    np.savez_compressed(path, twilight_data=np.array(events, dtype=[('type', 'U20'), ('mjd', 'f8')]))
    with pytest.raises(ValueError):
        utilities['load_twilight_windows'](path, verbose=False)


@pytest.mark.parametrize('verbose', [False, True])
def test_twilight_archive_writer_and_loader_round_trip_schema_and_complete_windows(utilities, tmp_path, capsys, verbose):
    """String metadata permits schema inspection without unpickling time objects."""
    stem = tmp_path / 'nested' / 'twilight'
    path = Path(f'{stem}.npz')
    records = [
        [('am twilight start', 'synthetic', 60000.2)],
        [('am twilight end', 'synthetic', 60000.25)],
        [('pm twilight start', 'synthetic', 60000.75)],
        [('pm twilight end', 'synthetic', 60000.8)],
    ]
    capsys.readouterr()

    assert utilities['save_twilight_windows_numpy']('start', 'end', records, stem) is None
    windows = utilities['load_twilight_windows'](path, verbose=verbose)

    output = capsys.readouterr().out
    assert f'Writing to {path.name}...' in output
    assert f'Twilight data saved to {path.name}' in output
    assert f'Reading from {path}...' in output
    assert ('am twilight start at MJD' in output) == verbose
    assert windows == [[60000.2, 60000.25], [60000.75, 60000.8]]
    print(f'Reading from {path}...')
    with np.load(path) as archive:
        assert set(archive.files) == {'twilight_data', 'metadata'}
        events = archive['twilight_data']
        assert events.dtype.names == ('type', 'mjd')
        assert events['type'].tolist() == [record[0][0] for record in records]
        np.testing.assert_array_equal(events['mjd'], [record[0][2] for record in records])
        np.testing.assert_array_equal(archive['metadata'], ['start', 'end'])


@pytest.fixture
def target_observations():
    dtype = [
        ('observationStartMJD', 'f8'), ('observationId', 'i8'), ('visitExposureTime', 'f8'),
        ('numExposures', 'i4'), ('filter', 'U1'), ('fieldRA', 'f8'), ('fieldDec', 'f8'),
        ('fiveSigmaDepth', 'f8'), ('seeingFwhmEff', 'f8'), ('skyBrightness', 'f8'),
    ]
    return np.array([(60000.25, 7, 15.0, 1, 'r', 12.0, -30.0, 24.0, 0.8, 21.0)], dtype=dtype)


@pytest.mark.parametrize('rows', [
    [[], ['SATELLITE-A', 'SATELLITE-B']],
    [['SATELLITE-A', 'SATELLITE-B'], []],
    [None, [], ['SATELLITE-A', 'SATELLITE-B']],
])
def test_fits_sanitizer_preserves_string_lists_after_empty_or_missing_rows(utilities, tmp_path, rows):
    values = np.empty(len(rows), dtype=object)
    values[:] = rows
    table = Table({'names': values})

    sanitized = utilities['fits_sanitize_table'](table)

    filename = tmp_path / 'string_lists.fits'
    print(f'Writing to {filename}...')
    sanitized.write(filename, format='fits')
    print(f'Reading from {filename}...')
    with fits.open(filename) as restored:
        assert restored[1].data['names'].tolist() == [','.join(value) if value is not None else '' for value in rows]
    assert table['names'].tolist() == rows


@pytest.mark.parametrize('rows', [
    [[1.0], [2.0]],
    [[], [2.0, 3.0]],
    [[1.0], [2.0, 3.0]],
    [None, [2.0, 3.0]],
    [[], []],
    [[1], [2.5, 3.5]],
    [[1.5], [2, 3]],
    [[], [2**60 + 1, 2**60 + 2]],
])
def test_fits_sanitizer_round_trips_numeric_list_rows_without_truncation(utilities, tmp_path, rows):
    values = np.empty(len(rows), dtype=object)
    values[:] = rows
    table = Table({'samples': values})

    sanitized = utilities['fits_sanitize_table'](table)

    filename = tmp_path / 'numeric_lists.fits'
    print(f'Writing to {filename}...')
    sanitized.write(filename, format='fits')
    print(f'Reading from {filename}...')
    with fits.open(filename) as restored:
        for actual, expected in zip(restored[1].data['samples'], rows):
            np.testing.assert_array_equal(actual, expected if expected is not None else [])
    assert table['samples'].tolist() == rows


@pytest.mark.parametrize(('column', 'value', 'field'), [
    ('observationStartMJD', 60001.25, 'mjd'),
    ('visitExposureTime', 30.0, 'exposure_time'),
    ('numExposures', 2, 'num_exposures'),
    ('filter', 'g', 'filter'),
    ('fieldRA', 13.0, 'field_ra'),
    ('fieldDec', -31.0, 'field_dec'),
    ('visitExposureTime', np.nan, 'exposure_time'),
    ('fieldRA', np.nan, 'field_ra'),
])
def test_visit_grouping_rejects_conflicting_shared_metadata(utilities, target_observations, column, value, field):
    inconsistent = target_observations.copy()
    inconsistent[column] = value
    bundle = SimpleNamespace(metric_values=[target_observations, inconsistent])
    targets = {'FIRST': {'ra': 13.0, 'dec': -29.0}, 'SECOND': {'ra': 14.0, 'dec': -28.0}}

    with pytest.raises(ValueError, match=rf'Conflicting metadata for visit .*: {field}'):
        utilities['get_targets_by_visit_exposure']([bundle], targets, 60000.0, 60002.0)


def test_visit_grouping_keeps_target_specific_values_and_consistent_unknown_metadata(utilities, target_observations):
    target_observations['fieldRA'] = np.nan
    second = target_observations.copy()
    second['seeingFwhmEff'] = 1.2
    second['skyBrightness'] = 22.0
    second['fiveSigmaDepth'] = 23.0
    bundle = SimpleNamespace(metric_values=[target_observations, second])
    targets = {'FIRST': {'ra': 13.0, 'dec': -29.0}, 'SECOND': {'ra': 14.0, 'dec': -28.0}}

    grouped = utilities['get_targets_by_visit_exposure']([bundle], targets, 60000.0, 60002.0)['7']

    assert np.isnan(grouped['field_ra'])
    assert [target['target_id'] for target in grouped['targets']] == ['FIRST', 'SECOND']
    assert [target['seeing'] for target in grouped['targets']] == [0.8, 1.2]
    assert [target['sky_brightness'] for target in grouped['targets']] == [21.0, 22.0]
    assert [target['five_sigma_depth'] for target in grouped['targets']] == [24.0, 23.0]


def test_consistent_grouped_visits_round_trip_to_fits_with_one_time_conversion(utilities, target_observations, tmp_path):
    second = target_observations.copy()
    second['seeingFwhmEff'] = 1.2
    bundle = SimpleNamespace(metric_values=[target_observations, second])
    targets = {'FIRST': {'ra': 13.0, 'dec': -29.0}, 'SECOND': {'ra': 14.0, 'dec': -28.0}}
    clock = Mock(wraps=Time)
    utilities['Time'] = clock

    groups = utilities['get_targets_by_visit_exposure']([bundle], targets, 60000.25, 60000.25)

    clock.assert_called_once_with(60000.25, format='mjd')
    utilities['pathToFiles'] = str(tmp_path)
    filename = utilities['create_exposure_windows'](list(groups.values()))
    print(f'Reading from {filename}...')
    with fits.open(filename) as restored:
        table = restored[1].data
        assert table['target_id'].tolist() == ['FIRST', 'SECOND']
        np.testing.assert_array_equal(table['mjd'], [60000.25, 60000.25])
        np.testing.assert_array_equal(table['target_ra'], [13.0, 14.0])
        np.testing.assert_array_equal(table['seeing'], np.array([0.8, 1.2], dtype=np.float32))
