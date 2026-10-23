import ast
from datetime import datetime, timedelta, timezone
import io
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock, Mock

from astropy.coordinates import EarthLocation, SkyCoord
from astropy.io import fits
from astropy.table import Column, Table
from astropy.time import Time
from astropy import units
import numpy as np
import pytest


@pytest.fixture
def exposure_module():
    """Use real class/Astropy code with synthetic Skyfield position interfaces."""
    path = Path(__file__).parents[1] / 'src' / 'tracking' / 'exposure_tracker.py'
    print(f'Reading from {path}...')
    tree = ast.parse(path.read_text())
    definition = next(node for node in tree.body if isinstance(node, ast.ClassDef) and node.name == 'ExposureTracker')
    def terrestrial_position(latitude, longitude, elevation_m=0.0):
        location = EarthLocation.from_geodetic(longitude * units.deg, latitude * units.deg, elevation_m * units.m)
        xyz = np.array([coordinate.to_value(units.m) for coordinate in location.to_geocentric()])
        return SimpleNamespace(itrs_xyz=SimpleNamespace(m=xyz))
    earth = MagicMock()
    earth.__add__.side_effect = lambda location: SimpleNamespace(location=location)
    namespace = {
        'np': np, 'u': units, 'EarthLocation': EarthLocation, 'Table': Table, 'fits': fits,
        'Column': Column,
        'wgs84': SimpleNamespace(latlon=Mock(side_effect=terrestrial_position)),
        'earth': earth, 'S': -1, 'W': -1,
    }
    utility_path = Path(__file__).parents[1] / 'src' / 'utils' / 'utils.py'
    print(f'Reading from {utility_path}...')
    utility_tree = ast.parse(utility_path.read_text())
    sanitizer = next(node for node in utility_tree.body if isinstance(node, ast.FunctionDef) and node.name == 'fits_sanitize_table')
    exec(compile(ast.Module(body=[sanitizer], type_ignores=[]), str(utility_path), 'exec'), namespace)
    exec(compile(ast.Module(body=[definition], type_ignores=[]), str(path), 'exec'), namespace)
    return namespace


@pytest.fixture
def exposure(exposure_module):
    """Initialize collection state without the constructor's site lookup."""
    tracker_class = exposure_module['ExposureTracker']
    tracker = tracker_class.__new__(tracker_class)
    for name in (
        'timestep_data', 'all_satellite_positions', 'all_intensities',
        'all_ab_magnitudes', 'all_satellite_names', 'all_times',
        'all_angular_velocities', 'all_contaminated_targets',
        'all_contaminated_names', 'all_separations',
    ):
        setattr(tracker, name, [])
    for name in (
        'total_events', 'trail_candidates', 'contaminations',
        'no_topo_skip', 'no_brightness_skip', 'no_targets_skip',
    ):
        setattr(tracker, name, 0)
    tracker.visit_info = None
    tracker.start_time = object()
    tracker.rubin_location_astrometric = Mock()
    return tracker


@pytest.mark.parametrize('visit', [None, {}, {'field_ra': 12.0}, {'visit_center_as_star': None}])
def test_missing_visit_center_returns_none_without_observation(exposure, visit):
    exposure.visit_info = visit

    assert exposure.get_visit_region_apparent() is None
    exposure.rubin_location_astrometric.at.assert_not_called()


def test_available_visit_center_preserves_observation_forwarding(exposure):
    center = object()
    exposure.visit_info = {'visit_center_as_star': center}
    expected = exposure.rubin_location_astrometric.at.return_value.observe.return_value.apparent.return_value

    assert exposure.get_visit_region_apparent() is expected
    exposure.rubin_location_astrometric.at.assert_called_once_with(exposure.start_time)
    exposure.rubin_location_astrometric.at.return_value.observe.assert_called_once_with(center)


@pytest.mark.parametrize('count', [0, 1, 2])
@pytest.mark.parametrize('as_array', [False, True])
def test_contamination_lists_and_arrays_keep_records_and_separations_aligned(exposure, count, as_array):
    targets = [object() for _ in range(count)]
    separations = [float(index + 1) for index in range(count)]
    contamination = {
        'contaminated_targets': np.array(targets, dtype=object) if as_array else targets,
        'contaminated_names': [f'TARGET-{index}'.encode() for index in range(count)],
        'intensities': np.arange(count, dtype=float) + 1,
        'ab_magnitudes': np.arange(count, dtype=float) + 20,
        'separations': np.array(separations) if as_array else separations,
    }
    instant = object()

    exposure.add_timestep_data(instant, {}, {}, contamination)

    assert len(exposure.timestep_data) == 1
    assert exposure.timestep_data[0]['time'] is instant
    assert len(exposure.all_contaminated_targets) == count
    assert exposure.all_contaminated_names == [f'TARGET-{index}' for index in range(count)]
    assert exposure.all_separations == separations
    for index, record in enumerate(exposure.all_contaminated_targets):
        assert record['position'] is targets[index]
        assert record['time'] is instant
        assert record['intensity'] == index + 1
        assert record['ab_magnitude'] == index + 20
        assert record['separation'] == index + 1


def test_none_contamination_fields_are_treated_as_unavailable(exposure):
    exposure.add_timestep_data(object(), {}, {}, {'contaminated_targets': None, 'separations': None})
    assert len(exposure.timestep_data) == 1
    assert exposure.all_contaminated_targets == []
    assert exposure.all_separations == []


@pytest.fixture
def location():
    return EarthLocation.from_geodetic(-70 * units.deg, -30 * units.deg, 2000 * units.m)


def test_constructor_honors_explicit_skyfield_location_without_site_lookup(exposure_module, location, monkeypatch):
    lookup = Mock(side_effect=AssertionError('Explicit location must not trigger site lookup'))
    monkeypatch.setattr(EarthLocation, 'of_site', lookup)
    supplied = SimpleNamespace(itrs_xyz=SimpleNamespace(m=np.array([
        coordinate.to_value(units.m) for coordinate in location.to_geocentric()
    ])))

    exposure = exposure_module['ExposureTracker'](1, object(), object(), observer_location=supplied)

    assert exposure.rubin_location_geocentric is supplied
    assert exposure.rubin_location_astrometric.location is supplied
    np.testing.assert_allclose(
        [coordinate.to_value(units.m) for coordinate in exposure.rubin_location_geocentric_astropy.to_geocentric()],
        supplied.itrs_xyz.m,
    )
    lookup.assert_not_called()


def test_constructor_accepts_astropy_location_and_preserves_positional_visit_info(exposure_module, location, monkeypatch):
    monkeypatch.setattr(EarthLocation, 'of_site', Mock(side_effect=AssertionError('Unexpected lookup')))
    visit = {'field_ra': 12.0}

    exposure = exposure_module['ExposureTracker'](1, object(), object(), visit, observer_location=location)

    assert exposure.visit_info is visit
    assert exposure.rubin_location_geocentric_astropy is location
    assert exposure.rubin_location_astrometric.location is exposure.rubin_location_geocentric
    np.testing.assert_allclose(
        exposure.rubin_location_geocentric.itrs_xyz.m,
        [coordinate.to_value(units.m) for coordinate in location.to_geocentric()],
        atol=1e-6,
    )


def test_constructor_preserves_default_rubin_lookup(exposure_module, location, monkeypatch):
    lookup = Mock(return_value=location)
    monkeypatch.setattr(EarthLocation, 'of_site', lookup)

    exposure = exposure_module['ExposureTracker'](1, object(), object())

    assert exposure.rubin_location_geocentric_astropy is location
    lookup.assert_called_once_with('rubin')
    assert exposure_module['wgs84'].latlon.call_count == 2


@pytest.mark.parametrize('filename', ['skyfield_sim_new _optimized.py', 'skyfield_sim_new _optimized_july.py'])
def test_current_pipeline_constructor_calls_match_the_public_interface(exposure_module, location, monkeypatch, filename):
    monkeypatch.setattr(EarthLocation, 'of_site', Mock(return_value=location))
    path = Path(__file__).parents[1] / 'scripts' / filename
    print(f'Reading from {path}...')
    tree = ast.parse(path.read_text())
    call = next(node for node in ast.walk(tree) if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id == 'ExposureTracker')
    start, end = object(), object()
    visit = {'visit_id': 7}
    namespace = {
        **exposure_module, 'exposure_counter': 3, 'exp': [start, end],
        'rubin_obs': location, 'visit': visit,
    }

    exposure = eval(compile(ast.Expression(call), str(path), 'eval'), namespace)

    assert exposure.start_time is start
    assert exposure.end_time is end
    assert exposure.exposure_id == (3 if filename == 'skyfield_sim_new _optimized.py' else 7)
    assert exposure.visit_info is (None if filename == 'skyfield_sim_new _optimized.py' else visit)


def test_real_constructor_collection_and_statistics_preserve_counts_and_values(exposure_module, location):
    start = datetime(2026, 1, 1, tzinfo=timezone.utc)
    initial = SimpleNamespace(utc_datetime=lambda: start)
    final = SimpleNamespace(utc_datetime=lambda: start + timedelta(seconds=15))
    exposure = exposure_module['ExposureTracker'](7, initial, final, observer_location=location)
    target = object()
    for instant, intensity, magnitude, separation in [
        (initial, 2.0, 20.0, 5.0), (final, 4.0, 19.0, 3.0),
    ]:
        satellites = {'SATELLITE': {
            'alt': 40.0, 'az': 120.0, 'height': 500000.0,
            'intensity': intensity, 'ab_magnitude': magnitude, 'angular_velocity': 2.0,
        }}
        contamination = {
            'contaminated_targets': np.array([target], dtype=object),
            'contaminated_names': np.array([b'TARGET']),
            'intensities': np.array([intensity]), 'ab_magnitudes': np.array([magnitude]),
            'separations': np.array([separation]),
        }
        exposure.add_timestep_data(instant, satellites, {}, contamination)
        exposure.increment_counters(total_events=1, trail_candidates=1, contaminations=1)

    statistics = exposure.get_statistics()

    assert statistics['duration_seconds'] == 15.0
    assert statistics['total_events'] == statistics['trail_candidates'] == statistics['contaminations'] == 2
    assert statistics['total_satellite_observations'] == 2
    assert statistics['total_contamination_events'] == 2
    assert statistics['unique_contaminated_targets'] == 1
    assert statistics['peak_intensity'] == 4.0
    assert statistics['average_intensity'] == 3.0
    assert statistics['peak_ab_magnitude'] == 19.0
    assert statistics['average_ab_magnitude'] == 19.5
    assert statistics['closest_approach_arcsec'] == 3.0
    assert statistics['average_separation_arcsec'] == 4.0
    assert exposure.all_satellite_names == ['SATELLITE', 'SATELLITE']
    assert exposure.all_times == [initial, final]


@pytest.mark.parametrize('ellipsoid', ['GRS80', 'WGS72'])
def test_astropy_location_conversion_preserves_itrs_position_across_ellipsoids(exposure_module, ellipsoid):
    original = EarthLocation.from_geodetic(
        -70 * units.deg, -30 * units.deg, 2000 * units.m, ellipsoid=ellipsoid,
    )

    exposure = exposure_module['ExposureTracker'](1, object(), object(), observer_location=original)

    expected = np.array([coordinate.to_value(units.m) for coordinate in original.to_geocentric()])
    np.testing.assert_allclose(exposure.rubin_location_geocentric.itrs_xyz.m, expected, rtol=0, atol=1e-6)


@pytest.fixture
def recorded_exposure(exposure_module, location):
    """Construct a synthetic exposure with one real-serializable time sample."""
    start = Time('2026-01-01T00:00:00')
    end = start + 15 * units.s
    def instant(value):
        return SimpleNamespace(
            utc_datetime=lambda: value.to_datetime(timezone=timezone.utc),
            utc_iso=lambda: value.isot + 'Z', to_astropy=lambda: value,
        )
    initial = instant(start)
    exposure = exposure_module['ExposureTracker'](7, initial, instant(end), observer_location=location)
    exposure.add_timestep_data(initial, {'SATELLITE': {
        'alt': 40.0, 'az': 120.0, 'height': 500000.0,
        'intensity': 2.0, 'ab_magnitude': 20.0, 'angular_velocity': 3.0,
    }}, {}, {})
    return exposure


@pytest.mark.parametrize('trail_count', [0, 1])
def test_recorded_satellite_fits_rows_do_not_depend_on_trail_counter(recorded_exposure, trail_count):
    recorded_exposure.trail_candidates = trail_count

    extensions = recorded_exposure.create_fits_extensions()

    assert [extension.name for extension in extensions] == ['SATELLITE_DATA']
    artifact = io.BytesIO()
    fits.HDUList([fits.PrimaryHDU()] + extensions).writeto(artifact)
    artifact.seek(0)
    with fits.open(artifact) as restored:
        table = restored['SATELLITE_DATA'].data
        assert len(table) == 1
        assert table['Satellite_Name'].tolist() == ['SATELLITE']
        assert table['Altitude'][0] == 40.0
        assert table['Azimuth'][0] == 120.0
        assert table['Height'][0] == 500000.0
        assert table['Intensity'][0] == 2.0
        assert table['AB_Magnitude'][0] == 20.0
        assert table['Angular_Velocity'][0] == 3.0
        native = Table.read(restored['SATELLITE_DATA'], format='fits', astropy_native=True)
        assert native['Time'][0] == Time('2026-01-01T00:00:00')
        assert restored['SATELLITE_DATA'].header['TRAILCAN'] == trail_count
        assert restored['SATELLITE_DATA'].header['NSAT'] == 1
        assert restored['SATELLITE_DATA'].header['PEAKI'] == 2.0


@pytest.mark.parametrize('unavailable', [np.nan, np.inf, -np.inf])
def test_nonfinite_measurements_keep_table_values_and_use_undefined_summary_headers(recorded_exposure, unavailable):
    recorded_exposure.all_satellite_positions[0]['ab_magnitude'] = unavailable
    recorded_exposure.all_ab_magnitudes = [unavailable]
    recorded_exposure.all_separations = [unavailable]

    extensions = recorded_exposure.create_fits_extensions()

    artifact = io.BytesIO()
    fits.HDUList([fits.PrimaryHDU()] + extensions).writeto(artifact)
    artifact.seek(0)
    with fits.open(artifact) as restored:
        satellite = restored['SATELLITE_DATA']
        assert satellite.header['PEAKM'] is None
        assert satellite.header['AVGM'] is None
        assert satellite.header['PEAKI'] == 2.0
        assert satellite.header['AVGI'] == 2.0
        np.testing.assert_equal(satellite.data['AB_Magnitude'][0], unavailable)
        separations = restored['SEPARATIONS']
        assert separations.header['MINSEP'] is None
        assert separations.header['AVGSEP'] is None
        np.testing.assert_equal(separations.data['Separation_arcsec'][0], unavailable)


@pytest.mark.parametrize('name', ['TARGET', 'TARGET-' + 'X' * 100])
def test_contaminated_target_identifiers_and_coordinates_round_trip_without_truncation(recorded_exposure, name):
    target = SkyCoord(alt=50 * units.deg, az=130 * units.deg, frame='altaz')
    recorded_exposure.add_timestep_data(recorded_exposure.start_time, {}, {}, {
        'contaminated_targets': [target], 'contaminated_names': [name], 'separations': [2.0],
    })

    extensions = recorded_exposure.create_fits_extensions()

    artifact = io.BytesIO()
    fits.HDUList([fits.PrimaryHDU()] + extensions).writeto(artifact)
    artifact.seek(0)
    with fits.open(artifact) as restored:
        contamination = restored['CONTAMINATED_TARGETS']
        assert len(contamination.data) == 1
        assert contamination.data['Target_Name'][0] == name
        assert contamination.data['Target_Alt'][0] == 50.0
        assert contamination.data['Target_Az'][0] == 130.0
        assert contamination.header['NCONTAM'] == contamination.header['NTARG'] == 1
        native = Table.read(contamination, format='fits', astropy_native=True)
        assert native['Time'][0] == Time('2026-01-01T00:00:00')


@pytest.mark.parametrize('period_count', [1, 3])
@pytest.mark.parametrize('intersects_visit', [False, True])
def test_july_pipeline_counts_each_satellite_duration_once(exposure, period_count, intersects_visit):
    """Use synthetic periods/visit intersections to isolate real counter placement."""
    path = Path(__file__).parents[1] / 'scripts' / 'skyfield_sim_new _optimized_july.py'
    print(f'Reading from {path}...')
    tree = ast.parse(path.read_text())
    loop = next(
        node for node in ast.walk(tree)
        if isinstance(node, ast.For)
        and any(
            isinstance(statement, ast.Assign)
            and any(isinstance(target, ast.Name) and target.id == 'satellite_visible_duration' for target in statement.targets)
            for statement in node.body
        )
    )
    satellites = [MagicMock(), MagicMock()]
    satellites[0].name, satellites[1].name = 'FIRST', 'SECOND'
    periods = [(2 * index, 2 * index + 1) for index in range(period_count)]
    position = Mock(return_value=object())
    accepted = []
    namespace = {
        'satellites_up': satellites,
        'satellite_trackers': {satellite: object() for satellite in satellites},
        'satellite_periods': {
            satellites[0]: {'combined_periods': periods, 'total_visible_duration': 12.7},
            satellites[1]: {'combined_periods': periods, 'total_visible_duration': 5.8},
        },
        'exposure': exposure, 'point_trail_event': position,
        'check_visit_region': Mock(return_value=intersects_visit),
        'rubin_obs': object(), 'visit_region': object(),
        'starlinks_path_over_visit': accepted, 'logger': Mock(),
    }

    exec(compile(ast.Module(body=[loop], type_ignores=[]), str(path), 'exec'), namespace)

    assert exposure.total_events == 17
    assert exposure.no_targets_skip == (0 if intersects_visit else 2 * period_count)
    assert accepted == (satellites if intersects_visit else [])
    assert position.call_count == 4 * period_count
