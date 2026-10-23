import ast
from datetime import datetime, timedelta, timezone
import io
import logging
from pathlib import Path
import subprocess
import sys
from types import ModuleType, SimpleNamespace
from unittest.mock import MagicMock, Mock, call

from astropy import units
from astropy.coordinates import SkyCoord, angular_separation
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
import numpy as np
import pytest

from satellite_photometry import calculate_observer_brightness


@pytest.fixture
def calculations():
    """Load the checked helper without importing optional satellite pipelines."""
    path = Path(__file__).parents[1] / 'src' / 'utils' / 'calculations.py'
    print(f'Reading from {path}...')
    tree = ast.parse(path.read_text())
    definitions = [node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name in {'target_check', 'trail_event', 'record_radec_trail', 'filter_targets', 'manual_angular_separation', 'brightness_calculation', 'point_brightness_calculation'}]
    namespace = {
        'time': SimpleNamespace(time=lambda: 0.0),
        'logger': logging.getLogger('rubinsforge.test.target_check'),
        'u': units,
        'np': np,
        'angular_separation': angular_separation,
        'calculate_observer_brightness': calculate_observer_brightness,
    }
    exec(compile(ast.Module(body=definitions, type_ignores=[]), str(path), 'exec'), namespace)
    return namespace


@pytest.mark.parametrize(('altitudes', 'azimuths', 'heights'), [
    ([30.0, 40.0], [100.0, 120.0], [500000.0, 550000.0]),
    (30.0, 100.0, 500000.0),
    ([], [], []),
])
def test_shared_photometry_preserves_grid_model_and_conversion(altitudes, azimuths, heights):
    """Check geometry plumbing with explicit synthetic observer and Lumos interfaces."""
    observer = Mock()
    observation_time = object()
    sun = object()
    surfaces = object()
    observer.at.return_value.observe.return_value.apparent.return_value.altaz.return_value = (
        SimpleNamespace(degrees=-12.0), SimpleNamespace(degrees=220.0), object(),
    )
    expected_altitudes, expected_azimuths = np.meshgrid(altitudes, azimuths)
    intensities = np.ones_like(expected_altitudes)
    magnitudes = np.full_like(expected_altitudes, 9.0)
    calculator = Mock(return_value=intensities)
    converter = Mock(return_value=magnitudes)

    result = calculate_observer_brightness(
        observer, observation_time, sun, altitudes, azimuths, heights,
        surfaces=surfaces, intensity_calculator=calculator, magnitude_converter=converter,
    )

    assert result[0] is intensities
    assert result[1] is magnitudes
    observer.at.assert_called_once_with(observation_time)
    observer.at.return_value.observe.assert_called_once_with(sun)
    actual_surfaces, actual_heights, actual_altitudes, actual_azimuths = calculator.call_args.args
    assert actual_surfaces is surfaces
    assert actual_heights is heights
    np.testing.assert_array_equal(actual_altitudes, expected_altitudes)
    np.testing.assert_array_equal(actual_azimuths, expected_azimuths)
    assert calculator.call_args.kwargs == {'sun_altitude': -12.0, 'sun_azimuth': 220.0, 'include_earthshine': False}
    converter.assert_called_once_with(intensities)


@pytest.mark.parametrize('point', [False, True])
@pytest.mark.parametrize('timing_guard', [False, True])
def test_maintained_photometry_adapters_preserve_surface_choice_shape_and_guard(calculations, point, timing_guard):
    """Preserve public adapters with explicit synthetic observer/Lumos interfaces."""
    observer = Mock()
    observation_time = object()
    sun = object()
    observer.at.return_value.observe.return_value.apparent.return_value.altaz.return_value = (
        SimpleNamespace(degrees=-12.0), SimpleNamespace(degrees=220.0), object(),
    )
    intensities = np.full((1, 1) if point else (2, 2), 2.0)
    magnitudes = np.full_like(intensities, 9.0)
    calculator = Mock(return_value=intensities)
    surfaces = object()
    calculations.update(
        lumos=SimpleNamespace(calculator=SimpleNamespace(get_intensity_observer_frame=calculator), conversions=SimpleNamespace(intensity_to_ab_mag=Mock(return_value=magnitudes))),
        starlink_sat=SimpleNamespace(SURFACES_LAB_BRDFS=surfaces),
        time=SimpleNamespace(time=Mock(side_effect=[0.0, np.nan if timing_guard else 0.0])),
    )

    if point:
        result = calculations['point_brightness_calculation'](observer, observation_time, sun, 30.0, 100.0, 500000.0)
    else:
        result = calculations['brightness_calculation'](observer, observation_time, sun, 1, 3, [0, 30, 40], [0, 100, 120], [0, 500000, 550000])

    if timing_guard:
        assert result is None
        calculator.assert_not_called()
        observer.at.assert_not_called()
    else:
        assert calculator.call_args.args[0] is surfaces
        np.testing.assert_array_equal(calculator.call_args.args[1], 500000.0 if point else [500000, 550000])
        if point:
            assert result == (2.0, 9.0)
        else:
            assert result[0] is intensities
            assert result[1] is magnitudes


@pytest.mark.parametrize('path', sorted(Path(__file__).parents[1].glob('skyfield_101224_positions_and_contamination_final_ver*.py')), ids=lambda path: path.name)
@pytest.mark.parametrize('timing_guard', [False, True])
def test_historical_photometry_adapters_keep_current_time_and_inferred_surfaces(path, timing_guard):
    print(f'Reading from {path}...')
    tree = ast.parse(path.read_text())
    imports = [node for node in tree.body if isinstance(node, ast.ImportFrom) and node.module == 'satellite_photometry']
    assert len(imports) == 1
    node = next(node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == 'brightness_calculation')
    observer = Mock()
    observer.at.return_value.observe.return_value.apparent.return_value.altaz.return_value = (
        SimpleNamespace(degrees=-12.0), SimpleNamespace(degrees=220.0), object(),
    )
    intensities = np.ones((2, 2))
    magnitudes = np.full((2, 2), 9.0)
    calculator = Mock(return_value=intensities)
    surfaces = object()
    timestamp = object()
    namespace = {
        'time': SimpleNamespace(time=Mock(side_effect=[0.0, np.nan]) if timing_guard else lambda: 0.0),
        'lumos': SimpleNamespace(calculator=SimpleNamespace(get_intensity_observer_frame=calculator), conversions=SimpleNamespace(intensity_to_ab_mag=Mock(return_value=magnitudes))),
        'starlink_sat': SimpleNamespace(SURFACES_INFER_BRDFS=surfaces), 'ti': timestamp,
    }
    exec(compile(ast.Module(body=imports + [node], type_ignores=[]), str(path), 'exec'), namespace)
    args = (observer, object(), 1, 3, [0, 30, 40], [0, 100, 120], [0, 500000, 550000])

    result = namespace['brightness_calculation'](*args)

    if timing_guard:
        assert result is None
        observer.at.assert_not_called()
        calculator.assert_not_called()
    else:
        assert result[0] is intensities and result[1] is magnitudes
        assert calculator.call_args.args[0] is surfaces
        np.testing.assert_array_equal(calculator.call_args.args[1], [500000, 550000])
        observer.at.assert_called_once_with(timestamp)
        next_time = object()
        namespace['ti'] = next_time
        namespace['brightness_calculation'](*args)
        assert observer.at.call_args_list == [call(timestamp), call(next_time)]


def test_shared_photometry_cold_import_does_not_load_optional_pipeline_dependencies():
    subprocess.run(
        [sys.executable, '-c', "import sys; import satellite_photometry; assert all(name not in sys.modules for name in ('skyfield', 'lumos', 'src'))"],
        cwd=Path(__file__).parents[1], capture_output=True, text=True, check=True,
    )


@pytest.mark.parametrize('debug_enabled', [False, True])
@pytest.mark.parametrize(('arcseconds', 'expected'), [(5.0, True), (10.0, True), (11.0, False)])
def test_target_check_uses_supplied_altaz_and_preserves_threshold_behavior(calculations, caplog, debug_enabled, arcseconds, expected):
    target = object()
    coordinates = [SimpleNamespace(degrees=40.0), SimpleNamespace(degrees=120.0), object()]
    topocentric = Mock()
    topocentric.separation_from.return_value = SimpleNamespace(arcseconds=lambda: arcseconds)
    separations = []
    level = logging.DEBUG if debug_enabled else logging.CRITICAL

    with caplog.at_level(level, logger=calculations['logger'].name):
        result = calculations['target_check'](target, coordinates, topocentric, separations, 10.0)

    assert result is expected
    assert separations == [arcseconds]
    topocentric.separation_from.assert_called_once_with(target)
    if debug_enabled and expected:
        assert any('40.0 , 120.0' in message for message in caplog.messages)


@pytest.mark.parametrize('filename', ['skyfield_sim.py', 'skyfield_sim_new.py'])
def test_pipeline_target_check_calls_match_the_shared_signature(calculations, filename):
    path = Path(__file__).parents[1] / 'scripts' / filename
    print(f'Reading from {path}...')
    tree = ast.parse(path.read_text())
    calls = [
        node for node in ast.walk(tree)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id == 'target_check'
    ]
    assert len(calls) == 1
    topocentric = Mock()
    topocentric.separation_from.return_value = SimpleNamespace(arcseconds=lambda: 5.0)
    separations = []
    namespace = {
        **calculations,
        'target_appr': object(),
        'target_altaz': [SimpleNamespace(degrees=40.0), SimpleNamespace(degrees=120.0), object()],
        'topocentric': topocentric, 'separation': separations, 'treshold': 10.0,
    }

    result = eval(compile(ast.Expression(calls[0]), str(path), 'eval'), namespace)

    assert result is True
    assert separations == [5.0]


@pytest.mark.parametrize('elapsed', [3.0, 3.1])
def test_target_check_preserves_timeout_return_without_computing_separation(calculations, elapsed):
    calculations['time'] = SimpleNamespace(time=Mock(side_effect=[0.0, elapsed]))
    topocentric = Mock()
    separations = []

    result = calculations['target_check'](object(), [], topocentric, separations, 10.0)

    assert result is None
    assert separations == []
    topocentric.separation_from.assert_not_called()


@pytest.fixture
def trail_case():
    """Vary synthetic range and height independently to test the API, not an orbit."""
    case = SimpleNamespace(
        starlink=MagicMock(), observer=object(), position=Mock(),
        midnight=datetime(2026, 1, 1, tzinfo=timezone.utc), zone=timezone.utc,
        distance=1500 * units.km, height=500 * units.km,
        buffers={name: [] for name in ('sat_ra', 'sat_dec', 'sat_dist', 'sat_t', 'sat_al', 'sat_az', 'sat_hght')},
    )
    case.instant = SimpleNamespace(
        astimezone=lambda zone: datetime(2026, 1, 1, 1, 30, tzinfo=timezone.utc),
        to_astropy=lambda: Time('2026-01-01T01:30:00'),
    )
    case.starlink.__sub__.return_value.at.return_value = case.position
    case.position.radec.side_effect = lambda: (SimpleNamespace(hours=2.5), SimpleNamespace(degrees=-15.0), case.distance)
    case.position.altaz.side_effect = lambda: (SimpleNamespace(degrees=60.0), SimpleNamespace(degrees=120.0), case.height)
    return case


@pytest.mark.parametrize(
    ('range_km', 'height_km', 'accepted'),
    [(1500.0, 2500.0, True), (1999.0, 500.0, True), (2000.0, 500.0, False), (2100.0, 500.0, False), (float('nan'), 500.0, False)],
)
def test_legacy_radec_recorder_preserves_range_filter_units_and_time_origin(calculations, trail_case, range_km, height_km, accepted):
    trail_case.distance = range_km * units.km
    trail_case.height = height_km * units.km

    result = calculations['record_radec_trail'](
        trail_case.starlink, trail_case.observer, trail_case.instant,
        trail_case.midnight, trail_case.zone, **trail_case.buffers,
    )

    if accepted:
        assert result is trail_case.position
        assert trail_case.buffers == {
            'sat_ra': [2.5], 'sat_dec': [-15.0], 'sat_dist': [range_km * 1000],
            'sat_t': [1.5], 'sat_al': [60.0], 'sat_az': [120.0],
            'sat_hght': [height_km * 1000],
        }
    else:
        assert result is None
        assert all(not values for values in trail_case.buffers.values())
    trail_case.starlink.__sub__.assert_called_once_with(trail_case.observer)
    trail_case.starlink.__sub__.return_value.at.assert_called_once_with(trail_case.instant)


def test_older_pipeline_calls_the_range_selected_shared_recorder(calculations, trail_case):
    path = Path(__file__).parents[1] / 'scripts' / 'skyfield_sim.py'
    print(f'Reading from {path}...')
    tree = ast.parse(path.read_text())
    call = next(node for node in ast.walk(tree) if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id == 'record_radec_trail')
    namespace = {
        **calculations, **trail_case.buffers,
        'starlink': trail_case.starlink, 'rubin_obs': trail_case.observer,
        'ti': trail_case.instant, 'midnight': trail_case.midnight, 'zone': trail_case.zone,
    }

    result = eval(compile(ast.Expression(call), str(path), 'eval'), namespace)

    assert result is trail_case.position
    assert trail_case.buffers['sat_t'] == [trail_case.instant]


@pytest.mark.parametrize('range_km', [1500.0, 1999.0, 2000.0, 2100.0])
def test_range_selected_recorder_matches_the_documented_parent(calculations, trail_case, range_km):
    path = Path(__file__).parents[1] / 'skyfield_101224_positions_and_contamination_final_ver4.2.py'
    print(f'Reading from {path}...')
    tree = ast.parse(path.read_text())
    definition = next(node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == 'trail_event')
    namespace = {
        'u': units, 'time': SimpleNamespace(time=lambda: 0.0),
        'midnight': trail_case.midnight,
    }
    exec(compile(ast.Module(body=[definition], type_ignores=[]), str(path), 'exec'), namespace)
    original_buffers = {name: [] for name in trail_case.buffers}
    trail_case.distance = range_km * units.km

    expected = namespace['trail_event'](
        trail_case.starlink, trail_case.observer, trail_case.instant,
        trail_case.zone, **original_buffers,
    )
    actual = calculations['record_radec_trail'](
        trail_case.starlink, trail_case.observer, trail_case.instant,
        trail_case.midnight, trail_case.zone, **trail_case.buffers,
    )

    assert actual is expected
    assert trail_case.buffers == original_buffers


@pytest.mark.parametrize('filename', ['skyfield_sim.py', 'skyfield_sim_new.py'])
@pytest.mark.parametrize('result_kind', ['timeout', 'empty_pair', 'success'])
def test_brightness_callers_check_failure_before_unpacking(filename, result_kind):
    path = Path(__file__).parents[1] / 'scripts' / filename
    print(f'Reading from {path}...')
    tree = ast.parse(path.read_text())
    call = next(node for node in ast.walk(tree) if isinstance(node, ast.Assign) and isinstance(node.value, ast.Call) and isinstance(node.value.func, ast.Name) and node.value.func.id == 'brightness_calculation')
    body = next(node.body for node in ast.walk(tree) if isinstance(node, ast.If) and call in node.body)
    index = body.index(call)
    statements = body[index:index + 3]
    completed = ast.Assign(targets=[ast.Name(id='completed', ctx=ast.Store())], value=ast.Constant(True))
    loop = ast.For(
        target=ast.Name(id='iteration', ctx=ast.Store()),
        iter=ast.List(elts=[ast.Constant(0)], ctx=ast.Load()),
        body=statements + [completed], orelse=[],
    )
    module = ast.fix_missing_locations(ast.Module(body=[loop], type_ignores=[]))
    intensity = np.array([[1.0, 2.0]])
    magnitude = np.array([[3.0, 4.0]])
    result = None if result_kind == 'timeout' else ((None, None) if result_kind == 'empty_pair' else (intensity, magnitude))
    namespace = {
        argument.id: object() for argument in call.value.args if isinstance(argument, ast.Name)
    }
    namespace.update(
        brightness_calculation=Mock(return_value=result),
        trailcan_counter=2, starlink=SimpleNamespace(name='synthetic'), completed=False,
    )

    exec(compile(module, str(path), 'exec'), namespace)

    namespace['brightness_calculation'].assert_called_once()
    assert namespace['completed'] == (result_kind == 'success')
    assert namespace['trailcan_counter'] == (2 if result_kind == 'success' else 1)
    if result_kind == 'success':
        assert namespace['intensities'] is intensity
        assert namespace['ab_magnitudes'] is magnitude


@pytest.mark.parametrize('filename', ['skyfield_sim.py', 'skyfield_sim_new.py'])
def test_real_trail_calls_produce_equal_length_fits_columns_and_times(calculations, trail_case, filename):
    path = Path(__file__).parents[1] / 'scripts' / filename
    print(f'Reading from {path}...')
    tree = ast.parse(path.read_text())
    helper_name = 'record_radec_trail' if filename == 'skyfield_sim.py' else 'trail_event'
    call = next(node for node in ast.walk(tree) if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id == helper_name)
    namespace = {
        **calculations, **trail_case.buffers,
        'starlink': trail_case.starlink, 'rubin_obs': trail_case.observer,
        'ti': trail_case.instant, 'midnight': trail_case.midnight, 'zone': trail_case.zone,
        'fits': fits, 'Table': Table,
    }
    bounds = {
        target.id: node
        for node in ast.walk(tree) if isinstance(node, ast.Assign)
        for target in node.targets if isinstance(target, ast.Name) and target.id in {'s', 'f'}
    }
    exec(compile(ast.Module(body=[bounds['s']], type_ignores=[]), str(path), 'exec'), namespace)

    result = eval(compile(ast.Expression(call), str(path), 'eval'), namespace)

    exec(compile(ast.Module(body=[bounds['f']], type_ignores=[]), str(path), 'exec'), namespace)
    assert result is trail_case.position
    assert namespace['f'] - namespace['s'] == 1
    recorded = trail_case.buffers if filename == 'skyfield_sim.py' else {
        name: trail_case.buffers[name] for name in ('sat_t', 'sat_al', 'sat_az', 'sat_hght')
    }
    assert all(len(values) == 1 for values in recorded.values())
    selections = [
        node for node in tree.body
        if isinstance(node, ast.Assign) and any(isinstance(target, ast.Name) and target.id in {'sat_astrpytime', 'position_hdu'} for target in node.targets)
    ]
    if not selections:
        selections = sorted([
            node for node in ast.walk(tree)
            if isinstance(node, ast.Assign) and any(isinstance(target, ast.Name) and target.id in {'sat_astrpytime', 'position_hdu'} for target in node.targets)
        ], key=lambda node: node.lineno)
    assert len(selections) == 2
    exec(compile(ast.Module(body=selections, type_ignores=[]), str(path), 'exec'), namespace)
    artifact = io.BytesIO()
    namespace['position_hdu'].writeto(artifact)
    artifact.seek(0)
    restored = Table.read(artifact, format='fits', astropy_native=True)
    assert len(restored) == 1
    if filename == 'skyfield_sim.py':
        assert restored['RA'][0] == 2.5
        assert restored['Dec'][0] == -15.0
        assert restored['Distance'][0] == 1500000.0
    else:
        assert restored.colnames == ['Time', 'Azimuth', 'Altitude', 'Height']
    assert restored['Altitude'][0] == 60.0
    assert restored['Azimuth'][0] == 120.0
    assert restored['Height'][0] == 500000.0
    assert Time(restored['Time'])[0] == Time('2026-01-01T01:30:00')


def satellite_position(altitude, azimuth):
    return SimpleNamespace(altaz=lambda: (
        SimpleNamespace(degrees=altitude), SimpleNamespace(degrees=azimuth), object(),
    ))


def test_target_filter_includes_either_endpoint_without_changing_single_epoch_calls(calculations):
    current = SkyCoord(alt=[60.0, 60.0] * units.deg, az=[40.0, 80.0] * units.deg, frame='altaz')
    previous = SkyCoord(alt=[60.0, 60.0] * units.deg, az=[0.0, 80.0] * units.deg, frame='altaz')
    start, end = satellite_position(60.0, 359.0), satellite_position(60.0, 1.0)

    current_only = calculations['filter_targets'](current, end, start)
    both_epochs = calculations['filter_targets'](current, end, start, previous_targets_altaz=previous)

    np.testing.assert_array_equal(current_only, [False, False])
    np.testing.assert_array_equal(both_epochs, [True, False])


@pytest.mark.parametrize(('start_az', 'end_az'), [(359.0, 1.0), (1.0, 359.0)])
def test_target_filter_preserves_wraparound_and_inclusive_boundaries(calculations, start_az, end_az):
    azimuths = np.array([359.0, 359.5, 0.0, 0.5, 1.0, 1.5, 180.0])
    targets = SkyCoord(alt=np.full(azimuths.size, 60.0) * units.deg, az=azimuths * units.deg, frame='altaz')

    selected = calculations['filter_targets'](
        targets, satellite_position(60.0, end_az), satellite_position(60.0, start_az), 0.0,
    )

    np.testing.assert_array_equal(selected, [True, True, True, True, True, False, False])


def test_target_filter_rejects_unaligned_endpoint_lengths(calculations):
    with pytest.raises(ValueError, match='same length'):
        calculations['filter_targets']([], object(), object(), previous_targets_altaz=[object()])


def test_optimized_pipeline_passes_endpoint_coordinates_explicitly(calculations):
    path = Path(__file__).parents[1] / 'scripts' / 'skyfield_sim_new _optimized.py'
    print(f'Reading from {path}...')
    tree = ast.parse(path.read_text())
    call = next(node for node in ast.walk(tree) if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id == 'filter_targets')
    current = SkyCoord(alt=[60.0, 60.0] * units.deg, az=[40.0, 80.0] * units.deg, frame='altaz')
    previous = SkyCoord(alt=[60.0, 60.0] * units.deg, az=[0.0, 80.0] * units.deg, frame='altaz')
    tracker = Mock()
    tracker.get_targets_at_time.side_effect = lambda instant: {'skycoords': {'start': previous, 'end': current}[instant]}
    namespace = {
        **calculations, 'target_tracker': tracker, 'satellite_period': [('start', 'end')],
        'period_end_pos': satellite_position(60.0, 1.0), 'period_start_pos': satellite_position(60.0, 359.0),
    }

    selected = eval(compile(ast.Expression(call), str(path), 'eval'), namespace)

    np.testing.assert_array_equal(selected, [True, False])


def test_target_filter_preserves_altitude_and_azimuth_padding(calculations):
    targets = SkyCoord(
        alt=[60.0, 60.005, 60.01, 60.0] * units.deg,
        az=[120.015, 120.0, 120.0, 120.02] * units.deg,
        frame='altaz',
    )
    position = satellite_position(60.0, 120.0)

    selected = calculations['filter_targets'](targets, position, position, 30.0)

    np.testing.assert_array_equal(selected, [True, True, False, False])


@pytest.mark.parametrize('with_previous', [False, True])
def test_target_filter_preserves_empty_boolean_masks(calculations, with_previous):
    position = satellite_position(60.0, 120.0)
    kwargs = {'previous_targets_altaz': []} if with_previous else {}

    selected = calculations['filter_targets']([], position, position, **kwargs)

    assert selected.shape == (0,)
    assert selected.dtype == bool


@pytest.mark.parametrize(
    ('first_longitude', 'first_latitude', 'second_longitude', 'second_latitude', 'expected_arcseconds'),
    [
        (0.0, 0.0, 0.0, 0.0, 0.0),
        (0.0, 0.0, 90.0, 0.0, 324000.0),
        (0.0, 0.0, 180.0, 0.0, 648000.0),
        (359.0, 0.0, 1.0, 0.0, 7200.0),
        (0.0, 89.99, 180.0, 89.99, 72.0),
    ],
)
def test_shared_angular_separation_preserves_degree_and_arcsecond_units(calculations, first_longitude, first_latitude, second_longitude, second_latitude, expected_arcseconds):
    actual = calculations['manual_angular_separation'](
        first_longitude, first_latitude, second_longitude, second_latitude,
    )
    assert actual == pytest.approx(expected_arcseconds, rel=1e-10, abs=1e-8)


def test_angular_separation_remains_finite_for_antipodes_and_preserves_array_shapes(calculations):
    latitudes = np.linspace(-89.0, 89.0, 1001)
    with np.errstate(invalid='raise'):
        actual = calculations['manual_angular_separation'](
            np.zeros_like(latitudes), latitudes, np.full_like(latitudes, 180.0), -latitudes,
        )
    np.testing.assert_allclose(actual, 648000.0, rtol=1e-13)
    assert actual.shape == latitudes.shape
    broadcast = calculations['manual_angular_separation'](
        np.array([[0.0, 90.0], [180.0, 270.0]]), 0.0, np.array([0.0, 90.0]), 0.0,
    )
    np.testing.assert_allclose(broadcast, [[0.0, 0.0], [648000.0, 648000.0]], atol=1e-8)


def test_satellite_velocity_separation_passes_azimuth_before_altitude(calculations):
    path = Path(__file__).parents[1] / 'src' / 'tracking' / 'satellite_tracker.py'
    print(f'Reading from {path}...')
    tree = ast.parse(path.read_text())
    call = next(node for node in ast.walk(tree) if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id == 'manual_angular_separation')
    state = SimpleNamespace(previous_alt=80.0, previous_az=0.0, current_alt=80.0, current_az=1.0)
    first = SkyCoord(alt=80.0 * units.deg, az=0.0 * units.deg, frame='altaz')
    second = SkyCoord(alt=80.0 * units.deg, az=1.0 * units.deg, frame='altaz')

    actual = eval(compile(ast.Expression(call), str(path), 'eval'), {'self': state, **calculations})

    assert actual == pytest.approx(first.separation(second).arcsecond, rel=1e-12)


@pytest.fixture
def angular_velocity(calculations, monkeypatch):
    """Execute the real velocity method with its real numerical dependency only."""
    path = Path(__file__).parents[1] / 'src' / 'tracking' / 'satellite_tracker.py'
    print(f'Reading from {path}...')
    tree = ast.parse(path.read_text())
    method = next(node for node in ast.walk(tree) if isinstance(node, ast.FunctionDef) and node.name == 'calculate_angular_velocity')
    dependency = ModuleType('_rubinsforge_velocity.utils.calculations')
    dependency.manual_angular_separation = calculations['manual_angular_separation']
    monkeypatch.setitem(sys.modules, dependency.__name__, dependency)
    namespace = {'__package__': '_rubinsforge_velocity.tracking'}
    exec(compile(ast.Module(body=[method], type_ignores=[]), str(path), 'exec'), namespace)
    return namespace['calculate_angular_velocity']


@pytest.mark.parametrize('seconds', [2.0, 0.0, -1.0])
def test_velocity_uses_arcseconds_per_second_and_rejects_nonpositive_intervals(angular_velocity, seconds):
    initial = datetime(2026, 1, 1, tzinfo=timezone.utc)
    state = SimpleNamespace(
        previous_alt=80.0, previous_az=0.0, current_alt=80.0, current_az=1.0,
        previous_time=SimpleNamespace(utc_datetime=lambda: initial),
        current_time=SimpleNamespace(utc_datetime=lambda: initial + timedelta(seconds=seconds)),
    )

    actual = angular_velocity(state)

    if seconds <= 0:
        assert actual is None
    else:
        first = SkyCoord(alt=80.0 * units.deg, az=0.0 * units.deg, frame='altaz')
        second = SkyCoord(alt=80.0 * units.deg, az=1.0 * units.deg, frame='altaz')
        assert actual == pytest.approx(first.separation(second).arcsecond / seconds, rel=1e-12)


def test_velocity_preserves_zero_motion_and_missing_position_returns(angular_velocity):
    initial = datetime(2026, 1, 1, tzinfo=timezone.utc)
    state = SimpleNamespace(
        previous_alt=80.0, previous_az=0.0, current_alt=80.0, current_az=0.0,
        previous_time=SimpleNamespace(utc_datetime=lambda: initial),
        current_time=SimpleNamespace(utc_datetime=lambda: initial.replace(second=1)),
    )
    assert angular_velocity(state) == 0.0
    state.previous_time = None
    assert angular_velocity(state) is None
    state.previous_alt = None
    assert angular_velocity(state) is None
