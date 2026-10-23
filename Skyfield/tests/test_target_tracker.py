from datetime import datetime, timedelta, timezone
from pathlib import Path
import runpy
from types import SimpleNamespace
from unittest.mock import Mock

from astropy.coordinates import EarthLocation
from astropy.time import Time
from astropy import units
import pytest


@pytest.fixture
def tracker_module():
    path = Path(__file__).parents[1] / 'src' / 'tracking' / 'target_tracker.py'
    print(f'Reading from {path}...')
    return runpy.run_path(str(path))


@pytest.fixture
def location():
    return EarthLocation.from_geodetic(-70 * units.deg, -30 * units.deg, 2000 * units.m)


def test_explicit_earth_location_is_preserved_without_site_lookup(tracker_module, location, monkeypatch):
    lookup = Mock(side_effect=AssertionError('Explicit location must not trigger a site lookup'))
    monkeypatch.setattr(EarthLocation, 'of_site', lookup)
    observer = object()

    tracker = tracker_module['TargetTracker']([], [], observer, location)

    assert tracker.rubin_location_geocentric is location
    assert tracker.rubin_location_astrometric is observer
    lookup.assert_not_called()


def test_default_location_uses_the_existing_site_lookup(tracker_module, location, monkeypatch):
    lookup = Mock(return_value=location)
    monkeypatch.setattr(EarthLocation, 'of_site', lookup)

    tracker = tracker_module['TargetTracker']([], [], object())

    assert tracker.rubin_location_geocentric is location
    lookup.assert_called_once_with('rubin')


class ExposureInstant:
    """A hashable time interface for synthetic cache tests, not an orbital clock."""

    def __init__(self, seconds):
        self.seconds = seconds

    def astimezone(self, zone):
        return (datetime(2026, 1, 1, tzinfo=timezone.utc) + timedelta(seconds=self.seconds)).astimezone(zone)

    def to_astropy(self):
        return Time(self.astimezone(timezone.utc))


@pytest.fixture
def moving_observer():
    """Return explicit synthetic positions that change with requested time."""
    def at(instant):
        def observe(target):
            altitude = 35.0 + instant.seconds if target == 'up' else 20.0 + 3 * instant.seconds
            apparent = SimpleNamespace(
                instant=instant, target=target,
                altaz=Mock(return_value=(
                    SimpleNamespace(degrees=altitude),
                    SimpleNamespace(degrees=120.0 + instant.seconds),
                    object(),
                )),
            )
            return SimpleNamespace(apparent=lambda: apparent)
        return SimpleNamespace(observe=observe)
    return SimpleNamespace(at=Mock(side_effect=at))


def test_exposure_cache_uses_each_timestamp_and_preserves_initial_target_order(tracker_module, location, moving_observer):
    tracker = tracker_module['TargetTracker'](['up', 'down'], ['UP', 'DOWN'], moving_observer, location)
    first, second = ExposureInstant(0), ExposureInstant(10)

    selected = tracker.calculate_targets_for_exposure([first, second], moving_observer, timezone.utc)

    assert selected == ['up']
    initial = tracker.get_targets_at_time(first)
    final = tracker.get_targets_at_time(second)
    assert initial is not final
    for instant, entry, altitude, azimuth in [
        (first, initial, 35.0, 120.0), (second, final, 45.0, 130.0),
    ]:
        assert entry['targets_up'] == ['up']
        assert entry['names_up'] == ['UP']
        assert entry['apparent_positions'][0].instant is instant
        entry['apparent_positions'][0].altaz.assert_called_once_with(temperature_C=12.0, pressure_mbar=750)
        assert entry['skycoords'][0].alt.degree == altitude
        assert entry['skycoords'][0].az.degree == azimuth
        assert entry['skycoords'][0].obstime == instant.to_astropy()
        assert entry['alt_az_tuples'][0][0].degrees == altitude
        assert entry['alt_az_tuples'][0][1].degrees == azimuth


def test_exposure_cache_replaces_previous_exposure_positions(tracker_module, location, moving_observer):
    tracker = tracker_module['TargetTracker'](['up'], ['UP'], moving_observer, location)
    initial, replacement = ExposureInstant(0), ExposureInstant(20)
    tracker.calculate_targets_for_exposure([initial], moving_observer, timezone.utc)

    tracker.calculate_targets_for_exposure([replacement], moving_observer, timezone.utc)

    assert not tracker.has_time_cached(initial)
    assert tracker.has_time_cached(replacement)
    assert tracker.get_skycoords(replacement)[0].alt.degree == 55.0


def test_no_visible_targets_keep_empty_per_timestamp_entries(tracker_module, location, moving_observer):
    tracker = tracker_module['TargetTracker'](['down'], ['DOWN'], moving_observer, location)
    instants = [ExposureInstant(0), ExposureInstant(10)]

    selected = tracker.calculate_targets_for_exposure(instants, moving_observer, timezone.utc)

    assert selected == []
    assert tracker.get_targets_at_time(instants[0]) is not tracker.get_targets_at_time(instants[1])
    for instant in instants:
        assert tracker.get_skycoords(instant) == []
        assert tracker.get_names_up(instant) == []


def test_default_name_getter_tracks_selected_targets_and_cache_reset(tracker_module, location, moving_observer):
    tracker = tracker_module['TargetTracker'](['down', 'up'], ['DOWN', 'UP'], moving_observer, location)
    instant = ExposureInstant(0)
    assert tracker.get_names_up() == ['DOWN', 'UP']

    tracker.calculate_targets_for_exposure([instant], moving_observer, timezone.utc)

    assert tracker.get_targets_up() == ['up']
    assert tracker.get_names_up() == ['UP']
    assert tracker.get_names_up(instant) == ['UP']
    tracker.clear_cache()
    assert tracker.get_targets_up() is None
    assert tracker.get_names_up() == ['DOWN', 'UP']
    assert tracker.get_cache_info() == {'cached_times': 0, 'exposure_duration': 0, 'targets_up': 0}


def test_empty_exposure_is_rejected_without_erasing_a_valid_cache(tracker_module, location, moving_observer):
    tracker = tracker_module['TargetTracker'](['up'], ['UP'], moving_observer, location)
    instant = ExposureInstant(0)
    tracker.calculate_targets_for_exposure([instant], moving_observer, timezone.utc)
    previous_entry = tracker.get_targets_at_time(instant)

    with pytest.raises(ValueError, match='At least one exposure time'):
        tracker.calculate_targets_for_exposure([], moving_observer, timezone.utc)

    assert tracker.get_targets_at_time(instant) is previous_entry
    assert tracker.exposure_times == [instant]
    assert tracker.get_targets_up() == ['up']
    assert tracker.get_names_up() == ['UP']
