"""Right-ascension tick formatting without loading satellite-analysis modules."""

from astropy.coordinates import Angle
from astropy import units


def ra_formatter(x, pos):
    """Format hour-valued ticks as HH:MM:SS rounded to the nearest second."""
    return Angle(x, unit=units.hourangle).to_string(
        unit=units.hourangle, sep=':', precision=0, pad=True
    )
