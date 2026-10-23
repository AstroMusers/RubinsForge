import numpy as np


def calculate_observer_brightness(
    observer, observation_time, sun, altitudes, azimuths, heights, *,
    surfaces, intensity_calculator, magnitude_converter,
):
    """Evaluate supplied surfaces using degree angles and heights in meters."""
    apparent = observer.at(observation_time).observe(sun).apparent()
    sun_altitude, sun_azimuth, _ = apparent.altaz()
    satellite_altitudes, satellite_azimuths = np.meshgrid(altitudes, azimuths)
    intensities = intensity_calculator(
        surfaces, heights, satellite_altitudes, satellite_azimuths,
        sun_altitude=sun_altitude.degrees, sun_azimuth=sun_azimuth.degrees,
        include_earthshine=False,
    )
    return intensities, magnitude_converter(intensities)
