from lsst.pipe.tasks.calexpCutout import CalexpCutoutTask
import lsst.geom
from lsst.geom import Point2D
from astropy.coordinates import SkyCoord
from astropy.visualization import astropy_mpl_style
from astropy.wcs import WCS
import numpy as np
import astropy.table as at
import astropy.units as u
import datetime as dt
from astropy.visualization.wcsaxes import Quadrangle

c = dt.datetime.now()
directoryf = "/data/a.saricaoglu/lsst_pipeline/Files/"  +  str(c.strftime("%m.%d"))
directoryp = "/data/a.saricaoglu/lsst_pipeline/Plots/"  +  str(c.strftime("%m.%d"))
def pixel_to_skycoord(wcs, reference_pixel, offsets):
    """
    Convert pixel offsets into astropy SkyCoord objects.

    Parameters:
    -----------
    wcs : astropy.wcs.WCS
        The WCS object associated with the image.
    reference_pixel : tuple
        The reference pixel (x, y) in image coordinates.
    offsets : list of tuple
        A list of (x_offset, y_offset) pixel offsets.

    Returns:
    --------
    sky_coords : list of SkyCoord
        List of sky coordinates corresponding to the input offsets.
    """
    # Convert reference pixel and offsets to absolute pixel positions
    absolute_pixels1 = [(reference_pixel[0] + dx + 50, reference_pixel[1] + dy + 50) for dx, dy in offsets]
    absolute_pixels2 = [(reference_pixel[0] + dx, reference_pixel[1] + dy) for dx, dy in offsets]
    # Convert pixel coordinates to celestial coordinates
    sky_coords1 = []
    for x, y in absolute_pixels1:
        sky_point = wcs.pixelToSky(Point2D(x, y))  # Convert pixel to sky
        ra = sky_point.getLongitude().asDegrees()  # RA in degrees
        dec = sky_point.getLatitude().asDegrees()  # Dec in degrees
        sky_coords1.append(SkyCoord(ra=ra, dec=dec, frame="icrs", unit="deg"))
        sky_coords2 = []
    for x, y in absolute_pixels2:
        sky_point = wcs.pixelToSky(Point2D(x, y))  # Convert pixel to sky
        ra = sky_point.getLongitude().asDegrees()  # RA in degrees
        dec = sky_point.getLatitude().asDegrees()  # Dec in degrees
        sky_coords2.append(SkyCoord(ra=ra, dec=dec, frame="icrs", unit="deg"))
    return sky_coords1, sky_coords2


def extract_cutouts_with_task(n, calexp, centers, size):
    """
    Extract cutouts using CalexpCutoutTask.

    Parameters:
    -----------
    calexp : lsst.afw.image.ExposureF
        The calibrated exposure (image) from which to extract cutouts.
    centers : list of tuple
        A list of (x, y) pixel coordinates for the cutout centers.
    size : int
        Size of the cutout in pixels (assumes square cutouts).

    Returns:
    --------
    cutouts : list of lsst.afw.image.ExposureF
        List of extracted cutout exposures.
    """
    # Initialize the CalexpCutoutTask
    config = CalexpCutoutTask.ConfigClass()
    cutout_task = CalexpCutoutTask(config=config)

    wcs = calexp.getWcs()
    reference_pixel = wcs.getPixelOrigin()
    # Convert to sky coordinates
    sky_coords1, sky_coords2 = pixel_to_skycoord(wcs, reference_pixel, centers)
    
    print(len(sky_coords2))
    print(sky_coords2[0])
    size = np.ones(len(sky_coords1))*size*u.pixel
    # # Print results
    # for i, coord in enumerate(sky_coords):
    #     print(f"Cutout {i + 1}: RA={coord.ra.deg}, Dec={coord.dec.deg}")
    calexp_sky_coords = calexp
    in_table = at.table.QTable([sky_coords2,size,size], names=('position','xspan','yspan'))
    # Run the task
    result = cutout_task.run(in_table, calexp)
    print('skipped ', len(result.skipped_positions))
    plot_cutouts(n, result.cutouts, calexp_sky_coords, sky_coords1)
    # Return the extracted cutouts
    return result.cutouts

# Example usage
# calexp = ...  # Your calibrated exposure
# cutout_centers = [(1000, 1500), (1200, 1700), (1400, 1900)]  # Replace with actual centers
# cutout_size = 100  # Size of the cutout in pixels

# # Extract the cutouts
# cutouts = extract_cutouts_with_task(calexp, cutout_centers, cutout_size)

# Access the cutout data


import matplotlib.pyplot as plt

def plot_cutouts(n, cutouts, calexpskycoord, skycoord):
    """
    Plot extracted cutouts.

    Parameters:
        cutouts: list of lsst.afw.image.ExposureF
            List of extracted cutout exposures.
    """
    
    cutouts = cutouts.getMaskedImages()
    full = calexpskycoord.getMaskedImage().getImage().getArray()
    f2, ax2 = plt.subplots(4, 6, figsize=(20, 13), sharex=True, sharey=True)
    #print(len(cutouts))
    m = 0
    k = 0
    p = 1
    for cutout in cutouts:
        # Get the image array from the cutout
        cut = cutout.getImage().getArray()
        #print('cutout shape', np.shape(cut))
        # Display the cutout
        im1 = ax2[k][m].imshow(cut, origin='lower', cmap='viridis')
        ax2[k][m].set_title(f'Cutout {p}')
        #plt.colorbar(im1, label='Pixel Value ', ax=ax2[k][m], shrink=0.9)
        m = m + 1
        p = p + 1
        if p == 7:
            k = 1
            m = 0
        if p == 13:
            k = 2
            m = 0
        if p == 19:
            k = 3
            m = 0
    plt.savefig(directoryp +  "/" + c.strftime('%H%M') + "/cutouts_" + str(n) + ".png",bbox_inches='tight')

    wcs = WCS(calexpskycoord.getWcs().getFitsMetadata())
    deltapix = calexpskycoord.getWcs().getPixelScale().asArcseconds()
    fig = plt.figure(figsize=(8, 6))
    #print(deltapix)
    # Using WCSAxes to create the plot with sky coordinates
    ax = fig.add_subplot(111, projection=wcs)
    ax.imshow(full, origin='lower', cmap='viridis')
    for sky in skycoord:
        # Get the image array from the cutout
        # Display the cutout
        # ax.scatter(sky.ra, sky.dec, transform=ax.get_transform('icrs'), s=100, marker='s',
        #    edgecolor='white', facecolor='none')
        q = Quadrangle((sky.ra,sky.dec )*u.deg, 100*u.arcsec*deltapix, 100*u.arcsec*deltapix,
                    label='Quadrangle', edgecolor='blue', facecolor='none',
                    transform=ax.get_transform('icrs'))
        ax.add_patch(q)
    plt.savefig(directoryp +  "/" + c.strftime('%H%M') + "/calexpsky_" + str(n) + ".png",bbox_inches='tight')


# Use the function

