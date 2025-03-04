# Parent version: satmodel_141024 copy 3
import numpy as np
from astropy import units as u
from astropy.io import fits
from astropy.coordinates import SkyCoord
import skyfield as sf
import os 
from skyfield.api import load, wgs84, EarthSatellite, N, W,S,E, Star, Angle
import datetime as dt
import time
from pytz import timezone
from skyfield import almanac
from skyfield.framelib import galactic_frame
from skyfield.data import hipparcos
import lumos.calculator
import numpy as np
import simple_sat
import lumos.plot
import lumos.brdf.library
from astropy.wcs import WCS
from astropy.wcs import Wcsprm
from astropy.io import fits
from astropy.wcs import utils
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import colors as mcolors
from astropy.visualization import quantity_support, time_support
quantity_support() 

start_time = time.time()
ts = load.timescale()
starlinks = load.tle_file('https://celestrak.org/NORAD/elements/gp.php?GROUP=starlink&FORMAT=tle')
event_names = 'rise above 30°', 'culminate', 'set below 30°'
eph = load('de421.bsp')
earth = eph['earth']
sun = eph['sun']
targets = []

c = dt.datetime.now()

wcs_rubin_list = WCS(naxis=2)
wcs_rubin_list.wcs.crpix = [512, 512]
wcs_rubin_list.wcs.crval = [84.52,  -47.59]
wcs_rubin_list.wcs.cunit = ["deg", "deg"]
wcs_rubin_list.wcs.ctype = ["RA---TAN", "DEC--TAN"]
wcs_rubin_list.wcs.cdelt = [-0.0002777777778, 0.0002777777778]
wcs_rubin_list.array_shape = [1024, 1024]
wcs = wcs_rubin_list

if not os.path.exists("/data/a.saricaoglu/Lumos-Sat/Plots/" +  str(c.strftime("%m.%d")) + "/" + c.strftime('%H%M')): 
    directoryp = "/data/a.saricaoglu/Lumos-Sat/Plots/"  +  str(c.strftime("%m.%d"))+ "/" + c.strftime('%H%M')
    os.makedirs(directoryp) 
directoryp = "/data/a.saricaoglu/Lumos-Sat/Plots/"  +  str(c.strftime("%m.%d"))

with fits.open('/grad/a.saricaoglu/lumos-sat/targetlist.fit') as hdul:
    hdul.info()
    targetlist = hdul[1].data
    ax = plt.subplot()

    print(targetlist[1][1])
    for target in targetlist:
        name_i =  target[0]
        ra_i = target[1]
        dec_i = target[2]
        target_i = Star(Angle(degrees=ra_i),Angle(degrees=dec_i))
        targets.append(target_i)
        t = ts.now()
        astrometric = earth.at(t).observe(target_i)
        ra, dec, distance = astrometric.radec()
        print(f'Ra {ra} vs ra_i {ra_i}')
        print(f'Dec {dec} vs dec_i {dec_i}')
        gal = SkyCoord(ra.hours, dec.degrees, frame='icrs', unit=u.deg)
        ax.scatter(gal.ra.wrap_at('180d').radian, gal.dec.radian, marker="*", c = 'b')
    
    targett = Star(Angle(degrees=84.52),Angle(degrees= -47.59))
    ts = load.timescale()
    t = ts.now()
    eph = load('de421.bsp')
    earth = eph['earth']
    astrometric = earth.at(t).observe(targett)
    rat, dect, distance = astrometric.radec()
    #ax.scatter(rat.hours, dect.degrees, marker="o",s=300, edgecolor='red', facecolor='none')
    legend1 = plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='w',markeredgecolor='r', markersize=10, label='WCS reference target')
    legend2 = plt.Line2D([0], [0], marker='*', color='w', markerfacecolor='b',markeredgecolor='k', markersize=10, label='Target')   
    ax.legend(handles=[legend1, legend2],  loc='lower center', bbox_to_anchor=(0.5, 1.1), ncol=3)
    plt.tight_layout()
    plt.title('Strong Lens Candidate Distrubition')
    plt.savefig(directoryp +  "/" + c.strftime('%H%M') + "/Target_positions1" + ".png")
    plt.close()
plt.show()

####
start_time = time.time()
ts = load.timescale()
starlinks = load.tle_file('https://celestrak.org/NORAD/elements/gp.php?GROUP=starlink&FORMAT=tle')
event_names = 'rise above 30°', 'culminate', 'set below 30°'
eph = load('de421.bsp')
earth = eph['earth']
sun = eph['sun']
targets = []

c = dt.datetime.now()

wcs_rubin_list = WCS(naxis=2)
wcs_rubin_list.wcs.crpix = [512, 512]
wcs_rubin_list.wcs.crval = [84.52,  -47.59]
wcs_rubin_list.wcs.cunit = ["deg", "deg"]
wcs_rubin_list.wcs.ctype = ["RA---TAN", "DEC--TAN"]
wcs_rubin_list.wcs.cdelt = [-0.0002777777778, 0.0002777777778]
wcs_rubin_list.array_shape = [1024, 1024]
wcs = wcs_rubin_list
target = Star(Angle(degrees=84.52),Angle(degrees= -47.59))
ts = load.timescale()
t = ts.now()
eph = load('de421.bsp')
earth = eph['earth']
astrometric = earth.at(t).observe(target)
rat, dect, distance = astrometric.radec()

if not os.path.exists("/data/a.saricaoglu/Lumos-Sat/Plots/" +  str(c.strftime("%m.%d")) + "/" + c.strftime('%H%M')): 
    directoryp = "/data/a.saricaoglu/Lumos-Sat/Plots/"  +  str(c.strftime("%m.%d"))+ "/" + c.strftime('%H%M')
    os.makedirs(directoryp) 
directoryp = "/data/a.saricaoglu/Lumos-Sat/Plots/"  +  str(c.strftime("%m.%d"))

with fits.open('/grad/a.saricaoglu/lumos-sat/targetlist.fit') as hdul:
    hdul.info()
    targetlist = hdul[1].data
    ax = plt.subplot(projection=wcs)
    pixel_scale = utils.proj_plane_pixel_scales(wcs)
    print(pixel_scale)
    print(targetlist[1][1])
    for target in targetlist:
        try:
            name_i =  target[0]
            ra_i = target[1]
            dec_i = target[2]
            target_i = Star(Angle(degrees=ra_i),Angle(degrees=dec_i))
            targets.append(target_i)
            t = ts.now()
            astrometric = earth.at(t).observe(target_i)
            ra, dec, distance = astrometric.radec()
            print(f'Ra {ra} vs ra_i {ra_i}')
            print(f'Dec {dec} vs dec_i {dec_i}')

            coord = SkyCoord(ra=ra_i*u.deg, dec=dec_i*u.deg) #overplotting with an object in the field of view
            target = utils.skycoord_to_pixel(coord, wcs)
            ax.scatter(target[0]*pixel_scale[0],target[1]*pixel_scale[1], facecolors='none', marker='*', edgecolors='b',transform=ax.get_transform('pixel'))
        except:
            print(f'Problem with {ra_i} , {dec_i}')
        # ax.scatter(ra.hours, dec.degrees, marker="*", c = 'b', trans)
    coord = SkyCoord(ra=84.52*u.deg, dec=-47.59*u.deg) #overplotting with an object in the field of view
    target = utils.skycoord_to_pixel(coord, wcs)
    ax.scatter(target[0]*pixel_scale[0],  target[1]*pixel_scale[1], marker="o",s=300, edgecolor='red', facecolor='none',transform=ax.get_transform('pixel'))
    # ax.scatter(target[0]*pixel_scale[0],  target[1]*pixel_scale[1], marker="o",s=300, edgecolor='green',facecolor='none')
  
    legend1 = plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='w',markeredgecolor='r', markersize=10, label='WCS reference target')
    legend2 = plt.Line2D([0], [0], marker='*', color='w', markerfacecolor='b',markeredgecolor='k', markersize=10, label='Target')   
    ax.legend(handles=[legend1, legend2],  loc='lower center', bbox_to_anchor=(0.5, 1.1), ncol=3)
    ax.coords.grid(color='black', alpha=0.5, linestyle='-')
    plt.tight_layout()
    plt.title('Strong Lens Candidate Distrubition')
    plt.savefig(directoryp +  "/" + c.strftime('%H%M') + "/Target_positions2" + ".png")
    plt.close()
plt.show()
