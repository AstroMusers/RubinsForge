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
from matplotlib.ticker import FuncFormatter

from astropy.wcs import WCS
from astropy.wcs import Wcsprm
from astropy.io import fits
from astropy.wcs import utils
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import colors as mcolors
from astropy.visualization import quantity_support, time_support
quantity_support() 

def ra_formatter(x, pos):
    hours = int(x)
    minutes = int((x - hours) * 60)
    seconds = int(((x - hours) * 60 - minutes) * 60)
    return f'{hours:02d}:{minutes:02d}:{seconds:02d}'
####
start_time = time.time()
ts = load.timescale()
starlinks = load.tle_file('https://celestrak.org/NORAD/elements/gp.php?GROUP=starlink&FORMAT=tle')
eph = load('de421.bsp')
earth = eph['earth']
sun = eph['sun']
targets = []

c = dt.datetime.now()

ts = load.timescale()
t = ts.now()
eph = load('de421.bsp')

c = dt.datetime.now()
if not os.path.exists("/data/a.saricaoglu/Lumos-Sat/Files/" +  str(c.strftime("%m.%d"))): 
    directoryf = "/data/a.saricaoglu/Lumos-Sat/Files/"  +  str(c.strftime("%m.%d"))
    os.makedirs(directoryf) 
if not os.path.exists("/data/a.saricaoglu/Lumos-Sat/Plots/" +  str(c.strftime("%m.%d")) + "/" + c.strftime('%H%M')): 
    directoryp = "/data/a.saricaoglu/Lumos-Sat/Plots/"  +  str(c.strftime("%m.%d"))+ "/" + c.strftime('%H%M')
    os.makedirs(directoryp) 
directoryf = "/data/a.saricaoglu/Lumos-Sat/Files/"  +  str(c.strftime("%m.%d"))
directoryp = "/data/a.saricaoglu/Lumos-Sat/Plots/"  +  str(c.strftime("%m.%d"))

with fits.open('/data/a.saricaoglu/lumos-sat/master.fits') as hdul:
    hdul.info()

    data = hdul[1].data
    ra_i = data['RAJ2000_deg']
    dec_i = data['DEJ2000_deg']
    print(ra_i, dec_i)
    # for i in range(0,len(data)):
    #     target_i = Star(Angle(degrees=ra_i[i]),Angle(degrees=dec_i[i]))
    #     targets.append(target_i)
    fig, ax = plt.subplots(figsize=(8, 8))
    for i in range(0,len(data)):
        target_i = sf.positionlib.position_of_radec(ra_i[i], dec_i[i])
        #targets.append(target_i.radec())
        # Observe the target to get accurate coordinates
        ra, dec, distance = target_i.radec()
        #print(ra, dec)
        # Plot the target's RA and Dec
        targets.append([ra,dec])
        ax.scatter(ra.hours, dec.degrees, marker="*", c='b', label='Target')

    # Add custom legend
    legend2 = plt.Line2D([0], [0], marker='*', color='w', markerfacecolor='b', markeredgecolor='k', markersize=3, label='Target')
    #ax.legend(handles=[legend2], loc='lower center', bbox_to_anchor=(0.5, 1.1), ncol=2)
    
    # Set axis labels
    ax.set_xlabel('RA (degrees)')
    ax.set_ylabel('Dec (degrees)')
        # Format RA axis to hh:mm:ss and Dec axis to degrees


    ax.xaxis.set_major_formatter(FuncFormatter(ra_formatter))
    ax.set_xlabel('RA (hh:mm:ss)')

    # Dec is already in degrees, so just label it
    ax.set_ylabel('Dec (degrees)')
    
    # Plot title and adjustments
    #plt.title('Strong Lens Candidate Distribution for ' + str(len(ra_i)) + ' Targets')
    plt.tight_layout()
    plt.savefig(directoryp +  "/" + c.strftime('%H%M') + "/Master_target_positions.png")
    plt.show()
end_time = time.time()
total_seconds = end_time - start_time

# Convert seconds to hours, minutes, and seconds
hours = int(total_seconds // 3600)
minutes = int((total_seconds % 3600) // 60)
seconds = total_seconds % 60
current_time = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime(time.time()))
# Print the formatted runtime
print(f"Current time: {current_time} Runtime: {hours}h {minutes}m {seconds:.2f}s")

plt.figure(figsize=(7, 4), constrained_layout=True)
ax = plt.subplot(111, projection='aitoff')

# add systems
plt.scatter(ra_i, dec_i, label='Candidate', )

plt.grid(True, alpha=0.4)
# plt.xlabel('Galactic Longitude')
# plt.ylabel('Galactic Latitude')
plt.title('Aitoff Projection no wrap', y=1.1)
plt.legend(loc='upper right')

plt.savefig(directoryp +  "/" + c.strftime('%H%M') + "/all_known_and_candidate1.png")
#plt.show()

plt.figure(figsize=(7, 4), constrained_layout=True)
ax = plt.subplot(111)

# add systems
plt.scatter(ra_i, dec_i, label='Candidate')

plt.grid(True, alpha=0.4)
# plt.xlabel('Galactic Longitude')
# plt.ylabel('Galactic Latitude')
plt.title('no Aitoff Projection no wrap', y=1.1)
plt.legend(loc='upper right')

plt.savefig(directoryp +  "/" + c.strftime('%H%M') + "/all_known_and_candidate2.png")
#plt.show()
print(targets)
ra_rad_candidates = [t[0].ra.wrap_at(180*u.deg).radian for t in targets]
dec_rad_candidates = [t[1].dec.radian for t in targets]
plt.figure(figsize=(7, 4), constrained_layout=True)
ax = plt.subplot(111, projection='aitoff')

# add systems
plt.scatter(ra_rad_candidates, dec_rad_candidates, label='Candidate')

plt.grid(True, alpha=0.4)
# plt.xlabel('Galactic Longitude')
# plt.ylabel('Galactic Latitude')
plt.title('Aitoff Projection wrap', y=1.1)
plt.legend(loc='upper right')

plt.savefig(directoryp +  "/" + c.strftime('%H%M') + "/all_known_and_candidate3.png")
#plt.show()

plt.figure(figsize=(7, 4), constrained_layout=True)
ax = plt.subplot(111)

# add systems
plt.scatter(ra_rad_candidates, dec_rad_candidates, label='Candidate')

plt.grid(True, alpha=0.4)
# plt.xlabel('Galactic Longitude')
# plt.ylabel('Galactic Latitude')
plt.title('no Aitoff Projection wrap', y=1.1)
plt.legend(loc='upper right')

plt.savefig(directoryp +  "/" + c.strftime('%H%M') + "/all_known_and_candidate4.png")
#plt.show()