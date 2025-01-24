# Parent : targets_2.py
# Creates a mock target distribution with a population of 10k sources above 30 degrees.
import numpy as np
import os 
import sys
from astropy.io import fits
from skyfield.api import load, Star, Angle
import datetime as dt
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from astropy.table import Table
from astropy import units as u
import matplotlib.ticker as ticker
import time

start_time = time.time()

c = dt.datetime.now()
if not os.path.exists("/data/a.saricaoglu/Lumos-Sat/Files/" +  str(c.strftime("%m.%d"))): 
    directoryf = "/data/a.saricaoglu/Lumos-Sat/Files/"  +  str(c.strftime("%m.%d"))
    os.makedirs(directoryf) 
if not os.path.exists("/data/a.saricaoglu/Lumos-Sat/Plots/" +  str(c.strftime("%m.%d")) + "/" + c.strftime('%H%M')): 
    directoryp = "/data/a.saricaoglu/Lumos-Sat/Plots/"  +  str(c.strftime("%m.%d"))+ "/" + c.strftime('%H%M')
    os.makedirs(directoryp) 
directoryf = "/data/a.saricaoglu/Lumos-Sat/Files/"  +  str(c.strftime("%m.%d"))
directoryp = "/data/a.saricaoglu/Lumos-Sat/Plots/"  +  str(c.strftime("%m.%d"))

# Load necessary data
ts = load.timescale()
eph = load('de421.bsp')
earth = eph['earth']



# Parameters for the RA and Dec distribution
num_sources = 10000  # Number of random sources you want to generate
ra_min, ra_max = 0, 360         # RA ranges from 0 to 360 degrees
dec_min, dec_max = -80,30       # Dec ranges from 30 to 90 degrees

# Generate random RA and Dec values
ra_values = np.random.uniform(ra_min, ra_max, num_sources)  # Random RA values
dec_values = np.random.uniform(dec_min, dec_max, num_sources)  # Random Dec values

# Create a table with the data
data_table = Table([ra_values, dec_values], names=('RA', 'Dec'))

# Define the FITS file and save the data
fits_filename = directoryf + "/mock_targets_list_" + str(num_sources) + ".fit"
data_table.write(fits_filename, format='fits', overwrite=True)

print(f"FITS file '{fits_filename}' created with {num_sources} random RA and Dec values.")


# Load the FITS file data (update the path as per your file location)
with fits.open(fits_filename) as hdul:
    data = hdul[1].data
    hdul.info()
    ra_values = data['RA']
    dec_values = data['DEC']
    
    fig, ax = plt.subplots(figsize=(8, 8))

    # Loop through each target and plot it
    for i in range(0,num_sources):
        ra_i = ra_values[i]
        dec_i = dec_values[i]
        target_i = Star(Angle(degrees=ra_i), Angle(degrees=dec_i))
        
        # Observe the target to get accurate coordinates
        t = ts.now()
        astrometric = earth.at(t).observe(target_i)
        ra, dec, distance = astrometric.radec()
        print(ra, dec)
        # Plot the target's RA and Dec
        ax.scatter(ra.hours, dec.degrees, marker="*", c='b', label='Target')

    # Add custom legend
    legend2 = plt.Line2D([0], [0], marker='*', color='w', markerfacecolor='b', markeredgecolor='k', markersize=10, label='Target')
    ax.legend(handles=[legend2], loc='lower center', bbox_to_anchor=(0.5, 1.1), ncol=2)
    
    # Set axis labels
    ax.set_xlabel('RA (degrees)')
    ax.set_ylabel('Dec (degrees)')
        # Format RA axis to hh:mm:ss and Dec axis to degrees
    def ra_formatter(x, pos):
        hours = int(x)
        minutes = int((x - hours) * 60)
        seconds = int(((x - hours) * 60 - minutes) * 60)
        return f'{hours:02d}:{minutes:02d}:{seconds:02d}'

    ax.xaxis.set_major_formatter(FuncFormatter(ra_formatter))
    ax.set_xlabel('RA (hh:mm:ss)')

    # Dec is already in degrees, so just label it
    ax.set_ylabel('Dec (degrees)')
    
    # Plot title and adjustments
    plt.title('Mock Strong Lens Candidate Distribution for ' + str(num_sources) + ' Targets')
    plt.tight_layout()
    plt.savefig(directoryp +  "/" + c.strftime('%H%M') + "/Random_target_positions.png")
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