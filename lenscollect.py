
import numpy as np
from matplotlib.ticker import FuncFormatter

from astropy.wcs import WCS
from astropy.wcs import Wcsprm
from astropy.io import fits
from astropy.wcs import utils
from astropy.coordinates import SkyCoord
from astropy import units as u
import logging
import os
import pandas as pd

import sys
import datetime as dt
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import colors as mcolors
from astropy.visualization import quantity_support, time_support
quantity_support() 

c = dt.datetime.now()

# Get the script name
script_name = os.path.basename(__file__)
# Configure logging
log_filename = f"/data/a.saricaoglu/repo/RubinsForge/LenSim/logs/{c.strftime('%m.%d')}/{c.strftime('%H%M')}/script.log"
os.makedirs(os.path.dirname(log_filename), exist_ok=True)
logging.basicConfig(filename=log_filename, level=logging.INFO, 
                    format='%(asctime)s - %(levelname)s - %(message)s')

# Redirect stdout and stderr to the log file
class StreamToLogger:
    def __init__(self, logger, log_level):
        self.logger = logger
        self.log_level = log_level
        self.linebuf = ''

    def write(self, buf):
        for line in buf.rstrip().splitlines():
            self.logger.log(self.log_level, line.rstrip())

    def flush(self):
        pass

# sys.stdout = StreamToLogger(logging.getLogger('STDOUT'), logging.INFO)
with fits.open('/grad/a.saricaoglu/Downloads/master.fits') as hdul:
    print(hdul[2].data)
    print(hdul[2].columns)

    data = hdul[2].data
    # Example input: lens catalog with RA/DEC

    radius_deg = 15.0  # ~2 arcseconds search radius

    sql_queries = []

    coord = SkyCoord(data['RAJ2000'], data['DEJ2000'], unit=(u.hourangle, u.hourangle))
    ra = coord.ra.deg
    dec = coord.dec.deg

    for rai, deci in zip(ra, dec):
        query = f"""SELECT \n   object_id, ra, dec, g_cmodel_flux, g_cmodel_fluxerr \n   FROM \n   pdr3_wide.forced \n   WHERE \n   isprimary \n   AND conesearch(coord, {rai}, {deci}, {radius_deg}) \n   AND NOT g_cmodel_flag;"""
        sql_queries.append(query)

    # Write all queries into a text file
    with open("hsc_lens_queries.sql", "w") as f:
        for q in sql_queries:
            f.write(q + "\n\n")

    print("SQL queries generated and saved as hsc_lens_queries.sql")

# with fits.open('/grad/a.saricaoglu/Downloads/773529.fits') as hdul:
#     print(hdul[1].data)
#     print(hdul[1].columns)
#     hdul.info()
import json

# Specify the path to the .json file
file_path = '/grad/a.saricaoglu/Downloads/lenses.json'

# Open and read the .json file
with open(file_path, 'r') as json_file:
    data = json.load(json_file)  # Load the JSON data into a Python dictionary

# Print the data to verify
print(data)

# Access specific keys in the JSON data
if isinstance(data, dict):  # Ensure it's a dictionary
    for key, value in data.items():
        if isinstance(value, dict) and value.get('instr') == 'Legacy Survey (DR10)':
            print(f"Lens ID: {key}, Details: {value}")