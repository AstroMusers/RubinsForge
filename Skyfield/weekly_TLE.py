# Load the weekly Starlink TLE data from Celestrak and save it to a file
from skyfield.api import load
import datetime as dt
import os
import glob

def download_new_tle(tle_dir):
    """
    Download new TLE data with year in filename
    """
    c = dt.datetime.now()
    date = str(c.strftime("%m.%d.%y"))  # Include year as YY
    name = os.path.join(tle_dir, f'starlinks_{date}.txt')
    
    print(f"Downloading new TLE data to: {name}")
    try:
        starlinks = load.download('https://celestrak.org/NORAD/elements/gp.php?GROUP=starlink&FORMAT=tle', filename=name)
        print(f"Successfully downloaded TLE data to {name}")
        return name
    except Exception as e:
        print(f"Error downloading TLE data: {e}")
        return None

if __name__ == "__main__":
    # Directory containing TLE files
    tle_dir = "/data/a.saricaoglu/repo/RubinsForge/Skyfield/tle_data/"

    
    print("\n=== Downloading new TLE data ===")
    download_new_tle(tle_dir)
    
 
