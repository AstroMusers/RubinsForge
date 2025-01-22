# Load the weekly Starlink TLE data from Celestrak and save it to a file
from skyfield.api import load
import datetime as dt

c = dt.datetime.now()
date = str(c.strftime("%m.%d"))
dir = "/data/a.saricaoglu/Lumos-Sat/Files/Weekly_Starlink_Archive/"
name = dir + 'starlinks_'+ date + '.txt'
starlinks = load.download('https://celestrak.org/NORAD/elements/gp.php?GROUP=starlink&FORMAT=tle', filename=name)
