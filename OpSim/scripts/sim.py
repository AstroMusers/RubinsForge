import numpy as np
import rubin_sim.maf as maf
import rubin_sim.sim_archive as arch
import rubin_sim.data as data
import os

 

print(os.getcwd())
opsdb = data.get_baseline()
print(data.get_data_dir())
print(opsdb)
row = maf.db.DisplayRow()
print(row)
# Pull columns you need
vis = opsdb.fetchMetricData(['observationStartMJD','fieldRA','fieldDec','filter',
                             'fiveSigmaDepth','seeingFwhmEff','airmass','rotSkyPos'])

# select visits within r degrees of (ra0, dec0)
def visits_near(vis, ra0, dec0, radius_deg=0.8):
    m = (np.abs(vis['fieldRA']-ra0) < radius_deg) & (np.abs(vis['fieldDec']-dec0) < radius_deg)
    return vis[m]

my_visits = visits_near(vis, ra0=150.0, dec0=2.0)
