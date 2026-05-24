import os
import sys
import numpy as np
import pandas as pd
import pylab as plt

import subprocess
from scipy.interpolate import interp1d
from astroquery.jplhorizons import Horizons
from astroquery.jplsbdb import SBDB

## retrieve JPL Horizons data
jpl_path = '/home/ellie/research/lsst/lsst_seti/nongrav/sorcha_ng'

## import the comparison module
jpl_dir = os.path.abspath(jpl_path)
sys.path.insert(0, jpl_dir)

from compare_sorcha_output import get_jpl_output as gjo

startday = '2025-06-20'
endday = '2025-7-24'
step = '1d'

mjd_jpl, ra_jpl, dec_jpl = gjo('3I', startday, endday, step) #'3I', startday, endday, step)

## read in the gravity-only data: 

fpath = '/home/ellie/research/lsst/sorcha_output/3iatlas/3iatlas.csv'

df = pd.read_csv(fpath) #sfpath+obj_id+'.csv')
df_grav = df[df['ObjID'] == 'ATLAS'] #'ATLAS']

ra_grav = df_grav['RA_deg']
mjd_grav = df_grav['fieldMJD_TAI']

df_ng = df[df['ObjID'] == 'ATLAS_NG'] #'ATLAS_NG']

ra_ng = df_ng['RA_deg']
mjd_ng = df_ng['fieldMJD_TAI']

df_ng_1 = df[df['ObjID'] == 'ATLAS_NG_1'] 
ra_ng_1 = df_ng_1['RA_deg']
mjd_ng_1 = df_ng_1['fieldMJD_TAI']

plt.plot(mjd_grav, ra_grav, label='grav only')
plt.plot(mjd_ng, ra_ng, label='nongravs')
plt.plot(mjd_jpl, ra_jpl, label='jpl')
plt.plot(mjd_ng_1, ra_ng_1, label='bigger nongravs')
plt.legend()
plt.show()
