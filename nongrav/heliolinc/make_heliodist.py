## create heliocentric hypothesis files 
## for HelioLinC

## MEW 1 Nov 2025

import numpy as np
import pandas as pd

G = 2.95912741*10**(-4) #G in au*(au/day)^2/M_sun units -- this is the one we use

def make_heliodist(fname):

    delta = 4 ## AU
    num_dists = 100
    num_vels = 10
    outfname = 'heliohypo_atlas.csv'
    acc = [0.00, 0.20, 0.40, 0.60, 0.80]

    df = pd.read_csv(fname)
    
    x = df['x']
    y = df['y']
    z = df['z']
    
    hdist = np.sqrt((x[0]**2) + (y[0]**2) + (z[0]**2))
    hdist_min = hdist - (delta/2)
    hdist_max = hdist + (delta/2)
    
    hdist_arr = np.linspace(hdist_min, hdist_max, num_dists)
    
    vx = df['xdot']
    vy = df['ydot']
    vz = df['zdot']
    
    vel_max = np.sqrt((vx[0]**2) + (vy[0]**2) + (vz[0]**2))
    #vel_max = np.sqrt(2*G/hdist)
    vel_min = -vel_max
    
    vel_arr = np.linspace(vel_min, vel_max, num_vels)
    
    hlist = []
    vlist = []
    alist = []
    norms = []
    
    for d in range(len(hdist_arr)):
        for v in range(len(vel_arr)):
            for a in acc:
                hlist.append(hdist_arr[d])
                vlist.append(vel_arr[v])
                alist.append(a)
                norms.append(v+1)
    
    df_new = pd.DataFrame()
    
    df_new['#r(AU)'] = hlist #hdist_arr
    df_new['rdot(AU/day)'] = vlist #[0.0]*num_dists
    df_new['norm'] = norms #[1]*len(df_new['#r(AU)'])
    df_new['mean_accel'] = alist
    
    df_new.to_csv(outfname, index=False, sep=' ')
    
def main():

    fname = 'c2017m4_atlas_orb_nongrav.csv'
    make_heliodist(fname)
    
if __name__ == '__main__':
    main()
