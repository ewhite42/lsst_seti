## convert Keplerian elements to Cartesian state vector
## using Sorcha's built-in orbit conversion utilities

## 29 Nov 2025, MEW

import numpy as np
import pandas as pd

from sorcha.ephemeris.orbit_conversion_utilities import universal_cartesian

def kep2cart(fname, outfname):

    df = pd.read_csv(fname)
    ObjID = df['ObjID']
    FORMAT = df['FORMAT']
    a = df['a']
    e = df['e']
    inc = np.radians(df['inc'])
    node = np.radians(df['node'])
    argPeri = np.radians(df['argPeri'])
    ma = np.radians(df['ma'])
    epochMJD_TDB = df['epochMJD_TDB']
    
    ## define the gravitational parameter mu
    mu = 2.95912741*10**(-4)

    ## calculate mean daily motion
    n = np.sqrt(mu / (a**3)) ## units - radians per day

    ## calculate the perihelion distance
    q = a*(1-e)
    
    ## calculate the epoch of periapsis, tp
    t0 = epochMJD_TDB
    tp = t0 - ma/n
    
    ## compute the corresponding Cartesian orbital elements
    ## for each object in the provided orbits file
    
    x = []
    y = []
    z = []
    vx = []
    vy = []
    vz = []
    
    for i in range(len(a)):
        xi,yi,zi,vxi,vyi,vzi = universal_cartesian(mu, q[i], e[i], inc[i], node[i], argPeri[i], tp[i], epochMJD_TDB[i])
        x.append(xi)
        y.append(yi)
        z.append(zi)
        
        vx.append(vxi)
        vy.append(vyi)
        vz.append(vzi)
        
    df_cart = pd.DataFrame()
    
    df_cart['ObjID'] = ObjID
    df_cart['FORMAT'] = 'CART'
    
    df_cart['x'] = x
    df_cart['y'] = y
    df_cart['z'] = z
    
    df_cart['xdot'] = vx
    df_cart['ydot'] = vy
    df_cart['zdot'] = vz
    
    df_cart['epochMJD_TDB'] = epochMJD_TDB
    
    df_cart.to_csv(outfname, index=False)
    print("New Cartesian orbits file {} written.".format(outfname))

def main():

    ## define the name of the orbits file in Keplerian format
    kep_fname = "pratchett_orb_kep.csv"
    
    ## define the name of the output file, where the Cartesian
    ## elements will be written
    cart_fname = "pratchett_orb_cart.csv"
    
    ## define the name of the new Cartesian orbits file
    ## you would like to create
    cart_fname = "pratchett_orb_cart.csv"
    
    ## generate the new Cartesian orbits file
    kep2cart(kep_fname,cart_fname)
    
if __name__ == '__main__':
    main()

