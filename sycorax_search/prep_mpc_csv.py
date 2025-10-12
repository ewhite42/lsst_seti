## the purpose of this script is to create a csv file with the 
## parameters we are interested in, calculated from the fields
## in the prepared MPC CSV file obtained using the code in 
## mpc_preprocessing.py and mpc_analysis.ipynb

## There is a lot of hard-coding here, so would need to go in and
## tweak this file if you change the query / add more fields. 
## Will try to keep other files as flexible as possible.

## MEW, 7 July 2025

import numpy as np
import pylab as plt
import pandas as pd

import time

from astropy import units as u
from astropy.units import cds
from astropy.coordinates import SkyCoord, GCRS, HeliocentricTrueEcliptic, ITRS, GeocentricTrueEcliptic, AltAz
from astropy.coordinates import CartesianRepresentation, CartesianDifferential
from astropy.coordinates import EarthLocation
from astropy.time import Time

from poliastro.bodies import Sun
from poliastro.twobody import Orbit
from poliastro.twobody.angles import M_to_E, E_to_nu, nu_to_E
from poliastro.frames.ecliptic import HeliocentricEclipticJ2000

from velocity_calcs import measured_velocity
from velocity_calcs import keplerian_velocity
from velocity_calcs import keplerian_velocity_xyz
from velocity_calcs import measured_velocity_spherical
from velocity_calcs import keplerian_velocity_spherical
from velocity_calcs import return_nu_ap_mpc

def make_csv(indata):
    # import original file
    # calculate quantities for new columns
    # create new csv file
    
    #infile_name = '/home/ellie/research/lsst/mpcorb_extended_complete.csv'
    #outfile_name = '/home/ellie/research/lsst/mpcorb_extended_ready.csv'
    
    ## create variables
    a = indata['a']
    e = indata['e']
    ap = a*(1+e)
    
    ## obtain inclination, node, and peri in radians
    node = indata['Node'] * np.pi/180 ## longitude of the ascending node
    peri = indata['Peri'] * np.pi/180 ## argument of the perihelion
    i = indata['i']* np.pi/180
    
    x = indata['X']
    y = indata['Y']
    z = indata['Z']
    
    ## add the aphelion column
    indata['ap'] = ap
    
    ## add the sin(i) column
    indata['sini'] = np.sin(indata['i'] * np.pi/180)
        
    ## xyz keplerian velocities
    # need to figure out how to calculate nu
    heliocentricDist = np.sqrt(x**2 + y**2 + z**2)
    indata['heliocentricDist'] = heliocentricDist
    
    ## xyz keplerian velocities
    # need to figure out how to calculate nu
    nu = indata.apply(return_nu_ap_mpc, axis=1)
    ecc = u.Quantity(e)
    
    E_list = []
    nu_list = []
    x_dotk = []
    y_dotk = []
    z_dotk = []
    
    poli_elapsed = 0
    
    #print(indata['E'])   
    eps = np.radians(23.43928)
    
    for ind, row in indata.iterrows():
    
        E = (nu_to_E(nu[ind], ecc[ind])) #+ 1.1*np.pi*u.rad
        E_list.append(E)
        #nu = E_to_nu(E, ecc[ind])
        #nu_list.append(nu)
        
        st_poli = time.time()
        
        obj_orbit = Orbit.from_classical(Sun, u.Quantity(a[ind], u.AU), ecc[ind], \
                                       u.Quantity(i[ind], u.rad), u.Quantity(node[ind], u.rad), \
                                       u.Quantity(peri[ind], u.rad), nu[ind], \
                                       u.Quantity(indata['Epoch'][ind], cds.MJD), plane='EARTH_ECLIPTIC')
                                       
                                       
        #orb_heliocentric = obj_orbit.transform_to(heliocentric_frame)
        vk_vector = obj_orbit.v
        vx = vk_vector.to(u.AU / u.day)[0]
        vy = vk_vector.to(u.AU / u.day)[1]
        vz = vk_vector.to(u.AU / u.day)[2]
        
        x_dotk_i = vx.value
        y_dotk_i = np.cos(eps)*vy.value - np.sin(eps)*vz.value
        z_dotk_i = np.sin(eps)*vy.value + np.cos(eps)*vz.value
        
        et_poli = time.time()
        poli_delta = et_poli - st_poli
        
        poli_elapsed += poli_delta
        
        r_vector = obj_orbit.r
        rx = r_vector.to(u.AU)[0]
        ry = r_vector.to(u.AU)[1]
        rz = r_vector.to(u.AU)[2]
        
        x_dotk.append(vx.value)
        y_dotk.append(np.cos(eps)*vy.value - np.sin(eps)*vz.value) #y_dotk.append(vy.value) 
        z_dotk.append(np.sin(eps)*vy.value + np.cos(eps)*vz.value)  #vz.value)
        #print(vk_vector.to(u.AU / u.day))
        #print(vk_vector.to(u.AU / u.day)[0])
        
    #indata['E'] = E_list
    #indata['nu'] = nu

    print(f"Elapsed time (poliastro): {poli_elapsed:.6f} seconds")
    
    st_mew = time.time()
            
    x_dotk, y_dotk, z_dotk = keplerian_velocity_xyz(indata['heliocentricDist'], a, e, i, node, peri, x, y, indata, mpc=True)
    
    et_mew = time.time()
    
    mew_elapsed = et_mew - st_mew
    
    print(f"Elapsed time (our calcs): {mew_elapsed:.6f} seconds")
    
    indata['x_dotk'] = x_dotk 
    indata['y_dotk'] = y_dotk 
    indata['z_dotk'] = z_dotk
    
    ## xyz velocity difference columns
    indata['delta_x_dot'] = indata['x_dotk']/indata["X'"]
    indata['delta_y_dot'] = indata['y_dotk']/indata["Y'"]
    indata['delta_z_dot'] = indata['z_dotk']/indata["Z'"]
    
    ## spherical measured and keplerian velocities 
    r_dot, az_dot, el_dot = measured_velocity_spherical(indata["X'"], indata["Y'"], indata["Z'"])
    r_dotk, az_dotk, el_dotk = keplerian_velocity_spherical(indata['heliocentricDist'], a, e, i, node, peri, x, y, indata, mpc=True)
    
    indata['r_dot'] = r_dot
    indata['az_dot'] = az_dot
    indata['el_dot'] = el_dot
    
    indata['r_dotk'] = r_dotk
    indata['az_dotk'] = az_dotk
    indata['el_dotk'] = el_dotk
    
    ## differences between keplerian and measured spherical velocities
    indata['delta_r_dot'] = indata['r_dotk']/indata['r_dot']
    indata['delta_az_dot'] = indata['az_dotk']/indata['az_dot']
    indata['delta_el_dot'] = indata['el_dotk']/indata['el_dot']

    return indata, [mew_elapsed, poli_elapsed]
    
def efficiency_comparison():
    infile_name = '/home/ellie/research/lsst/poliastro_comparison.csv' #mpcorb_extended_complete.csv' #
    #outfile_name = '/home/ellie/research/lsst/poliastro_comparison.csv' #
    
    #df = pd.read_csv(infile_name)

    df_indata = pd.read_csv(infile_name)
    '''effs_mew = []
    effs_poli = []
    nsamp = []
    
    df_1 = df_indata[:1]
    df_results1, effsi = make_csv(df_1)
    
    for i in range(6):
        end = 10**i + 1
        df_i = df_indata[:end]
        df_resultsi, effsi = make_csv(df_i)
        
        effs_mew.append(effsi[0])
        effs_poli.append(effsi[1])
        nsamp.append(10**i)
    
    df = pd.DataFrame()
    
    df['nsamp'] = nsamp
    df['mew_elapsed'] = effs_mew
    df['poli_elapsed'] = effs_poli
    
    df['mew_rate'] = df['mew_elapsed']/df['nsamp']
    df['poli_rate'] = df['poli_elapsed']/df['nsamp']
    
    df.to_csv(outfile_name)'''
    
    nsamp = df_indata['nsamp']
    effs_mew = df_indata['mew_elapsed']
    effs_poli = df_indata['poli_elapsed']
    
    plt.semilogy(nsamp, effs_mew, marker='.', label='velocity_calcs')
    plt.semilogy(nsamp, effs_poli, marker='.', label='poliastro')
    plt.title('Velocity Calculation Efficiency Comparison')
    plt.legend()

    plt.xlabel('Number of samples')
    plt.ylabel('Time taken to compute velocity vector (s)')
    plt.tight_layout()
    plt.savefig('/home/ellie/research/lsst/yvelocity_eff_comparison.png')
    plt.show() 
    
    plt.plot(df['nsamp'], df['mew_rate'])
    plt.plot(df['nsamp'], df['poli_elapsed'])
    plt.show()
    
def main():
    infile_name = '/home/ellie/research/lsst/mpcorb_ceres.csv' #mpcorb_extended_complete.csv' #
    outfile_name = '/home/ellie/research/lsst/mpcorb_ceres_vel_22aug_poli.csv' #mpcorb_extended_complete_vel_2sept.csv' #

    df_indata = pd.read_csv(infile_name)
 
    df_results = make_csv(df_indata)
    df_results.to_csv(outfile_name)

    '''num_rows = len(df_indata)
    iter_size = 1000000
    num_iter = int(num_rows/iter_size)+1
    start_idx = 0
    
    for i in range(num_iter):
        end_idx = start_idx+iter_size
        if end_idx > num_rows:
            df = df_indata.loc[start_idx:num_rows-1].copy()
            df_result_i = main(df)
            df_results = pd.concat([df_results, df_result_i], ignore_index=True)
            start_idx = end_idx
            print("Reached the end")
            break
    
        df = df_indata.loc[start_idx:end_idx].copy()
        df_result_i = main(df)
        df_results = pd.concat([df_results, df_result_i], ignore_index=True)
    
        start_idx = end_idx
        print("just finished the {}th iteration".format(i))

    ## write to a new csv file
    df_results.to_csv(outfile_name)'''
    
if __name__ == '__main__':
    #main()
    efficiency_comparison()
