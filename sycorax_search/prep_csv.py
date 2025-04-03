## the purpose of this script is to create a csv file with the 
## parameters we are interested in, calculated from the fields
## in the downloaded Rubin CSV file obtained using the query in
## the file, DP03_sycorax_query1.txt

## There is a lot of hard-coding here, so would need to go in and
## tweak this file if you change the query / add more fields. 
## Will try to keep other files as flexible as possible.

## M.E. White, 15 Jan 2025

import numpy as np
import pandas as pd

from astropy import units as u
from astropy.units import cds
from astropy.coordinates import SkyCoord, GCRS, HeliocentricTrueEcliptic

from poliastro.bodies import Sun
from poliastro.twobody import Orbit
from poliastro.twobody.angles import M_to_E, E_to_nu
from poliastro.frames.ecliptic import HeliocentricEclipticJ2000

from velocity_calcs import measured_velocity
from velocity_calcs import keplerian_velocity
from velocity_calcs import keplerian_velocity_xyz
from velocity_calcs import measured_velocity_spherical
from velocity_calcs import keplerian_velocity_spherical
from velocity_calcs import return_M_ap

def main():
    # import original file
    # calculate quantities for new columns
    # create new csv file
    
    #infile_name = '/home/ellie/research/lsst/table_dp03_catalogs_10yr.MPCORB-AS-mpc-JOIN-dp03_c.csv'
    #outfile_name = '/home/ellie/research/lsst/LSST_sim.csv'
    
    infile_name = '/home/ellie/research/lsst/s1003Hmna_data.csv' #
    outfile_name = '/home/ellie/research/lsst/s1003Hmna_data_vel.csv' #
    
    #infile_name = '/home/ellie/research/lsst/S100a6n8a_data.csv' #s1003Hmna_data.csv' #
    #outfile_name = '/home/ellie/research/lsst/S100a6n8a_data_vel.csv' #s1003Hmna_data_vel.csv' #
    
    #infile_name = '/home/ellie/research/lsst/2015_VA53_data.csv' #s1003Hmna_data.csv' #
    #outfile_name = '/home/ellie/research/lsst/2015_VA53_data_vel.csv'
    
    #infile_name = '/home/ellie/research/lsst/LSST_500k_objects.csv'
    #outfile_name = '/home/ellie/research/lsst/LSST_500k_objects_ready.csv'
    
    indata = pd.read_csv(infile_name)
    
    ## add the semi-major axis column
    q = indata['q']  
    e = indata['e']  ## eccentricity
    a = q/(1-e)  ## semi major axis
    #n = indata['n']
    node = indata['node'] * np.pi/180 ## longitude of the ascending node
    peri = indata['peri'] * np.pi/180 ## argument of the perihelion
    i = indata['incl']* np.pi/180
    x = indata['heliocentricX']
    y = indata['heliocentricY']
    
    ## set all z velocities to zero
    #indata['heliocentricVZ'] = 0.00
    
    ## add the semi major axis column
    indata['a'] = a
    
    ## add the sin(i) column
    indata['sini'] = np.sin(indata['incl'] * np.pi/180)
    
    ## add the i - z color column
    indata['i-z'] = indata['i_H'] - indata['z_H']
    
    ## add the a* composite color column
    indata['a*'] = 0.89*(indata['g_H'] - indata['r_H']) + 0.45*(indata['r_H'] - indata['i_H']) - 0.57
    
    ## add the perihelion distance column
    indata['r'] = indata['heliocentricDist']
    
    ## xyz keplerian velocities
    # need to figure out how to calculate nu
    M = indata.apply(return_M_ap, axis=1)
    ecc = u.Quantity(e)
    
    #heliocentric_frame = HeliocentricEclipticJ2000()
    
    E_list = []
    nu_list = []
    x_dotk = []
    y_dotk = []
    z_dotk = []
    
    #print(indata['E'])   
    
    '''for ind, row in indata.iterrows():
        E = (M_to_E(M[ind], ecc[ind]))+ 1.1*np.pi*u.rad
        E_list.append(E)
        nu = E_to_nu(E, ecc[ind])
        nu_list.append(nu)
        
        obj_orbit = Orbit.from_classical(Sun, u.Quantity(a[ind], u.AU), ecc[ind], \
                                       u.Quantity(i[ind], u.rad), u.Quantity(node[ind], u.rad), \
                                       u.Quantity(peri[ind], u.rad), nu, \
                                       u.Quantity(indata['epoch'][ind], cds.MJD))
                                       
        #orb_heliocentric = obj_orbit.transform_to(heliocentric_frame)
        vk_vector = obj_orbit.v
        vx = vk_vector.to(u.AU / u.day)[0]
        vy = vk_vector.to(u.AU / u.day)[1]
        vz = vk_vector.to(u.AU / u.day)[2]
        
        r_vector = obj_orbit.r
        rx = r_vector.to(u.AU)[0]
        ry = r_vector.to(u.AU)[1]
        rz = r_vector.to(u.AU)[2]
        
        x_dotk.append(vx.value)
        y_dotk.append(vy.value)
        z_dotk.append(vz.value)
        #print(vk_vector.to(u.AU / u.day))
        #print(vk_vector.to(u.AU / u.day)[0])
        
    indata['E'] = E_list
    indata['nu'] = nu_list'''
        
    x_dotk, y_dotk, z_dotk = keplerian_velocity_xyz(indata['heliocentricDist'], a, e, i, node, peri, x, y, indata)
    
    indata['x_dotk'] = x_dotk
    indata['y_dotk'] = y_dotk
    indata['z_dotk'] = z_dotk
    
    ## xyz velocity difference columns
    '''indata['delta_x_dot'] = indata['heliocentricVX'] - indata['x_dotk']
    indata['delta_y_dot'] = indata['heliocentricVY'] - indata['y_dotk']
    indata['delta_z_dot'] = indata['heliocentricVZ'] - indata['z_dotk']'''
    
    ## spherical measured and keplerian velocities 
    r_dot, az_dot, el_dot = measured_velocity_spherical(indata['heliocentricVX'], indata['heliocentricVY'], indata['heliocentricVZ'])
    r_dotk, az_dotk, el_dotk = keplerian_velocity_spherical(indata['heliocentricDist'], a, e, i, node, peri, x, y, indata)
    
    indata['r_dot'] = r_dot
    indata['az_dot'] = az_dot
    indata['el_dot'] = el_dot
    
    indata['r_dotk'] = r_dotk
    indata['az_dotk'] = az_dotk
    indata['el_dotk'] = el_dotk
    
    ## differences between keplerian and measured spherical velocities
    indata['delta_r_dot'] = indata['r_dot'] - indata['r_dotk']
    indata['delta_az_dot'] = indata['az_dot'] - indata['az_dotk']
    indata['delta_el_dot'] = indata['el_dot'] - indata['el_dotk']
    
    ## add the 1-dimensional measured velocities column
    vel = measured_velocity(indata['heliocentricVX'], indata['heliocentricVY'], indata['heliocentricVZ'])
    indata['v'] = vel
    
    ## add the 1-dimensional keplerian velocities column
    vk = keplerian_velocity(indata['heliocentricDist'], indata['a'])
    indata['vk'] = vk
    
    ## add the 1-dimensional difference-in-velocities column
    indata['v-vk'] = indata['v'] - indata['vk']
    
    ## add the 1-dimensional velocity difference versus perihelion distance column
    indata['diffv/r'] = indata['v-vk']/indata['heliocentricDist']
    
    ## write to a new csv file
    indata.to_csv(outfile_name) #, index=False) 
    
if __name__ == '__main__':
    main()
