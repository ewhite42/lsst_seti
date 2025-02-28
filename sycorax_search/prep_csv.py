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
from velocity_calcs import measured_velocity
from velocity_calcs import keplerian_velocity
from velocity_calcs import keplerian_velocity_xyz
from velocity_calcs import measured_velocity_spherical
from velocity_calcs import keplerian_velocity_spherical

def main():
    # import original file
    # calculate quantities for new columns
    # create new csv file
    
    infile_name = '/home/ellie/research/lsst/table_dp03_catalogs_10yr.MPCORB-AS-mpc-JOIN-dp03_c.csv'
    outfile_name = '/home/ellie/research/lsst/LSST_sim.csv'
    
    indata = pd.read_csv(infile_name)
    
    ## add the semi-major axis column
    q = indata['q']  
    e = indata['e']  ## eccentricity
    a = q/(1-e)  ## semi major axis
    n = indata['n']
    
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
    x_dotk, y_dotk, z_dotk = keplerian_velocity_xyz(indata['heliocentricDist'], a, e, n, i)
    
    indata['x_dotk'] = x_dotk
    indata['y_dotk'] = y_dotk
    indata['z_dotk'] = z_dotk
    
    ## xyz velocity difference columns
    indata['delta_x_dot'] = indata['heliocentricVX'] - indata['x_dotk']
    indata['delta_y_dot'] = indata['heliocentricVY'] - indata['y_dotk']
    indata['delta_z_dot'] = indata['heliocentricVZ'] - indata['z_dotk']
    
    ## spherical measured and keplerian velocities 
    r_dot, az_dot, el_dot = measured_velocity_spherical(indata['heliocentricVX'], indata['heliocentricVY'], indata['heliocentricVZ'])
    r_dotk, az_dotk, el_dotk = keplerian_velocity_spherical(indata['heliocentricDist'], a, e, n, i)
    
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
