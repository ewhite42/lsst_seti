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

def main():
    # import original file
    # calculate quantities for new columns
    # create new csv file
    
    infile_name = 'table_dp03_catalogs_10yr.MPCORB-data.csv'
    indata = pd.read_csv(infile_name)
    
    ## add the semi-major axis column
    q = indata['q']
    e = indata['e']
    a = q/(1-e)
    indata['a'] = a
    
    ## add the sin(i) column
    indata['sini'] = np.sin(indata['incl'] * np.pi/180)
    
    ## add the i - z color column
    indata['i-z'] = indata['i_H'] - indata['z_H']
    
    ## add the a* composite color column
    indata['a*'] = 0.89*(indata['g_H'] - indata['r_H']) + 0.45*(indata['r_H'] - indata['i_H']) - 0.57

    ## add the measured velocities column
    vel = measured_velocity(indata['heliocentricVX'], indata['heliocentricVY'], indata['heliocentricVZ'])
    indata['v'] = vel
    
    ## add the keplerian velocities column
    vk = keplerian_velocity(indata['heliocentricDist'], indata['a'])
    indata['vk'] = vk
    
    ## add the difference-in-velocities column
    indata['v-vk'] = indata['v'] - indata['vk']
    
    ## add the velocity difference versus perihelion distance column
    indata['diffv/r'] = indata['v-vk']/indata['heliocentricDist']
    
    ## write to a new csv file
    indata.to_csv('LSST_sim.csv', index=False) 
    
if __name__ == '__main__':
    main()
