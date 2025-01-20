## the functions here will allow us to calculate the 
## measured velocity and predicted Keplerian velocity
## of a given object based on parameters in the LSST
## simulation data

## M. E. White, 15 Jan 2025

import numpy as np

## define our constants
G_au_yr = 39.422888 #G in au*(au/yr)^2/M_sun units
G_km_s = 887.1287008 #G in au*(km/s)^2/M_sun units
G_au_s = 3.96402313*10**(-14) #G in au*(au/s)^2/M_sun units
G_au_d = 2.95912741*10**(-4) #G in au*(au/day)^2/M_sun units -- this is the one we use

def measured_velocity(vx, vy, vz):
    # this function calculates the total velocity from the x,y,z 
    # heliocentric velocities from the SSSource table
    
    # not totally sure what units these are in -- data product
    # document says that the heliocentric velocities are in "AU" 
    # whatever that means... I am thinking they are in AU/day, based on 
    # comparing Keplerian velocities calculated for various units, bc
    # AU/day agreed the best.
    
    return np.sqrt(vx**2 + vy**2 + vz**2)
    
def keplerian_velocity(r, a):
    # this function allows us to calculate the predicted,
    # Keplerian velocity of an object using: 
    # v = sqrt(GM(2/r - 1/a))
    
    # note r = distance from Sun to object in AU
    # and a = semi-major axis in AU
    
    # note also that we are assuming M = M_sun + M_object ~ M_sun, 
    # and since M_sun is the unit of mass in the G definition
    # we are omitting it here. 

    # this returns the Keplerian velocity in AU/day
    return np.sqrt((G_au_d)*((2/r)-(1/a)))
