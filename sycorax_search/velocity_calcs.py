## the functions here will allow us to calculate the 
## measured velocity and predicted Keplerian velocity
## of a given object based on parameters in the LSST
## simulation data

## M. E. White, 15 Jan 2025

import numpy as np
from scipy.optimize import fsolve

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
    

def measured_velocity_spherical(vx, vy, vz):
    # this function calculates the total velocity from the x,y,z 
    # heliocentric velocities from the SSSource table
    
    # not totally sure what units these are in -- data product
    # document says that the heliocentric velocities are in "AU" 
    # whatever that means... I am thinking they are in AU/day, based on 
    # comparing Keplerian velocities calculated for various units, bc
    # AU/day agreed the best.
    
    r_dot = np.sqrt(vx**2 + vy**2 + vz**2)
    az_dot = np.arctan(vy/vx)
    el_dot = np.arccos(vz/r_dot)
    
    return r_dot, az_dot, el_dot
    
def kepler_equation(E, M, e):
    ## define Kepler's equation to find
    ## the eccentric anomaly, E
    
    return (M-E+e*np.sin(E))
    
def return_M(row):
    #return np.arccos(row['heliocentricX']/row['heliocentricDist']) + ((1.5*np.pi) if row['heliocentricY'] < 0 else 0)
    if row['heliocentricY'] < 0:
        return np.arccos(row['heliocentricX']/row['heliocentricDist']) 
        
    else: 
        return (2*np.pi) - np.arccos(row['heliocentricX']/row['heliocentricDist'])

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
    

def keplerian_velocity_xyz(r, a, e, i, node, peri, x, y, df): 

    # r = heliocentric distance (~ distance from 1 focus)
    # a = semi-major axis
    # e = eccentricity
    # n = mean daily motion
    # i = inclination
    
    # n = mean motion
    # mean motion is the average angular distance traveled per day
    
    n = np.sqrt(G_au_d / (a**3))
    
    ## mean anomaly
    M = df.apply(return_M, axis=1)

    #M = np.arccos(x/r)
    
    ## calculate eccentric anomaly
    guess_E = M
    solution_E = fsolve(kepler_equation, M, args=(M,e))
    #print(solution_E)
    E = solution_E
    
    #print((a*e + x)/a)
    #E = np.arccos((a*e + x)/a)#np.arccos((1 - r/a)/e) ## calculate eccentric anomaly

    '''x_dot = -(a*n*np.sin(E)*np.cos(i))/(1-e*np.cos(E))
    y_dot = (a*n*np.sqrt(1-(e**2))*np.cos(E))/(1-(e*np.cos(E)))
    z_dot = (a*n*np.sin(E)*np.sin(i))/(1-(e*np.cos(E)))'''
    
    ## unrotated coordinates and vector
    X_dot = -(a*n*np.sin(E))/(1-e*np.cos(E))
    Y_dot = (a*n*np.sqrt(1-(e**2))*np.cos(E))/(1-(e*np.cos(E)))
    
    v_XYZ = [X_dot, Y_dot, 0]

    ## rotation matrix
    P11 = np.cos(node)*np.cos(peri)-np.sin(node)*np.cos(i)*np.sin(peri)
    P12 = -np.cos(node)*np.sin(peri)-np.sin(node)*np.cos(i)*np.cos(peri)
    P13 = np.sin(node)*np.sin(i)
    
    P21 = np.sin(node)*np.cos(peri)+np.cos(node)*np.cos(i)*np.sin(peri)
    P22 = -np.sin(node)*np.sin(peri)+np.cos(node)*np.cos(i)*np.cos(peri)
    P23 = -np.cos(node)*np.sin(i)
    
    P31 = np.sin(i)*np.sin(peri)
    P32 = np.sin(i)*np.cos(peri)
    P33 = np.cos(i)
    
    P = np.array([[P11, P12, P13], [P21, P22, P23], [P31, P32, P33]])
    
    ## numpy matrix multiplication throwing errors, too tired to debug it, 
    ## so writing it out by hand for now...
    #v_xyz = np.dot(P, v_XYZ)
    
    v_xyz = [P11*v_XYZ[0]+P12*v_XYZ[1]+P13*v_XYZ[2], P21*v_XYZ[0]+P22*v_XYZ[1]+P23*v_XYZ[2], P31*v_XYZ[0]+P32*v_XYZ[1]+P33*v_XYZ[2]]
    
    x_dot = v_xyz[0]
    y_dot = v_xyz[1]
    z_dot = v_xyz[2]
    
    return x_dot, y_dot, z_dot
    

def keplerian_velocity_spherical(r, a, e, i, node, peri, x, y, df): 
    
    # n = mean motion
    # mean motion is the average angular distance traveled per day
    
    n = np.sqrt(G_au_d / (a**3))
    
    ## mean anomaly
    #M = np.arccos(x/r)
        ## mean anomaly
    M = df.apply(return_M, axis=1)
    
    ## calculate eccentric anomaly
    guess_E = M
    solution_E = fsolve(kepler_equation, M, args=(M,e))
    #print(solution_E)
    E = solution_E
    #print((1 - r/a)/e)
    #E = np.arccos((a*e + x)/a) 
    #np.arccos((1 - r/a)/e) 
    #print(E)
    
    '''x_dot = -(a*n*np.sin(E)*np.cos(i))/(1-e*np.cos(E))
    y_dot = (a*n*np.sqrt(1-(e**2))*np.cos(E))/(1-(e*np.cos(E)))
    z_dot = (a*n*np.sin(E)*np.sin(i))/(1-(e*np.cos(E)))'''
    
    ## unrotated coordinates and vector
    X_dot = -(a*n*np.sin(E))/(1-e*np.cos(E))
    Y_dot = (a*n*np.sqrt(1-(e**2))*np.cos(E))/(1-(e*np.cos(E)))
    
    v_XYZ = [X_dot, Y_dot, 0]

    ## rotation matrix
    P11 = np.cos(node)*np.cos(peri)-np.sin(node)*np.cos(i)*np.sin(peri)
    P12 = -np.cos(node)*np.sin(peri)-np.sin(node)*np.cos(i)*np.cos(peri)
    P13 = np.sin(node)*np.sin(i)
    
    P21 = np.sin(node)*np.cos(peri)+np.cos(node)*np.cos(i)*np.sin(peri)
    P22 = -np.sin(node)*np.sin(peri)+np.cos(node)*np.cos(i)*np.cos(peri)
    P23 = -np.cos(node)*np.sin(i)
    
    P31 = np.sin(i)*np.sin(peri)
    P32 = np.sin(i)*np.cos(peri)
    P33 = np.cos(i)
    
    P = np.array([[P11, P12, P13], [P21, P22, P23], [P31, P32, P33]])
    
    ## numpy matrix multiplication throwing errors, too tired to debug it, 
    ## so writing it out by hand for now...
    #v_xyz = np.dot(P, v_XYZ)
    
    v_xyz = [P11*v_XYZ[0]+P12*v_XYZ[1]+P13*v_XYZ[2], P21*v_XYZ[0]+P22*v_XYZ[1]+P23*v_XYZ[2], P31*v_XYZ[0]+P32*v_XYZ[1]+P33*v_XYZ[2]]
    
    x_dot = v_xyz[0]
    y_dot = v_xyz[1]
    z_dot = v_xyz[2]
    
    r_dot = np.sqrt((x_dot**2) + (y_dot**2) + (z_dot**2))
    az_dot = np.arctan(y_dot/x_dot)
    el_dot = np.arccos(z_dot/r_dot)    
    
    return r_dot, az_dot, el_dot
