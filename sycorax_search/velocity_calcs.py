## the functions here will allow us to calculate the 
## measured velocity and predicted Keplerian velocity
## of a given object based on parameters in the LSST
## simulation data

## M. E. White, 15 Jan 2025

import numpy as np

from scipy.optimize import fsolve
from astropy import units as u

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
    
def return_nu(row):
    #return np.arccos(row['heliocentricX']/row['heliocentricDist']) + ((1.5*np.pi) if row['heliocentricY'] < 0 else 0)
    
    node = row['node'] * np.pi/180 
    peri = row['peri'] * np.pi/180 ## argument of the perihelion
    i = row['incl']* np.pi/180
    
    P11 = np.cos(node)*np.cos(peri)-np.sin(node)*np.cos(i)*np.sin(peri)
    P21 = np.sin(node)*np.cos(peri)+np.cos(node)*np.cos(i)*np.sin(peri)
    P31 = np.sin(i)*np.sin(peri)
    
    x = row['heliocentricX']
    y = row['heliocentricY']
    z = row['heliocentricZ']
    
    X = P11*x + P21*y + P31*z
    
    if row['heliocentricY'] > 0:
        return (2*np.pi) - np.arccos(X/row['heliocentricDist'])
        
    else: 
        return np.arccos(X/row['heliocentricDist'])
        
def return_nu_mpc(row):
    
    node = row['Node'] * np.pi/180 
    peri = row['Peri'] * np.pi/180 ## argument of the perihelion
    i = row['i']* np.pi/180
    
    P11 = np.cos(node)*np.cos(peri)-np.sin(node)*np.cos(i)*np.sin(peri)
    P21 = np.sin(node)*np.cos(peri)+np.cos(node)*np.cos(i)*np.sin(peri)
    P31 = np.sin(i)*np.sin(peri)
    
    P12 = -np.cos(node)*np.sin(peri)-np.sin(node)*np.cos(i)*np.cos(peri)
    P22 = -np.sin(node)*np.sin(peri)+np.cos(node)*np.cos(i)*np.cos(peri)
    P32 = np.sin(i)*np.cos(peri)
    
    P13 = np.sin(node)*np.sin(i)
    P23 = -np.cos(node)*np.sin(i)
    P33 = np.cos(i)
    
    x = row['X']
    y = row['Y']
    z = row['Z']
    
    X = P11*x + P21*y + P31*z
    Y = P12*x + P22*y + P32*z
    Z = P13*x + P23*y + P33*z
    
    ## obliquity
    eps = np.radians(23.43928)
    
    Y_eq = np.cos(eps)*Y - np.sin(eps)*Z
    
    ## note the switch to the < sign below (instead of
    ## the > sign used for the Rubin data, above). 
    ## I am guessing this means that MPC uses a different
    ## nu convention than Rubin.
    
    if Y_eq > 0:
        return np.arccos(X/row['heliocentricDist'])        
    else: 
        return (2*np.pi) - np.arccos(X/row['heliocentricDist'])
    #return np.arccos(X/row['heliocentricDist'])
        
def return_nu_ap(row):

    node = row['node'] * np.pi/180 
    peri = row['peri'] * np.pi/180 ## argument of the perihelion
    i = row['incl']* np.pi/180
    
    P11 = np.cos(node)*np.cos(peri)-np.sin(node)*np.cos(i)*np.sin(peri)
    P21 = np.sin(node)*np.cos(peri)+np.cos(node)*np.cos(i)*np.sin(peri)
    P31 = np.sin(i)*np.sin(peri)
    
    x = row['heliocentricX']
    y = row['heliocentricY']
    z = row['heliocentricZ']
    
    X = P11*x + P21*y + P31*z

    #return np.arccos(row['heliocentricX']/row['heliocentricDist']) + ((1.5*np.pi) if row['heliocentricY'] < 0 else 0)
    if row['heliocentricY'] > 0:
        return u.Quantity(2*np.pi - np.arccos(X/row['heliocentricDist']), u.rad)
        
    else: 
        return u.Quantity(np.arccos(X/row['heliocentricDist']), u.rad)
        
def return_nu_ap_mpc(row):

    node = row['Node'] * np.pi/180 
    peri = row['Peri'] * np.pi/180 ## argument of the perihelion
    i = row['i']* np.pi/180
    
    P11 = np.cos(node)*np.cos(peri)-np.sin(node)*np.cos(i)*np.sin(peri)
    P21 = np.sin(node)*np.cos(peri)+np.cos(node)*np.cos(i)*np.sin(peri)
    P31 = np.sin(i)*np.sin(peri)
    
    x = row['X']
    y = row['Y']
    z = row['Z']
    
    X = P11*x + P21*y + P31*z

    #return np.arccos(row['heliocentricX']/row['heliocentricDist']) + ((1.5*np.pi) if row['heliocentricY'] < 0 else 0)
    if row['Y'] > 0:
        return u.Quantity(2*np.pi - np.arccos(X/row['heliocentricDist']), u.rad)
        
    else: 
        return u.Quantity(np.arccos(X/row['heliocentricDist']), u.rad)

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
    

def keplerian_velocity_xyz(r, a, e, i, node, peri, x, y, df, mpc=False): 

    # r = heliocentric distance (~ distance from 1 focus)
    # a = semi-major axis
    # e = eccentricity
    # n = mean daily motion
    # i = inclination
    
    # n = mean motion
    # mean motion is the average angular distance traveled per day
    n = np.sqrt(G_au_d / (a**3))
    
    ## obliquity
    eps = np.radians(23.43928)
    
    ## true anomaly
    if mpc:
        nu = df.apply(return_nu_mpc, axis=1)
    else:
        nu = df.apply(return_nu, axis=1)
    #M = np.arccos(x/r)
    
    ## calculate eccentric anomaly
    #guess_E = M
    #solution_E = fsolve(kepler_equation, M, args=(M,e))
    
    sin_E = (np.sin(nu) * np.sqrt(1 - e**2)) / (1 + e * np.cos(nu)) 
    cos_E = (e + np.cos(nu)) / (1 + e * np.cos(nu)) 
    E = np.arctan2(sin_E, cos_E) #+ np.pi 

    #E = solution_E + np.pi
    
    #print((a*e + x)/a)
    #E = np.arccos((a*e + x)/a)#np.arccos((1 - r/a)/e) ## calculate eccentric anomaly
    
    #node = 0 #2*np.pi - node # node+np.pi
    #peri = 0 #2*np.pi - peri #peri+np.pi #5*np.pi/180
    #i = i #*180/np.pi #0.5*np.pi - i
    
    ## unrotated coordinates and vector
    X_dot = (-(a*n*np.sin(E))/(1-e*np.cos(E)))
    Y_dot = ((a*n*np.sqrt(1-(e**2))*np.cos(E))/(1-(e*np.cos(E))))
    
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
    
    #i = 0.5*np.pi - i
    
    ## try a different rotation matrix...
    '''P11 = np.cos(i)
    P12 = 0
    P13 = np.sin(i)
    
    P21 = 0
    P22 = 1
    P23 = 0
    
    P31 = -np.sin(i)
    P32 = 0
    P33 = np.cos(i)'''
    
    #P31 = np.sin(i)*np.sin(peri)
    #P32 = np.sin(i)*np.cos(peri)
    #P33 = np.cos(i)
    
    v_xyz = [P11*v_XYZ[0]+P12*v_XYZ[1]+P13*v_XYZ[2], P21*v_XYZ[0]+P22*v_XYZ[1]+P23*v_XYZ[2], P31*v_XYZ[0]+P32*v_XYZ[1]+P33*v_XYZ[2]]
    
    ## convert from ecliptic plane to equatorial J2000 plane
    
    x_dot = v_xyz[0] #X_dot #v_xyz[0] #X_dot #
    y_dot = np.cos(eps)*v_xyz[1] - np.sin(eps)*v_xyz[2] #Y_dot*np.cos(i) #v_xyz[1] #Y_dot #
    z_dot = np.sin(eps)*v_xyz[1] + np.cos(eps)*v_xyz[2] #-Y_dot*np.sin(i) #v_xyz[2]
    
    return x_dot, y_dot, z_dot

def keplerian_velocity_spherical(r, a, e, i, node, peri, x, y, df, mpc=False): 
    
    # n = mean motion
    # mean motion is the average angular distance traveled per day
    
    v_xyz = keplerian_velocity_xyz(r, a, e, i, node, peri, x, y, df, mpc)
    
    '''n = np.sqrt(G_au_d / (a**3))
    
    nu = df.apply(return_nu, axis=1)
    #M = np.arccos(x/r)
    
    ## calculate eccentric anomaly
    #guess_E = M
    #solution_E = fsolve(kepler_equation, M, args=(M,e))
    
    sin_E = (np.sin(nu) * np.sqrt(1 - e**2)) / (1 + e * np.cos(nu)) 
    cos_E = (e + np.cos(nu)) / (1 + e * np.cos(nu)) 
    E = np.arctan2(sin_E, cos_E) #+ np.pi 
    
    #print((1 - r/a)/e)
    #E = np.arccos((a*e + x)/a) 
    #np.arccos((1 - r/a)/e) 
    #print(E)
    
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
    
    v_xyz = [P11*v_XYZ[0]+P12*v_XYZ[1]+P13*v_XYZ[2], P21*v_XYZ[0]+P22*v_XYZ[1]+P23*v_XYZ[2], P31*v_XYZ[0]+P32*v_XYZ[1]+P33*v_XYZ[2]]'''
    
    x_dot = v_xyz[0]
    y_dot = v_xyz[1]
    z_dot = v_xyz[2]
    
    r_dot = np.sqrt((x_dot**2) + (y_dot**2) + (z_dot**2)) #X_dot**2 + Y_dot**2) #np.sqrt((x_dot**2) + (y_dot**2) + (z_dot**2))
    az_dot = np.arctan2(y_dot, x_dot)
    el_dot = np.arccos(z_dot/r_dot)    
    
    return r_dot, az_dot, el_dot

def new_keplerian_xyz(r, a, e, i, node, peri, x, y, df, mpc=False): 

    # r = heliocentric distance (~ distance from 1 focus)
    # a = semi-major axis
    # e = eccentricity
    # n = mean daily motion
    # i = inclination
    
    # n = mean motion
    # mean motion is the average angular distance traveled per day
    n = np.sqrt(G_au_d / (a**3))
    
    # obliquity: 
    eps = 23.43928
    
    ## true anomaly
    if mpc:
        nu = df.apply(return_nu_mpc, axis=1)
    else:
        nu = df.apply(return_nu, axis=1)
    #M = np.arccos(x/r)
    
    ## calculate eccentric anomaly
    #guess_E = M
    #solution_E = fsolve(kepler_equation, M, args=(M,e))
    
    sin_E = (np.sin(nu) * np.sqrt(1 - e**2)) / (1 + e * np.cos(nu)) 
    cos_E = (e + np.cos(nu)) / (1 + e * np.cos(nu)) 
    E = np.arctan2(sin_E, cos_E) #+ np.pi 
    
    ## compute object's heliocentric coordinates in the orbital plane
    x_prime = a*(np.cos(E) - e)
    y_prime = a*np.sqrt(1 - (e**2))*np.sin(E)
    
    ## coordinates in the J2000 ecliptic plane
    xecl = (np.cos(peri)*np.cos(node) - np.sin(peri)*np.sin(node)*np.cos(i))*x_prime + \
           (-np.sin(peri)*np.cos(node) - np.cos(peri)*np.sin(node)*np.cos(i))*y_prime
           
    yecl = (np.cos(peri)*np.sine(node) + np.sin(peri)*np.cos(node)*np.cos(i))*x_prime + \
           (-np.sin(peri)*np.sin(node) + np.cos(peri)*np.cos(node)*np.cos(i))*y_prime
           
    zecl = np.sin(peri)*np.sin(i)*x_prime + np.cos(peri)*np.sin(i)*y_prime
    
    return x_dot, y_dot, z_dot   
