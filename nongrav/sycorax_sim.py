## This program will allow us to generate simulated
## Solar System objects with unusual velocities and
## accelerations using ASSIST

## MEW, 24 Aug 2025

import numpy as np
import pylab as plt
import random
from datetime import datetime
import rebound
import assist

def create_3d_vector(magnitude):

    ## this function will allow us to generate
    ## the components of 3d vectors with the 
    ## magnitude provided
    
    theta = np.random.uniform(0, math.pi)
    phi = np.random.uniform(0, 2*math.pi)
    
    x = magnitude * math.sin(theta) * math.cos(phi)
    y = magnitude * math.sin(theta) * math.sin(phi)
    z = magnitude * math.cos(theta)

    return [x, y, z]
    
def gen_starting_values():

    ## this function will allow us to generate
    ## locations, velocities, and acceleration 
    ## parameters that we can use to initialize
    ## the ASSIST integrations
    
    ## how many sycoraxes (sycori? sycoraxi?) do you want to generate?
    n = 100
    
    ## set limits on locations
    min_dist = 0.00 ## AU
    max_dist = 20.00 ## AU
    
    ## set limits on velocities
    min_v = 0.00 ## AU / day
    max_v = 35.00 ## AU / day - note this is 20% the speed of light - Starshot's target speed
    
    ## set limits on accelerations
    min_a = 0.00 ## AU / day**2
    max_a = 4990.00 ## AU / day**2 - note this is Starshot's target acceleration
    
    ## generate arrays of random values for
    ## location, velocity, and acceleration
    dist_arr = np.random.uniform(low=min_dist, high=max_dist, size=n)
    v_arr = np.random.uniform(low=min_v, high=max_v, size=n)
    a_arr = np.random.uniform(low=min_a, high=max_a, size=n)
    
    ## transform each magnitude of distance, velocity, 
    ## and acceleration into 3d vectors, choosing the
    ## components at random
    
    dist_vectors = []
    v_vectors = []
    a_vectors = []
    
    for i in range(n):
        dist_vectors.append(create_3d_vector(dist_arr[i]))
        v_vectors.append(create_3d_vector(v_arr[i]))
        a_vectors.append(create_3d_vector(a_arr[i]))
        
    return dist_vectors, v_vectors, a_vectors
    
def gen_random_date(startyear, startmonth, startday, endyear, endmonth, endday):
    start_date = datetime(startyear, startmonth, startday)
    end_date = datetime(endyear, endmonth, endday)
    
    julian_startdate = start_date.toordinal() + 1721425
    julian_enddate = end_date.toordinal() + 1721425
    
    return np.random.uniform(low=julian_startdate, high=julian_enddate)
    
def sim_ephemerides():

    ## use arrays of 3d vectors to simulate 
    ## ephemerides of weird objects with ASSIST
    
    ## define possible range of integration times
    min_int = 10 ## days
    max_int = 3650 ## days
    
    ## get the arrays of vectors
    dist_vectors, v_vectors, a_vectors = gen_starting_values()
    
    ## simulate each ephemeris
    
    ## note I will probably need to adjust the path to the
    ## files below
    ephem = assist.Ephem("../data/linux_p1550p2650.440", "../data/sb441-n16.bsp")
    
    ## https://assist.readthedocs.io/en/latest/jupyter_examples/GettingStarted/ 
    
    for i in range(len(dist_vectors)):
    
        ## create the rebound simulation and attach ASSIST
        sim = rebound.Simulation()
        ex = assist.Extras(sim, ephem)
    
        ## add each of our generated objects to the simulation
        sycorax_initial_helio = sim.add(x = dist_vectors[0], y=dist_vectors[1], z=dist_vectors[2], \
                                        vx = v_vectors[0], vy=v_vectors[1], vz=v_vectors[2])  
                
        ## assign a start date to each object  
        start_date = gen_random_date(2026, 1, 1, 2030, 12, 31)
        t_initial = start_date - ephem.jd_ref
        
        ## convert from heliocentric to barycentric
        sun_initial = ephem.get_particle("sun", t_initial)
        sycorax_initial = sycorax_initial_helio + sun_initial

        ##        
        sim.ri_ias15.min_dt = 0.001
        extras = assist.Extras(sim, ephem)
        
        ## set up the acceleration parameters
        extras.particle_params = np.array([a_vectors[0], a_vectors[1], a_vectors[2]])
        
        ## create the ephemeris
        int_days = random.randint(min_int, max_int)
        t_final = t_initial + int_days  - ephem.jd_ref
        
        N_samples = int_days*3 ## this is like saying you observe it 3x per day
        times = np.linspace(t_initial, t_final, N_samples, endpoint=True)
        
        sycorax_pos = np.zeros((N_samples, 3))
        earth_pos = np.zeros((N_samples, 3))
        
        for i, t in enumerate(times):
            extras.integrate_or_interpolate(t)
            sycorax_pos[i] = sim.particles[0].xyz
            earth_pos[i] = ephem.get_particle("earth", t).xyz
            
        plt.xlabel("x [AU]")
        plt.ylabel("y [AU]")
        plt.plot(sycorax_pos[:,0], sycorax_pos[:,1], label="Sycorax {}".format(i))
        plt.plot(earth_pos[:,0],earth_pos[:,1], label="Earth")
        plt.legend()    
        plt.show()
        
def main():
    sim_ephemerides()
    
if __name__ == '__main__':
    main()
