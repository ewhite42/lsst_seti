## This program will allow us to generate a simulated
## Solar System object with nongravitational
## accelerations using ASSIST, and plot its
## trajectories

## MEW, 25 Jan 2026

import numpy as np
import pylab as plt
import math
import random

from datetime import datetime
import rebound
import assist

## need to heavily modify the below function, 
## including plotting the orbit of known Solar
## System bodies for reference (possibly)

def sim_ephemerides():

    ## use arrays of 3d vectors to simulate 
    ## ephemerides of weird objects with ASSIST
    
    ## define possible range of integration times
    min_int = 10 ## days
    max_int = 3650 ## days
    
    ## define how many sycoraxes you want to generate
    n = 10
    
    ## get the arrays of vectors
    dist_vectors, v_vectors, a_vectors = gen_starting_values(n)
    
    ## simulate each ephemeris
    ephem = assist.Ephem("/home/ellie/src/assist/data/linux_p1550p2650.440", "/home/ellie/src/assist/data/sb441-n16.bsp")
    
    for i in range(len(dist_vectors)):
    
        print('simulating object {}'.format(i))
        ## create the rebound simulation and attach ASSIST
        sim = rebound.Simulation()
        ex = assist.Extras(sim, ephem)
    
        ## add each of our generated objects to the simulation
        sycorax_initial_helio = rebound.Particle(x = dist_vectors[i][0], \
                                                 y=dist_vectors[i][1], \
                                                 z=dist_vectors[i][2], \
                                                 vx = v_vectors[i][0], \
                                                 vy=v_vectors[i][1], \
                                                 vz=v_vectors[i][2])  
                
        ## assign a start date to each object  
        start_date = gen_random_date(2025, 1, 1, 2030, 1, 1)
        t_initial = start_date - ephem.jd_ref
        
        ## convert from heliocentric to barycentric
        sun_initial = ephem.get_particle("sun", t_initial)
        sycorax_initial = sycorax_initial_helio + sun_initial

        ##      
        sim = rebound.Simulation()
        sim.add(sycorax_initial)
        sim.t = t_initial   
        sim.ri_ias15.min_dt = 0.3
        extras = assist.Extras(sim, ephem)
        
        ## set up the acceleration parameters
        extras.particle_params = np.array([a_vectors[i][0], a_vectors[i][1], a_vectors[i][2]])
        
        ## create the ephemeris
        int_days = random.randint(min_int, max_int)
        t_final = t_initial + int_days
        
        N_samples = int_days*3 ## this is like saying you observe it 3x per day
        times = np.linspace(t_initial, t_final, N_samples, endpoint=True)
        
        sycorax_pos = np.zeros((N_samples, 3))
        earth_pos = np.zeros((N_samples, 3))
        
        for idx, t in enumerate(times):
            #print('integrating object {}, iteration {}/{}, time {}'.format(i, idx, N_samples, t))
            extras.integrate_or_interpolate(t)
            sycorax_pos[idx] = sim.particles[0].xyz
            earth_pos[idx] = ephem.get_particle("earth", t).xyz
            
        plt.plot(sycorax_pos[:,0], sycorax_pos[:,1], label="Sycorax {}".format(i))
    plt.plot(earth_pos[:,0],earth_pos[:,1], label="Earth")
    
    plt.xlabel("x [AU]")
    plt.ylabel("y [AU]")
    plt.legend()    
    plt.show()
