## a quick program to calculate expected and measured velocities 
## of Solar System objects from Rubin data

## 1 Sept 2024 -MEW

import numpy as np

## define our constants
G_au_yr = 39.422888 #G in au*(au/yr)^2/M_sun units
G_km_s = 887.1287008 #G in au*(km/s)^2/M_sun units
G_au_s = 3.96402313*10**(-14) #G in au*(au/s)^2/M_sun units

def keplerian_velocity(r, a):
    # note r = distance from Sun to object in AU
    # and a = semi-major axis in AU
    
    # note also that we are assuming M = M_sun + M_object ~ M_sun, 
    # and since M_sun is the unit of mass in the G definition
    # we are omitting it here. 
    
    return np.sqrt((G_au_yr)*((2/r)-(1/a))), np.sqrt((G_km_s)*((2/r)-(1/a))), np.sqrt((G_au_s)*((2/r)-(1/a)))

def measured_velocity(vx, vy, vz):
    # not totally sure what units these are in -- data product
    # document says that the heliocentric velocities are in "AU" 
    # which is nonsense...
    
    return np.sqrt(vx**2 + vy**2 + vz**2)
    
    
def main():
    # run our calculations...
    
    # first, import an LSST simulation table
    fname = 'table_dp03_catalogs_10yr.MPCORB-data.lsst.cloud.csv'
    infile = open(fname, 'r')
    
    rows = []
    
    # this is a quick and rough way to load in the data -- remember to port it over 
    # to a nicer, pandas version (similar to Brian's data loader) soon - this
    # is just to get the ideas down
    for line in infile.readlines()[1:]:
        l = line.split(',')
        rows.append(l)
        
    print(rows[0])
    
    velocities = []
    
    for cols in rows:
        r = float(cols[7])
        q = float(cols[2])
        e = float(cols[1])
        a = q/(1-e)
        vk_au, vk_km, vk_au_s = keplerian_velocity(r, a)
        
        vx = float(cols[8])
        vy = float(cols[9])
        vz = float(cols[10])
        v_actual = measured_velocity(vx, vy, vz)
        
        v = [vk_au, vk_km, vk_au_s, v_actual]
        velocities.append(v)
        print(v)       
    
    
    '''# get r from heliocentricDist in SSSource table
    r = 2.8450608
    
    # get a by calculating it from q and e from MPCORB table
    q = 2.169822915818627 #in AU
    e = 0.14689253115098788

    a = q/(1-e)
    
    vk_au, vk_km, vk_au_s = keplerian_velocity(r, a)
    
    # get heliocentric vx, vy, and vz from SSSource table
    vx = 0.0065061725
    vy = -0.0066827396
    vz = -0.0021914262
    
    v_actual = measured_velocity(vx, vy, vz)
    
    print("keplerian velocity (au/yr): "+str(vk_au))
    print("keplerian velocity (km/s): "+str(vk_km))
    print("keplerian velocity (au/s): "+str(vk_au_s))
    print("measured velocity (au/?): "+str(v_actual))'''
    
if __name__ == '__main__':
    main()
