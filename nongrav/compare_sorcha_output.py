## code from Google

import pylab as plt
import pandas as pd
from astroquery.jplhorizons import Horizons
import astropy.units as u

def read_sorcha_output(fname):
    df = pd.read_csv(fname)
    mjd = df['fieldMJD_TAI']
    ra = df['RA_deg']
    dec = df['Dec_deg']
    
    return mjd, ra, dec
    
def get_jpl_output(target_id, startday, endday, step, obscode='X05'):

    # use 'yyyy-mm-dd' format for startday and endday
    # make sure to indicate unit in step, e.g.: '1d' for 1 day

    observer_location = {'body': '399'} # Earth's ID

    # Define the time range for the ephemerides
    epochs = {
        'start': startday,
        'stop': endday,
        'step': step
    }

    # Create a Horizons object
    obj = Horizons(id=target_id, location=obscode, epochs=epochs)

    # Retrieve ephemeris data
    ephemerides = obj.ephemerides() # RA/DEC, position vectors, apparent magnitude

    # Convert the Astropy Table object to a Pandas DataFrame
    df = ephemerides.to_pandas()

    ra_jpl = df['RA']
    dec_jpl = df['DEC']
    mjd_jpl = df['datetime_jd']-2400000.5
    
    return mjd_jpl, ra_jpl, dec_jpl
    
def main(): 

    fname_kep = "pratchett.csv"
    fname_cart = "pratchett_cart.csv"

    mjd1,ra1,dec1 = read_sorcha_output(fname_cart)
    #mjd2,ra2,dec2 = read_sorcha_output(fname_cart)
     
    ## Here is code to use when comparing Sorcha output to JPL
    ## Horizons output (instead of comparing 2 Sorcha outputs): 
        
    target_id = '127005' ## for Pratchett
    startday = '2024-02-10'
    endday = '2024-08-08'
    step = '1d'
   
    mjd2,ra2,dec2 = get_jpl_output(target_id, startday, endday, step)     
    
    plt.plot(mjd1, ra1, label="Sorcha")
    plt.plot(mjd2, ra2, label="JPL Horizons")
    plt.legend()
    plt.xlabel("MJD TDB")
    plt.ylabel("RA")
    plt.show()
    
    plt.plot(mjd1, dec1, label="Sorcha")
    plt.plot(mjd2, dec2, label="JPL Horizons")
    plt.legend()
    plt.xlabel("MJD TDB")
    plt.ylabel("Dec")
    plt.show()
    
if __name__ == '__main__':
    main()
