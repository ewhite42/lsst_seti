## In this program, we will walk through steps of implementing 
## nongravitational accelerations in Sorcha, plotting the results
## after various steps for verification

## MEW, 21 October 2025

import numpy as np
import pylab as plt
import pandas as pd
from scipy.interpolate import interp1d

def compare_nongrav_with_grav(orb_fname, sorcha_fname):

    df_ids = pd.read_csv(orb_fname)
    df = pd.read_csv(sorcha_fname)
    
    id_mask = df_ids['ObjID'].str.contains("ONLY_GRAV")
    df_ids_masked = df_ids[~id_mask]
    ids = df_ids_masked['ObjID']
    
    hdist = np.sqrt(df_ids_masked['x']**2 + df_ids_masked['y']**2 + df_ids_masked['z']**2)
    
    mjd_nongrav_lists = []
    ra_nongrav_lists = []
    dec_nongrav_lists = []
    
    mjd_grav_lists = []
    ra_grav_lists = []
    dec_grav_lists = []
    
    for i in range(len(ids)):
    
        id = ids[i]
        d = hdist[i]
    
        mask1 = df['ObjID'] == id
        df_masked = df[mask1]
        
        mask2 = df['ObjID'] == str(id)+"_ONLY_GRAV"
        df_masked_g = df[mask2]
        
        mjd = df_masked['fieldMJD_TAI']
        ra = df_masked['RA_deg']
        dec = df_masked['Dec_deg']
        
        mjd_nongrav_lists.append(mjd)
        ra_nongrav_lists.append(ra)
        dec_nongrav_lists.append(dec)
        
        mjd_g = df_masked_g['fieldMJD_TAI']
        ra_g = df_masked_g['RA_deg']
        dec_g = df_masked_g['Dec_deg']
        
        mjd_grav_lists.append(mjd_g)
        ra_grav_lists.append(ra_g)
        dec_grav_lists.append(dec_g)
        
        if len(ra) != 0:
        
            mjd = np.array(mjd)
            mjd_g = np.array(mjd_g)
            
            fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(12, 8))
        
            axes[0,0].plot(mjd, ra, label=id)
            axes[0,0].plot(mjd_g, ra_g, label=str(id)+"_ONLY_GRAV")
            axes[0,0].set_xlabel("MJD TAI")
            axes[0,0].set_ylabel("RA (deg)")
            axes[0,0].set_title("RA vs MJD")
        
            axes[0,1].plot(mjd, dec, label=id)
            axes[0,1].plot(mjd_g, dec_g, label=str(id)+"_ONLY_GRAV")
            axes[0,1].set_xlabel("MJD TAI")
            axes[0,1].set_ylabel("Dec (deg)")
            axes[0,1].set_title("Dec vs MJD")
            
            ## determine endpoints of interpolation:
            if mjd[0] < mjd_g[0]:
                mjd_min = mjd_g[0]
            else:
                mjd_min = mjd[0]
                
            if mjd[-1] > mjd_g[-1]:
                mjd_max = mjd_g[-1]
            else:
                mjd_max = mjd[-1]
                
            ## determine how many points to interpolate
            ## for, and also create the interpolation
            ## x-dataset
            
            n_interp = len(mjd)
            mjd_new = np.linspace(mjd_min, mjd_max, 100) #n_interp)
            
            ## interpolate RA for both nongrav and only-grav
            ## cases, so we can take the difference
            
            ra_cubic = interp1d(mjd, ra, kind='cubic')
            ra_interp = ra_cubic(mjd_g)
            
            ra_g_cubic = interp1d(mjd_g, ra_g, kind='cubic')
            ra_g_interp = ra_g_cubic(mjd_new)
            
            ## interpolate Dec for both nongrav and only-grav
            ## cases, so we can take the difference
            
            dec_cubic = interp1d(mjd, dec, kind='cubic')
            dec_interp = dec_cubic(mjd_g)#_new)
            
            dec_g_cubic = interp1d(mjd_g, dec_g, kind='cubic')
            dec_g_interp = dec_g_cubic(mjd_new)
            
            axes[1,0].plot(mjd_g, ra_interp-ra_g)#_interp)
            axes[1,0].set_xlabel("MJD TAI")
            axes[1,0].set_ylabel("Delta RA (deg)")
            axes[1,0].set_title("Difference in RA for nongrav vs grav-only")
        
            axes[1,1].plot(mjd_g, dec_interp-dec_g)#_interp)
            axes[1,1].set_xlabel("MJD TAI")
            axes[1,1].set_ylabel("Dec (deg)")
            axes[1,1].set_title("Dec vs MJD")
        
            axes[0,0].legend()
            axes[0,1].legend()
            
            plt.tight_layout()
            #plt.title("{0}, heliocentric dist. = {1} AU".format(id, d))
            plt.show()

def read_c2017m4atlas_data(): 
    ## read in the csv file resulting from running Sorcha on the
    ## C2017 M4 ATLAS data with only gravitational acceleration

    df = pd.read_csv('/home/ellie/research/lsst/sorcha_output/c2017m4_atlas_2.csv')

    mask_grav = df["ObjID"] == "C2017_M4_ATLAS"
    df_g = df[mask_grav]
    x = df_g['fieldMJD_TAI'] 
    y1 = df_g['RA_deg']
    y2 = df_g['Dec_deg']

    mask_nongrav = df["ObjID"] == "C2017_M4_ATLAS_NG"
    df_nongrav = df[mask_nongrav]
    x_ng = df_nongrav['fieldMJD_TAI']
    y1_ng = df_nongrav['RA_deg']
    y2_ng = df_nongrav['Dec_deg']

    mask_comet_model = df["ObjID"] == "C2017_M4_ATLAS_NG_a"
    df_cm = df[mask_comet_model]
    x_cm = df_cm['fieldMJD_TAI']
    y1_cm = df_cm['RA_deg']
    y2_cm = df_cm['Dec_deg']
    
    print(len(x))
    print(len(x_ng))
    
    f_cubic = interp1d(x, y1, kind='cubic')
    y1_interp = f_cubic(x_ng)
    
    #plt.plot(x,y1, label="Gravity Only")
    plt.plot(x_ng, y1_ng-y1_interp, label="Asteroid Model")
    #plt.plot(x_cm, y1_cm, label="Comet Model")
    plt.xlabel('Field MJD')
    plt.ylabel('RA (deg)')
    plt.title('Delta-RA')
    plt.legend()
    plt.show()

    plt.plot(x,y1, label="Gravity Only")
    plt.plot(x_ng, y1_ng, label="Asteroid Model")
    #plt.plot(x_cm, y1_cm, label="Comet Model")
    plt.xlabel('Field MJD')
    plt.ylabel('RA (deg)')
    plt.legend()
    plt.show()
    
    f2_cubic = interp1d(x, y2)
    y2_interp = f2_cubic(x_ng)
    
    #plt.plot(x,y1, label="Gravity Only")
    plt.plot(x_ng, y2_ng-y2_interp, label="Asteroid Model")
    #plt.plot(x_cm, y1_cm, label="Comet Model")
    plt.xlabel('Field MJD')
    plt.ylabel('RA (deg)')
    plt.title('Delta-Dec')
    plt.legend()
    plt.show()

    plt.plot(x,y2, label="Gravity Only")
    plt.plot(x_ng, y2_ng, label="Asteroid Model")
    #plt.plot(x_cm, y2_cm, label="Comet Model")
    plt.xlabel('Field MJD')
    plt.ylabel('Dec (deg)')
    plt.legend()
    plt.show()
    
def main():

    ofname = 'nongrav_top20_orb.csv'
    sfname = '/home/ellie/research/lsst/sorcha_output/nongrav_top20.csv'
    compare_nongrav_with_grav(ofname, sfname)
    
if __name__ == '__main__':
    main()
