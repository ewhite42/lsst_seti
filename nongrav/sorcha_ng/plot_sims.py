## In this program, we will walk through steps of implementing 
## nongravitational accelerations in Sorcha, plotting the results
## after various steps for verification

## MEW, 21 October 2025

import pylab as plt
import pandas as pd
from scipy.interpolate import interp1d

def main(): 
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
    
if __name__ == '__main__':
    main()
