## script to query objects from JPL Horizons
## with nongravitational acceleration parameters

## MEW 1 Nov 2025

import pandas as pd
from astroquery.jplhorizons import Horizons

def make_orbits(fname):

    df = pd.read_csv(fname)
    
    df_obj = pd.DataFrame()
    
    ids = []
    
    x = []
    y = []
    z = []
    
    vx = []
    vy = []
    vz = []
    
    a1 = []
    a2 = []
    a3 = []
    
    mjd_tdb = []
    
    num_objs = len(df['spkid'])
    
    for i,row in df.iterrows():
        spkid = row['spkid']
        id_str = "DES= {};".format(spkid)
        
        if spkid == 54427459:
            print("Skipping {}, as an error occurred.".format(spkid))
            continue
            
        else:
            obj = Horizons(id=id_str, location='X05', epochs={'start':'2025-6-15', 'stop':'2025-7-15','step':'1d'})
             
        try:	                       
            vec = obj.vectors()
        except Exception as e:
            last_row = str(e).split('\n')[-2]
            new_id = last_row.split()[0]
            obj = Horizons(id=new_id, location='X05', epochs={'start':'2025-6-15', 'stop':'2025-7-15','step':'1d'})
            vec = obj.vectors()
            #continue
        
        target = vec['targetname'][0].lstrip()
        target = target.replace(' ', '_')
        target = target.replace(')', '_')
        target = target.replace('(', '_')
        target = target.replace('__', '_')
        target = target.lstrip('_')
        target = target.rstrip('_')
        
        ids.append(target)

        x.append(vec['x'][0])
        y.append(vec['y'][0])
        z.append(vec['z'][0])
        
        vx.append(vec['vx'][0])
        vy.append(vec['vy'][0])
        vz.append(vec['vz'][0])
        
        a1.append(row['A1'])
        a2.append(row['A2'])
        a3.append(row['A3'])
        
        mjd_tdb.append(vec['datetime_jd'][0]-2400000.5)
        
        print("{0} out of {1} objects written".format(i+1, num_objs))
        
    df_obj['ObjID'] = ids
    df_obj['FORMAT'] = ['NONGRAV']*len(ids)
    
    df_obj['x'] = x
    df_obj['y'] = y
    df_obj['z'] = z
    
    df_obj['xdot'] = vx
    df_obj['ydot'] = vy
    df_obj['zdot'] = vz
    
    df_obj['a1'] = a1
    df_obj['a2'] = a2
    df_obj['a3'] = a3
    
    df_obj['model'] = ['ASTEROID']*len(ids)
    df_obj['epochMJD_TDB'] = mjd_tdb
    
    df_sorted = df_obj.sort_values(by='a1', ascending=False)
    #df_top_20 = df_sorted.head(20)
    
    df_sorted.to_csv('nongrav_sorted_orb.csv', index=False)
    print("Orbits file created.")
    
def make_physical_parameters_file(orb_fname):
    df = pd.read_csv(orb_fname)
    ids = df['ObjID']
    
    df_phy = pd.DataFrame()
    
    df_phy['ObjID'] = ids
    df_phy['H_r'] = [5.63]*len(df_phy['ObjID'])
    df_phy['u-r'] = [2.55]*len(df_phy['ObjID'])
    df_phy['g-r'] = [0.92]*len(df_phy['ObjID'])
    df_phy['i-r'] = [-0.38]*len(df_phy['ObjID'])
    df_phy['z-r'] = [-0.59]*len(df_phy['ObjID'])
    df_phy['y-r'] = [-0.70]*len(df_phy['ObjID'])
    df_phy['GS'] = [0.15]*len(df_phy['ObjID'])
    
    df_phy.to_csv('nongrav_sorted_phy.csv', index=False)
    
def insert_duplicates_no_nongrav(orbfname, phyfname):

    '''This function inserts copies of rows for each
       object with the Ai parameters set to zero, so
       that you can easily compare the version of the
       objects with nongrav acceleration and those
       without. 
       
       New versions of both the orbits file and the
       physical parameters files are generated with 
       both the original rows and the new duplicates.'''

    df_orb = pd.read_csv(orbfname)
    df_phy = pd.read_csv(phyfname)
    
    df_orb_copy = df_orb.copy()
    df_phy_copy = df_phy.copy()
    
    for i,row in df_orb_copy.iterrows():
        df_orb_copy.at[i, 'ObjID'] = str(row['ObjID'])+'_ONLY_GRAV'
        df_orb_copy.at[i, 'a1'] = 0.00
        df_orb_copy.at[i, 'a2'] = 0.00
        df_orb_copy.at[i, 'a3'] = 0.00
        
    for j, row2 in df_phy_copy.iterrows():
        df_phy_copy.at[j, 'ObjID'] = str(row2['ObjID'])+'_ONLY_GRAV'
        
    df_phy_final = pd.concat([df_phy, df_phy_copy], axis=0)
    df_phy_final.to_csv('nongrav_top20_phy.csv', index=False)
    
    df_orb_final = pd.concat([df_orb, df_orb_copy], axis=0)
    df_orb_final.to_csv('nongrav_top20_orb.csv', index=False)
        
def main():
    #fname = '/home/ellie/research/lsst/sbdb_nongrav_data.csv'
    #make_orbits(fname)
    #make_physical_parameters_file(orb_fname)
    
    orbfname = '/home/ellie/research/lsst/lsst_seti/nongrav/sorcha_ng/nongrav_sorted_orb.csv'
    phyfname = '/home/ellie/research/lsst/lsst_seti/nongrav/sorcha_ng/nongrav_sorted_phy.csv'
    insert_duplicates_no_nongrav(orbfname, phyfname)
    
if __name__ == '__main__':
    main()
        
