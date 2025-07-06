''' The purpose of this program is to read in data from the
    extended MPC ORB table, downloaded from minorplanetcenter.org.
    
    ~ MEW, 6 July 2025 
'''

import numpy as np
import bigjson
import pandas as pd

def main():
    mpc_fname = '/home/ellie/Downloads/mpcorb_extended.json'
    out_fname = '/home/ellie/research/lsst/mpcorb_extended_full.csv'
    
    start_idx = 0
    num_rows = 1000
    
    desigs = []
    epochs = []
    a = []
    e = []
    incl = []
    node = []
    peri = []
    M = []
    n = []
    last_obs = []
    orbit_types = [] 
    
    with open(mpc_fname, 'rb') as f:
        j = bigjson.load(f)
        
        for i in range(start_idx, start_idx+num_rows):
            element = j[i]    
            
            desigs.append(element['Principal_desig'])
            epochs.append(element['Epoch'])
            a.append(element['a'])
            e.append(element['e'])
            incl.append(element['i'])
            node.append(element['Node'])
            peri.append(element['Peri'])
            M.append(element['M'])
            n.append(element['n'])
            last_obs.append(element['Last_obs'])
            orbit_types.append(element['Orbit_type'])
    
    df = pd.DataFrame()
    df['Principle_desig'] = desigs
    df['Epoch'] = epochs
    df['a'] = a
    df['e'] = e
    df['i'] = incl
    df['Node'] = node
    df['Peri'] = peri
    df['M'] = M
    df['n'] = n
    df['Last_obs'] = last_obs
    df['Orbit_type'] = orbit_types
    
    df.to_csv(out_fname)
    print(df)

if __name__ == '__main__':
    main()
