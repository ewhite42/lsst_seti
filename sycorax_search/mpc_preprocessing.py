''' The purpose of this program is to read in data from the
    extended MPC ORB table, downloaded from minorplanetcenter.org.
    
    ~ MEW, 6 July 2025 
'''

import numpy as np
import bigjson
import pandas as pd

from astroquery.mpc import MPC

def read_mpcorb(mpc_fname, numrows=10, random=False, exclude=None):

    ''' this function reads an MPCORB json file and returns
        a pandas dataframe with relevant fields. It also saves
        the dataframe to a csv file. 
        
        mpc_fname is the filepath for the extended MPCORB json file
        
        random - if random is False, objects are read in order. if random 
                 is True, objects are selected randomly. 
        
        exclude is a list of indices of entries in the json file to exclude.
        
        By default, exclude is set to None, so no objects are excluded. 
        The purpose of exclude is to ensure that if you are randomly 
        selecting objects from the json (instead of doing it in order), 
        you will have the choice to leave out a list of objects. '''
    
    start_idx = 0
    
    ## for some reason, random doesn't work just yet
    
    if random:
        idx_list = np.random.randint(0, 800000, numrows)
    else:
        idx_list = range(start_idx, start_idx+numrows)
    
    numbers = []
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
        
        for i in idx_list:
            element = j[i]    
            
            numbers.append(element['Number'])
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
    df['Number'] = numbers
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

    return df
    
def read_MPCephem(nums):
    
    ephs = []

    for i in nums: 
    
        try: 
            eph = MPC.get_ephemeris(str(i), eph_type='heliocentric', number=1)
            eph = eph.to_pandas()
            ephs.append(eph)
            
        except RuntimeError:
            print('No data found')
        except ValueError:
            print('Object {} not found'.format(i))
        except Exception as e:
            print('some other type of error occurred')
    
    ephs_df = pd.concat(ephs)
    ephs_df['Number'] = nums
    
    ## fix weird issue with Y' and Z' having a number, space, then the correct value
    ## so that now the columns just have the correct values
    
    ephs_df[["junk y'","Y'"]] = ephs_df["Y'"].str.split(' ', n=1, expand=True)
    ephs_df = ephs_df.drop("junk y'", axis=1)
    
    ephs_df[["junk z'","Z'"]] = ephs_df["Z'"].str.split(' ', n=1, expand=True)
    ephs_df = ephs_df.drop("junk z'", axis=1)
    
    return ephs_df    

if __name__ == '__main__':
    mpc_fname = '/home/ellie/Downloads/mpcorb_extended.json'
    out_fname = '/home/ellie/research/lsst/mpcorb_extended_pt1.csv'
    
    num_rows = 10
    
    mpcorb_df = read_mpcorb(mpc_fname, numrows=num_rows, random=False)
    nums = mpcorb_df['Number'].tolist()
    
    ephemerides_df = read_MPCephem(nums)
    
    df = pd.merge(mpcorb_df, ephemerides_df, on='Number')
    df.to_csv(out_fname)
    print(df.head)
    
