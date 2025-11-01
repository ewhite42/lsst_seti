## this program adds any needed columns to the
## Sorcha output file so that it will be compatible 
## with make_tracklets and hence HelioLinC2

## MEW 31 Oct 2025

import numpy as np
import pandas as pd

def main():
    fname = 'c2017m4_atlas.csv'
    df = pd.read_csv(fname)
    nrows = len(df['ObjID'])
    df['ObsCode'] = np.full((nrows, 1), 'I11')
    
    outfilename = 'c2017m4_atlas_helio.csv'
    df.to_csv(outfilename, index=False)

if __name__ == '__main__':
    main()
