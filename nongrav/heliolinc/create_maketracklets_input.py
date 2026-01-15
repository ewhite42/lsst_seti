## this program adds any needed columns to the
## Sorcha output file so that it will be compatible 
## with make_tracklets and hence HelioLinC2

## MEW 31 Oct 2025

import numpy as np
import pandas as pd

def create_maketracklets_input(fname, ofname):
    df = pd.read_csv(fname)
    nrows = len(df['ObjID'])
    df['ObsCode'] = np.full((nrows, 1), 'I11')

    df.to_csv(ofname, index=False)
    
def main():
    fname = '/home/ellie/research/lsst/sorcha_output/c2012v1_panstarrs/c2012v1_panstarrs_ng_1mo.csv'
    ofname = '/home/ellie/research/lsst/sorcha_output/c2012v1_panstarrs/c2012v1_panstarrs_ng_1mo_helio.csv'
    create_maketracklets_input(fname,ofname)

if __name__ == '__main__':
    main()
