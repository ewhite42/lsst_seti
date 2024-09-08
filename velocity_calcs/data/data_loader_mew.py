## this program reads in data from a CSV file generated via
## a query to the LSST DP03 simulation database. This program 
## was originally written by Brian Rogers, and in the following 
## version tweaks have been added by Ellie White to accomodate 
## reading in different data columns. 

## 1 Sept 2024 -MEW

import pandas as pd
import numpy as np

from sklearn.preprocessing import MinMaxScaler
import os

# define gravitational constant in units: au*(au/day)^2/M_sun units
G_au_d = 2.95912741*10**(-4)

# define subspaces
subspace = ["a*", "i-z", "a", "sini", "e", "v", "vk-v"]
subspace_err = []
colorspace = ["a*", "i-z"]

def loader(
        fname = "LSST_DP03_data_velocities_1Sept2024.csv", # file path
        size = None, # Number of rows
        cols = ["a*", "i-z", "a", "sini", "e","v", "vk-v" ] # ensure values are present
):
 
    df = pd.read_csv(f"{os.path.dirname(os.path.realpath(__file__))}/"+fname, index_col=0)

    df["u-g"] = df["u_H"] - df["g_H"]
    df["u-r"] = df["u_H"] - df["r_H"]
    df["u-i"] = df["u_H"] - df["i_H"]
    df["u-z"] = df["u_H"] - df["z_H"]
    df["u-y"] = df["u_H"] - df["y_H"]


    df["g-r"] = df["g_H"] - df["r_H"]
    df["g-i"] = df["g_H"] - df["i_H"]
    df["g-z"] = df["g_H"] - df["z_H"]
    df["g-y"] = df["g_H"] - df["y_H"]



    df["r-i"] = df["r_H"] - df["i_H"]
    df["r-z"] = df["r_H"] - df["z_H"]
    df["r-y"] = df["r_H"] - df["y_H"]


    df["i-z"] = df["i_H"] - df["z_H"]
    df["i-y"] = df["i_H"] - df["z_H"]

    df["z-y"] = df["z_H"] - df["y_H"]

    # Need to calculate errors for colors?


    df["a*"] = .89*(df["g-r"]) + .45*(df["r-i"]) -.57

    df["sini"] = np.sin(df["incl"] * np.pi/180)
    df["a"] = df["q"]/(1-df["e"])
    
    ## add columns for the measured velocity, the predicted Keplerian 
    # velocity, and the difference between the two
    
    df["v"] = np.sqrt(df["heliocentricVX"]**2 + df["heliocentricVY"]**2 + df["heliocentricVZ"]**2)
    df["vk"] = np.sqrt((G_au_d)*((2/df["heliocentricDist"])-(1/df["a"])))
    df["vk-v"] = df["vk"] - df["v"]

    df = df.dropna(subset=cols)
    if size:
        return df.iloc[:size, :]

    else:
        return df.dropna(subset=subspace)

def preprocessed_subspace():
    return loader()[subspace]

def normalised_subspace():
    return scaler.fit_transform(preprocessed_subspace().to_numpy())

scaler = MinMaxScaler()
scaler.fit(preprocessed_subspace().to_numpy())

