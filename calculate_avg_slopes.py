import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
import os
from datetime import date

input_filepath ='C:\\Users\\bondh\Box\\HaileyBond\\phd\\arch cape\\surveys\\4_analysis\\cobble_slopes\\'

files = []
for file in glob(input_filepath + '*.csv'):
    files.append(file)
filenames = pd.Series(files)

data_out = []
for file in filenames:
    data_in = pd.read_csv(file, names = ['date', 'transect_name', 'slope'], header=1)
    transects = data_in['transect_name'].unique()
    for transect in transects:
        data = data_in.loc[data_in['transect_name'] == transect]
        avg_slope = data['slope'].mean()
        std_dev = data['slope'].std()
        data_out.append([transect, avg_slope, std_dev])
        