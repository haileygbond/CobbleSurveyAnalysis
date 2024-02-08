import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
import os
from datetime import date

%matplotlib qt5

input_filepath ='C:\\Users\\bondh\Box\\HaileyBond\\phd\\arch cape\\surveys\\4_analysis\\'

files = []
for file in glob(survey_filepath + '*.csv'):
    files.append(file)
filenames = pd.Series(files)

for file in filenames
    pd.read_csv(file, columns = ['date', 'transect_name', 'slope'])
    transects = data_in['feature'].unique()
    for transect in transects:
        data = data = data_in.loc[data_in['transect_name'] == transect]