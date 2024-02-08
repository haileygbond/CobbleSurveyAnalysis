#%%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
import os
from datetime import date

survey_filepath = 'C:\\Users\\bondh\\Box\\Dynamic_Revetments_OSG_2022\\SurveyData\\3Processed\\Westport\\'
output_filepath = 'C:\\Users\\bondh\\Box\\HaileyBond\\phd\\arch cape\\surveys\\4_analysis\\'
output_filename = 'Westport_cobble_slopes.csv'
dates = [20231219]
transects = ['WORM100', '85', '84.8']
use_all_dates = True
use_all_transects = False

get_beach_slope = True
beach_slope_method = 'feature_code'
feature_code = 'rev_toe'

#%%
if use_all_dates:
    files = []
    for file in glob(survey_filepath + '*.csv'):
        files.append(file)
    survey_filenames = pd.Series([os.path.basename(x) for x in files])
elif len(dates)>1:
    files = []
    for date in dates:
        file = glob(survey_filepath + '*' + str(date) + '.csv')
        files.append(file[0])
    survey_filenames = pd.Series([os.path.basename(x) for x in files])
else:
    file = glob(survey_filepath + '*' + str(dates[0]) + '.csv')
    survey_filenames = pd.Series(os.path.basename(file[0])) # this might not work, check!
#%%
data_out = []
for survey_filename in survey_filenames:
    data_in = pd.read_csv(survey_filepath + survey_filename, header = 0,  names = ['x','y', 'z', 'dist', 'offset', 'feature', 'transect_name'])
    print(data_in)
    if get_beach_slope:
        if use_all_transects:
            transects = data_in['feature'].unique()
        for transect in transects:
            data = data_in.loc[data_in['transect_name'] == transect]
            if data['feature'].str.contains(feature_code, na=False).any():
                feature_dist = data[data['feature'].str.contains(feature_code, na=False)]['dist'].iloc[0]
                if len(data[data['dist']<feature_dist])>0:
                    beach = data[data['dist']>feature_dist]
                    beach = beach.reset_index()
                    beach.drop([0], axis=0, inplace=True) # temp fix for now to get rid of start point
                    if len(beach)>1:
                        beach_best_fit = np.polyfit(beach['dist'], beach['z'], 1)
                        slope = beach_best_fit[0]
                        surveydate = survey_filename[len(survey_filename)-12:len(survey_filename) - 4]
                        data_out.append([surveydate, transect, slope])
                    else:
                        continue
                else:
                    continue
            else:
                continue
            
# %%
data_out = pd.DataFrame(data_out, columns = ['date', 'transect_name', 'beach_slope'])
data_out.to_csv(output_filepath + output_filename, index=False)