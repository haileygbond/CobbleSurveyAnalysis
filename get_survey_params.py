#%%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
import os
from datetime import date

site = 'Westport'

survey_filepath = 'C:\\Users\\bondh\\Box\\Dynamic_Revetments_OSG_2022\\SurveyData\\3Processed\\' + site + '\\'
output_filepath = 'C:\\Users\\bondh\\Box\\HaileyBond\\phd\\analysis_field\\survey_parameters\\'

dates = [20231219]
transects = ['FC1', 'FC2', 'FC3']
use_all_dates = True
use_all_transects = True

get_beach_slope = True
beach_slope_method = 'feature_code'
feature_code = 'rev_toe'

get_cobble_slope = True
cobble_slope_method = 'feature_code'
feature_code = 'rev_toe'

get_cobble_toe_elevation = True
cobble_slope_method = 'feature_code'
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
#%% beach slope calculation
if get_beach_slope:
    data_out = []
    for survey_filename in survey_filenames:
        data_in = pd.read_csv(survey_filepath + survey_filename, header = 0,  names = ['x','y', 'z', 'dist', 'offset', 'feature', 'transect_name'])
        if use_all_transects:
            transects = data_in['transect_name'].unique()
        for transect in transects:
            data = data_in.loc[data_in['transect_name'] == transect]
            if data['feature'].isna().all():
                continue
            if data['feature'].str.contains(feature_code, na=False).any():
                feature_dist = data[data['feature'].str.contains(feature_code, na=False)]['dist'].iloc[0]
                if len(data[data['dist']<feature_dist])>0:
                    beach = data[data['dist']<feature_dist]
                    beach = beach.reset_index()
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
    data_out = pd.DataFrame(data_out, columns = ['date', 'transect_name', 'beach_slope'])
    data_out.to_csv(output_filepath + site + '_beach_slopes.csv', index=False)
#%% cobble slope slope calculation
if get_cobble_slope:
    data_out = []
    for survey_filename in survey_filenames:
        data_in = pd.read_csv(survey_filepath + survey_filename, header = 0,  names = ['x','y', 'z', 'dist', 'offset', 'feature', 'transect_name'])
        if use_all_transects:
            transects = data_in['feature'].unique()
        for transect in transects:
            data = data_in.loc[data_in['transect_name'] == transect]
            if data['feature'].isna().all():
                continue
            if data['feature'].str.contains(feature_code, na=False).any():
                feature_dist = data[data['feature'].str.contains(feature_code, na=False)]['dist'].iloc[0]
                if len(data[data['dist']<feature_dist])>0:
                    cobble = data[data['dist']>feature_dist]
                    cobble = cobble.reset_index()
                    cobble = cobble.dropna(subset = ['dist', 'z'])
                    if len(cobble)>1:
                        cobble_best_fit = np.polyfit(cobble['dist'], cobble['z'], 1)
                        slope = cobble_best_fit[0]
                        surveydate = survey_filename[len(survey_filename)-12:len(survey_filename) - 4]
                        data_out.append([surveydate, transect, slope])
                    else:
                        continue
                else:
                    continue
            else:
                continue
    data_out = pd.DataFrame(data_out, columns = ['date', 'transect_name', 'cobble_slope'])
    data_out.to_csv(output_filepath + site + '_cobble_slopes.csv', index=False)
    
#%% cobble slope slope calculation
if get_cobble_toe_elevation:
    data_out = []
    for survey_filename in survey_filenames:
        data_in = pd.read_csv(survey_filepath + survey_filename, header = 0,  names = ['x','y', 'z', 'dist', 'offset', 'feature', 'transect_name'])
        if use_all_transects:
            transects = data_in['feature'].unique()
        for transect in transects:
            data = data_in.loc[data_in['transect_name'] == transect]
            if data['feature'].isna().all():
                continue
            if data['feature'].str.contains(feature_code, na=False).any():
                feature_elev = data[data['feature'].str.contains(feature_code, na=False)]['z'].iloc[0]
                surveydate = survey_filename[len(survey_filename)-12:len(survey_filename) - 4]
                data_out.append([surveydate, transect, feature_elev])
    data_out = pd.DataFrame(data_out, columns = ['date', 'transect_name', 'cobble_toe_elevation'])
    data_out.to_csv(output_filepath + site + '_cobble_toe_elevations.csv', index=False)