# process survey data to get transects
# data should be in a csv file with columns: easting, northing, elevation, pointID, point code
# any surveying mistakes should be removed from data (i.e., fix incorrect point codes, delete duplicate points)
# data should be stored at file path input below

# code translated from Dept of Ecology matlab code 1/2024 by H Bond

#%%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
import os
from datetime import date

#%% INPUTS - ONLY CHANGE THIS CELL

# survey information - must be csv file(s)
# survey must have columns named 'easting', 'northing', 'elevation','pointid','pointcode'
# there can be more columns and the mandatory ones can be in any order
# change all single backslashes to double backslashes, make sure there are backslashes at the end
survey_filepath = 'C:\\Users\\bondh\\Box\\Dynamic_Revetments_OSG_2022\\SurveyData\\2Cleaned\\Westport\\'
# if you want to process all surveys in folder, change the following to true
# all surveys in folder must be from same site
process_all_surveys = True
survey_filename = '' + '.csv' # can leave blank if process_all_surveys = True
survey_columns = ['easting', 'northing', 'elevation', 'pointid', 'pointcode']

# transect information - must be csv file. 
# must have columns 'easting_start', 'northing_start', 'easting_end','northing_end', 'transect_name'
# there can be more columns and the mandatory columns can be in any order
# change all single backslashes to double backslashes, make sure there are backslashes at the end
transects_filepath = 'C:\\Users\\bondh\Box\\Dynamic_Revetments_OSG_2022\\SurveyData\\survey_transects\\'
transects_filename = 'Westport_lines' + '.csv'
transects_columns = ['easting_start', 'northing_start', 'easting_end', 'northing_end', 'transect_name']

# output directory for processed data
# change all single backslashes to double backslashes, make sure there are backslashes at the end
output_filepath = 'C:\\Users\\bondh\\Box\\Dynamic_Revetments_OSG_2022\\SurveyData\\3Processed\\Westport\\'

# processing options
max_offset = 3 #any points beyond this distance from line will be thrown out
max_gap = 5 # if gap between points is bigger than this, insert nans
interval = .5 # separates transect into bins of this size, and picks point with smallest offset within bin to save
include_feature_codes = True
featurecodes = ['cobble_toe','rev_toe', 'photo']

# processing figures
# change to true if you want to create figures
processing_figures_on = True
# output directory for quality control figures
# change all single backslashes to double backslashes, make sure there are backslashes at the end
fig_output_filepath = 'C:\\Users\\bondh\\Box\\Dynamic_Revetments_OSG_2022\\SurveyData\\4ProcessingFigs\\Westport\\'
# each transect will create 3 subplots in final figure, use below variable to control their location
# look up subplotmosaic to understand how to set this up
profile_fig_mosaic ="""
    AB..
    ACGH
    DEGI
    DF..
    """
# change number of letters to 3 * number of transects at your site
profile_axes_list = ['A','B','C','D','E','F','G','H','I']
# Prep visualisation
S = 14
M = 18
L = 20

plt.rc('font', size=S)          # controls default text sizes
plt.rc('axes', titlesize=S)     # fontsize of the axes title
plt.rc('axes', labelsize=M)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=S)    # fontsize of the tick labels
plt.rc('ytick', labelsize=S)    # fontsize of the tick labels
plt.rc('legend', fontsize=M)    # legend fontsize
plt.rc('figure', titlesize=L)  # fontsize of the figure title
plt.rc('axes', titlesize=L)  # fontsize of the subplot title
plt.rc('figure', figsize=[35, 18])  # fontsize of the subplot title



#%% read in data
if process_all_surveys:
    files = []
    for file in glob(survey_filepath + '*.csv'):
        files.append(file)
    survey_filename = pd.Series([os.path.basename(x) for x in files])
else:
    survey_filename = pd.Series(survey_filename)
for name in survey_filename:
    file = survey_filepath + name
    data = pd.read_csv(file, header = None, names = survey_columns, index_col=False)

    #%% read in transect lines
    transect_filename = transects_filepath + transects_filename
    transect_lines = pd.read_csv(transect_filename, names = transects_columns)

    #%% rotate data to align with transect lines and save points with smallest offsets

    # initialize list
    ls = []
    # loop throught transects
    for index, row in transect_lines.iterrows():
        # calculate angle of transect
        m = (row['northing_end'] - row['northing_start'])/(row['easting_end'] - row['easting_start'])
        theta = np.arctan(m)
        # rotate transect points
        xt = data['easting']*np.cos(theta) + data['northing'] * np.sin(theta)
        yt = -data['easting']*np.sin(theta) + data['northing'] * np.cos(theta)
        # rotate collected data points
        xl = [row['easting_start']*np.cos(theta) + row['northing_start'] * np.sin(theta),
            row['easting_end']*np.cos(theta) + row['northing_end'] * np.sin(theta)]
        yl = [-row['easting_start']*np.sin(theta) + row['northing_start'] * np.cos(theta),
            -row['easting_end']*np.sin(theta) + row['northing_end'] * np.cos(theta)]
        # calculate distance along line (x) and offset from line (y)
        dist = xt - xl[1] 
        offset = yt - yl[1]
        # find points that are less than input max offset
        oidx = np.abs(offset)<max_offset
        # create bins along transect based on input interval
        xi = pd.Series(np.arange(min(dist), max(dist),interval))
        # loop through bins
        for idx, rw in xi.items():
            # prevent out of range error
            if idx + 1 > len(xi) - 1:
                continue
            else:
                # find points within bin
                didx = (dist> xi[idx]) & (dist < xi[idx + 1])
                # find points that are within bin and have offset less than max offset
                # check first if there is a point with a desired feature code that meets criteria
                if include_feature_codes:
                    fidx = data['pointcode'].isin(featurecodes)
                    intersection = np.logical_and(oidx, didx, fidx)
                    # if there isn't a feature code point, use a regular topo point
                    if not intersection.any():
                        intersection = np.logical_and(oidx, didx) 
                else:
                    intersection = np.logical_and(oidx, didx)
                # prevent error from empty series            
                if not intersection.any():
                    continue
                else:
                    # find index of point within bin that has smallest offset
                    tidx = np.abs(offset[intersection]).idxmin()
                    # use that index to get all other desired columns
                    dout_x = data['easting'][tidx]
                    dout_y = data['northing'][tidx]
                    dout_z = data['elevation'][tidx]
                    dout_dist = dist[tidx]
                    dout_offset = offset[tidx]
                    dout_feature = data['pointcode'][tidx]
                    dout_transect_name = row['transect_name']
                    # if there is a gap bigger than max_gap, insert row of nans
                    if len(ls)>0:
                        gap = dout_dist - ls[-1][3]
                        if gap > max_gap:
                            ls.append(np.full(7, np.nan).tolist())
                    # add points to a list
                    
                    ls.append([dout_x, dout_y, dout_z, dout_dist, dout_offset, dout_feature, dout_transect_name])
    # create dataframe from list
           
    data_out = pd.DataFrame(data = ls, columns = ['x','y', 'z', 'dist', 'offset', 'feature', 'transect_name'])
    data_out.to_csv(output_filepath + name, index=False)
    
##% processing figures
    fig, axd = plt.subplot_mosaic(profile_fig_mosaic)
    ii = 0
    while ii < len(profile_axes_list)-1:
        for index, row in transect_lines.iterrows():
            df = data_out.loc[data_out['transect_name'] == row['transect_name']]
            
            x = [row['easting_start'], row['easting_end']]
            y = [row['northing_start'], row['northing_end']]
            axd[profile_axes_list[ii]].plot(x,y)
            x = df['x']
            y =  df['y']
            axd[profile_axes_list[ii]].scatter(x,y)
            
            axd[profile_axes_list[ii]].set_xlabel('Easting (m)')
            axd[profile_axes_list[ii]].set_ylabel('Northing (m)')
            axd[profile_axes_list[ii]].set_title(row['transect_name'])
            ii = ii + 1
            
            x = df['dist']
            y = df['z']
            axd[profile_axes_list[ii]].plot(x,y)
            axd[profile_axes_list[ii]].set_xlabel('Distance (m)')
            axd[profile_axes_list[ii]].set_ylabel('Elevation (m NAVD88)')
            axd[profile_axes_list[ii]].set_title(row['transect_name'])
            ii = ii + 1
            
            if df.empty:
                x = [0,1]
            else:
                x = df['dist'].iloc[[0, -1]]
            y = [0,0]
            axd[profile_axes_list[ii]].plot(x,y)
            
            x = df['dist']
            y = df['offset']
            axd[profile_axes_list[ii]].plot(x,y)
            axd[profile_axes_list[ii]].plot(x,y)
            axd[profile_axes_list[ii]].set_xlabel('Distance (m)')
            axd[profile_axes_list[ii]].set_ylabel('Offset from line (m)')
            axd[profile_axes_list[ii]].set_title(row['transect_name'])
            ii = ii + 1
            
            fig.tight_layout()
            fig_filename = fig_output_filepath + name.replace('.csv', '') + '_profiles'
            fig.savefig(fig_filename)
            plt.close()       
            
    fig, ax = plt.subplots()
    x = data['easting']
    y = data['northing']
    ax.scatter(x,y)
    
    x = data_out['x']
    y = data_out['y']
    ax.scatter(x,y)
    
    ax.set_xlabel('Easting (m)')
    ax.set_ylabel('Northing (m)')
    
    txt = 'max offset = ' + str(max_offset) +'\nmax gap = ' + str(max_gap) + '\ninterval = ' + str(interval) +  '\nfeature codes = ' + str(featurecodes) + '\nprocessed on ' + str(date.today())
    
    ax.text(0.1, 0.1, txt, transform=ax.transAxes)
    
    fig_filename = fig_output_filepath + name.replace('.csv', '') + '_map'
    fig.savefig(fig_filename)
    plt.close()       