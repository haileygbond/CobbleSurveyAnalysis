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
survey_filepath = 'C:\\Users\\bondh\\Box\\Dynamic_Revetments_OSG_2022\\SurveyData\\2Cleaned\\FC\\'
# if you want to process all surveys in folder, change the following to true
# all surveys in folder must be from same site
process_all_surveys = False
survey_filename = 'FC_20240105' + '.csv' # can leave blank if process_all_surveys = True
survey_columns = ['easting', 'northing', 'elevation', 'pointid', 'pointcode']

# transect information - must be csv file. 
# must have columns 'easting_start', 'northing_start', 'easting_end','northing_end', 'transect_name'
# there can be more columns and the mandatory columns can be in any order
# change all single backslashes to double backslashes, make sure there are backslashes at the end
transects_filepath = 'C:\\Users\\bondh\Box\\Dynamic_Revetments_OSG_2022\\SurveyData\\survey_transects\\'
transects_filename = 'FC_lines' + '.csv'
transects_columns = ['easting_start', 'northing_start', 'easting_end', 'northing_end', 'transect_name']

# output directory for processed data
# change all single backslashes to double backslashes, make sure there are backslashes at the end
output_filepath = 'C:\\Users\\bondh\\Box\\Dynamic_Revetments_OSG_2022\\SurveyData\\3Processed\\FC\\'

# processing options
max_offset = 3 #any points beyond this distance from line will be thrown out
interval = .5 # separates transect into bins of this size, and picks point with smallest offset within bin to save
min_points = 20 # transect must have this number of points to be included in file
include_feature_codes = True
featurecodes = ['cobble_toe','rev_toe', 'photo']

# processing figures
# change to true if you want to create figures
processing_figures_on = True
# output directory for quality control figures
# change all single backslashes to double backslashes, make sure there are backslashes at the end
fig_output_filepath = 'C:\\Users\\bondh\\Box\\Dynamic_Revetments_OSG_2022\\SurveyData\\4ProcessingFigs\\FC\\'
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
                # if no points meet criteria, add row of nans          
                if not intersection.any():
                    if include_feature_codes:
                        continue
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
                    if include_feature_codes:
                        dout_feature = data['pointcode'][tidx]
                    dout_transect_name = row['transect_name']
                    # add points to a list
                    if include_feature_codes:
                        ls.append([dout_x, dout_y, dout_z, dout_dist, dout_offset, dout_feature, dout_transect_name])
                    else:
                        ls.append([dout_x, dout_y, dout_z, dout_dist, dout_offset, dout_transect_name])
        num_points = len([x for x in ls if row['transect_name'] in x])
        if num_points < min_points:
            ls = [x for x in ls if row['transect_name'] not in x]
    # create dataframe from list
    if include_feature_codes:
        data_out = pd.DataFrame(data = ls, columns = ['x','y', 'z', 'dist', 'offset', 'feature', 'transect_name'])
    else:
        data_out = pd.DataFrame(data = ls, columns = ['x','y', 'z', 'dist', 'offset', 'transect_name'])
    data_out.to_csv(output_filepath + name, index=False)
    
##% processing figures
    fig1, ax1 = plt.subplot_mosaic(profile_fig_mosaic)
    ii = 0
    while ii < len(profile_axes_list)-1:
        for index, row in transect_lines.iterrows():
            df = data_out.loc[data_out['transect_name'] == row['transect_name']]
            
            x = [row['easting_start'], row['easting_end']]
            y = [row['northing_start'], row['northing_end']]
            ax1[profile_axes_list[ii]].plot(x,y)
            x = df['x']
            y =  df['y']
            ax1[profile_axes_list[ii]].scatter(x,y)
            
            ax1[profile_axes_list[ii]].set_xlabel('Easting (m)')
            ax1[profile_axes_list[ii]].set_ylabel('Northing (m)')
            ax1[profile_axes_list[ii]].set_title(row['transect_name'])
            ii = ii + 1
            
            x = df['dist']
            y = df['z']
            ax1[profile_axes_list[ii]].plot(x,y)
            ax1[profile_axes_list[ii]].set_xlabel('Distance (m)')
            ax1[profile_axes_list[ii]].set_ylabel('Elevation (m NAVD88)')
            ax1[profile_axes_list[ii]].set_title(row['transect_name'])
            ii = ii + 1
            
            if df.empty:
                x = [0,1]
            else:
                x = df['dist'].iloc[[0, -1]]
            y = [0,0]
            ax1[profile_axes_list[ii]].plot(x,y)
            
            x = df['dist']
            y = df['offset']
            ax1[profile_axes_list[ii]].plot(x,y)
            ax1[profile_axes_list[ii]].plot(x,y)
            ax1[profile_axes_list[ii]].set_xlabel('Distance (m)')
            ax1[profile_axes_list[ii]].set_ylabel('Offset from line (m)')
            ax1[profile_axes_list[ii]].set_title(row['transect_name'])
            ii = ii + 1
            
            fig1.tight_layout()
            fig_filename = fig_output_filepath + name.replace('.csv', '') + '_profiles'
            fig1.savefig(fig_filename)
            plt.close()       
            
    fig2, ax2 = plt.subplots()
    for index, row in transect_lines.iterrows():
        x = [row['easting_start'], row['easting_end']]
        y = [row['northing_start'], row['northing_end']]
        ax2.plot(x,y, zorder=1,color= '#1f77b4')
    
    x = data['easting']
    y = data['northing']
    ax2.scatter(x,y, zorder=2)
    
    x = data_out['x']
    y = data_out['y']
    ax2.scatter(x,y, zorder=3)
    
    xmin = np.min(data['easting']) - 5
    xmax = np.max(data['easting']) + 5
    ymin = np.min(data['northing']) - 5
    ymax = np.max(data['northing']) + 5
    
    ax2.set_xlim([xmin,xmax])
    ax2.set_ylim([ymin,ymax])
    
    ax2.set_xlabel('Easting (m)')
    ax2.set_ylabel('Northing (m)')
    
    txt = 'max offset = ' + str(max_offset) +'\nmin point = ' + str(min_points) + '\ninterval = ' + str(interval) +  '\nfeature codes = ' + str(featurecodes) + '\nprocessed on ' + str(date.today())
    
    ax2.text(0.1, 0.1, txt, transform=ax2.transAxes)
    
    fig_filename = fig_output_filepath + name.replace('.csv', '') + '_map'
    fig2.savefig(fig_filename)
    plt.close()       