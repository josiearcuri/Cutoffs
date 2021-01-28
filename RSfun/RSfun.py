import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def merge_dfs(cutoff_file, reach_file):
    """
    merges cutoff and reach data into single pandas dataframe
    
    INPUTS
    cutoff_file - column headings are: "reach	time	distance_downstream	bump"

    reaches_file - column headings are "reach	length	total_time	sinuosity	width"
    
    OUTPUT
    cutoffs_df - cutoff dataframe appended with matching values from reach dataframe based on reach key"""
    #load .csvs
    cutoffs_df = pd.read_csv(cutoff_file,sep = ',', header = 0,index_col=False)
    reaches_df = pd.read_csv(reach_file,sep = ',', header = 0,index_col=False)

    #count columns in reaches
    col_names = list(reaches_df.columns[1:])
    num_cols = len(col_names)

    #add columns to dataframe
    cutoffs_df[col_names] = 0

    #fill in with matching values
    for count, value in enumerate(reaches_df['reach'].values):
        rows = cutoffs_df[cutoffs_df['reach']==value].index
        for row in rows:
            cutoffs_df.loc[row, col_names]=reaches_df.loc[count]
    
    return cutoffs_df

def near_neigh(cutoffs_df):

    cutoffs_df['near_neigh_dist'] = cutoffs_df['length'].values
    cutoffs_df['near_neigh_time'] = cutoffs_df['total_time'].values
    for reach in np.unique(cutoffs_df['reach'].values):
        
        subset = cutoffs_df[cutoffs_df['reach']==reach]
        loc = subset['distance_downstream'].values
        time = subset['time'].values

        npts = len(subset.index)
        nn_space = np.zeros(shape=npts, dtype=np.double)
        nn_time = np.zeros(shape=npts, dtype=np.double)
        
        for i in range(npts):
            delta_time = abs(time[i] - time)
            delta_space = abs(loc[i] - loc)
            
            nn_space[i] = np.min(delta_space[np.nonzero(delta_space)])
            nn_time[i] = np.min(delta_time[np.nonzero(delta_time)])

        for count, row in enumerate(subset.index):
    
            cutoffs_df.loc[row,'near_neigh_dist'] = nn_space[count]
            cutoffs_df.loc[row,'near_neigh_time'] = nn_time[count]
    return cutoffs_df


def plot_nn_v1(cutoffs_df):
    x_s = np.asarray(cutoffs_df['near_neigh_dist']/cutoffs_df['length'])
    x_t = np.asarray(cutoffs_df['near_neigh_time']/cutoffs_df['total_time'])

    y = np.asarray(1/cutoffs_df['ncuts'])
    c = np.asarray(cutoffs_df['bump'])

    fig, (ax1, ax2) = plt.subplots(1,2, figsize = (4,8), sharey = True)

    ax1.scatter(x_s,y, c = c, s = .5)
    ax1.set_xlabel('distance to nearest neighbor / average reach length')
    ax1.set_ylabel('1 / number of neighboring cutoffs')
    ax1.set_xlim([0,1])
    ax1.set_ylim([0,1])
    ax1.plot([0,1], [0,1], 'k--', linewidth = .5)
    ax1.set_title('Observed Cutoff Spacing')

    ax2.scatter(x_t,y, c = c, s = .5)
    ax2.set_xlabel('time to nearest neighbor / observed time interval')
    ax2.plot([0,1], [0,1], 'k--', linewidth = .5)
    ax2.set_xlim([0,1])
    ax2.set_ylim([0,1])
    ax2.set_title('Observed Cutoff Timing')
    return fig

def plot_nn_v2(cutoffs_df):
    x_s = np.asarray(cutoffs_df['near_neigh_dist']/cutoffs_df['width'])
    x_t = np.asarray(cutoffs_df['near_neigh_time'])

    y_s = np.asarray(cutoffs_df['length']/(cutoffs_df['ncuts']*cutoffs_df['width']))
    y_t = np.asarray(cutoffs_df['total_time']/cutoffs_df['ncuts'])
    c = np.asarray(cutoffs_df['bump'])

    fig, (ax1, ax2) = plt.subplots(1,2, figsize = (4,8))

    ax1.scatter(x_s,y_s, c = c, s = .5)
    ax1.set_xlabel('distance to nearest neighbor/ch-w')
    ax1.set_ylabel('L /(ch-w*ncuts)')
    ax1.set_xlim([0,100])
    ax1.set_ylim([0,100])
    ax1.plot([0,100], [0,100], 'k--', linewidth = .5)
    ax1.set_title('Observed Cutoff Spacing')

    ax2.scatter(x_t,y_t, c = c, s = .5)
    ax2.set_xlabel('time to nearest neighbor')
    ax2.set_ylabel('T / ncuts')
    ax2.plot([0,35], [0,35], 'k--', linewidth = .5)
    ax2.set_xlim([0, 35])
    ax2.set_ylim([0,35])
    ax2.set_title('Observed Cutoff Timing')
    return fig