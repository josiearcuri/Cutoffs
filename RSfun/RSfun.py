import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
from scipy.stats import norm, linregress
from scipy import stats

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

def probs(cutoffs_df):
   
    n_bins_bump = 4
    n_bins_time = 10
    bad_idx = (cutoffs_df['duration']==0)|(cutoffs_df['bump']==0)|(cutoffs_df['cluster_flag']==0)
    # Generate a normal distribution, center at x=0 and y=5
    time = cutoffs_df['duration'].values[~bad_idx]#
    bump = cutoffs_df['bump'].values[~bad_idx]
    N_points = len(time)
    median_time = np.median(time)
    std_time = np.std(time)
 
    median_bump =np.median(bump)
    std_bump = np.std(bump)
    print(median_time)
    print(std_time)
    print(median_bump)
    print(std_bump)
    fig, axs = plt.subplots(1, 2, sharey=True, tight_layout=True)
    time_sig_points = np.arange(median_time-(2*std_time), median_time+(2*std_time), 1)
    bump_sig_points = np.arange(median_bump-(2*std_bump), median_time+(2*std_bump), 1)
    N, bins_time, patches_time = axs[0].hist(time, range = (2,12),bins=n_bins_time, density=True, label = 'n = '+str(N_points), histtype = 'bar', facecolor = 'lightgrey', edgecolor = 'k')
    N, bins_bump, patches_bump =axs[1].hist(bump, range = (1,5), bins=n_bins_bump, density = True,  label = 'n = '+str(N_points), histtype = 'bar', facecolor = 'lightgrey', edgecolor = 'k')
    # We can set the number of bins with the `bins` kwarg
   # N, bins_time, patches_time = axs[0].hist(time, range = (0,10),bins=n_bins_time, density=True, label = "95% CI", histtype = 'bar', facecolor = 'r', edgecolor = 'k')
    #N, bins_bump, patches_bump =axs[1].hist(bump, range = (0,5), bins=n_bins_bump, density = True,  label = "95% CI", histtype = 'bar', facecolor = 'r', edgecolor = 'k')

    #axs[0].legend(prop={'size': 10}, frameon = False)
    #axs[1].legend(prop={'size': 10}, frameon = False)
   # 
    #for p,b in zip(patches_time,bins_time[:-1]):
    #    if b-2<(median_time-(2*std_time)) or b>=(median_time+(2*std_time)):
    #        p.set_facecolor('lightgrey')
    #for p,b in zip(patches_bump,bins_bump[:-1]):
    #    if b-2<(median_bump-(2*std_bump)) or b>=(median_bump+(2*std_bump)):
    #        p.set_facecolor('lightgrey')
            


    axs[0].legend()
    axs[0].yaxis.set_major_formatter(PercentFormatter(xmax=1))
    
    axs[0].set_xlabel('duration [years]')
    axs[1].set_xlabel('bump [max(MR)/BR]')

    axs[0].set_title('P(duration)')
    axs[1].set_title('P(bump)')
    

    
    axs[0].spines['top'].set_visible(False)
    axs[0].spines['right'].set_visible(False)
    
    axs[1].spines['top'].set_visible(False)
    axs[1].spines['right'].set_visible(False)
    
    return fig

def plot_nn_v1(cutoffs_df):
    x_s = np.asarray(cutoffs_df['near_neigh_dist']/cutoffs_df['width'])
    x_t = np.asarray(cutoffs_df['near_neigh_time'])
    
    y = np.asarray(1/cutoffs_df['ncuts'])
    c = np.asarray(cutoffs_df['cluster_flag'])

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
    x_t = np.asarray(cutoffs_df['near_neigh_time']/cutoffs_df['BR'])

    y_s = np.asarray(cutoffs_df['length']/cutoffs_df['ncuts']/cutoffs_df['width']) 
    y_t = np.asarray(cutoffs_df['total_time']/cutoffs_df['ncuts']/cutoffs_df['BR']) - x_t
    c = np.asarray(cutoffs_df['cluster_flag'])
    
    slope_s, intercept_s, r_s, p_s, se_s = linregress(x_s[cutoffs_df['cluster_flag']==0], y_s[cutoffs_df['cluster_flag']==0])
    print(slope_s)
    print(intercept_s)
    print(r_s)
    fig, (ax1, ax2) = plt.subplots(1,2, figsize = (4,8))

    ax1.scatter(x_s,y_s, c = c, s = .5, label = '_nolegend_')
    ax1.set_xlabel('distance to nearest neighbor/ch-w')
    #ax1.set_ylabel('') #    ax1.set_ylabel('L /(ch-w*ncuts)')

    #ax1.set_xlim([0,100])
    #ax1.set_ylim([0,100])
    ax1.plot([0,1000], [0,intercept_s+(1000*slope_s)], 'r--', linewidth = .5, label  ='linear regression with r^2 ='+str(r_s**2))
    ax1.set_title('Observed Cutoff Spacing')
    ax1.legend()
    ax2.scatter(x_t,y_t, c = c, s = .5)
    ax2.set_xlabel('time to nearest neighbor')
    #ax2.set_ylabel('T / ncuts')
    ax2.plot([0,35], [0,35], 'k--', linewidth = .5)
    ax2.set_xlim([0, 35])
   # ax2.set_ylim([0,35])
    ax2.set_title('Observed Cutoff Timing')
    return fig
def nn_probs(cutoffs_df):
    N_points = len(cutoffs_df.index)

    n_bins_t = 20

    # Generate a normal distribution, center at x=0 and y=5
    x_t_clust = np.asarray((cutoffs_df['near_neigh_time'])[cutoffs_df['cluster_flag']==1]) 
    #-(cutoffs_df['total_time']/cutoffs_df['ncuts'])
    x_t_notclust = np.asarray((cutoffs_df['near_neigh_time'])[cutoffs_df['cluster_flag']==0])
    #-(cutoffs_df['total_time']/cutoffs_df['ncuts'])
    median_t = np.median(x_t_clust)
    std_t = np.std(x_t_clust)
 
    median_t_not =np.median(x_t_notclust)
    std_t_not = np.std(x_t_notclust)
    

    t,p = stats.ttest_ind_from_stats(mean1=median_t, std1=std_t, nobs1=len(x_t_clust),
                     mean2=median_t_not, std2=std_t_not, nobs2=len(x_t_notclust))
    #print(stats.ks_2samp(x_t_clust, x_t_notclust))
    print(t)
    print(p)
    print(len(x_t_clust))
    print(len(x_t_notclust))
    fig, axs = plt.subplots(1, 1, figsize = (8,4), tight_layout=True)
    axs.hist(x_t_clust,range=(0,30), bins=n_bins_t, color= 'r', histtype='step', cumulative=True, density = True, label = "clustered")
    axs.hist(x_t_notclust, range=(0,30),bins=n_bins_t, color= 'k', histtype='step', cumulative=True, density = True, label = "not clustered")
    


    axs.yaxis.set_major_formatter(PercentFormatter(xmax=1))
    axs.set_xlabel('years to nearest neighbor')
    axs.legend(loc = "lower center", frameon = False)
    axs.set_title('cumulative P(nn_t)')
    
    axs.spines['top'].set_visible(False)
    axs.spines['right'].set_visible(False)
    
    
    return fig
def plot_corr(cutoffs_df):
    bad_idx = (cutoffs_df['duration']==0)|(cutoffs_df['bump']==0)#|(cutoffs_df['cluster_flag']==0)
    # Generate a normal distribution, center at x=0 and y=5
    x = cutoffs_df['duration'].values[~bad_idx]#
    y = cutoffs_df['near_neigh_time'][~bad_idx]
    
    fig, axs = plt.subplots(1, 1, figsize = (8,8), tight_layout=True)
    axs.scatter(x,y, c=cutoffs_df['cluster_flag'][~bad_idx])
    return fig

    