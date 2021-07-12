import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
import cartopy as cp
#import geopandas as gpd
from scipy.stats import norm, linregress
from scipy import stats
from scipy.spatial import distance
#import cartopy.crs as ccrs
#import contextily
#import rasterio
#from rasterio.plot import show as rioshow
#from rasterio.mask import mask
import math
from scipy.interpolate import InterpolatedUnivariateSpline, UnivariateSpline
#from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from matplotlib.colors import ListedColormap, LinearSegmentedColormap, Normalize
from matplotlib import cm, transforms, patches
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.font_manager as fm
#from shapely.geometry import Polygon, box
from scipy.ndimage.interpolation import rotate
from matplotlib.ticker import AutoMinorLocator
import os

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
   
    n_bins_bump = 8
    n_bins_time = 10
    bad_idx_dur = (cutoffs_df['duration']==0)
    bad_idx_mag = (cutoffs_df['bump']==0)#|(cutoffs_df['cluster_flag']==0)
    # Generate a normal distribution, center at x=0 and y=5
    dur = cutoffs_df['duration'].values[~bad_idx_dur]#
    mag = cutoffs_df['bump'].values[~bad_idx_mag]-1
    #N_points = len(dur)
    median_time = np.median(dur)
    std_time = np.std(dur)
 
    median_bump =np.median(mag)
    std_bump = np.std(mag)
    print(np.percentile(dur, 50))
    print(np.percentile(dur, 75))
    print(np.max(dur))
    print(np.percentile(mag, 50))
    print(np.percentile(mag, 75))
    print(np.max(mag))
    fig, axs = plt.subplots(1, 2, figsize = (6,4), tight_layout=True)

    #N, bins_time, patches_time = axs[0].hist(dur, range = (2,12),bins=n_bins_time, label = 'n = '+str(N_points), histtype = 'bar', facecolor = 'lightgrey', edgecolor = 'k')
    hist = axs[0].boxplot(dur, vert=False)
    hist = axs[1].boxplot(mag, vert=False)
    #N, bins_bump, patches_bump =axs[1].hist(mag, range = (1,5), bins=n_bins_bump,  label = 'n = '+str(N_points), histtype = 'bar', facecolor = 'lightgrey', edgecolor = 'k')
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
            


   
    #axs[0].yaxis.set_major_formatter(PercentFormatter(xmax=1))
    
    axs[0].set_xlabel('years')
    axs[1].set_xlabel('max(% diff)')

    axs[0].set_title('NE Duration'+' n = '+str(len(dur)))
    axs[1].set_title('NE Magnitude'+' n = '+str(len(mag)))
    

    axs[0].set_yticks([])
    axs[1].set_yticks([])
    axs[0].spines['top'].set_visible(False)
    axs[0].spines['right'].set_visible(False)
    #axs[0].spines['left'].set_visible(False)                      
    
    axs[1].spines['top'].set_visible(False)
    axs[1].spines['right'].set_visible(False)
    #axs[1].spines['left'].set_visible(False)                      
    
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
    plt.rcParams["font.family"] = "Arial"
    font = {'family' : 'Arial',
        'size'   : 8}

    plt.rc('font', **font)
    N_points = len(cutoffs_df.index)

    n_bins_t = 30

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
    print("t = " + str(t))
    print("p = " + str(p))
    print("n_inside = " + str(len(x_t_clust)))
    print("n_outside = " + str(len(x_t_notclust)))
    print("mean not clustered = " + str(median_t_not))
    print("mean clustered = " + str(median_t))
    fig, axs = plt.subplots(1, 1, figsize = (5, 4), tight_layout = True)
    axs.hist(x_t_clust,range=(0,30), bins=range(0,30,2),color= 'r', linewidth = 1.5, histtype='stepfilled',label = "inside n="+str(len(x_t_clust))) #,cumulative=True, density = True, 
    axs.hist(x_t_notclust, range=(0,30),bins=range(0,30,2), color= 'grey', linewidth = 1.5, histtype='stepfilled',  label = "outside n="+str(len(x_t_notclust))) #cumulative=True, density = True,
    axs.set_title("Years Between Cutoffs", fontsize = 10)
    #axs.set_yticks([.5, 1])
    #axs.yaxis.set_major_formatter(PercentFormatter(xmax=1))
    axs.set_xlabel('t [years]', fontsize = 9)
    axs.legend(bbox_to_anchor=(1, 1),loc = "upper right", title = "proximity to \nnonlocal effects",  edgecolor = 'w', fontsize = 9, title_fontsize = 10)
    axs.set_ylabel('count (nn_t < t)', fontsize = 9)
    axs.xaxis.set_minor_locator(AutoMinorLocator())
    axs.yaxis.set_minor_locator(AutoMinorLocator())
    axs.tick_params(which='minor', color='grey')
    axs.spines['top'].set_visible(False)
    axs.spines['right'].set_visible(False)
    axs.set_xlim([0, 30.1])
    axs.set_ylim([0, 7])
    plt.grid(True, which = 'major')
    
    return fig
def plot_corr(cutoffs_df):
    bad_idx = (cutoffs_df['duration']==0)|(cutoffs_df['bump']==0)#|(cutoffs_df['cluster_flag']==0)
    # Generate a normal distribution, center at x=0 and y=5
    x = cutoffs_df['duration'].values[~bad_idx]#
    y = cutoffs_df['near_neigh_time'][~bad_idx]
    
    fig, axs = plt.subplots(1, 1, figsize = (8,8), tight_layout=True)
    axs.scatter(x,y, c=cutoffs_df['cluster_flag'][~bad_idx])
    return fig
def plot_sample_locations(gdf, topimage):
    plt.rcParams["font.family"] = "Arial"
    
    fig = plt.figure(figsize = (3,3))
    ax = fig.add_subplot(projection=ccrs.PlateCarree())
    #ax_upper = fig.add_subplot(gs[1,0], projection=ccrs.PlateCarree(), frameon = False)
    #world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))

    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    gdf.plot(ax=ax, marker = '*', facecolors='w', edgecolors='k', linewidths = .5, markersize = 100, label = 'samples')

    basin = gpd.read_file("D:/hydrosheds/sa_bas_30s_beta/sa_bas_30s_beta.shp")
    
    ARB = basin[basin.AREA_SQKM == np.max(basin.AREA_SQKM)]   
    #ARB.boundary.plot(ax = ax,linewidth = .5, alpha = .5,color = 'lightgrey', facecolor = 'None', edgecolor = 'darkgrey')
    ax.legend(bbox_to_anchor=(.1, .1), handletextpad = .1, fontsize = 10, loc='lower left', frameon = False, borderpad = 0, borderaxespad=0)
    #extent = ARB.geometry.bounds
    #ax.add_image(raster)
   # ax.set_extent((extent.minx, extent.miny, extent.maxx, extent.maxy))

    
#bcbddc
#756bb1
    current_cmap = cm.get_cmap('cividis')
    #current_cmap = get_continuous_cmap(['#7fcdbb', '#2c7fb8'])
    newcolors = current_cmap(np.linspace(0, 1, 101))

    back = newcolors[0, :]
    newcolors[1:25, :] = back
    newcolors[0, :] = np.array([0,0,0,0])
    current_cmap = ListedColormap(newcolors)
    

    with rasterio.open("sample_data/RSCutoffData/ARB_occurrence_km_unmasked.tif") as r:
        out_image, out_transform = mask(r, ARB.geometry, crop=True, nodata = -1)
        rioshow(out_image, transform=out_transform, ax=ax, cmap = current_cmap)

    ax.invert_yaxis()
   
    x_ticks = ax.get_xticks()
    ax.set_xticks(x_ticks[0::2])
    y_ticks = ax.get_yticks()
    ax.set_yticks(y_ticks[1::2])
    

    lon_formatter = LongitudeFormatter(degree_symbol='',
                                       dateline_direction_label=True)
    lat_formatter = LatitudeFormatter(degree_symbol='')

    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    
    fontprops = fm.FontProperties(size=6)

    scalebar = AnchoredSizeBar(ax.transData, 
                           500/111, '500 km', 'upper left', 
                           pad=2.5,
                           color='dimgrey',
                           frameon=False, 
                           size_vertical=.1,
                           fontproperties=fontprops)
    ax.add_artist(scalebar)
    #y
    ax.tick_params(axis='y', which = 'both', direction= 'in', pad = -15, colors='dimgrey', right=False, left = True, labelright=False, labelleft = True)
    #x
    ax.tick_params(axis='x', which = 'both', direction= 'in', pad = -10, colors='dimgrey',  top=True, bottom = False, labeltop=True, labelbottom = False)
    
    ax.tick_params(which='minor', color='lightgrey')
    x, y, arrow_length = 0.15, 0.35, .1
    ax.annotate('N', xy=(x, y), xytext=(x, y-arrow_length),
            arrowprops=dict(facecolor='dimgrey',edgecolor='None', arrowstyle='simple'),
            ha='center', va='center', fontsize=6,color = 'dimgrey',
            xycoords='axes fraction')
    ax.set_extent((-85, -50, -24, 8), crs=ccrs.PlateCarree())
    plt.setp(ax.get_xticklabels(), fontsize=6)
    plt.setp(ax.get_yticklabels(), fontsize=6)
    
    ax.outline_patch.set_visible(False)
    ax.spines['left'].set_visible(True)  
    ax.spines['top'].set_visible(True) 
    ax.spines['left'].set_linewidth(1)  
    ax.spines['top'].set_linewidth(1) 
    ax.spines['left'].set_color('grey')  
    ax.spines['top'].set_color('grey') 

        
    return fig


def get_continuous_cmap(hex_list, float_list=None):
    ''' creates and returns a color map that can be used in heat map figures.
        If float_list is not provided, colour map graduates linearly between each color in hex_list.
        If float_list is provided, each color in hex_list is mapped to the respective location in float_list. 
        
        Parameters
        ----------
        hex_list: list of hex code strings
        float_list: list of floats between 0 and 1, same length as hex_list. Must start with 0 and end with 1.
        
        Returns
        ----------
        colour map'''
    rgb_list = [rgb_to_dec(hex_to_rgb(i)) for i in hex_list]
    if float_list:
        pass
    else:
        float_list = list(np.linspace(0,1,len(rgb_list)))
        
    cdict = dict()
    for num, col in enumerate(['red', 'green', 'blue']):
        col_list = [[float_list[i], rgb_list[i][num], rgb_list[i][num]] for i in range(len(float_list))]
        cdict[col] = col_list
    cmp = LinearSegmentedColormap('my_cmp', segmentdata=cdict, N=256)
    return cmp
def hex_to_rgb(value):
    '''
    Converts hex to rgb colours
    value: string of 6 characters representing a hex colour.
    Returns: list length 3 of RGB values'''
    value = value.strip("#") # removes hash symbol if present
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))


def rgb_to_dec(value):
    '''
    Converts rgb to decimal colours (i.e. divides each value by 256)
    value: list (length 3) of RGB values
    Returns: list (length 3) of decimal values'''
    return [v/256 for v in value]
def plot_sample(topimage):
    plt.rcParams["font.family"] = "Arial"

    fig = plt.figure(figsize = (8,2))
    ax_upper = fig.add_subplot(projection=ccrs.PlateCarree(), frameon = False)
   

    
    #top_cmap = get_continuous_cmap(['#9e9ac8', '#54278f', '000000'])
    top_cmap = current_cmap = cm.get_cmap('cividis')#get_continuous_cmap(['#edf8b1','#7fcdbb', '#2c7fb8'])

    newcolors = top_cmap(np.linspace(0, 1, 101))
    white = np.array([0, 0, 0, 0])
    newcolors[0, :] = white
    top_cmap = ListedColormap(newcolors)
    
    #box_up = Polygon([(r.bounds[0], r.bounds[1]),(r.bounds[2], r.bounds[1]), (r.bounds[2], r.bounds[3]), (r.bounds[1], r.bounds[3])])
    fontprops = fm.FontProperties(size=8)

 
    with rasterio.open(topimage, 'r+') as r:
        image_hidden = ax_upper.imshow(r.read(1), 
                         cmap=top_cmap, 
                         vmin=0, 
                         vmax=100)
        img=rioshow(rotate(r.read(1), angle = -48), transform=r.transform, ax=ax_upper, vmin = 0, vmax = 100, cmap = top_cmap) 
    scalebar = AnchoredSizeBar(ax_upper.transData,
                        5/111, '5 km', 'lower center', 
                        pad=.1,
                        color='k',
                        frameon=False,
                        size_vertical=.5/555,
                        fontproperties=fontprops)
    #scalebar = patches.Rectangle((-65.2,-16.45),5/111, .1, color = 'k')

    #scalebar.set_transform(t_end)
    ax_upper.add_artist(scalebar)
    #scalebar.set_transform(t)
    cbar = fig.colorbar(image_hidden, ax = ax_upper, shrink = .6)
    cbar.ax.set_ylabel("% Occurrence", rotation=-90, va="bottom")
    cbar.set_ticks([0, 100])
    ax_upper.set_extent((-65.45, -65.12, -16.48, -16.40), crs=ccrs.PlateCarree())
    #x_ticks = ax_upper.get_xticks()
    #ax_upper.set_xticks(x_ticks[1:-1])
    #y_ticks = ax_upper.get_yticks()
    #ax_upper.set_yticks(y_ticks[1:-1])
    return fig
def compute_curvature(x,y):
    """function for computing first derivatives and curvature of a curve (centerline)
    x,y are cartesian coodinates of the curve
    outputs:
    dx - first derivative of x coordinate
    dy - first derivative of y coordinate
    ds - distances between consecutive points along the curve
    s - cumulative distance along the curve
    curvature - curvature of the curve (in 1/units of x and y)"""
    dx = np.gradient(x) # first derivatives
    dy = np.gradient(y)      
    ds = np.sqrt(dx**2+dy**2)
    s = np.hstack((0,np.cumsum(ds[1:])))
    ddx = np.gradient(dx) # second derivatives 
    ddy = np.gradient(dy) 
    curvature = (dx*ddy-dy*ddx)/((dx**2+dy**2)**1.5)
    return curvature,s
def load_data(dirname):
    """
    simplified from curvaturepy @zoltansylvester
    """
   
    fnames = os.listdir(dirname)
    fnames.sort()
    clxs = []
    clys = []
    ages = []
    downstream = []
    curvatures = []
    for i in range(0,len(fnames)):
        df = pd.read_csv(dirname+fnames[i])
        x = np.array(df.iloc[:,0])*30
        y = -1*np.array(df.iloc[:,1])*30
        points = np.vstack([x, y]).T
        curv, s = compute_curvature(x,y)
        # Build a list of the spline function, one for each dimension:
        #splines = [UnivariateSpline(s, coords, k=3, s=5) for coords in points.T]

        # Compute the spline for the smoothed(sampled) distances:
        #points_fitted = np.vstack([spl(np.linspace(0, s[-1],int(s[-1]/(25)))) for spl in splines])
        #curv, s = compute_curvature(x,y)
        clxs.append(x)
        clys.append(y)
        curvatures.append(curv)
        downstream.append(s)
        ages.append(int(fnames[i][:4]))##-8:-4]))
        
    current_cmap = cm.get_cmap('cividis')
    #current_cmap = get_continuous_cmap(['#7fcdbb', '#2c7fb8'])
    
    newcolors = current_cmap(np.linspace(0,1, ages[-1]-ages[0]))
    current_cmap = ListedColormap(newcolors)
    
    fig, ax = plt.subplots(1,1, figsize=(8,8))
    
    

    for i in range(len(clxs)):
        ax.plot(clxs[i],clys[i],linewidth = .5,alpha = .75, c=current_cmap(i))
    plt.axis('equal')
    cbar = fig.colorbar(cm.ScalarMappable(norm=Normalize(ages[0], ages[-1]), cmap=current_cmap), ax=ax, ticks = ages[::4], shrink = .75)
    #cbar = plt.colorbar(lines, cax=ax,orientation="vertical")
#cbar.ax.set_ylabel("year", rotation=-90, va="bottom")
    cbar.ax.set_title('year', fontsize = 10)
    cbar.ax.tick_params(labelsize=10) 
    return fnames,clxs,clys,curvatures,ages, downstream

def get_migr_rate(x1,y1,x2,y2,age1,age2):
    """use dynamic time warping to correlate centerlines
    inputs:
    x1, y1 - coordinates of first centerline
    x2, y2 - coordinates of second centerline
    years - time between the two centerlines, in years
    s - distance downsream, in m
    
    outputs:
    migr_rate - migration rate (in m/years)
    migr_sign - migration sign
    p - indices of correlation in second centerline
    q - indices of correlation in first centerline"""
    p,q,sm = correlate_clines(x1,x2,y1,y2)
    
    years = int(age2-age1)
    #to account for two nodes migrateing into one.
    qn = np.delete(np.array(q),np.where(np.diff(q)==0)[0]+1)
    pn = np.delete(np.array(p),np.where(np.diff(q)==0)[0]+1)
    
    xa = x1[qn][:-1]
    xb = x1[qn][1:]
    
    ya = y1[qn][:-1]
    yb = y1[qn][1:]
    
    x = x2[pn][1:]
    y = y2[pn][1:]
    
    migr_sign = np.sign((x-xa)*(yb-ya) - (y-ya)*(xb-xa))
    migr_rate = migr_sign*sm[pn,qn][1:]/years
    migr_rate = np.hstack((0,migr_rate))
    return migr_rate, pn, qn

def correlate_clines(x1,x2,y1,y2):
    """correlate points from new centerlines to their nearest neighbor in the previous centerline from a cost matrix or euclidian distances
    x1, y1 - coordinates of first centerline
    x2, y2 - coordinates of second centerline

    outputs:
    p - indices of correlation in second centerline
    q - indices of correlation in first centerline
    sm - distance matrix"""
    c = len(x1)
    r = len(x2)
    
    cl1 = np.asarray([(x1[i],y1[i]) for i in range(c)])
    
    cl2 = np.asarray([(x2[j],y2[j]) for j in range(r)])
    #sm = distance.cdist(cl2, cl1, metric =  'euclidean')*30
    sm = np.zeros((r,c))
    for i in range(0,r):
        sm[i,:] = ((x1-x2[i])**2 + (y1-y2[i])**2)**0.5
    ###PPPPP#### #array of row indices#
    p = np.arange(r)
    ###QQQQ#### #array of column indices#
    q = [np.where(abs(sm[i,:]) == np.min(abs(sm[i, :])))[0][0] for i in range(r)]

    return p,q,sm

def bendbybend(xs1, ys1, xs2, ys2, ages):
        #where curvature changes direction, =1
    curv
    nodes = np.array([(curv[i-2]*curv[i])< 0 and (curv[i+2]*curv[i])> 0 for i in range(2,len(curv)-2)])
    idx = np.where(nodes==1)[0]
    idx = [idx[i] for i in range(len(idx)-1) if idx[i+1]-idx[i] >10]
     
    
    MR = [np.mean(np.abs(R1[idx[i]:idx[i+1]])) for i in range(len(idx)-1)]
def segmented_MR(curv, R1, s, n=100):
    """
    approximate a bend-by-bend nth percentile lateral migration rate. 
    
    Inputs:
    curve: array of curvature for every node along the centerline.
    R1: array of already computed migration distances for every node along the centerline.
    ds: array of cumulative distance downstream between nodes.
    n: percentile
    
    Output:
    MR: array of nth percentile migration rate for each segment.  
    upstream: distance downstream of each segment start
    downstream: distance downstream of each segment end
    """
    R1 = np.array(R1)

    #where curvature changes direction, =1
    nodes = np.array([(curv[i-2]*curv[i])< 0 and (curv[i+2]*curv[i])> 0 for i in range(2,len(curv)-2)])
    idx = np.where(nodes==1)[0]
    idx = [idx[i] for i in range(len(idx)-1) if idx[i+1]-idx[i] >5]
    dists =  s[idx]
    
    MR = [np.nanpercentile(np.abs(R1[idx[i]:idx[i+1]]), n) for i in range(len(idx)-1)]

    return MR, dists
def moving_average(matrix, window):
    """
    averages migration rates on each bend over a set time window upwind. 
    
    inputs
    matrix: n years by m bend migration rate numpy array
    window: how many indices to average over, as years in the past
    
    output
    mid: moving average migration rate over window years for each bend
    """
 
    years = len(matrix[:,0])
    bends = np.min([np.count_nonzero(~np.isnan(matrix[i, :])) for i in range(years)])

    mid = np.zeros(shape = (years, bends))
    mid[0,:bends] = matrix[0,:bends]
    for year in range(1,years):
        if year < window:
            mid[year, :] = np.nanmean(matrix[:year, :bends], axis = 0)
        else:
            mid[year, :] = np.nanmean(matrix[year-window:(year+1), :bends], axis = 0)
     
    return mid