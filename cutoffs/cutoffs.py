import pandas as pd
from scipy.interpolate import InterpolatedUnivariateSpline
import numpy as np
import matplotlib.pyplot as plt
from astropy.stats import RipleysKEstimator
import os
import math
import meanderpyalt as mp



def generate_channel_from_file(filelist, slope = .01, D_in= 10, smooth_factor=.25, matlab_corr= -1):
    """function for creating a MeanderPy Channel object from an externally-sourced centerline in .csv file format.
        inputs:
        filelist - filelist must be a list of filepaths.  thie first should be a csv containing x and y values for each point on a centerline.  The second should be the widths of each point along the centerline
        slope - average channel slope (default = .01)
        D_in - channel depth (m, default = 10)
        smooth_factor - fraction of centerline points to sample from spline (default = 1/4)
        matlab_corr - 1 if y-axis need not be flipped, -1 for centerlines exported from matlab and need flipping
        
        outputs:
        ch - MeanderPy object of channel centerline
        x - uninterpolated array of x coordinates
        y - uninterpolated array of y coordinates
        z - array of z coordinates
        cl_len - length of  centerline(m) """
    #use pandas to load x,y coordinates and widths of centerline from csv
    varlist = [pd.read_csv(file, sep = ',', header=None).values for file in filelist]
    x = varlist[0][:,0]*30  ## x-dim array in Sm
    y = varlist[0][:,1]*matlab_corr*30 ##southern hemisphere y-dim array in m
    
    # shift coordinates so all are positive
    if sum(n < 0 for n in y) > 0:
        y = y-min(y)## y-dim earray in m 
    if sum(n < 0 for n in x) > 0:
        x = x-min(x)
        
    #average over widths to get a reach-constant width scalar
    W = np.mean(varlist[1][:,0])*30
  
    ## water depth scalar#
    D = D_in  
    # Linear length along the line, add a zero for first point:
    points = np.vstack([x, y]).T
    distance = np.cumsum( np.sqrt(np.sum( np.diff(points, axis=0)**2, axis=1 )) )
    distance = np.insert(distance, 0, 0)

    # Build a list of the spline function, one for each dimension:
    splines = [InterpolatedUnivariateSpline(distance, coords) for coords in points.T]

    # Compute the spline for the smoothed(sampled) distances:
    points_fitted = np.vstack([spl(np.linspace(0, distance[-1],round(len(x)*smooth_factor))) for spl in splines])
    
    ## z-dim array, interpolated with constant slope along points of centerline.  assumes centerline points are equidistantly placed along original centerline. 
    z = np.interp(np.asarray(range(len(points_fitted[0]))), [1, len(points_fitted[0])], [(slope*distance[-1]), 0]) 
    deltas = round(distance[-1]/(len(points_fitted[0])-1)) 
    return [mp.Channel(points_fitted[0],points_fitted[1],z,W,D), x, y, z, distance[-1], deltas]


def cutoff_distributions(cutoffs, year, filepath):
    #pull cutoff data from channel belt object and export csv, return dataframe for plotting
    distances = [i.dist for i in cutoffs]
    times = [i.time for i in cutoffs]
    
    cuts = pd.DataFrame({'downstream_distance': distances, 'time': times})

    newcuts = cuts.to_csv(filepath+str(year)+"years_cutoff_distributions.csv", index_label = "Cutoff")
    return cuts
    
def plot_cutoff_distributions(cuts, year, filepath):

    f = plt.figure()
    
    sc = plt.scatter(cuts['time'], cuts['downstream_distance'], c = 'black', s = 1, edgecolor = 'black')
    
    plt.title("cutoff spatiotemporal distribution")
    plt.xlabel("time (years)")
    plt.ylabel("relative distance downstream (m)")
    
    plt.savefig(filepath+str(year) + "yrs_timevsspace.jpg")    

def mc_envelope(cutoffs, year, resultdir, nit = 100, mode = ' modeled'): 
     
    #load estimator
    Kest = RipleysKEstimator(area=math.ceil(np.max(cutoffs.time)), x_max=math.ceil(np.max(cutoffs.time)), y_max=1, x_min=0, y_min=0)
    #load sample
    data = cutoffs[['downstream_distance', 'time']].to_numpy() 
   
    
    #generate random distibutions in same space + time ranges as data
    num_samples = len(cutoffs.time)
    r = np.linspace(0, math.ceil((np.max(cutoffs.time)/2)**.5))
    r =r[np.where(r<=50)]
    mc = np.zeros((len(r), nit))
    z = np.zeros((num_samples, 2))
    for i in range(nit):
        z[:,0] = np.random.random(size = num_samples)
        z[:,1] = np.random.random(size = num_samples)*math.ceil(np.max(cutoffs.time))
        mc[:,i] = Kest(data=z, radii=r, mode='ohser')

   
    ###check if normal distribution shows up as clustered 
    #data[:,0] = np.random.normal(.5, .01,num_samples)
    #data[:,1] = np.random.normal(.5, .01, num_samples)
    
    # compute bounds of CSR envelope, transform y values for plotting
    upper = np.subtract(np.sqrt(np.divide(np.ma.max(mc, axis = 1),math.pi)), r)
    lower = np.subtract(np.sqrt(np.divide(np.ma.min(mc, axis = 1),math.pi)), r)
    middle = np.subtract(np.sqrt(np.divide(np.ma.mean(mc, axis = 1),math.pi)), r)
    data = np.subtract(np.sqrt(np.divide(Kest(data=data, radii=r, mode='ohser'),math.pi)), r)
    
    #check for significant nonrandomness
    clustered = data[np.where(data>upper)]
    r_clus = r[np.where(data>upper)]
    regular = data[np.where(data<lower)]
    r_reg = r[np.where(data<lower)]
    
    fig = plt.figure()
    
    #plot CSR envelope
    plt.plot(r, upper, color='red', ls=':', label='_nolegend_')
    plt.plot(r, lower, color='red', ls=':', label='_nolegend_')
    plt.plot(r, middle, color='red', ls=':', label='CSR')
    
    #plot data
    plt.plot(r, data, color='black',label=str(num_samples)+mode +" events")

    #plot non random k values, flag if they exist
    cluster_flag = 0
    regular_flag = 0
    if len(clustered)>0:
        plt.scatter(r_clus, clustered, c='purple', s=30, marker = '*', label='clustered', alpha = .5)
        cluster_flag = 1
    if len(regular)>0:
        regular_flag = 1
        plt.scatter(r_reg, regular, c='purple', s=30, marker = '*', label='regularly spaced', alpha = .5)
    
    #plot specs
    plt.title("Monte Carlo CSR Envelope with ohser edge correction")
    plt.legend(loc = 'upper left')
    plt.xlabel("search radius (years)")
    plt.ylabel("sqrt(Ripley's K/pi) - r")
    plt.savefig(resultdir+str(year) + "yrs_clustertest.jpg")  
    return cluster_flag, regular_flag
