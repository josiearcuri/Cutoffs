import pandas as pd
import scipy.interpolate
from scipy.interpolate import InterpolatedUnivariateSpline
import numpy as np
import meanderpy as mp
import matplotlib.pyplot as plt
from astropy.stats import RipleysKEstimator
import os
from pointpats import ripley
import math
import random

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


def calibrate_migrate(self,nit,saved_ts,deltas,pad,crdist,Cf,kl,kv,dt,dens,t1,t2,t3,aggr_factor,cutoff_limit,*D):
    """function for computing migration rates along channel centerlines and moving the centerlines accordingly
    inputs:
    nit - number of iterations
    saved_ts - which time steps will be saved
    deltas - distance between nodes on centerline
    pad - padding (number of nodepoints along centerline)
    crdist - threshold distance at which cutoffs occur
    Cf - dimensionless Chezy friction factor
    kl - migration rate constant (m/s)
    kv - vertical slope-dependent erosion rate constant (m/s)
    dt - time step (s)
    dens - density of fluid (kg/m3)
    t1 - time step when incision starts
    t2 - time step when lateral migration starts
    t3 - time step when aggradation starts
    aggr_factor - aggradation factor
    D - channel depth (m)
    cutoff_limit - stop migrating once channel edperiences this umber of cutoffs"""
    channel = self.channels[-1] # first channel is the same as last channel of input
    x = channel.x; y = channel.y; z = channel.z
    W = channel.W;
        
    if len(D)==0: 
        D = channel.D
    else:
        D = D[0]
    k = 1.0 # constant in HK equation
    xc = [] # initialize cutoff coordinates
    # determine age of last channel:
    if len(self.cl_times)>0:
        last_cl_time = self.cl_times[-1]
    else:
        last_cl_time = 0
    dx, dy, dz, ds, s = compute_derivatives(x,y,z)
    slope = np.gradient(z)/ds
    # padding at the beginning can be shorter than padding at the downstream end:
    pad1 = int(pad/10.0)
    if pad1<5:
        pad1 = 5
    omega = -1.0 # constant in curvature calculation (Howard and Knutson, 1984)
    gamma = 2.5 # from Ikeda et al., 1981 and Howard and Knutson, 1984
    itn = -1
    while len(self.cutoff_times) < cutoff_limit:
        itn += 1
        update_progress(len(self.cutoff_times)/cutoff_limit)
        x, y = migrate_one_step(x,y,z,W,kl,dt,k,Cf,D,pad,pad1,omega,gamma)
        x,y,z,xc,yc,zc,cl_dist = cut_off_cutoffs(x,y,z,s,crdist,deltas)
        print(cl_dist)# find and execute cutoffs
        x,y,z,dx,dy,dz,ds,s = resample_centerline(x,y,z,deltas) # resample centerline
        slope = np.gradient(z)/ds
         #for itn<t1, z is unchanged
   
        if (itn>t1) & (itn<=t2): # incision
            if np.min(np.abs(slope))!=0:
                z = z + kv*dens*9.81*D*slope*dt 
            else:
                z = z - kv*dens*9.81*D*dt*0.00001
        if (itn>t2) & (itn<=t3): # lateral migration
            if np.min(np.abs(slope))!=0:
                z = z + kv*dens*9.81*D*slope*dt - kv*dens*9.81*D*np.median(slope)*dt
            else:
                z = z # no change in z
        if (itn>t3): # aggradation
            if np.min(np.abs(slope))!=0:
                z = z + kv*dens*9.81*D*slope*dt - aggr_factor*kv*dens*9.81*D*np.mean(slope)*dt 
            else:
                z = z + aggr_factor*dt
        if len(xc)>0: # save cutoff data
            self.cutoff_times.append(last_cl_time+(itn+1)*dt/(365*24*60*60.0))
            self.cutoff_dists.append(cl_dist)
            cutoff = Cutoff(xc,yc,zc,W,D) # create cutoff object
            self.cutoffs.append(cutoff)
            # saving centerlines:
        if np.mod(itn,saved_ts)==0:
            self.cl_times.append(last_cl_time+(itn)*dt/(365*24*60*60.0))
            channel = Channel(x,y,z,W,D) # create channel object
            self.channels.append(channel)
def cutoff_distributions(cutoffs, year, filepath):
    #pull cutoff data from channel belt object and export csv, return dataframe for plotting
    distances = [i.dist for i in cutoffs]
    times = [i.time for i in cutoffs]
    
    cuts = pd.DataFrame({'downstream_distance': distances, 'time': times})

    newcuts = cuts.to_csv(filepath+str(year)+"years_cutoff_distributions.csv", index_label = "Cutoff")
    return cuts
    
def plot_cutoff_distributions(cuts, year, filepath):
    #histograms of time and space
    #fig, axes = plt.subplots(1, 2)
    #for k in range(0, 2):
    #    axes[k].set_title(str(cuts.columns[k]))
    #    axes[k].hist(cuts[str(cuts.columns[k])]) 
    #axes[0].set_xlim(0,1)
    #axes[1].set_xlim(0,year)
    #axes[2].set_xlim(0,.5)
    #plt.savefig(filepath+str(year) + "yrs_cutoffdist.jpg")
    
    #colorscale = cuts['lengths']/max(cuts['lengths'])
    #scatter plot of cutoff time vs. rlative distance downstream
    f = plt.figure()
    
    sc = plt.scatter(cuts['time'], cuts['downstream_distance'], c = 'black', s = 1, edgecolor = 'black')
    
    plt.title("cutoff spatiotemporal distribution")
    plt.xlabel("time (years)")
    plt.ylabel("relative distance downstream (m)")
    
    plt.savefig(filepath+str(year) + "yrs_timevsspace.jpg")    
def Ripleys_K(cutoffs, filepath, year):
    
    y_max = 1
    y_min = 0
    x_max = year
    x_min = 0
    area = x_max*y_max
    maxr = round((x_max/2)**.5)
    
    z = cutoffs[['downstream_distance', 'time']].to_numpy()
    #z[:,0] = (100*z[:,0]).round()
    #z[:,1] = z[:,1].round()
    Kest = RipleysKEstimator(area=area, x_max=x_max, y_max=y_max, x_min=x_min, y_min=y_min)
    r = np.linspace(0, maxr, 100)
    fig = plt.figure()
    plt.plot(r, Kest.poisson(r), color='green', ls=':', label='poisson')
    plt.plot(r, Kest(data=z, radii=r, mode='none'), color='red', ls='--',
             label='none')
    #plt.plot(r, Kest(data=z, radii=r, mode='translation'), color='black',
     #        label='translation')
    plt.plot(r, Kest(data=z, radii=r, mode='ohser'), color='blue', ls='-.',
             label='ohser')
    plt.plot(r, Kest(data=z, radii=r, mode='var-width'), color='green',
             label=r'var_width')
    #plt.plot(r, Kest(data=z, radii=r, mode='ripley'), color='yellow',
     #        label=r'ripley')
    plt.title("Ripleys K estimator with various edge correction modes")
    plt.legend(loc = 'upper left')
    #plt.ylim((0,10))
    plt.xlabel("search radius, max = " + str(maxr))
    plt.ylabel("K")
    plt.savefig(filepath+str(year) + "yrs_RipleysK.jpg")   
      ##old
def k_envelope(cutoffs, filepath, year):
    k_test = ripley.k_test(cutoffs[['downstream_distance', 'time']].to_numpy(),distances = np.asarray(range(0,30,5)), keep_simulations=True)
    plt.plot(k_test.support, k_test.simulations.T, color='k', alpha=.01)
    plt.plot(k_test.support, k_test.statistic, color='orangered')

    plt.scatter(k_test.support, k_test.statistic, 
                cmap='viridis', c=k_test.pvalue < .05,
                zorder=4 # make sure they plot on top
               )

    plt.xlabel('Distance')
    plt.ylabel('K Function')
    plt.title('K Function Plot')
    plt.show()

     ##old
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
