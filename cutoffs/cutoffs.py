import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from SpaceTime import RipleysKEstimator_spacetime
import os
import math

def cutoff_distributions(cutoffs, year, filepath):
    """pull cutoff data from channel belt object and export csv, return dataframe for plotting
    """
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

def mc_envelope(cutoffs, year, spacing,resultdir, nit = 99, d_max = 1, mode = ' modeled'): 
    #load estimator
    Kest = RipleysKEstimator_spacetime(t_max=year, d_max=d_max, t_min=0, d_min=0)
    #load sample
    data = cutoffs[['downstream_distance', 'time']].to_numpy() 
    data[:,0] = data[:,0]*d_max
    #generate random distibutions in same space + time ranges as data
    num_samples = len(cutoffs.time)
    r_time = np.linspace(1, 50, spacing)
    r_space = np.linspace(.000001,50, spacing)
    mc_dt = np.zeros((len(r_space),len(r_time), nit))
    mc_d = np.zeros((len(r_space), nit))
    mc_t = np.zeros((len(r_time), nit))
    z = np.zeros((num_samples, 2))
    for i in range(nit):
        z[:,0] = np.random.random(size = num_samples)*d_max
        z[:,1] = np.random.random(size = num_samples)*math.ceil(np.max(cutoffs.time))
        k_dt, k_d, k_t = Kest(data=z, dist_time=r_time, dist_space=r_space) 
        mc_dt[:,:,i] =  k_dt
        mc_d[:,i] = k_d
        mc_t[:,i] =k_t
    print("Monte Carlo Simulation Complete")

    #data2 = np.zeros((num_samples,2)) 
    ###check if normal distribution shows up as clustered 
    #data2[:,0] = np.random.random(size = num_samples)
    #data2[:,1] = np.random.normal(size = num_samples)*year
    #plt.plot(data2[:,0], data2[:,1])
    #plt.show()
    # compute bounds of CSR envelope, transform y values for plotting
    upper_dt = np.ma.max(mc_dt, axis = 2)
    upper_d = np.ma.max(mc_d, axis = 1)
    upper_t = np.ma.max(mc_t, axis = 1)
    lower_dt = np.ma.min(mc_dt, axis = 2)
    lower_d = np.ma.min(mc_d, axis = 1)
    lower_t = np.ma.min(mc_t, axis = 1)
    middle_dt = np.ma.mean(mc_dt, axis = 2)
    middle_d = np.ma.mean(mc_d, axis = 1)
    middle_t = np.ma.mean(mc_t, axis = 1)
    
    
    K_dt, K_d, K_t = Kest(data=data, dist_time=r_time, dist_space=r_space)
   
    #check for significant nonrandomness
   # clustered = K_dt[np.where(K_dt>upper)]
    #r_clus_ = r[np.where(K_dt>upper)]
    #regular = K_dt[np.where(K_dt<lower)]
    #r_reg = r[np.where(K_dt<lower)]
    fig = plt.figure()
     #plot CSR envelope
    plt.plot(r_space, upper_d, color='red', ls=':', label='_nolegend_')
    plt.plot(r_space, lower_d, color='red', ls=':', label='_nolegend_')
    plt.plot(r_space, middle_d, color='red', ls=':', label='CSR')
    plt.plot(r_space, K_d, color = "black", label = 'spatial K')
    plt.legend(loc = 'upper left')
    plt.xlabel("d")
    plt.ylabel("Ripley's K/2 - d")
    plt.title("Homegrown 1D space Ripley's K")
    plt.show()
    fig2 = plt.figure()
     #plot CSR envelope
    plt.plot(r_time, upper_t, color='red', ls=':', label='_nolegend_')
    plt.plot(r_time, lower_t, color='red', ls=':', label='_nolegend_')
    plt.plot(r_time, middle_t, color='red', ls=':', label='CSR')
    plt.plot(r_time, K_t, color = "black", label = 'temporal K')
    plt.legend(loc = 'upper left')
    plt.xlabel("t")
    plt.ylabel("Ripley's K/2 - t")
    plt.title("Homegrown 1D time Ripley's K")
    plt.show()
    fig3 = plt.figure()
    ax = fig3.gca(projection='3d')
    #plot CSR envelope
    #plt.plot(r, upper, color='red', ls=':', label='_nolegend_')
    #plt.plot(r, lower, color='red', ls=':', label='_nolegend_')
    #plt.plot(r, middle, color='red', ls=':', label='CSR')
    
    #plot data
  
    X, Y = np.meshgrid(r_space, r_time)
    ax.plot_surface(X,Y,K_dt,cmap=cm.coolwarm)
    ax.plot_surface(X,Y,middle_dt, cmap=cm.coolwarm)
    #ax.plot_surface(X,Y,lower_dt, cmap=cm.coolwarm)
    #ax.plot_surface(X,Y,upper_dt, cmap=cm.coolwarm)
    #im = ax.matshow(K_dt)
    #fig.colorbar(im)
    #plot non random k values, flag if they exist
    cluster_flag = 0
    regular_flag = 0
    #if len(clustered)>0:
    #    plt.scatter(r_clus, clustered, c='purple', s=30, marker = '*', label='clustered', alpha = .5)
    #    cluster_flag = 1
    #if len(regular)>0:
    #    regular_flag = 1
    #    plt.scatter(r_reg, regular, c='green', s=30, marker = '*', label='regularly spaced', alpha = .5)
    
    #plot specs
    plt.title("Homegrown space time Ripley's K")
    #plt.legend(loc = 'upper left')
    #plt.xlabel("search radius (years)")
    #plt.ylabel("sqrt(Ripley's K/pi) - r")
    plt.show()
    return cluster_flag, regular_flag
