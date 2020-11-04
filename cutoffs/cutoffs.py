<<<<<<< HEAD
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from astropy.stats import RipleysKEstimator
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

def mc_envelope(cutoffs, year, resultdir, nit = 99, mode = ' modeled'): 
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
        plt.scatter(r_clus, clustered, c='purple', s=40, marker = '*', label='clustered', alpha = .5)
        cluster_flag = 1
    if len(regular)>0:
        regular_flag = 1
        plt.scatter(r_reg, regular, c='green', s=40, marker = '*', label='regularly spaced', alpha = .5)
    
    #plot specs
    plt.title("Monte Carlo CSR Envelope with ohser edge correction")
    plt.legend(loc = 'upper left')
    plt.xlabel("search radius (years)")
    plt.ylabel("sqrt(Ripley's K/pi) - r")
    plt.savefig(resultdir+str(year) + "yrs_clustertest.jpg")  
    return cluster_flag, regular_flag
=======
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
    x_location = [i.x for i in cutoffs]
    y_location = [i.y for i in cutoffs]
    cuts = pd.DataFrame({'downstream_distance': distances, 'time': times, 'x': x_location, 'y': y_location})
    newcuts = cuts.to_csv(filepath+str(year)+"years_cutoff_distributions.csv", index_label = "Cutoff")
    return cuts
    
def plot_cutoff_distributions(cuts, year, filepath):
    """save scatter plot of cutoff locations downstream and occurrence in model time
    """
    f = plt.figure()
    
    sc = plt.scatter(cuts['time'], cuts['downstream_distance'], c = 'black', s = 1, edgecolor = 'black')
    
    plt.title("cutoff spatiotemporal distribution")
    plt.xlabel("time (years)")
    plt.ylabel("relative distance downstream (m)")
    
    plt.savefig(filepath+str(year) + "yrs_timevsspace.jpg", dpi = 500)    

def mc_envelope(cutoffs, year, spacing,resultdir, nit = 99, d_max = 1000, mode = ' modeled'): 
    """pull cutoff data from model output, test distribution for nonrandomness in space and time using 1D Ripley's K test with monte carlo simulations to test statistical significance against complete spatial randomness.  
    """
    #load estimator
    Kest = RipleysKEstimator_spacetime(t_max=year, d_max=d_max, t_min=0, d_min=0)
    #load sample
    data = cutoffs[['downstream_distance', 'time']].to_numpy() 
    data[:,0] = data[:,0]*d_max
    #generate random distibutions in same space + time ranges as data
    num_samples = len(cutoffs.time)
    r_time = np.linspace(1, 50, spacing)
    r_space = np.linspace(1,.01*d_max, spacing)
    
    mc_d = np.zeros((len(r_space), nit))
    mc_t = np.zeros((len(r_time), nit))
    z = np.zeros((num_samples, 2))
    for i in range(nit):
        z[:,0] = np.random.random(size = num_samples)*d_max
        z[:,1] = np.random.random(size = num_samples)*math.ceil(np.max(cutoffs.time))
        k_d, k_t = Kest(data=z, dist_time=r_time, dist_space=r_space) 

        mc_d[:,i] = k_d
        mc_t[:,i] =k_t
    

    #data2 = np.zeros((num_samples,2)) 
    ###check if random distribution shows up as clustered 
    #data2[:,0] = np.random.random(size = num_samples)
    #data2[:,1] = np.random.normal(size = num_samples)*year
    #plt.plot(data2[:,0], data2[:,1])
    #plt.show()
    # compute bounds of CSR envelope, transform y values for plotting

    upper_d = np.ma.max(mc_d, axis = 1)
    upper_t = np.ma.max(mc_t, axis = 1)
 
    lower_d = np.ma.min(mc_d, axis = 1)
    lower_t = np.ma.min(mc_t, axis = 1)
   
    middle_d = np.ma.mean(mc_d, axis = 1)
    middle_t = np.ma.mean(mc_t, axis = 1)
    
    
    K_d, K_t = Kest(data=data, dist_time=r_time, dist_space=r_space)

    fig = plt.figure()
     #plot CSR envelope
    plt.plot(r_space/d_max, upper_d, color='red', ls=':', label='_nolegend_')
    plt.plot(r_space/d_max, lower_d, color='red', ls=':', label='_nolegend_')
    plt.plot(r_space/d_max, middle_d, color='red', ls=':', label='CSR')
    plt.plot(r_space/d_max, K_d, color = "black", label = str(num_samples)+ ' cutoffs')
    plt.legend(loc = 'upper left')
    plt.xlabel("d in relative distance along centerline")
    plt.ylabel("Ripley's K - 2d")
    plt.title("Homegrown 1D space Ripley's K with " + mode)
    plt.savefig(resultdir + str(year)+"yrs_Space_Ripley_"+mode+".jpg", dpi = 500)
    fig2 = plt.figure()
     #plot CSR envelope
    plt.plot(r_time, upper_t, color='red', ls=':', label='_nolegend_')
    plt.plot(r_time, lower_t, color='red', ls=':', label='_nolegend_')
    plt.plot(r_time, middle_t, color='red', ls=':', label='CSR')
    plt.plot(r_time, K_t, color = "black", label =str(num_samples)+ ' cutoffs')
    plt.legend(loc = 'upper left')
    plt.xlabel("t in years")
    plt.ylabel("Ripley's K - 2t")
    plt.title("Homegrown 1D time Ripley's K with " + mode)
    plt.savefig(resultdir + str(year)+"yrs_Time_Ripley_"+mode+".jpg", dpi = 500)

    return 
>>>>>>> master
