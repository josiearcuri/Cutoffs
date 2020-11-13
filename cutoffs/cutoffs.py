import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from SpaceTime import RipleysKEstimator_spacetime
import os
import math
import imageio
import glob

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
    
def plot_cutoff_distributions(cuts, year, filepath, mode):
    """save scatter plot of cutoff locations downstream and occurrence in model time
    """
    f = plt.figure()
    
    sc = plt.scatter(cuts['time'], cuts['downstream_distance'], c = 'black', s = 1, edgecolor = 'black')
    ncuts = len(cuts['time'])
    plt.title(str(ncuts)+" cutoffs with "+ mode)
    plt.xlabel("time (years)")
    plt.ylabel("relative distance downstream")
    
    plt.savefig(filepath+str(year) + mode+"_yrs_timevsspace.jpg", dpi = 500)    

def mc_envelope(cutoffs, year,resultdir, nit = 99, d_max = 1000, mode = ' modeled'): 
    """pull cutoff data from model output, test distribution for nonrandomness in space and time using 1D Ripley's K test with monte carlo simulations to test statistical significance against complete spatial randomness.  
    """
    #load estimator
    Kest = RipleysKEstimator_spacetime(t_max=year, d_max=d_max, t_min=0, d_min=0)
    #load sample
    data = cutoffs[['downstream_distance', 'time']].to_numpy() 
    data[:,0] = data[:,0]*d_max
    #generate random distibutions in same space + time ranges as data
    num_samples = len(cutoffs.time)
    r_time = np.linspace(5, 100, 20)
    r_space = np.linspace(10,d_max/5, 20)
    
    mc_d = np.zeros((len(r_space), nit))
    mc_t = np.zeros((len(r_time), nit))
    z = np.zeros((num_samples, 2))
    for i in range(nit):
        z[:,0] = np.random.random(size = num_samples)*d_max
        z[:,1] = np.random.random(size = num_samples)*year
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
    plt.legend(loc = 'lower left')
    plt.xlabel("d in relative distance along centerline")
    plt.ylabel("Ripley's K - 2d")
    plt.title("Homegrown 1D space Ripley's K with " + mode)
    plt.savefig(resultdir + str(year)+"yrs_Space_Ripley_"+mode+".jpg", dpi = 500)
    plt.show()
    fig2 = plt.figure()
     #plot CSR envelope
    plt.plot(r_time, upper_t, color='red', ls=':', label='_nolegend_')
    plt.plot(r_time, lower_t, color='red', ls=':', label='_nolegend_')
    plt.plot(r_time, middle_t, color='red', ls=':', label='CSR')
    plt.plot(r_time, K_t, color = "black", label =str(num_samples)+ ' cutoffs')
    plt.legend(loc = 'lower left')
    plt.xlabel("t in years")
    plt.ylabel("Ripley's K - 2t")
    plt.title("Homegrown 1D time Ripley's K with " + mode)
    plt.show()
    plt.savefig(resultdir + str(year)+"yrs_Time_Ripley_"+mode+".jpg", dpi = 500)

    return
def save_animations(source, name,sp,  filt):
    images = []
    original_files=list(glob.glob(source))
    original_files.sort(reverse=False)
    
    for file_ in original_files[11:121:filt]:
        images.append(imageio.imread(file_))
    imageio.mimsave(name, images, duration=1/sp, subrectangles=True)
    print(name)
def mc_envelope_comp(cutoffs_curv,cutoffs_ne, year,resultdir, nit = 99, d_max = 100): 
    """pull cutoff data from model output, test distribution for nonrandomness in space and time using 1D Ripley's K test with monte carlo simulations to test statistical significance against complete spatial randomness.  
    """
    #load estimator
    Kest = RipleysKEstimator_spacetime(t_max=year, d_max=d_max, t_min=0, d_min=0)
    #load sample
    data_curv = cutoffs_curv[['downstream_distance', 'time']].to_numpy() 
    data_ne = cutoffs_ne[['downstream_distance', 'time']].to_numpy() 
    data_curv[:,0] = data_curv[:,0]*d_max
    data_ne[:,0] = data_ne[:,0]*d_max
    #generate random distibutions in same space + time ranges as data
    num_samples = max(len(cutoffs_curv.time), len(cutoffs_ne.time))
    r_time = np.linspace(1, 101, 10)
    r_space = np.linspace(1,11, 20)
    mc_d = np.zeros((len(r_space), nit))
    mc_t = np.zeros((len(r_time), nit))
    z = np.zeros((num_samples, 2))
    for i in range(nit):
        z[:,0] = np.random.random(size = num_samples)*d_max
        z[:,1] = np.random.random(size = num_samples)*year
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
    
    
    K_d_curv, K_t_curv = Kest(data=data_curv, dist_time=r_time, dist_space=r_space)
    K_d_ne, K_t_ne = Kest(data=data_ne, dist_time=r_time, dist_space=r_space)

    fig = plt.figure()
     #plot CSR envelope
    plt.plot(r_space/d_max, upper_d, color='red', ls=':', label='_nolegend_')
    plt.plot(r_space/d_max, lower_d, color='red', ls=':', label='_nolegend_')
    plt.plot(r_space/d_max, middle_d, color='red', ls=':', label='CSR')
    plt.plot(r_space/d_max, K_d_curv, color = "black", label = str(len(cutoffs_curv.time))+ ' cutoffs from only curvature')
    plt.plot(r_space/d_max, K_d_ne, color = "purple", label =str(len(cutoffs_ne.time))+ ' cutoffs with from nonlocal effects')
    plt.legend(loc = 'lower left')
    plt.xlabel("d in relative distance along centerline")
    plt.ylabel("Ripley's K - 2d")
    plt.title("Homegrown 1D space Ripley's K")
    plt.savefig(resultdir + str(year)+"yrs_Space_Ripley_comp+.jpg", dpi = 500)
    plt.close()
    fig2 = plt.figure()
     #plot CSR envelope
    plt.plot(r_time, upper_t, color='red', ls=':', label='_nolegend_')
    plt.plot(r_time, lower_t, color='red', ls=':', label='_nolegend_')
    plt.plot(r_time, middle_t, color='red', ls=':', label='CSR')
    plt.plot(r_time, K_t_curv, color = "black", label =str(len(cutoffs_curv.time))+ ' cutoffs from only curvature')
    plt.plot(r_time, K_t_ne, color = "purple", label =str(len(cutoffs_ne.time))+ ' cutoffs from nonlocal effects')
    plt.legend(loc = 'lower left')
    plt.xlabel("t in years")
    plt.ylabel("Ripley's K - 2t")
    plt.title("Homegrown 1D time Ripley's K")
    plt.savefig(resultdir + str(year)+"yrs_Time_Ripley_comp.jpg", dpi = 500)
    plt.close()
    return