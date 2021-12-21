"""
This script uses the HKplus functions to generate a a bunch of centerlines, tracking their sinuosity as they develop with different nonlocal effects, and plot change in sinuosity through time. 

@JA 2021
"""

import math
import matplotlib.pyplot as plt
import HKplus as hkp
import numpy as np
import pandas as pd


#output file location
filepath = "~/Desktop/Cutoffs/results/sinuositytest/"
 #Set Variables for centerline and curvature calculation

#values to test
Magnitude_range = np.linspace(1,3,2)
Tau_range = np.linspace(2,10,3)
 
#array of paramaters with a non-nonlocal effect test included
params = np.asarray([(int(M),int(tau))for M in Magnitude_range for tau in Tau_range])
params = np.concatenate(([[0,0]],params), axis = 0)
print(params)
for k in [len(params)-1]:
    #create centerline

    #Choose name for centerline
    name= "sin_test_2_"+str(k)

    #Initiate Channel Object
    ch = hkp.generate_initial_channel(W=100,D=3.4,deltas=50,pad=100,Length=50000)

    #Set Variables fror Cutoff nonlocal efects
    tau = params[k][1]#e-folding timescale for nonlocal effects in years
    #this is the half-life on nonlocal effects, in units of seconds
    if tau ==0:
        decay_rate = 0
    else:
        decay_rate = .25/(tau)
    bump_scale = params[k][0]             #this is the magntiude of nonlocal effects in relative difference 
    print(decay_rate)
    print(bump_scale)
    #Initiate Channel Belt for migrating channel object
    chb = hkp.ChannelBelt(channels=[ch],cutoffs=[],cl_times=[0.0],cutoff_times=[], cutoff_dists = [], decay_rate = decay_rate, bump_scale = bump_scale, cut_thresh = 300, sinuosity=[1])
 
    #migrate and track sinuosity
    chb.migrate_cuts(saved_ts=20,deltas=50,pad=100,crdist=200,Cf=.005,kl=25/(365*24*60*60.0)  ,dt=.25*(365*24*60*60.0)) 
    chb.plot_channels()
    plt.show()

    #Save Sinuosity time series
    times = chb.cl_times
    sins = chb.sinuosity
    #Save sinuosity throught time in csv
    #uncomment to save sinuosity series
    sinseries = pd.DataFrame({'time':times, 'sinuosity': sins})
    sinseries.to_csv(filepath+name+".csv", header = True, index = False)
    print(str(k+1)+" run(s) complete")