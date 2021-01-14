"""
This script uses the HKplus functions to load an existing centerline, run it for a certain number of years, then save the resulting cutoff distribution.  

"""
import math
import matplotlib.pyplot as plt
import HKplus as hkp
import numpy as np
import pandas as pd
import os


#Set Variables for centerline and curvature calculation
D = 10;   
W = 150                     #constant width
deltas = 50;            #spacing of nodes along centerline
Cf = 0.022              # dimensionless Chezy friction factor
kl = 10/(365*24*60*60.0) # migration rate constant (m/s)
dt = .5*365*24*60*60.0     # time step (s)
pad= 100                     # dont change
saved_ts = 400               # which time steps centerline will be saved at
crdist = W                    # how close  banks get before cutoff in m

#Set Variables for nonlocal efects
decay_rate = dt/(10*(365*24*60*60.0));   #ranges between 1/3 to 1/10, to be developed
bump_scale = 4              #to multiple kl by,amplitude of ne bump, range between 1 and 4, set to 0 for no nonlocal effects
cut_thresh = 500            #how many cutoffs to simulate

#Set mode for titles
mode = "NonlocalEffects"

#Set Result Directory
result_dir = "sample_results/10mpyr_small/" 

#Load existing Centerline
filepath ="sample_data/InitialChannel/InitialCL_10mpyr_small.csv"
ch= hkp.load_initial_channel(filepath, W, D, deltas)

#Ititalize Channel Belt for migration
chb = hkp.ChannelBelt(channels=[ch],cutoffs=[],cl_times=[0.0],cutoff_times=[], cutoff_dists = [], decay_rate = decay_rate, bump_scale = bump_scale, cut_thresh= cut_thresh)

#Migrate
chb.migrate_cuts(saved_ts,deltas,pad,crdist,Cf,kl,dt) 

#Plot resulting Centerline
chb.plot_channels()
plt.title(str(int(chb.cutoff_times[-1]))+ " years at "+ str(kl*(365*24*60*60.0))+ "m/yr")
plt.show()

# Save Cutoff Distributions for Clustering Tests #
chb.cutoff_distributions(int(chb.cutoff_times[-1]), result_dir, mode)

