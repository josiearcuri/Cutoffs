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
D = 3.27;   
W = 100                     #constant width
deltas = W//2;            #spacing of nodes along centerline
Cf = 0.005              # dimensionless Chezy friction factor
kl = 5/(365*24*60*60.0) # migration rate constant (m/s)
dt = .5*365*24*60*60.0     # time step (s)
pad= 100                     # dont change
saved_ts = 200               # which time steps centerline will be saved at
crdist = 2*W                    # how close  banks get before cutoff in m
nit = 1000
#Set Variables for nonlocal efects
decay_rate = dt/(10*(365*24*60*60.0));   #ranges between 1/3 to 1/10, to be developed
bump_scale = 4              #to multiple kl by,amplitude of ne bump, range between 1 and 4, set to 0 for no nonlocal effects
cut_thresh = 10            #how many cutoffs to simulate

#Set mode for titles
mode = "example"

#Set Result Directory
result_dir = "sample_results/" 

#Load existing Centerline
filepath ="sample_data/InitialChannel/InitialCL_init.csv"
ch= hkp.load_initial_channel(filepath, W, D, deltas)

#Ititalize Channel Belt for migration
chb = hkp.ChannelBelt(channels=[ch],cutoffs=[],cl_times=[0.0],cutoff_times=[], cutoff_dists = [], decay_rate = decay_rate, bump_scale = bump_scale, cut_thresh= cut_thresh, sinuosity = [])

chb.plot_channels()

plt.show()

#Migrate
chb.migrate_cuts(saved_ts,deltas,pad,crdist,Cf,kl,dt) 

#Plot resulting Centerline
chb.plot_channels()

plt.show()

# Save Cutoff Distributions for Clustering Tests #
chb.cutoff_distributions(int(chb.cutoff_times[-1]), result_dir, mode)
plt.title(str(len(chb.cutoff_times))+" cutoffs")
plt.savefig(result_dir + mode+str(len(chb.cutoff_times))+"_cutoffs_timevsspace.png",bbox_inches='tight')
plt.close()

#Save Resulting Centerline
xes = chb.channels[-1].x
yes = chb.channels[-1].y
cl = pd.DataFrame({'x': xes, 'y': yes});
cl.to_csv(result_dir+"InitialCL_init_result.csv", header = False, index = False)
