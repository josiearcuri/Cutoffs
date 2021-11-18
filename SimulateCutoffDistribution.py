"""
This script uses the HKplus functions to load an existing centerline, run it for a certain number of years, then save the resulting cutoff distribution and bend-by-bend migration rate measurements.  

"""
import math
import matplotlib.pyplot as plt
import HKplus as hkp
import numpy as np
import pandas as pd
import os


#Set Variables for centerline and curvature calculation

D=3.4                       #constant width-average channel depth (m)
W = 100                     #constant channel width (m)
deltas = W//2;              #spacing of nodes along centerline
Cf = 0.005                  #dimensionless Chezy friction factor
kl = 25/(365*24*60*60.0)    #migration rate constant (m/s)
dt = .5*365*24*60*60.0     #time step (s)
pad= 100                    #number of nodes for periodic boundary
saved_ts = 100               #timeteps between saving centerlines
crdist = 4*W                #how close banks get before cutoff in m
cut_thresh = 500            #how many cutoffs to simulate, arbitrary if running for time

#Set Variables fror Cutoff nonlocal efects
tau = 1#e-folding timescale for nonlocal effects in years
decay_rate = dt/(tau*(365*24*60*60.0));   #this is the half-life on nonlocal effects, in units of seconds
bump_scale = 0             #this is the magntiude of nonlocal effects in relative difference 



#Set mode for titles
mode = "7"

#Set Result Directory
result_dir = "results/experiments/"+str(mode)+"/" 

#Load existing Centerline
filepath ="data/InitialChannel/InitialCL_"+str(mode)+".csv"

ch= hkp.load_initial_channel(filepath, W, D, deltas)

#Initalize Channel Belt for migration
chb = hkp.ChannelBelt(channels=[ch],cutoffs=[],cl_times=[0],cutoff_times=[], cutoff_dists = [], decay_rate = decay_rate, bump_scale = bump_scale, cut_thresh= cut_thresh, sinuosity = [])

#Plot initial centerline
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
plt.savefig(result_dir + mode+"_"+str(len(chb.cutoff_times))+"_timevsspace.png",bbox_inches='tight')
plt.show()


# Save Resulting Centerline
xes = chb.channels[-1].x
yes = chb.channels[-1].y
cl = pd.DataFrame({'x': xes, 'y': yes});

cl.to_csv(result_dir+"InitialChannel_result.csv", header = False, index = False)



