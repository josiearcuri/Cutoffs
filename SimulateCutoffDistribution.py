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
D=3.4
W = 100                   #constant width (m)
deltas = W//2;            #spacing of nodes along centerline (m)
nit = 2*(21**2)             # number of iterations
Cf = .005            # dimensionless Chezy friction factor
kl = 45/(365*24*60*60.0) # migration rate constant (m/s) 
dt = .1*365*24*60*60.0     # time step (s)
pad= 200                     # dont change
saved_ts = int(365*24*60*60.0/dt)               # which time steps centerline will be saved at
crdist = 2*W                    # how close  banks get before cutoff in m
#Set Variables for nonlocal efects
decay_rate = dt/(5*(365*24*60*60.0));   #ranges between 1/3 to 1/10, to be developed
bump_scale = 2              #to multiple kl by,amplitude of ne bump, range between 1 and 4, set to 1 for no nonlocal effects
cut_thresh = 300           #how many cutoffs to simulate

#Set mode for titles
mode = "30"

#Set Result Directory
result_dir = "sample_results/case4/"+mode+"/" 

#Load existing Centerline
filepath =result_dir+"InitialChannel_"+mode+".csv"
#filepath = "sample_results/case2/24/" +"InitialCL_result.csv"

ch= hkp.load_initial_channel(filepath, W, D, deltas)

#Ititalize Channel Belt for migration
chb = hkp.ChannelBelt(channels=[ch],cutoffs=[],cl_times=[0],cutoff_times=[], cutoff_dists = [], decay_rate = decay_rate, bump_scale = bump_scale, cut_thresh= cut_thresh, sinuosity = [3])


chb.plot_channels()

plt.show()

#Migrate
chb.migrate_bendtracking(saved_ts,deltas,pad,crdist,Cf,kl,dt) 

#Plot resulting Centerline
chb.plot_channels()

plt.show()

# Save Cutoff Distributions for Clustering Tests #
chb.cutoff_distributions(int(chb.cutoff_times[-1]), result_dir, mode)
plt.title(str(len(chb.cutoff_times))+" cutoffs")
plt.savefig(result_dir + mode+"_"+str(len(chb.cutoff_times))+"_cutoffs_timevsspace.png",bbox_inches='tight')
plt.close()


hkp.plot_sinuosity(chb.cl_times, chb.sinuosity)
plt.show()

#Save Resulting Centerline
xes = chb.channels[-1].x
yes = chb.channels[-1].y
cl = pd.DataFrame({'x': xes, 'y': yes});

cl.to_csv(filepath, header = False, index = False)

chb.MR_time(result_dir+"bendbybendmr.csv")
plt.title(str(bump_scale)+ " NE magnitude with "+ str(round(kl*(365*24*60*60.0)/100, 3))+ "ch-w/yr MR constant")
plt.savefig(result_dir+"bendbybend.png", dpi = 500)
plt.show()


