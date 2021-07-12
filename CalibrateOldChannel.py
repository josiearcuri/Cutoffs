"""
This script uses the HKplus functions to load an existing centerline, run it for a certain number of years, then save the resulting centerline for later use.  

"""
import math
import matplotlib.pyplot as plt
import HKplus as hkp
import numpy as np
import pandas as pd

#Set Variables for centerline and curvature calculation

Cf = 0.005              # dimensionless Chezy friction factor
kl =45/(365*24*60*60.0) # migration rate constant (m/s)
dt = .25*365*24*60*60.0 
nit = 20# time step (s)
pad= 200                     # dont change
saved_ts = 1# int(365*24*60*60.0/dt)  
W = 100# which time steps centerline will be saved at
deltas = W//2
crdist = 2*W
D = 3.4
dur = 4
#Set Variables fro nonlocal efects
decay_rate = dt/((dur/2)*(365*24*60*60.0));   #ranges between 1/3 to 1/10, to be developed 
bump_scale = .5              #to multiple kl by,amplitude of ne bump, range between 1 and 4, set to 0 for no nonlocal effects
cut_thresh = 5            #how many cutoffs to simulate, arbitrary if running for time

#Set Result Directory
result_dir = "sample_results/case4/30/InitialChannel_31.csv" ##change this to wherever you want to save your results

#Load Existing Centerline

filepath = "sample_data/InitialChannel/InitialCL_starter.csv"

ch = hkp.load_initial_channel(filepath, W, D, deltas)

#Initialize Channel Belt for Migration
chb = hkp.ChannelBelt(channels=[ch],cutoffs=[],cl_times=[0.0],cutoff_times=[], cutoff_dists = [], decay_rate = decay_rate, bump_scale = bump_scale, cut_thresh = cut_thresh, sinuosity = [3.5])

#Plot initial channel
chb.plot_channels()
plt.title("Initial Centerline")
plt.show()
 
#Migrate
chb.migrate_cuts(saved_ts,deltas,pad,crdist,Cf,kl,dt) 

#Plot resulting centerline
chb.plot_channels2(2,"/~/")
#plt.savefig("figure2_check.png", dpi = 800)
plt.show()


#Save Resulting Centerline
xes = chb.channels[-1].x
yes = chb.channels[-1].y
cl = pd.DataFrame({'x': xes, 'y': yes});

#cl.to_csv(result_dir, header = False, index = False)


