"""
This script uses the HKplus functions to load an existing centerline, run it for a certain number of years, then save the resulting centerline for later use.  

"""
import math
import matplotlib.pyplot as plt
import HKplus as hkp
import numpy as np
import pandas as pd

#Set Variables for centerline and curvature calculation
D=10
W = 100                   #constant width (m)
deltas = W//2;            #spacing of nodes along centerline (m)
nit = 200              # number of iterations
Cf = 0.022              # dimensionless Chezy friction factor
kl = 5/(365*24*60*60.0) # migration rate constant (m/s)
dt = 1*365*24*60*60.0     # time step (s)
pad= 100                     # dont change
saved_ts = 20               # which time steps centerline will be saved at
crdist = 2*W                    # how close  banks get before cutoff in m

#Set Variables fro nonlocal efects
decay_rate = dt/(15*(365*24*60*60.0));   #ranges between 1/3 to 1/10, to be developed
bump_scale = 0              #to multiple kl by,amplitude of ne bump, range between 1 and 4, set to 0 for no nonlocal effects
cut_thresh = 100            #how many cutoffs to simulate, arbitrary if running for time

#Set Result Directory
result_dir = "sample_data/InitialChannel/" ##change this to wherever you want to save your results

#Load Existing Centerline
filepath ="sample_data/InitialChannel/InitialCL_init.csv"
ch= hkp.load_initial_channel(filepath, W, D, deltas)

#Initialize Channel Belt for Migration
chb = hkp.ChannelBelt(channels=[ch],cutoffs=[],cl_times=[0.0],cutoff_times=[], cutoff_dists = [], decay_rate = decay_rate, bump_scale = bump_scale, cut_thresh = cut_thresh, sinuosity = [])

#Plot initial channel
chb.plot_channels()
plt.title("Initial Centerline")
plt.show()
 
#Migrate
chb.migrate_years(nit,saved_ts,deltas,pad,crdist,Cf,kl,dt) 

#Plot resulting centerline
chb.plot_channels()
plt.title(str(int(nit*dt/(365*24*60*60.0)))+ " years at "+ str(kl*(365*24*60*60.0))+ "m/yr")
plt.show()

#Save Resulting Centerline
xes = chb.channels[-1].x
yes = chb.channels[-1].y
cl = pd.DataFrame({'x': xes, 'y': yes});
cl.to_csv("sample_data/InitialChannel/InitialCL_init.csv", header = False, index = False)