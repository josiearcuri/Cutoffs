"""
This script uses the HKplus functions to load an existing centerline, run it for a certain number of years, then analyize how migration rate changes through time on a bend-by-bend resultion  

"""
import math
import matplotlib.pyplot as plt
import HKplus as hkp
import numpy as np
import pandas as pd

#Set Variables for centerline and curvature calculation

D=3.4
W = 100                     #constant width
deltas = W//2;            #spacing of nodes along centerline
nit = 1000              # number of iterations
Cf = 0.005              # dimensionless Chezy friction factor
kl =22.7/(365*24*60*60.0) # migration rate constant (m/s)
dt = .5*365*24*60*60.0     # time step (s)
pad= 200                     # dont change
saved_ts = 2              # which time steps centerline will be saved at
crdist = 4*W                    # how close  banks get before cutoff in m


#Set Variables fro nonlocal efects
decay_rate = dt/(5*(365*24*60*60.0));   #ranges between 1/3 to 1/10, to be developed 
bump_scale = 0              #to multiple kl by,amplitude of ne bump, range between 1 and 4, set to 0 for no nonlocal effects
cut_thresh = 10            #how many cutoffs to simulate, arbitrary if running for time

#Set Result Directory
result_dir = "sample_results/case2/7/" ##change this to wherever you want to save your results

#Load Existing Centerline

filepath =result_dir + "InitialCL_result.csv"
ch= hkp.load_initial_channel(filepath, W, D, deltas)

#Initialize Channel Belt for Migration
chb = hkp.ChannelBelt(channels=[ch],cutoffs=[],cl_times=[0.0],cutoff_times=[], cutoff_dists = [], decay_rate = decay_rate, bump_scale = bump_scale, cut_thresh = cut_thresh, sinuosity = [3])

#Plot initial channel
chb.plot_channels()
plt.title("Initial Centerline")
plt.show()

#Migrate
chb.migrate_bendtracking(saved_ts,deltas,pad,crdist,Cf,kl,dt) 

#Plot resulting centerline
chb.plot_channels()
plt.title(str(int(nit*dt/(365*24*60*60.0)))+ " years at "+ str(kl*(365*24*60*60.0))+ "m/yr")
plt.savefig(result_dir+"channelfrombendbybend.png", dpi = 500)
plt.close()

chb.MR_time(result_dir+"bendbybendmr.csv")
plt.show()
