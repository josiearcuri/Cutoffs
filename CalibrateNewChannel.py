
"""
This script uses the HKplus functions to generate a channel centerline, run it for a certain number of years, then save the resulting centerline for later use.  

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
nit = 2000              # number of iterations
Cf = 0.005              # dimensionless Chezy friction factor
kl = 20/(365*24*60*60.0) # migration rate constant (m/s)
dt = 2*365*24*60*60.0     # time step (s)
pad= 200                     # dont change
saved_ts = 200              # which time steps centerline will be saved at
crdist = 2*W                    # how close  banks get before cutoff in m


#Set Variables fro nonlocal efects
decay_rate = dt/(15*(365*24*60*60.0));   #ranges between 1/3 to 1/10, to be developed
bump_scale = 0              #to multiple kl by,amplitude of ne bump, range between 1 and 4, set to 0 for no nonlocal effects
cut_thresh = 100            #how many cutoffs to simulate, arbitrary if running for time


#Set Result Directory
result_dir = "sample_data/InitialChannel/" ##change this to wherever you want to save your results

name= "100width"


#Initiate Channel Object
ch = hkp.generate_initial_channel(W,D,deltas,pad)

#Initiate Channel Belt for migrating channel object

chb = hkp.ChannelBelt(channels=[ch],cutoffs=[],cl_times=[0.0],cutoff_times=[], cutoff_dists = [], decay_rate = decay_rate, bump_scale = bump_scale, cut_thresh = cut_thresh, sinuosity=[1])


#Plot initial centerline
chb.plot_channels()
plt.title(str(int(nit*dt/(365*24*60*60.0)))+ " years at "+ str(kl*(365*24*60*60.0))+ "m/yr")
plt.show()
 
#Migrate Centerline
chb.migrate_years(nit,saved_ts,deltas,pad,crdist,Cf,kl,dt) 

#Plot Resulting Centerline
chb.plot_channels()
plt.title(str(int(nit*dt/(365*24*60*60.0)))+ " years at "+ str(kl*(365*24*60*60.0))+ "m/yr")
plt.show()
#Save Sinuosity time series
times = chb.cl_times
sins = chb.sinuosity
#sinseries = pd.DataFrame({'time':times, 'sinuosity': sins});
#sinseries.to_csv(result_dir+"InitialCL_"+name+"sinseries"+".csv", header = True, index = False)
hkp.plot_sinuosity(times, sins)
plt.show()
#Save Resulting Centerline
xes = chb.channels[-1].x
yes = chb.channels[-1].y
cl = pd.DataFrame({'x': xes, 'y': yes});
cl.to_csv(result_dir+"InitialCL_"+name+".csv", header = False, index = False)