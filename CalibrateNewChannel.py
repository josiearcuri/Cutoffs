
"""

This script uses the HKplus functions to generate a channel centerline, 
run it for a certain number of years, 
then save the resulting centerline for later use.  

@JA 2021
"""
import math
import matplotlib.pyplot as plt
import HKplus as hkp
import numpy as np
import pandas as pd

#Set Variables for centerline and curvature calculation

D=3.4                       #constant width-average channel depth (m)
W = 100                     #constant channel width (m)
deltas = W//2;              #spacing of nodes along centerline
nit = 3500                   #number of iterations/how many timesteps to migrate the centerline
Cf = 0.005                  #dimensionless Chezy friction factor
kl = 25/(365*24*60*60.0)    #migration rate constant (m/s)
dt = .25*365*24*60*60.0       #time step (s)
pad= 20                     #number of nodes for periodic boundary
saved_ts = 100               #timesteps between saving centerlines
crdist = 4*W                #how close banks get before cutoff in m
Length = ((50*W)*10)
 
#Set Variables fror Cutoff nonlocal efects
tau = 1#e-folding timescale for nonlocal effects in years
decay_rate = dt/(tau*(365*24*60*60.0));   #this is the half-life on nonlocal effects, in units of seconds
bump_scale = 0             #this is the magntiude of nonlocal effects in relative difference 


#Set Result Directory
result_dir = "data/InitialChannel/" ##change this to wherever you want to save your results

#Choose name for centerline
name= "7"

#Initiate Channel Object
ch = hkp.generate_initial_channel(W,D,deltas,pad,Length)

#Initiate Channel Belt for migrating channel object

chb = hkp.ChannelBelt(channels=[ch],cutoffs=[],cl_times=[0.0],cutoff_times=[], cutoff_dists = [], decay_rate = decay_rate, bump_scale = bump_scale, cut_thresh = 1000, sinuosity=[1])

#Plot initial centerline
#chb.plot_channels()
#plt.title(str(int(nit*dt/(365*24*60*60.0)))+ " years at "+ str(kl*(365*24*60*60.0))+ "m/yr")
#plt.show()
 
#Migrate Centerline for nit 
chb.migrate_years(nit,saved_ts,deltas,pad,crdist,Cf,kl,dt) 

#Plot Resulting Centerline
chb.plot_channels()
plt.title(str(int(nit*dt/(365*24*60*60.0)))+ " years at "+ str(kl*(365*24*60*60.0))+ "m/yr")
plt.show()

#Save Sinuosity time series
times = chb.cl_times
sins = chb.sinuosity

#uncomment to save sinuosity series
#sinseries = pd.DataFrame({'time':times, 'sinuosity': sins});
#sinseries.to_csv(result_dir+"InitialCL_"+name+"sinseries"+".csv", header = True, index = False)

hkp.plot_sinuosity(times, sins)
plt.show()

#Save Resulting Centerline
xes = chb.channels[-1].x
yes = chb.channels[-1].y
cl = pd.DataFrame({'x': xes, 'y': yes});
cl.to_csv(result_dir+"InitialCL_"+name+".csv", header = False, index = False)