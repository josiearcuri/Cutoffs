import math
import matplotlib.pyplot as plt
import HKplus as hkp
import numpy as np
import pandas as pd

#set variables
D = 10;   
deltas = 50;
nit = 1000              # number of iterations, run for at least 5000 years to see clustering
Cf = 0.022                # dimensionless Chezy friction factor
kl = 10/(365*24*60*60.0)   # migration rate constant (m/s)
kv =  1.0E-11              # vertical slope-dependent erosion rate constant (m/s)
dt = .5*365*24*60*60.0     # time step (s)
dens = 100                  # density of water (kg/m3)# which time steps will be saved every year
decay_rate = dt/(10*(365*24*60*60.0));   #ranges between 1/3 to 1/10, eventually this will not be a constant,
W = 150
bump_scale = 0           #to multiple kl by, range between 1 and 3, set to 0 for no nonlocal effects
pad= 100                     #depends on sample
saved_ts = 200               # which time steps will be saved
                 # approximate number of bends you want to model
Sl = 0.0001 


result_dir = "sample_results/InitialChannel/" ##change this to wherevery you want to save your results
#ch = mp.generate_initial_channel(W,D,Sl,deltas,100,100)

filepath ="sample_results/InitialChannel/InitialCL_4mpyr.csv"
ch= hkp.load_initial_channel(filepath, W, D, Sl, deltas)
mode = "OnlyCurvature"
crdist = ch.W 

chb = hkp.ChannelBelt(channels=[ch],cutoffs=[],cl_times=[0.0],cutoff_times=[], cutoff_dists = [], decay_rate = decay_rate, bump_scale = bump_scale, cut_thresh = 200)
chb.plot('strat',20,60,chb.cl_times[-1], len(chb.channels))
plt.title(str(int(nit*dt/(365*24*60*60.0)))+ " years at "+ str(kl*(365*24*60*60.0))+ "m/yr")
plt.show()
 

chb.migrate(nit,saved_ts,deltas,pad,crdist,Cf,kl,kv,dt,dens) 
chb.plot('strat',20,60,chb.cl_times[-1], len(chb.channels))
plt.title(str(int(nit*dt/(365*24*60*60.0)))+ " years at "+ str(kl*(365*24*60*60.0))+ "m/yr")
plt.show()


xes = chb.channels[-1].x
yes = chb.channels[-1].y
cl = pd.DataFrame({'x': xes, 'y': yes});
#cl.to_csv(result_dir+"InitialCL_4mpyr.csv", header = False, index = False)