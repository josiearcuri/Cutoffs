import math
import matplotlib.pyplot as plt
import HKplus as hkp
import numpy as np
import pandas as pd

#set variables
D = 10;   
deltas = 50;
nit = 6000              # number of iterations, run for at least 5000 years to see clustering
Cf = 0.022                # dimensionless Chezy friction factor
kl = 5/(365*24*60*60.0)   # migration rate constant (m/s)

dt = .5*365*24*60*60.0     # time step (s)
decay_rate = dt/(10*(365*24*60*60.0));   #ranges between 1/3 to 1/10, eventually this will not be a constant,
W = 150
bump_scale = 0           #to multiple kl by, range between 1 and 3, set to 0 for no nonlocal effects
pad= 100                     #depends on sample
saved_ts = 200               # which time steps will be saved every year



result_dir = "sample_results/" ##change this to wherevery you want to save your results

filepath ="sample_data/InitialChannel/InitialCL_5mpyr.csv"
ch= mp.load_initial_channel(filepath, W, D, deltas)
mode = "OnlyCurvature"
crdist = ch.W

chb = mp.ChannelBelt(channels=[ch],cutoffs=[],cl_times=[0.0],cutoff_times=[], cutoff_dists = [], decay_rate = decay_rate, bump_scale = bump_scale)
chb.plot('strat',20,60,chb.cl_times[-1], len(chb.channels))
plt.title(str(int(nit*dt/(365*24*60*60.0)))+ " years at "+ str(kl*(365*24*60*60.0))+ "m/yr")
plt.show()
 

chb.migrate(nit,saved_ts,deltas,pad,crdist,Cf,kl,dt,dens) 
chb.plot('strat',20,60,chb.cl_times[-1], len(chb.channels))
plt.title(str(int(nit*dt/(365*24*60*60.0)))+ " years at "+ str(kl*(365*24*60*60.0))+ "m/yr")
plt.show()


xes = chb.channels[-1].x
yes = chb.channels[-1].y
cl = pd.DataFrame({'x': xes, 'y': yes});
cl.to_csv(result_dir+"InitialCL_5mpyr.csv", header = False, index = False)