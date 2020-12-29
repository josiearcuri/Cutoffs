
import math
import matplotlib.pyplot as plt
import HKplus as hkp
import numpy as np
import pandas as pd

#set variables
D = 10;   
deltas = 50;
nit = 100              # number of iterations
Cf = 0.022                # dimensionless Chezy friction factor
kl = 5/(365*24*60*60.0)   # migration rate constant (m/s)

dt = 365*24*60*60.0     # time step (s)
decay_rate = dt/(10*(365*24*60*60.0));   #ranges between 1/3 to 1/10, eventually this will not be a constant,
W = 150
bump_scale = 4           #to multiple kl by, range between 1 and 3, set to 0 for no nonlocal effects
pad= 100                     #depends on sample
saved_ts = 10               # which time steps will be saved
                 # approximate number of bends you want to model

cut_thresh = 20
mode = "NonlocalEffects"

result_dir = "sample_results/5mpyr/" 
##change this to wherevery you want to save your 
filepath ="sample_data/InitialChannel/InitialCL_5mpyr.csv"
ch= hkp.load_initial_channel(filepath, W, D, deltas)

crdist = ch.W 

chb = hkp.ChannelBelt(channels=[ch],cutoffs=[],cl_times=[0.0],cutoff_times=[], cutoff_dists = [], decay_rate = decay_rate, bump_scale = bump_scale, cut_thresh= cut_thresh)

chb.migrate(nit,saved_ts,deltas,pad,crdist,Cf,kl,dt) 
chb.plot('strat',20,60,chb.cl_times[-1], len(chb.channels))
plt.title(str(int(nit*dt/(365*24*60*60.0)))+ " years at "+ str(kl*(365*24*60*60.0))+ "m/yr")
plt.show()
## Save Cutoff Distributions for Clustering Tests ##
chb.cutoff_distributions(int(nit*dt/(365*24*60*60.0)), result_dir, mode)

