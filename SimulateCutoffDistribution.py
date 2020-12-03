
import math
import matplotlib.pyplot as plt
import cutoffs as co
import meanderpyalt as mp
import numpy as np
import pandas as pd

#set variables
D = 10;   
deltas = 50;
nit = 2000              # number of iterations, run for at least 5000 years to see clustering
Cf = 0.022                # dimensionless Chezy friction factor
kl = 5/(365*24*60*60.0)   # migration rate constant (m/s)
kv =  1.0E-11              # vertical slope-dependent erosion rate constant (m/s)
dt = .5*365*24*60*60.0     # time step (s)
dens = 100                  # density of water (kg/m3)# which time steps will be saved every year
decay_rate = dt/(10*(365*24*60*60.0));   #ranges between 1/3 to 1/10, eventually this will not be a constant,
W = 150
bump_scale = 4           #to multiple kl by, range between 1 and 3, set to 0 for no nonlocal effects
pad= 100                     #depends on sample
saved_ts = 200               # which time steps will be saved
                 # approximate number of bends you want to model
Sl = 0.0001 
cut_thresh = 100
mode = "NonlocalEffects"

result_dir = "C:/Users/Josie/Desktop/Cutoffs/sample_results/5mpyr/NonlocalEffects/" 
##change this to wherevery you want to save your 
filepath ="sample_results/InitialChannel/InitialCL_5mpyr.csv"
ch= mp.load_initial_channel(filepath, W, D, Sl, deltas)

crdist = ch.W 

chb = mp.ChannelBelt(channels=[ch],cutoffs=[],cl_times=[0.0],cutoff_times=[], cutoff_dists = [], decay_rate = decay_rate, bump_scale = bump_scale, cut_thresh= cut_thresh)

chb.migrate(nit,saved_ts,deltas,pad,crdist,Cf,kl,kv,dt,dens) 
chb.plot('strat',20,60,chb.cl_times[-1], len(chb.channels))
plt.title(str(int(nit*dt/(365*24*60*60.0)))+ " years at "+ str(kl*(365*24*60*60.0))+ "m/yr")
plt.show()
## Save Cutoff Distributions for Clustering Tests ##
cuts = co.cutoff_distributions(chb.cutoffs, int(nit*dt/(365*24*60*60.0)), result_dir, mode)
co.plot_cutoff_distributions(cuts, int(nit*dt/(365*24*60*60.0)), result_dir, mode)

