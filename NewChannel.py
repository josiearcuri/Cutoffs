import math
import matplotlib.pyplot as plt
import cutoffs as co
import meanderpyalt as mp
import numpy as np
import pandas as pd

#set variables
D = 10;   
deltas = 50;
nit = 100                   # number of iterations, run for at least 5000 years to see clustering
Cf = 0.022                # dimensionless Chezy friction factor
kl = 5/(365*24*60*60.0)   # migration rate constant (m/s)
kv =  1.0E-11              # vertical slope-dependent erosion rate constant (m/s)
dt = 1*365*24*60*60.0     # time step (s)
dens = 100                  # density of water (kg/m3)# which time steps will be saved every year
decay_rate = dt/(10*(365*24*60*60.0));   #ranges between 1/3 to 1/10, eventually this will not be a constant,
W = 200
bump_scale = 0           #to multiple kl by, range between 1 and 3, set to 0 for no nonlocal effects
               # initial slope (matters more for submarine channels than rivers)
pad= 100                     #depends on sample
saved_ts = 10               # which time steps will be saved
                 # approximate number of bends you want to model
Sl = 0.0001 
mode = "old"

result_dir = "sample_results/InitialChannel/OnlyCurvature_new/" ##change this to wherevery you want to save your results

filepath = "C:/Users/Josie/Desktop/Cutoffs/sample_results/InitialChannel/OnlyCurvature/cl_onlycurvature_6000yrs.csv"
if mode == "new":
    ch = mp.generate_initial_channel(W,D,Sl,deltas,pad,100)
else:
    ch= mp.load_initial_channel(filepath, W, D, Sl, deltas)

mode = "OnlyCurvature"
crdist = ch.W 

chb = mp.ChannelBelt(channels=[ch],cutoffs=[],cl_times=[0.0],cutoff_times=[], cutoff_dists = [], decay_rate = decay_rate, bump_scale = bump_scale)

 

chb.migrate(nit,saved_ts,deltas,pad,crdist,Cf,kl,kv,dt,dens) 
chb.plot('strat',20,60,chb.cl_times[-1], len(chb.channels))
plt.title(str(int(nit*dt/(365*24*60*60.0)))+ " years at "+ str(kl*(365*24*60*60.0))+ "m/yr")
plt.show()


xes = chb.channels[-1].x
yes = chb.channels[-1].y
cl = pd.DataFrame({'x': xes, 'y': yes});
cl.to_csv(result_dir+"InitialCL_after3000yearsofonlycurvature.csv", header = False, index = False)
#chb.create_movie(np.min(ch.x), np.max(ch.x), "strat", "channel", "./NE_movie/channels/", 50, 50, 1, [chb.cl_times])#, chb.cutoff_times])
#co.save_animations("./NE_movie/effects/*.png", "./NE_movie/Effectmovie.gif")
#co.save_animations("./NE_movie/channels/*.png", "./NE_movie/Channelmovie.gif")

## Statistically Test Cutoff Distributions for Clustering ##
cuts = co.cutoff_distributions(chb.cutoffs, int(nit*dt/(365*24*60*60.0)), result_dir)
co.plot_cutoff_distributions(cuts, int(nit*dt/(365*24*60*60.0)), result_dir, mode = mode)