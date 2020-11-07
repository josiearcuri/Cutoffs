import math
import matplotlib.pyplot as plt
import cutoffs as co
import meanderpyalt as mp
import numpy as np



#set variables
D = 10;   
deltas = 50;
nit = 100                   # number of iterations, run for at least 5000 years to see clustering
Cf = 0.022                # dimensionless Chezy friction factor
kl = 50/(365*24*60*60.0)   # migration rate constant (m/s)
kv =  1.0E-11              # vertical slope-dependent erosion rate constant (m/s)
dt = .5*365*24*60*60.0     # time step (s)
dens = 100                  # density of water (kg/m3)# which time steps will be saved every year
decay_rate = dt/(10*(365*24*60*60.0));   #ranges between 1/3 to 1/10, eventually this will not be a constant,
W = 200
bump_scale = 0           #to multiple kl by, range between 1 and 3, set to 0 for no nonlocal effects
               # initial slope (matters more for submarine channels than rivers)
pad= 100                     #depends on sample
saved_ts = 20                # which time steps will be saved
                 # approximate number of bends you want to model
Sl = 0.0 
mode = "NonlocalEffects"
if bump_scale == 0:
    mode = "OnlyCurvature"

result_dir = "sample_results/InitialChannel/" +mode+"/"##change this to wherevery you want to save your results
filelist = ['sample_data/firstcl_reach6_1984.csv','sample_data/width_reach6_1984.csv']

#Simulate migration on real centerline, keeoing track of cutoff locationa nd times#
#initialize first channel and channel belt 
#[ch, x, y, z, cl_len, deltas] = mp.generate_channel_from_file(filelist, smooth_factor = .25)
ch = mp.generate_initial_channel(W,D,Sl,deltas,pad,100)

crdist = 2*ch.W 

chb = mp.ChannelBelt(channels=[ch],cutoffs=[],cl_times=[0.0],cutoff_times=[], cutoff_dists = [], decay_rate = decay_rate, bump_scale = bump_scale)
chb.migrate(nit,saved_ts,deltas,pad,crdist,Cf,kl,kv,dt,dens) 

chb.plot('strat',20,60)
plt.title(str(int(nit*dt/(365*24*60*60.0)))+ " years at "+ str(kl*(365*24*60*60.0))+ "m/yr")
plt.savefig(result_dir +mode+ "channelafter"+str(int(nit*dt/(365*24*60*60.0))) +"years.jpg", dpi = 800)

## Statistically Test Cutoff Distributions for Clustering ##
#cuts = co.cutoff_distributions(chb.cutoffs, int(nit*dt/(365*24*60*60.0)), result_dir)
#co.plot_cutoff_distributions(cuts, int(nit*dt/(365*24*60*60.0)), result_dir)

#co.mc_envelope(cutoffs = cuts, year = int(nit*dt/(365*24*60*60.0)), d_max = 1000, spacing = 25, resultdir=result_dir,mode = mode)
#chb.create_movie(0, np.max(chb.channels[-1].x), "strat", "35yrs_dt1_reachC_", result_dir+"movie/", 20, #60, 1, chb.cl_times[::5]+chb.cutoff_times)
#co.save_animations(result_dir+"movie/*.png", result_dir+"35yrs_dt1_reachC_anim.gif")
