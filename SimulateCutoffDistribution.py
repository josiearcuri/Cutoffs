import math
import matplotlib.pyplot as plt
import cutoffs as co
import meanderpyalt as mp




#set variables
D = 10;                       
nit = 10000                   # number of iterations, run for at least 5000 years to see clustering
Cf = 0.022                # dimensionless Chezy friction factor
kl = 5/(365*24*60*60.0)   # migration rate constant (m/s)
kv =  1.0E-11              # vertical slope-dependent erosion rate constant (m/s)
dt = .5*365*24*60*60.0     # time step (s)
dens = 100                  # density of water (kg/m3)
saved_ts = 2               # which time steps will be saved every year
decay_rate = dt/(10*(365*24*60*60.0));   #ranges between 1/3 to 1/10, eventually this will not be a constant,
bump_scale = 0           #to multiple kl by, range between 1 and 3, set to 0 for no nonlocal effects
Sl = 0.01                    # initial slope (matters more for submarine channels than rivers)
pad= 30                     #depends on sample

mode = "NonlocalEffects"
if bump_scale == 0:
    mode = "OnlyCurvature"

result_dir = "sample_results/ReachC/" +mode+"/"##change this to wherevery you want to save your results
filelist = ['sample_data/firstcl_reachC_1984.csv','sample_data/width_reachC_1984.csv']

#Simulate migration on real centerline, keeoing track of cutoff locationa nd times#
#initialize first channel and channel belt 
[ch, x, y, z, cl_len, deltas] = mp.generate_channel_from_file(filelist, smooth_factor = .25)

crdist = 2.0*ch.W 

chb = mp.ChannelBelt(channels=[ch],cutoffs=[],cl_times=[0.0],cutoff_times=[], cutoff_dists = [], decay_rate = decay_rate, bump_scale = bump_scale)
chb.migrate(nit,saved_ts,deltas,pad,crdist,Cf,kl,kv,dt,dens) 

chb.plot('strat',20,60)
plt.title(str(int(nit*dt/(365*24*60*60.0)))+ " years at "+ str(kl*(365*24*60*60.0))+ "m/yr")
plt.savefig(result_dir +mode+ "channelafter"+str(int(nit*dt/(365*24*60*60.0))) +"years.jpg", dpi = 500)

## Statistically Test Cutoff Distributions for Clustering ##
cuts = co.cutoff_distributions(chb.cutoffs, int(nit*dt/(365*24*60*60.0)), result_dir)
co.plot_cutoff_distributions(cuts, int(nit*dt/(365*24*60*60.0)), result_dir)

co.mc_envelope(cutoffs = cuts, year = int(nit*dt/(365*24*60*60.0)), d_max = 1000, spacing = 25, resultdir=result_dir,mode = mode)
