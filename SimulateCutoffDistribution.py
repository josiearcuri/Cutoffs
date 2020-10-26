
import math
import matplotlib.pyplot as plt

import cutoffs as co
import meanderpyalt as mp




#set variables
D = 10;                       
nit = 200# number of iterations
Cf = 0.022                # dimensionless Chezy friction factor
kl = 5/(365*24*60*60.0)   # migration rate constant (m/s)
kv =  1.0E-11              # vertical slope-dependent erosion rate constant (m/s)
dt = .5*365*24*60*60.0     # time step (s)
dens = 100                  # density of water (kg/m3)
saved_ts = 2               # which time steps will be saved every year
decay_rate = 1/10;         #ranges between 1/3 to 1/10, eventually this will not be a constant
bump_scale = 1.5            #to multiple kl by, range between 1 and 3
Sl = 0.01                    # initial slope (matters more for submarine channels than rivers)
pad= 20                     #depends on sample

result_dir = "results/"
filelist = ['data/Reach6CL1984.csv','data/Reach6CL_widths1984.csv']


#initialize first channel and channel belt 
[ch, x, y, z, cl_len, deltas] = co.generate_channel_from_file(filelist, smooth_factor = .5)

crdist = 2.0*ch.W 

chb = mp.ChannelBelt(channels=[ch],cutoffs=[],cl_times=[0.0],cutoff_times=[], cutoff_dists = [], decay_rate = decay_rate, bump_scale = bump_scale)
chb.migrate(nit,saved_ts,deltas,pad,crdist,Cf,kl,kv,dt,dens) 

chb.plot('strat',20,60)
plt.title(str(int(nit*dt/(365*24*60*60.0)))+ " years at "+ str(kl*(365*24*60*60.0))+ "m/yr")
plt.savefig(result_dir + "channelafter"+str(int(nit*dt/(365*24*60*60.0))) +"years.png")
cuts = co.cutoff_distributions(chb.cutoffs, int(nit*dt/(365*24*60*60.0)), result_dir)
co.plot_cutoff_distributions(cuts, int(nit*dt/(365*24*60*60.0)), result_dir)
co.mc_envelope(cuts, year=int(nit*dt/(365*24*60*60.0)), resultdir = result_dir, nit = 99, mode = ' modeled') 


