
import meanderpyalt as mp
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math
import cutoffs as co
import os

#set variables
D = 10;                       
nit = 200                   # number of iterations, run for at least 5000 years to see clustering
Cf = 0.022                # dimensionless Chezy friction factor
kl = 5/(365*24*60*60.0)   # migration rate constant (m/s)
kv =  1.0E-11              # vertical slope-dependent erosion rate constant (m/s)
dt = .5*365*24*60*60.0     # time step (s)
dens = 100                  # density of water (kg/m3)
saved_ts = 2               # which time steps will be saved every year
decay_rate = 1;         #ranges between 1/3 to 1/10, eventually this will not be a constant
bump_scale = 0           #to multiple kl by, range between 1 and 3, set to 0 for no nonlocal effects
Sl = 0.01                    # initial slope (matters more for submarine channels than rivers)
pad= 100                     #depends on sample

result_dir = "C:/Users/Josie/Desktop/dump_justcurvature/" ##change this to wherevery you want to save your results
filelist = ['sample_data/Reach6CL1984.csv','sample_data/Reach6CL_widths1984.csv']

#Simulate migration on real centerline, keeoing track of cutoff locationa nd times#
#initialize first channel and channel belt 
[ch, x, y, z, cl_len, deltas] = mp.generate_channel_from_file(filelist, smooth_factor = .5)

crdist = 2.0*ch.W 

chb = mp.ChannelBelt(channels=[ch],cutoffs=[],cl_times=[0.0],cutoff_times=[], cutoff_dists = [], decay_rate = decay_rate, bump_scale = bump_scale)
years = []
clusters_array = []
regulars_array= []
n_cuts_array = []
#migrate in batch, keep track of when clustering or regular spacing is detected
for j in range(0,100):
    chb.migrate(nit,saved_ts,deltas,pad,crdist,Cf,kl,kv,dt,dens) 

    if np.mod(j+1,10)==0:
        cuts = co.cutoff_distributions(chb.cutoffs, year*(j+1), result_dir)
        n_cuts_array.append(len(cuts.time))
        co.plot_cutoff_distributions(cuts, (j+1)*year, result_dir)
        cluster_flag, regular_flag = co.mc_envelope(cuts, year=(1+j)*year, resultdir = result_dir, nit = 200) 
        years.append((j+1)*year)
        clusters_array.append(cluster_flag)
        regulars_array.append(regular_flag)
        chb.plot('strat',20,60)
        plt.title(str(year)+ " years at "+ str(kl*(365*24*60*60.0))+ "m/yr")
        plt.savefig(result_dir+str((j+1)*year) + "yrs.png")

nonrandoms = pd.DataFrame({'years': years, 'n_cuts': n_cuts_array, 'clusters': clusters_array, 'regulars': regulars_array})
nonrandoms.to_csv(result_dir+"batch_results.csv", index = False)
batch_specs = pd.DataFrame({'nit': nit, 'background migration rate': kl*(365*24*60*60.0), 'dt': dt/(365*24*60*60.0), 'pad': pad, "cutoff critical distance": crdist, "decay rate":-decay_rate, "foreground migrations rate": kl*(365*24*60*60.0)*bump_scale}, index = [0])
batch_specs.to_csv(result_dir+'batch_specs.csv', index = False)