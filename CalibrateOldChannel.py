"""
This script uses the HKplus functions to load an existing centerline, allow that to migrate until a certain number of cutoffs have occurred, then save the resulting centerline for later use.  

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
nit = 2000                   #number of iterations/how many timesteps to migrate the centerline
Cf = 0.005                  #dimensionless Chezy friction factor
kl = 25/(365*24*60*60.0)    #migration rate constant (m/s)
dt = .5*365*24*60*60.0      #time step (s)
pad= 20                    #number of nodes for periodic boundary
saved_ts = 50               #timeteps between saving centerlines
crdist = 4*W                #how close banks get before cutoff in m
cut_thresh = 0
#Set Variables fror Cutoff nonlocal efects
tau = 1#e-folding timescale for nonlocal effects in years
decay_rate = dt/(tau*(365*24*60*60.0));   #this is the half-life on nonlocal effects, in units of seconds
bump_scale = 0             #this is the magntiude of nonlocal effects in relative difference 


#Load Existing Centerline

filepath = "data/InitialChannel/InitialCL_7.csv"

#suffix = "cal"
#Set Resulting file's suffix
result_dir = filepath#[:-4]+str(suffix)+".csv" ##change this to wherever you want to save your results


ch = hkp.load_initial_channel(filepath, W, D, deltas)

#Initialize Channel Belt for Migration
chb = hkp.ChannelBelt(channels=[ch],cutoffs=[],cl_times=[0.0],cutoff_times=[], cutoff_dists = [], decay_rate = decay_rate, bump_scale = bump_scale, cut_thresh = cut_thresh, sinuosity = [3])

#Plot initial channel
chb.plot_channels()
plt.title("Initial Centerline")
plt.show()
 
#Migrate Cneterline for cut_thresh cutoffs
#chb.migrate_cuts(saved_ts,deltas,pad,crdist,Cf,kl,dt) 

#Migrate Centerline for nit 
chb.migrate_years(nit,saved_ts,deltas,pad,crdist,Cf,kl,dt) 

#Plot resulting centerline
chb.plot_channels()
#plt.savefig("figure2_check.png", dpi = 800)
plt.show()


#Save Resulting Centerline
xes = chb.channels[-1].x
yes = chb.channels[-1].y
cl = pd.DataFrame({'x': xes, 'y': yes});

cl.to_csv(result_dir, header = False, index = False)


