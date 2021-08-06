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
nit = 600                   #number of iterations/how many timesteps to migrate the centerline
Cf = 0.005                  #dimensionless Chezy friction factor
kl = 20/(365*24*60*60.0)    #migration rate constant (m/s)
dt = .25*365*24*60*60.0      #time step (s)
pad= 200                    #number of nodes for periodic boundary
saved_ts = 20               #timeteps between saving centerlines
crdist = 4*W                #how close banks get before cutoff in m

#Set Variables for nonlocal efects
decay_rate = dt/((5)*(365*24*60*60.0));   #ranges between 1/3 to 1/10, to be developed 
bump_scale = 0              #to multiple kl by,amplitude of ne bump, range between 1 and 4, set to 0 for no nonlocal effects
cut_thresh = 5              #how many cutoffs to simulate, arbitrary if running for time


#Load Existing Centerline

filepath = "sample_data/InitialChannel/InitialCL_experiment007.csv"

#Set Resulting file's suffix
result_dir = filepath[:-4]+str(suffix)+".csv" ##change this to wherever you want to save your results


ch = hkp.load_initial_channel(filepath, W, D, deltas)

#Initialize Channel Belt for Migration
chb = hkp.ChannelBelt(channels=[ch],cutoffs=[],cl_times=[0.0],cutoff_times=[], cutoff_dists = [], decay_rate = decay_rate, bump_scale = bump_scale, cut_thresh = cut_thresh, sinuosity = [3.5])

#Plot initial channel
chb.plot_channels()
plt.title("Initial Centerline")
plt.show()
 
#Migrate
chb.migrate_cuts(saved_ts,deltas,pad,crdist,Cf,kl,dt) 

#Plot resulting centerline
chb.plot_channels2(2,"/~/")
#plt.savefig("figure2_check.png", dpi = 800)
plt.show()


#Save Resulting Centerline
xes = chb.channels[-1].x
yes = chb.channels[-1].y
cl = pd.DataFrame({'x': xes, 'y': yes});

cl.to_csv(result_dir, header = False, index = False)


