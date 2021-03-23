
"""This script runs statistical tests on model-output cutoff distributions - determines where cutoffs are more clustered or regularly spaced than randomly generated point patterns 
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import SpaceTime 
import os

#Set location of cutoff distribution to tests
file = "sample_results/case2/5/Case2_Run5_500_cutoffs_distribution.csv"


result_dir = "sample_results/case2/5/"

#read point pattern as events in space(distance downstream) and time (model years),
cutoffs = pd.read_csv(file, sep = ',')
year = int(np.max(cutoffs['time']))
length = int(np.max(cutoffs['downstream_distance']))


W = int((1.19*np.mean(cutoffs['cutlen'])))

dt = 2*int(year/500)
max_search_d = 1#int((1.19*np.mean(cutoffs['cutlen']))/W)
max_search_t = 1

print(length)
print(year)

#Check if data fits test criteria/boundary effects small enough to ignore
##print("Last cutoff occured at " + str(year)+ " years")
if year >= (max_search_t*dt)/4:
    print("time span is sufficently large for statistical tests")
else:
    print("Warning: model run long enough to only search over "+str(int(4*ymax_search_t*dt))+" year windows")
        

if int(length) >=(max_search_d*W)/4:
    print("centerline is sufficently long enough for statistical tests")
else: 
    print("Warning: centerline only long enough to search over " +str(int(4*max_search_d*W))+" ch-w windows")
#input()   
#Initialize Ripley's K test for 2d space-time
Kest = SpaceTime.RipleysKEstimator_spacetime(t_max=year, d_max=length, t_min=0, d_min=0, width = W, dt= dt)

#Run tests that output figures when complete
[D, D_sig]=Kest(cutoffs= cutoffs, mode = 'K_st', max_search_d= max_search_d, max_search_t = max_search_t, plotornot=1)

plt.savefig(result_dir +"Case2_Run5_"+str(len(cutoffs['time']))+"_cutoffs_Dplot.png", transparent=False, bbox_inches = 'tight')
plt.show()

print("D = " + str(D))
print("D_sig = " + str(D_sig))

