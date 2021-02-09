
"""This script runs statistical tests on model-output cutoff distributions - determines where cutoffs are more clustered or regularly spaced than randomly generated point patterns 
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import SpaceTime 
import os

#Set location of cutoff distribution to test
file = "sample_results/25/Run25_300_cutoffs_distribution.csv"
W = 100
result_dir = "sample_results/25/"
#read point pattern as events in space(distance downstream) and time (model years),
cutoffs = pd.read_csv(file, sep = ',')
year = int(np.max(cutoffs['time']))
length = int(np.max(cutoffs['downstream_distance']))

#Check if data fits test criteria/boundary effects small enough to ignore
##print("Last cutoff occured at " + str(year)+ " years")
if year >= 50**2:
    print("time span is sufficently large for statistical tests")
else:
    print("Warning: model run long enough to only search over "+str(int(np.sqrt(year)))+" year windows")
        

if int(length/W) >= 50**2:
    print("centerline is sufficently long enough for statistical tests")
else: 
    print("Warning: centerline only long enough to search over " +str(int(np.sqrt(length/W)))+" ch-w windows")
input()   
#Initialize Ripley's K test for 2d space-time
Kest = SpaceTime.RipleysKEstimator_spacetime(t_max=year, d_max=length, t_min=0, d_min=0, width = W)

#Run tests that output figures when complete
[D, count]=Kest(cutoffs= cutoffs, mode = 'K_st', max_search_d= 50,max_search_t = 30)
plt.savefig(result_dir +"Run25_"+str(len(cutoffs['time']))+"_cutoffs_spacetimektest.png", transparent=True, bbox_inches = 'tight')
plt.close()
print("D = " + str(D))
print("count = " + str(count))
