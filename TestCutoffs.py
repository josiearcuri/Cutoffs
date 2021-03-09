
"""This script runs statistical tests on model-output cutoff distributions - determines where cutoffs are more clustered or regularly spaced than randomly generated point patterns 
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import SpaceTime 
import os

#Set location of cutoff distribution to tests
file = "sample_results/case2/18/Case2_Run18_500_cutoffs_distribution.csv"


result_dir = "sample_results/case2/18/"

#read point pattern as events in space(distance downstream) and time (model years),
cutoffs = pd.read_csv(file, sep = ',')
year = int(np.max(cutoffs['time']))
length = int(np.max(cutoffs['downstream_distance']))


W = 100#int((1.19*np.mean(cutoffs['cutlen']))/100)

dt = 1
max_search_d = 50#int((1.19*np.mean(cutoffs['cutlen']))/W)
max_search_t = 50

print(length)
print(year)

#Check if data fits test criteria/boundary effects small enough to ignore
##print("Last cutoff occured at " + str(year)+ " years")
if year >= (max_search_t)**2:
    print("time span is sufficently large for statistical tests")
else:
    print("Warning: model run long enough to only search over "+str(int(np.sqrt(year)))+" year windows")
        

if int(length) >=(max_search_d*W)**2:
    print("centerline is sufficently long enough for statistical tests")
else: 
    print("Warning: centerline only long enough to search over " +str(int(np.sqrt(length/W)))+" ch-w windows")
input()   
#Initialize Ripley's K test for 2d space-time
Kest = SpaceTime.RipleysKEstimator_spacetime(t_max=year, d_max=length, t_min=0, d_min=0, width = W, dt= dt)

#Run tests that output figures when complete
[sumclust,sumreg, count]=Kest(cutoffs= cutoffs, mode = 'K_st', max_search_d= max_search_d, max_search_t = max_search_t)

plt.savefig(result_dir +"Case2_Run18_"+str(len(cutoffs['time']))+"_cutoffs_spacetimektest.png", transparent=False, bbox_inches = 'tight')
plt.show()
#plt.close()
print("sumclust = " + str(sumclust))
print("sumreg = " + str(sumreg))
print("count = " + str(count))
