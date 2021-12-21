
"""This script runs statistical tests on model-output cutoff distributions - determines where cutoffs are more clustered or regularly spaced than randomly generated point patterns
@JA 2021
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import SpaceTime
import os

num = 7
#Set location of cutoff distribution to tests
file = "results/experiments/"+str(num)+"/"+str(num)+"_500_cutoffs_distribution.csv"#"sample_results/case2/"+str(num)+"/Case2_Run"+str(num)+"_500_cutoffs_distribution.csv"

result_dir = "sample_results/case2/"+str(num)+"/"

#read point pattern as events in space(distance downstream) and time (model years),
cutoffs = pd.read_csv(file, sep = ',')

year = int(np.max(cutoffs['time']))
length = int(np.max(cutoffs['downstream_distance']))
#bends = pd.read_csv(bends, sep = ',', header=None, index_col=False).values
#mean_n_bends = np.mean([np.count_nonzero(~np.isnan(bends[i, :])) for i in range(len(bends[:,1]))])
#print(mean_n_bends)
#print(10*length/mean_n_bends)

#bendlength = length/mean_n_bends
W = int(1.19*(np.mean(cutoffs['cutlen'])))

dt = 2*int(year/500)
max_search_d = 1
max_search_t = 1
#print(dt/2)
#input()
#print(length)
#print(year)

#Check if data fits test criteria/boundary effects small enough to ignore
##print("Last cutoff occured at " + str(year)+ " years")
if year >= (max_search_t*dt)/10:
    print("time span is sufficently large for statistical tests")
else:
    print("Warning: model run long enough to only search over "+str(int(4*max_search_t*dt))+" year windows")


if int(length) >=(max_search_d*W)/10:
    print("centerline is sufficently long enough for statistical tests")
else:
    print("Warning: centerline only long enough to search over " +str(int(4*max_search_d*W))+" ch-w windows")
#input()


#Initialize Ripley's K test for 2d space-time
Kest = SpaceTime.RipleysKEstimator_spacetime(t_max=year, d_max=length, t_min=0, d_min=0, width = W, dt= dt)

#Run tests that output figures when complete
[D,D_null, D_upper, D_lower]=Kest(cutoffs= cutoffs, mode = 'K_st', max_search_d= max_search_d, max_search_t = max_search_t, plotornot=1)
normed = (2*D/(D_upper-D_lower))
plt.savefig("results/experiments/"+str(num)+"/"+str(len(cutoffs['time']))+"_cutoffs_Dplot.png", transparent=False, bbox_inches = 'tight')

#ax = plt.gca()

plt.close()
print("D = " + str(D))

print("normed = " +str(normed))
