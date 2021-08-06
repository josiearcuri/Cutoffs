"""
This script plots cumulative cutoffs through time from simulated migration outputs.  Lines are colored by their Ripley's K value.  
"""
import math
import matplotlib.pyplot as plt
import HKplus as hkp
import numpy as np
import pandas as pd
from matplotlib import cm
from matplotlib.colors import ListedColormap
from scipy import stats


results = pd.read_csv("sample_results/junemodelruns.csv", sep = ',', index_col = False)

mask = (results['ID'].values>=0)
Dvals = results['relative'].values*mask
IDs = results['ID'].values*mask
obmr = results['NE Duration [yr]'].values*mask
starter = "sample_results/case4/"

paths = np.array(["26", "28", "32", "13", "35", "07", "08"])#, "33"])#, "24"])
#paths = np.array(["27", "26", "25"])# "23", "34", "30"])#, "24", "32", "34", "33", "26", "27"])
#paths = np.array([ "32", "33", "34"])
fig, axs = plt.subplots(len(paths),1, figsize=(5,8), sharex=True, sharey = True)
D = np.array([Dvals[IDs == int(paths[i])][0] for i in range(len(paths))])
idx = (D).argsort()
print(D)
print(idx)

cmap = cm.get_cmap('cividis')#get_continuous_cmap(['#edf8b1','#7fcdbb', '#2c7fb8'])

newcolors = cmap(np.linspace(0.1, 1, len(paths)))

cmap = ListedColormap(newcolors)
    
for i in range(len(paths)):

    cutofffile = pd.read_csv(starter+str(paths[idx[i]])+"/"+paths[idx[i]]+"300_cutoffs_distribution.csv")
    
    downstream = cutofffile["downstream_distance"].values
    #cllen = cutofffile["cllen"].values[:int(len(cutofffile["cllen"])/2)]
    times = cutofffile["time"].values[:500]#[:int(len(cutofffile["cllen"])/2)]#[downstream<=np.percentile(downstream, 100)]
    ncuts = np.arange(1, len(times)+1)
    dcuts = ncuts
    #dcuts[:]

    frq = 300/times[-1]
    dt = (times[1:]-times[:-1])
    grad = np.gradient(ncuts, times)
    slope, intercept, r, p, se = stats.linregress(times, ncuts)
    #plt.plot(times/times[-1], intercept + (slope*times), c = cmap(i), label='r = '+str(r))
    #axs.plot(times[:]/times[-1],ncuts - (slope*times+intercept), c = cmap(idx[i]/len(idx)),linewidth= .5,  alpha = .8, label = 'D = '+str(D[i]))
    excess = (ncuts)-(slope*times +intercept)
    axs[i].plot(times[:]/times[-1], excess, c = cmap(i),linewidth= 1,  alpha = .8, label = "r2 = "+str(round(r**2,4)))
    #axs.scatter(times[:]/times[-1],excess, color = cmap(idx[i]/len(idx)),linewidth= .5,  alpha = .8, label = 'D = '+str(D[i]) + " se = " + str(se))
    axs[i].plot([0,1],[0,0],linestyle = ":", c ='k',linewidth = .5,  alpha = .8, label = "y = "+str(round(slope, 2))+"x + "+str(round(intercept, 2)))#+", r**2 = "+str(round(r**2, 2)))
    #slope, intercept, r, p, se = stats.linregress(times, 1/cllen)
    #axs[1].plot(times/times[-1], cllen , ls = '--', c= 'k', linewidth = .5, alpha = .8, label = '_nolegend_')
    #fig, ax, line = hkp.plot_segmented_MR(MR.to_numpy(), 1/obmr[IDs==int(paths[i])], fig, axs[i], 'b', "D = "+str(Dvals[IDs == int(paths[i])]))
    axs[i].legend(loc= 'upper right', fontsize = 6, frameon = False)
    axs[i].set_title("D = "+str(D[idx][i]))
    
    
axs[-1].set_xlabel("time/ total time [*]")

axs[1].set_ylabel("modeled - expected [# cutoffs]")

plt.savefig("Figure3_cutoffrate.png", dpi = 500)
plt.show()