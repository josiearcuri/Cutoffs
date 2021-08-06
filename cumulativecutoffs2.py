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
Dvals = results['D'].values#['D(d=1.19*Cutoff_Length, t=2*T/n)'].values
IDs = results['ID'].values
obmr = results['NE Duration [yr]'].values
starter = "sample_results/case4/"

paths = np.array(["13"])#, "33"])#, "24"])
#paths = np.array(["27", "26", "25"])# "23", "34", "30"])#, "24", "32", "34", "33", "26", "27"])
#paths = np.array([ "32", "33", "34"])
fig, axs = plt.subplots(1,2, figsize=(6,2))
D = np.array([Dvals[IDs == int(paths[i])][0] for i in range(len(paths))])


cmap = cm.get_cmap('cividis')#get_continuous_cmap(['#edf8b1','#7fcdbb', '#2c7fb8'])

newcolors = cmap(np.linspace(0.1, 1, len(paths)))

cmap = ListedColormap(newcolors)
    
cutofffile = pd.read_csv(starter+str(paths[0])+"/"+paths[0]+"300_cutoffs_distribution.csv")
    
downstream = cutofffile["downstream_distance"].values#[:100]
cllen = cutofffile["cllen"].values/100#[:int(len(cutofffile["cllen"])/2)]
times = cutofffile["time"].values#[:100]#[:int(len(cutofffile["cllen"])/2)]#[downstream<=np.percentile(downstream, 100)]
cutlen = cutofffile["cutlen"].values/100
ncuts = np.arange(1, len(times)+1)
dcuts = ncuts
    #dcuts[:]

frq = 300/times[-1]
dt = (times[1:]-times[:-1])#/np.mean((times[1:]-times[:-1]))
dcl = (cllen[1:] - cllen[:-1])#/np.mean(cllen[1:] - cllen[:-1])
grad = np.gradient(ncuts, times)
slope, intercept, r, p, se = stats.linregress(cutlen[:-1], dt)
    #plt.plot(times/times[-1], intercept + (slope*times), c = cmap(i), label='r = '+str(r))
    #axs.plot(times[:]/times[-1],ncuts - (slope*times+intercept), c = cmap(idx[i]/len(idx)),linewidth= .5,  alpha = .8, label = 'D = '+str(D[i]))
excess = dcl/dt-(slope*(cutlen[:-1]+cutlen[1:])/dt+intercept)
#axs[0].scatter((cutlen[:-1]), dt, c = dcl, cmap = "cividis",  alpha = .8, label = "r2 = "+str(round(r**2, 4)))
#axs[0].hist(dt)
#axs[0].plot(times, (slope*times+intercept), c = 'r',linewidth= 1,alpha = .5, label = "D = "+str(D))
    #axs.scatter(times[:]/times[-1],excess, color = cmap(idx[i]/len(idx)),linewidth= .5,  alpha = .8, label = 'D = '+str(D[i]) + " se = " + str(se))
slope, intercept, r, p, se = stats.linregress(ncuts, cllen)
axs[1].hist(cutlen, bins = np.linspace(0,round(np.max(cutlen))+1,num=int((round(np.max(cutlen))+1)/10)))
axs[0].hist(dt, bins = np.linspace(0,round(np.max(dt))+1, num=int((round(np.max(dt))+1))))
#axs[1].plot(times, cllen, c= 'b', linewidth = 1, alpha = .8, label = '_nolegend_')                                                                
#axs[1].plot(times, (slope*times+intercept), c= 'r', linewidth = 1, alpha = .5, label = '_nolegend_')
    #fig, ax, line = hkp.plot_segmented_MR(MR.to_numpy(), 1/obmr[IDs==int(paths[i])], fig, axs[i], 'b', "D = "+str(Dvals[IDs == int(paths[i])]))
axs[0].legend(loc= 'upper right', fontsize = 6, frameon = False)
    
#axs[0].set_xlim([0, 60000])    
axs[1].set_title("length removed by cutoff [ch-w]")
#axs[0].set_ylabel("centerline length change")

axs[0].set_title("cutoff period [years]")


plt.savefig(starter+str(paths[0])+"/"+paths[0]+"cutoffhists.png", dpi = 500)
plt.show()