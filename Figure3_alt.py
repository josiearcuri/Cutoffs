
"""
This script reproduces Figure 3 from 'Meander Bend Cutoffs Cluster from Self-induced Migration'

"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as tick
from scipy.stats import linregress as lr

def get_cutoff_period(ID):
    filepath = "C:/Users/Josie/Desktop/Cutoffs/results/experiments/"+str(ID)+"/"+str(ID)+"_500_cutoffs_distribution.csv"
    cutoffs = pd.read_csv(filepath, sep = ',', index_col = None)
    period = cutoffs['time'].values[-1]/500
    return period

#specify location of csv with model results
file = "C:/Users/Josie/Desktop/Cutoffs/results/experiments/results.csv"

#specify where you want to save the figure
result_dir = "results/figures/"

#read model results as a pandas dataframe
cutoffs = pd.read_csv(file, sep = ',', index_col=None).fillna(0)

#mask any rows which you do not want to plot
mask = (cutoffs['normed']!= 0)&(cutoffs['k1 [m/yr]'] >=20)#&(cutoffs['k1 [m/yr]'] >=20)


Dstar = cutoffs['normed'].values[mask] #for D values, use 'D' in column identifier, for D* values, use 'normed'
#only choose model runs with results

#specify channel width
W = 100
#extract arrays of main model parameters
ID = cutoffs['ID'].values[mask]
k1 = cutoffs['k1 [m/yr]'].values[mask]/W #in channel-widths
tau = cutoffs['NE timescale [yr]'].values[mask] #nonlocal effect half-life in years
M = cutoffs['NE magnitude [relative difference]'].values[mask]  #nonlocal effect magnitude in relative difference
period = np.asarray([get_cutoff_period(i) for i in ID])
freq = 1/period
NE = M*tau/period
tau = tau
#linear regression between NE and D(or D*)
results_NE = lr(NE, Dstar)

#linear regression between duration and D(or D*)
results_tau = lr(tau, Dstar)
#linear regression between duration and D(or D*)
results_freq = lr(freq, Dstar)

#linear regression between duration and D(or D*)
results_M = lr(M, Dstar)

#plotting set-up
fig, axes = plt.subplots(3,1,figsize=(9, 9))
[ax1, ax2, ax3] = axes
plt.rcParams.update({'font.size': 10})

#Upper axes

#ax1.fill([-5, -5, 11, 11], [-1,1,1,-1], color = 'lightgrey', edgecolor = None, alpha = .5, label = "Monte Carlo simulation envelope", zorder = 0)
ax3.plot([0, int(np.max(NE))], [results_NE.intercept,results_NE.intercept+int(np.max(NE))*results_NE.slope], ls = '--', linewidth = 1, c = 'k', zorder =1, label = "y = "+str(round(results_NE.slope, 2))+"x + "+str(round(results_NE.intercept, 2))+", $r^{ 2}$ = "+str(round(results_NE.rvalue**2,3)))
sc3=ax3.scatter(NE, Dstar, s=20, c=freq, edgecolor = 'k', cmap = 'cividis', linewidths = .5, alpha = .8, vmin = 0.1, vmax = .3, label = '_nolegend_', zorder = 2)
#ax2.plot([0, 1], [results.intercept,results.intercept+3*results.slope], ls = '--', linewidth = .5, c= 'k', zorder =0)



ax3.set_xlim([-1,10])
ax3.set_ylim([-2,2.5])
ax3.set_xlabel(r'$\tau$ * $M$ * $f$')

#Middle axes

#ax2.fill([-5, -5, 11, 11], [-1,1,1,-1], color = 'lightgrey', edgecolor = None, alpha = .5, label = "Monte Carlo simulation envelope", zorder = 0)
ax2.plot([0, int(np.max(tau))], [results_tau.intercept,results_tau.intercept+int(np.max(tau))*results_tau.slope], ls = '--', linewidth = 1, c = 'k', zorder =1, label = "y = "+str(round(results_tau.slope, 2))+"x + "+str(round(results_tau.intercept, 2))+", $r^{ 2}$ = "+str(round(results_tau.rvalue**2,3)))
sc2=ax2.scatter(tau, Dstar, s=20, c=freq, edgecolor = 'k', cmap = 'cividis', linewidths = .5, alpha = .8,  vmin = 0.1, vmax = .3, label = '_nolegend_', zorder = 2)


ax2.set_ylabel('D *', rotation = 'horizontal')
ax2.set_xlim([-.5,11])
ax2.set_ylim([-2,2.5])

ax2.set_xlabel(r'$\tau$')



#Lower axes

#ax2.fill([-5, -5, 11, 11], [-1,1,1,-1], color = 'lightgrey', edgecolor = None, alpha = .5, label = "Monte Carlo simulation envelope", zorder = 0)
ax1.plot([0, int(np.max(M))], [results_M.intercept,results_M.intercept+int(np.max(M))*results_M.slope], ls = '--', linewidth = 1, c = 'k', zorder =1, label = "y = "+str(round(results_M.slope, 2))+"x + "+str(round(results_M.intercept, 2))+", $r^{ 2}$ = "+str(round(results_M.rvalue**2,3)))
sc1=ax1.scatter(M, Dstar, s=20, c=freq, edgecolor = 'k', cmap = 'cividis', linewidths = .5, alpha = .8,  vmin = .1, vmax = .3, label = '_nolegend_', zorder = 2)


ax1.set_ylabel('D *', rotation = 'horizontal')
ax1.set_xlim([-.5,3.5])
ax1.set_ylim([-2,2.5])

ax1.set_xlabel('$M$')

#Format both axes
ax1.fill([-5, -5, 11, 11], [-1,1,1,-1], color = 'lightgrey', edgecolor = None, alpha = .5, label = "Monte Carlo simulation envelope", zorder = 0)
[ax.fill([-5, -5, 11, 11], [-1,1,1,-1], color = 'lightgrey', edgecolor = None, alpha = .5, label = "_nolegend", zorder = 0) for ax in axes[1:]]
ax1.axhline(y=0, xmin = 0,c='k',linewidth = 1, xmax = max(M)+10, zorder =1, label = "theoretical randomness")
[ax.axhline(y=0, xmin = 0,c='k',linewidth = 1, xmax = max(M)+10, zorder =1, label = "_nolegend_") for ax in axes[1:]]
[ax.legend(fontsize = 9, frameon = False, loc = "upper left") for ax in axes]
[ax.tick_params(axis='y', labelsize = 9) for ax in axes[1:]]
[ax.ticklabel_format(style='sci',scilimits=(0,0),axis='y') for ax in axes]

[ax.set_ylabel('D *', rotation = 'horizontal') for ax in axes]
[ax.tick_params(axis='x', labelsize = 9) for ax in axes]
#Format colorbar for both axes
c1 = fig.colorbar(sc1, ax = axes, shrink = .7)
c1.ax.set_title('$f$\n[1/year]\n', fontsize = 9, pad = .1)
c1.ax.tick_params(labelsize=8)

#Save figure
plt.savefig(result_dir+"Figure3.png", dpi =800)
plt.show()