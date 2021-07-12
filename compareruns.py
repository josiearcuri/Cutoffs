import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as tick
from scipy.stats import linregress as lr
import RSfun

file = "C:/Users/Josie/Desktop/Cutoffs/sample_results/modelrunsD.csv"
result_dir = "sample_results/"



#read point pattern as events in space(distance downstream) and time (model years),
cutoffs = pd.read_csv(file, sep = ',', index_col=None)

mask =(cutoffs['kl [m/yr]']>= 12) #& (cutoffs['NE Duration [yr]']==8)

obmr = cutoffs['obmr ch-1/yr'].values*mask
kl = cutoffs['kl [m/yr]'].values*mask/100
bump = cutoffs['NE Bump [*][(FR)/BR]'].values*mask -1
sumD = cutoffs['normed'].values*mask
duration = cutoffs['NE Duration [yr]'].values*mask
radii = cutoffs['mean radius of curvature [m]'].values*mask
length = cutoffs['mean cutoff length [m]'].values*mask

betweencuts = cutoffs['years between cutoffs [yr]']*mask
NE = bump

newmask = mask & ~np.isnan(sumD)
results = lr(NE[newmask>0], sumD[newmask>0])
print(sumD[newmask>0])
fig, ax1 = plt.subplots(1,figsize=(8, 4))
plt.rcParams.update({'font.size': 10})

ax1.axhline(y=0, xmin = 0,c='k',linewidth = 1, xmax = max(bump)+10, zorder =1, label = "randomness")
ax1.fill([-5, -5, 11, 11], [-1,1,1,-1], color = 'lightgrey', edgecolor = None, alpha = .5, label = "95% mc envelope", zorder = 0)
ax1.plot([0, int(np.max(NE[newmask>0]))], [results.intercept,results.intercept+int(np.max(NE[newmask>0]))*results.slope], ls = '--', linewidth = 1, c = 'k', zorder =1, label = "y = "+str(round(results.slope, 2))+"x + "+str(round(results.intercept, 2))+"\nr^2 = "+str(round(results.rvalue**2,3)))
sc1=ax1.scatter(NE[newmask>0], sumD[newmask>0], s=30, c=kl[newmask>0], edgecolor = 'k', cmap = 'cividis', linewidths = .5, alpha = .8, vmin = 0, label = '_nolegend_', zorder = 2)
#ax2.plot([0, 1], [results.intercept,results.intercept+3*results.slope], ls = '--', linewidth = .5, c= 'k', zorder =0)

ax1.legend(fontsize = 8, frameon = False, loc = "upper left")
#y_fmt = tick.FormatStrFormatter('%1.1e')
#ax1.yaxis.set_major_formatter(y_fmt)
#ax2.yaxis.set_major_formatter(y_fmt)
ax1.ticklabel_format(style='sci',scilimits=(0,0),axis='y')

ax1.set_ylabel('D *', rotation = 'horizontal')
ax1.set_xlim([-.5,3.5])
ax1.set_ylim([-2,2])
plt.xlabel('magnitude [*]')
ax1.tick_params(axis='x', labelsize = 8)

ax1.tick_params(axis='y', labelsize = 8)

c1 = fig.colorbar(sc1, ax = ax1, shrink = .6, pad = .1)
c1.ax.set_title('k1 [ch-w/yr]', fontsize = 9)

plt.savefig(result_dir+"Figure2_nondimed.png", dpi = 1000)
plt.show()