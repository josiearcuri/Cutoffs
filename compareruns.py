import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

file = "sample_results/modelruns.csv"
result_dir = "sample_results/"

#read point pattern as events in space(distance downstream) and time (model years),
cutoffs = pd.read_csv(file, sep = ',')

mask = cutoffs['cutoff critical distance [m]'] == 400

kl = cutoffs['kl [m/yr]']*mask
bump = cutoffs['NE Bump [*]']*mask
sumD = cutoffs['all']*mask
duration = cutoffs['NE Duration [yr]']*mask
radii = cutoffs['mean radii [m]']*mask
length = cutoffs['mean cutoff length [m]']*mask

betweencuts = cutoffs['years between cutoffs [yr]']
NE = bump
   
fig, (ax1, ax2) = plt.subplots(2,figsize=(12, 12))
plt.rcParams.update({'font.size': 10})
ax1.axhline(y=0, xmin = 0,c='k', xmax = max(bump)+10, zorder =0)
ax2.axhline(y=0, xmin = 0,c='k', xmax = max(bump)+10, zorder = 0)
sc1=ax1.scatter(NE, sumD, c=kl, cmap = 'viridis')
sc2=ax2.scatter(NE, sumD, c=duration/betweencuts, cmap = 'magma')

ax1.set_ylabel('K-null')
ax2.set_ylabel('K-null')
plt.xlabel('Bump Amplitude [*]')

c1 = fig.colorbar(sc1, ax = ax1)
c1.set_label('migration rate')

c2 = fig.colorbar(sc2, ax = ax2)
c2.set_label('duration/time between cuts')
plt.show()