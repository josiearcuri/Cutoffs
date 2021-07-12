
import math
import matplotlib.pyplot as plt
import HKplus as hkp
import numpy as np
import pandas as pd


results = pd.read_csv("sample_results/modelrunsD.csv", sep = ',', index_col = False)
Dvals = results['D(d=1.19*Cutoff_Length, t=2*T/n)'].values
IDs = results['ID'].values
obmr = results['NE Duration [yr]'].values
starter = "sample_results/case2/"

path1 = "32" 
path2 = "33" 
path3 = "34"
path4 = "24" 
path5 = "25" 
path6 = "23"

paths = np.array(["23","24","25", "27", "30","32", "33", "34"])
#path4 = "sample_results/case2/23/"
#most regular = 31, 25, 24
#most clustered = 32, 33, 22
#random = 34, 30, 29
name = "bendbybendmr.csv"
D = np.array([Dvals[IDs == int(paths[i])][0] for i in range(len(paths))])
print(D)
input()
fig, axs = plt.subplots(len(paths),1, figsize=(6,6), sharey = True, sharex = True)
idx = D.argsort()
print(idx)
input()
for i in range(len(paths)):
    
    MR = pd.read_csv(starter+paths[idx[i]]+"/"+name, sep = ',', index_col = 0)
    fig, ax, line = hkp.plot_segmented_MR(MR.to_numpy(), int(500/len(MR.to_numpy()[:,0])), fig, axs[i], 'k', "D = "+str(round(D[idx[i]], 3)))
    #fig, ax, line = hkp.plot_segmented_MR(MR.to_numpy(), 1/obmr[IDs==int(paths[i])], fig, axs[i], 'b', "D = "+str(Dvals[IDs == int(paths[i])]))

    axs[i].set_xlim((0,1))
#1/obmr[IDs==int(paths[i])]
#MR1 = pd.read_csv(path1+name, sep = ',', index_col = 0)
#fig, axs = hkp.plot_segmented_MR(MR1.to_numpy(), fig, axs, 'r', 'regular')
axs[0].set_title('Outer Bank Migration Periodograms')
axs[-1].set_xlabel('Signal Frequency / Cutoff Frequency [*]')
axs[int(len(axs)/2)].set_ylabel(' log spectral power\n[m^2/s^2]')
#axs[0].set_xlim((0,100))
#axs[1].set_xlim((0,100))
#axs[2].set_xlim((0,100))

#plt.savefig('figure3_periodograms.png', dpi = 500)
plt.show()
