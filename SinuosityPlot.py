"""
plot results from sinuosity test and save figure

@JA 2021
"""

import math
import matplotlib.pyplot as plt
import HKplus as hkp
import numpy as np
import pandas as pd

def plot_sin(sinseries, name):
    time = sinseries['time'].values
    sin = sinseries['sinuosity'].values
    dsdt = np.gradient(sin, time[1]-time[0])

    fig, ax = plt.subplots(2,1, sharex = True)
    line1 = ax[0].plot(time*4, sin, 'k')
    line2 = ax[1].plot(time*4, dsdt, 'k')
    ax[1].set_xlabel(r"time (yr)")
    ax[0].set_ylabel("S \n[$L_{c}$/$L_{v}$]\n")
    ax[1].set_ylabel("dS/dt\n[$\Delta$$L_{c}$/$L_{v}$/year]\n")
    fig.tight_layout()
    plt.savefig("./results/sinuositytest/sinuosityplot_"+str(name)+".png", dpi = 800)

    return


filepath = "~/Desktop/Cutoffs/results/sinuositytest/sin_test_0.csv"
df = pd.read_csv(filepath)

plot_sin(df, name = "0")
