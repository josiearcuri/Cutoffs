"""
This script reproduces Figure 3 from 'Meander Bend Cutoffs Cluster from Self-induced Migration'

"""
import numpy as np
import pandas as pd
import matplotlib as plt

import matplotlib.pyplot as plt
import matplotlib.ticker as tick
from scipy.stats import linregress as lr

#
# def get_cutoff_period(ID):
#     filepath = (
#         "/results/experiments/"
#         + str(ID)
#         + "/"
#         + str(ID)
#         + "_500_cutoffs_distribution.csv"
#     )
#     cutoffs = pd.read_csv(filepath, sep=",", index_col=None)
#     period = cutoffs["time"].values[-1] / 500
#     return period


def get_cutoff_length(ID):
    ID = int(ID)
    filepath = (
        "results/experiments/"
        + str(ID)
        + "/"
        + str(ID)
        + "_500_cutoffs_distribution.csv"
    )
    cutoffs = pd.read_csv(filepath, sep=",", index_col=None)
    cutofflength = np.mean(cutoffs["cutlen"].values)
    return cutofflength


def get_cutoff_radius(ID):
    ID = int(ID)

    filepath = (
        "results/experiments/"
        + str(ID)
        + "/"
        + str(ID)
        + "_500_cutoffs_distribution.csv"
    )
    cutoffs = pd.read_csv(filepath, sep=",", index_col=None)
    cutoffrad = np.mean(abs(cutoffs["radius"].values))
    return cutoffrad


def get_cutoff_period(ID):
    ID = int(ID)

    filepath = (
        "results/experiments/"
        + str(ID)
        + "/"
        + str(ID)
        + "_500_cutoffs_distribution.csv"
    )
    cutoffs = pd.read_csv(filepath, sep=",", index_col=None)
    cutoffperiod = cutoffs["cutlen"].values[-1] / 500
    return cutoffperiod


def get_cutoff_length(ID):
    ID = int(ID)

    filepath = (
        "results/experiments/"
        + str(ID)
        + "/"
        + str(ID)
        + "_500_cutoffs_distribution.csv"
    )
    cutoffs = pd.read_csv(filepath, sep=",", index_col=None)
    cutofflength = np.mean(cutoffs["cutlen"].values)
    return cutofflength


cmap = plt.cm.plasma
# specify location of csv with model results
file = "results/Ktestresults_V2_2_90th.csv"

# specify where you want to save the figure
result_dir = "results/figures/"

# read model results as a pandas dataframe
cutoffs = pd.read_csv(file, sep=",", index_col=None, header=2).fillna(0)
print(cutoffs.head())
# # mask any rows which you do not want to plot
mask = cutoffs["kl (m/yr)"] > 0  # 40


D = cutoffs["D"].values[mask]  # [
Dupper = cutoffs["D95"].values[mask]
Dlower = cutoffs["D5"].values[mask]
Dnull = cutoffs["Dnull"].values[mask]

#    mask
# for D values, use 'D' in column identifier, for D* values, use 'normed'
# only choose model runs with results

# specify channel width
W = 100
# extract arrays of main model parameters
ID = cutoffs["ID"].values[mask].astype("int")
print(ID)
kl = cutoffs["kl (m/yr)"].values[mask]  # / W  # in channel-widths
tau = cutoffs["tau (yr)"].values[mask]  # nonlocal effect half-life in years
M = cutoffs["M (relative difference ratio)"].values[
    mask
]  # nonlocal effect magnitude in relative difference
l = np.asarray([get_cutoff_length(i) for i in ID])
p = np.asarray([get_cutoff_period(i) for i in ID])
rad = np.asarray([get_cutoff_radius(i) for i in ID])
# cutoffs["R (m)"][mask] = rad
NE = kl * tau * M / (rad)
# cutoffs["phi"][mask] = NE

# NE = 1.19 * l /() tau * kl  # l)  # ** 2  # ** 0.5
# linear regression between NE and D(or D*)
results_NE = lr(NE, D)
cutoffs["PHI"] = np.nan
cutoffs["PHI"][mask] = NE
# # linear regression between duration and D(or D*)
# results_tau = lr(tau, Dstar)
# # linear regression between duration and D(or D*)
# results_freq = lr(freq, Dstar)
#
# # linear regression between duration and D(or D*)
# results_M = lr(M, Dstar)

# plotting set-up
fig, ax = plt.subplots(1, 1, figsize=(5, 2.5), tight_layout=True)
plt.rcParams.update({"font.size": 9})

# Upper axes

# ax1.fill([-5, -5, 11, 11], [-1,1,1,-1], color = 'lightgrey', edgecolor = None, alpha = .5, label = "Monte Carlo simulation envelope", zorder = 0)
# ax.plot(
#     [0, int(np.max(NE))],
#     [results_NE.intercept, results_NE.intercept + int(np.max(NE)) * results_NE.slope],
#     ls="--",
#     linewidth=1,
#     c="k",
#     zorder=1,
#     label="y = "
#     + str(round(results_NE.slope, 2))
#     + "x + "
#     + str(round(results_NE.intercept, 2))
#     + ", $r^{ 2}$ = "
#     + str(round(results_NE.rvalue**2, 3)),
# )
# ax.scatter(NE, D - Dnull)
# ax.hist(
#     NE[np.where((D > Dlower) & (D < Dupper))[0][:]],
#     color=cmap(1 / 3),
#     bins=10,
#     range=(0, 2.5),
#     # edgecolor="k",
#     alpha=0.6,
#     histtype="step",
#     cumulative=True,
#     #    vmin=0.1,
#     #    vmax=0.3,
#     label="random",
#     # zorder=2,
# )
nclust = np.sum((D > Dupper))
ndisp = np.sum((D < Dlower))
nrand = np.sum((D >= Dlower) & (D <= Dupper))
ax.hist(
    x=[
        NE[np.where((D < Dlower))[0][:]],
        NE[np.where((D > Dupper))[0][:]],
        NE[np.where((D >= Dlower) & (D <= Dupper))[0][:]],
    ],
    color=[cmap(2 / 3), cmap(1 / 3), cmap(3 / 3)],
    bins=8,
    range=(0, 3),
    edgecolor="k",
    # alpha=0.6,
    histtype="bar",
    rwidth=1,
    orientation="vertical",
    #    density=True,
    stacked=False,
    #    cumulative=True,
    #    vmin=0.1,
    #    vmax=0.3,
    label=[
        "Dispersed, n = " + str(int(ndisp)),
        "Clustered, n = " + str(int(nclust)),
        "Random, n = " + str(int(nrand)),
    ]
    # zorder=2,
)
# ax.hist(
#     NE[np.where((D <= Dlower))[0][:]],
#     color=cmap(1),
#     range=(0, 2.5),
#     bins=10,
#     # edgecolor="k",
#     histtype="step",
#     cumulative=True,
#     alpha=0.8,
#     #    vmin=0.1,
#     #    vmax=0.3,
#     label="dispersed",
# )
# zorder=2,
# ax2.plot([0, 1], [results.intercept,results.intercept+3*results.slope], ls = '--', linewidth = .5, c= 'k', zorder =0)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
# ax.spines['bottom'].set_visible(False)
# ax.spines['left'].set_visible(False)

# ax.set_xlim([-1, 10])
# ax.set_ylim([-2, 2.5])
ax.set_xlabel(r"$\phi$")
# ax.set_xlim((0, 3.01))
# ax.set_ylim((0, 15))

# Middle axes

ax.legend(fontsize=9, frameon=False, loc="upper right")
ax.tick_params(axis="y", labelsize=12)
# ax.ticklabel_format(style="sci", scilimits=(0, 0), axis="y")

ax.set_ylabel("count", labelpad=10)
ax.tick_params(axis="x", labelsize=12)
# Format colorbar for both axes
# c1 = fig.colorbar(sc, ax=ax, shrink=0.7)
# c1.ax.set_title(r"$k_l  (m/yr)$", fontsize=9, pad=1)
# c1.ax.tick_params(labelsize=8)
plt.title("Point Patterns Generated by Nonlocal Effects")
# Save figure
plt.savefig(result_dir + "Figure3_427.png", dpi=800)
print(cutoffs.columns)
plt.close()
cutoffs.to_csv("results\supp_table_427.csv")
