import cutoffs as co
import pandas as pd
import numpy as np
import SpaceTime as st

file = "sample_results/5000years_cutoff_distributions.csv"
file2 = "sample_results/5000years_cutoff_distributions.csv"
cutoffs = pd.read_csv(file, sep = ',')
cutoffs2 = pd.read_csv(file2, sep = ',')

year = 5000
resultdir = "sample_results/"
co.mc_envelope(cutoffs, year, spacing = 50, resultdir=resultdir, nit = 99, d_max = 1000,mode = ' nonlocal effects')
co.mc_envelope(cutoffs2, year, spacing = 50, resultdir=resultdir, nit = 99, d_max = 1000,mode = ' only curvature')