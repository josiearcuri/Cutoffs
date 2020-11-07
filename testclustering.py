import cutoffs as co
import pandas as pd
import numpy as np
import SpaceTime as st

file = "C:/Users/Josie/Desktop/Cutoffs/sample_results/InitialChannel/OnlyCurvature/2000.csv"
#file2 = "sample_results/5000years_cutoff_distributions.csv"
cutoffs = pd.read_csv(file, sep = ',')

co.plot_cutoff_distributions(cutoffs, 2000, "sample_results/InitialChannel/")
year = 2000
resultdir = "sample_results/InitialChannel/"
co.mc_envelope(cutoffs, year, resultdir=resultdir, nit = 99, d_max = 1,mode = 'NonlocalEffects')
