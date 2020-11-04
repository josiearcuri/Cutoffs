import cutoffs as co
import pandas as pd
import numpy as np
import SpaceTime as st

file = "sample_results/NonlocalEffects/5000years_cutoff_distributions.csv"
#file2 = "sample_results/5000years_cutoff_distributions.csv"
cutoffs = pd.read_csv(file, sep = ',')


year = 5000
resultdir = "sample_results/NonlocalEffects/"
co.mc_envelope(cutoffs, year, spacing = 25, resultdir=resultdir, nit = 99, d_max = 10000,mode = 'NonlocalEffects')
