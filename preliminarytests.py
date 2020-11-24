import cutoffs as co
import pandas as pd
import numpy as np
import SpaceTime as st

file = "C:/Users/Josie/Desktop/Cutoffs/sample_results/InitialChannel/OnlyCurvature/2000years_cutoff_distributions.csv"
cutoffs = pd.read_csv(file, sep = ',')
year = 2000
co.mc_envelope(cutoffs, yea, resultdir=resultdir, nit = 99, d_max = 1000, mode = 'OnlyCurvature')