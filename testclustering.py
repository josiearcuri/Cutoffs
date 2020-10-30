import cutoffs as co
import pandas as pd
import numpy as np

file = "C:/Users/Josie/Desktop/dump_nonlocaleffects_pad20/5000years_cutoff_distributions.csv"
cutoffs = pd.read_csv(file, sep = ',')

year = 5000
resultdir = "C:/Users/Josie/Desktop/"
co.mc_envelope(cutoffs, year, resultdir, nit = 99, mode = ' with nonlocal effects')