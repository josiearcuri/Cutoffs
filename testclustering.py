import cutoffs as co
import pandas as pd
import numpy as np
import SpaceTime as st

file = "C:/Users/Josie/Desktop/dump_justcurvature_pad20/5000years_cutoff_distributions.csv"
cutoffs = pd.read_csv(file, sep = ',')

year = 5000
resultdir = "C:/Users/Josie/Desktop/"
co.mc_envelope(cutoffs, year, spacing = 10, resultdir=resultdir, nit = 99, d_max = 100,mode = ' with nonlocal effects')