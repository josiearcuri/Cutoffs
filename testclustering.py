import cutoffs as co
import pandas as pd
import numpy as np
import SpaceTime as st

file = "C:/Users/Josie/Desktop/Cutoffs/sample_results/InitialChannel/OnlyCurvature/2000years_cutoff_distributions.csv"
#file2 = "sample_results/5000years_cutoff_distributions.csv"
cutoffs_curv = pd.read_csv(file, sep = ',')

#co.plot_cutoff_distributions(cutoffs_curv, 2000, "C:/Users/Josie/Desktop/Cutoffs/sample_results/InitialChannel/OnlyCurvature/", mode = 'OnlyCurvature')
year = 2000
#resultdir ="C:/Users/Josie/Desktop/Cutoffs/sample_results/InitialChannel/OnlyCurvature/"
#co.mc_envelope(cutoffs, year, resultdir=resultdir, nit = 99, d_max = 1000,mode = 'OnlyCurvature')


file = "C:/Users/Josie/Desktop/Cutoffs/sample_results/InitialChannel/NonlocalEffects_case1/2000years_cutoff_distributions.csv"
#file2 = "sample_results/5000years_cutoff_distributions.csv"
cutoffs_ne = pd.read_csv(file, sep = ',')

#co.plot_cutoff_distributions(cutoffs_ne, 2000, "C:/Users/Josie/Desktop/Cutoffs/sample_results/InitialChannel/NonlocalEffects_case2/", mode = 'NonlocalEffects_case2')
resultdir ="C:/Users/Josie/Desktop/Cutoffs/sample_results/InitialChannel/NonlocalEffects_case1/"
co.mc_envelope_comp(cutoffs_curv,cutoffs_ne, year, resultdir=resultdir, nit = 99, d_max = 100)