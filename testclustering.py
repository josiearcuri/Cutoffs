import cutoffs as co
import pandas as pd

file = "C:/Users/Josie/Desktop/dump_nonlocaleffects_pad20/9000years_cutoff_distributions.csv"
cutoffs = pd.read_csv(file, sep = ',')
print(cutoffs.head())
year = np.max(cutoffs.time)
resultdir = "C:/Users/Josie/Desktop/"
co.mc_envelope(cutoffs, year, resultdir, nit = 99, mode = ' modeled')