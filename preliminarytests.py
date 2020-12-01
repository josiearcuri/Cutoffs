import cutoffs as co
import pandas as pd
import numpy as np
from SpaceTime import PreliminaryTesting

file = "C:/Users/Josie/Desktop/Cutoffs/sample_results/InitialChannel/OnlyCurvature_new/2000years_cutoff_distributions.csv"
cutoffs = pd.read_csv(file, sep = ',')
year = 2000
print(np.max(cutoffs['downstream_distance']))
Hest = PreliminaryTesting(t_max=year, d_max=np.max(cutoffs['downstream_distance']), t_min=0, d_min=0)
Hest(cutoffs= cutoffs, mode = 'K')