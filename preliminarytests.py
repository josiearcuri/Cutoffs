import cutoffs as co
import pandas as pd
import numpy as np
import SpaceTime 

file = "C:/Users/Josie/Desktop/Cutoffs/sample_results/5mpyr/OnlyCurvature/OnlyCurvature9_cutoffs_distribution.csv"
cutoffs = pd.read_csv(file, sep = ',')
year = int(np.max(cutoffs['time']))
print(int(np.max(cutoffs['time'])))
Kest = SpaceTime.RipleysKEstimator_spacetime(t_max=year, d_max=np.max(cutoffs['downstream_distance']), t_min=0, d_min=0)
Kest(cutoffs= cutoffs, mode = 'G')