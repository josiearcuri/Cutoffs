import cutoffs as co
import pandas as pd
import numpy as np
import SpaceTime 

file = "sample_results/2mpyr/NonlocalEffects200_cutoffs_distribution.csv"

#file = "sample_results/5mpyr/OnlyCurvature/OnlyCurvature400_cutoffs_distribution.csv"
cutoffs = pd.read_csv(file, sep = ',')
year = int(np.max(cutoffs['time'])+1)
print(int(np.max(cutoffs['time'])))
Kest = SpaceTime.RipleysKEstimator_spacetime(t_max=year, d_max=int(np.max(cutoffs['downstream_distance'])+1), t_min=0, d_min=0, width = 150)
Kest(cutoffs= cutoffs, mode = 'K_st')