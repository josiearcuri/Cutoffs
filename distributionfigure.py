import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import HKplus as hkp


file = "sample_results/case1/1/Case1_Run1_500_cutoffs_distribution.csv"
W = 100
result_dir = "sample_results/case1/1/"

#read point pattern as events in space(distance downstream) and time (model years),
cutoffs = pd.read_csv(file, sep = ',')

hkp.plot_distribution(cutoffs, W, result_dir)
plt.savefig(result_dir+"Case1_Run11_cutoffdist.png", bbox_inches='tight')
plt.close()