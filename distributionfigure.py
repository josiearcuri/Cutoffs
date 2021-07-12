import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import HKplus as hkp


file = "sample_results/case2/24/Case2_Run24_500_cutoffs_distribution.csv"
W = 100
result_dir = "sample_results/case2/24/"

#read point pattern as events in space(distance downstream) and time (model years),
cutoffs = pd.read_csv(file, sep = ',').iloc[range(60, 110)]

hkp.plot_distribution(cutoffs, W, result_dir)
plt.savefig(result_dir+"Case2_Run24_cutoffdist.png", dpi = 800)#, bbox_inches='tight')
plt.show()