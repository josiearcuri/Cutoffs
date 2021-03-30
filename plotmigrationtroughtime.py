
import math
import matplotlib.pyplot as plt
import HKplus as hkp
import numpy as np
import pandas as pd



path = "sample_results/case2/17/"
name = "bendbybendmr.csv"
MR = pd.read_csv(path+name, sep = ',', index_col = 0)

hkp.plot_segmented_MR(MR.to_numpy())
plt.savefig(path+"bendbybend.png")
plt.show()