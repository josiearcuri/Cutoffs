import matplotlib.pyplot as plt
import RSfun as rs

cutoff_file = "sample_data/RSCutoffData/fakecutoffs.csv"
reach_file = "sample_data/RSCutoffData/fakereaches.csv"

cutoffs_df = rs.merge_dfs(cutoff_file, reach_file)
cutoffs_df = rs.near_neigh(cutoffs_df)


rs.probs(cutoffs_df)
plt.show()