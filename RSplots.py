import matplotlib.pyplot as plt
import RSfun as rs


cutoff_file = "D:/MeanderBends/cutoffresults.csv"
reach_file = "D:/MeanderBends/reachresults.csv"

cutoffs_df = rs.merge_dfs(cutoff_file, reach_file)
cutoffs_df = rs.near_neigh(cutoffs_df)

print(cutoffs_df.columns)
rs.probs(cutoffs_df)
plt.savefig("sample_results/bumphist.png", transparent = True, dpi = 1000)
plt.show()