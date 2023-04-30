"""This script runs statistical tests on model-output cutoff distributions - determines where cutoffs are more clustered or regularly spaced than randomly generated point patterns
@JA 2021
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import SpaceTime
import os

# # Set location of cutoff distribution to tests
# file = (
#     "results/experiments/" + str(num) + "/" + str(num) + "_500_cutoffs_distribution.csv"
# )  # "sample_results/case2/"+str(num)+"/Case2_Run"+str(num)+"_500_cutoffs_distribution.csv"
#
#
# # read point pattern as events in space(distance downstream) and time (model years),
# cutoffs = pd.read_csv(file, sep=",")
#
# year = int(np.max(cutoffs["time"]))
# length = int(np.max(cutoffs["downstream_distance"]))
# # bends = pd.read_csv(bends, sep = ',', header=None, index_col=False).values
# # mean_n_bends = np.mean([np.count_nonzero(~np.isnan(bends[i, :])) for i in range(len(bends[:,1]))])
# # print(mean_n_bends)
# # print(10*length/mean_n_bends)
df = pd.DataFrame(
    data={
        "ID": [],
        "D": [],
        "null_D": [],
        "D05": [],
        "D90": [],
    },
    # index=0,
)


def get_cutoff_info(ID):
    filepath = (
        "results/experiments/"
        + str(ID)
        + "/"
        + str(ID)
        + "_500_cutoffs_distribution.csv"
    )
    resultdir = "results/experiments/" + str(ID) + "/Dtest_dt3.csv"

    cutoffs = pd.read_csv(filepath, sep=",", index_col=None)
    year = int(np.max(cutoffs["time"]))
    length = int(np.max(cutoffs["downstream_distance"]))
    W = int(1.19 * (np.mean(cutoffs["cutlen"])))
    dt = 1 * int(year / 500)

    return cutoffs, dt, year, length, W, resultdir


idlist = np.concatenate((range(4, 10), range(11, 67)))

for ID in idlist:
    cutoffs, dt, year, length, W, resultdir = get_cutoff_info(ID)
    # bendlength = length/mean_n_bends
    # W = int(1.19 * (np.mean(cutoffs["cutlen"])))
    #
    # dt = 2 * int(year / 500)
    max_search_d = 1
    max_search_t = 1
    # print(dt/2)
    # input()
    # print(length)
    # print(year)

    # Check if data fits test criteria/boundary effects small enough to ignore
    ##print("Last cutoff occured at " + str(year)+ " years")
    if year < (max_search_t * dt) / 10:
        # print("time span is sufficently large for statistical tests")
        #    else:
        print(
            "Warning: model run long enough to only search over "
            + str(int(4 * max_search_t * dt))
            + " year windows"
        )

    if int(length) < (max_search_d * W) / 10:
        #    print("centerline is sufficently long for statistical tests")
        # else:
        print(
            "Warning: centerline only long enough to search over "
            + str(int(4 * max_search_d * W))
            + " ch-w windows"
        )
    # input()

    # Initialize Ripley's K test for 2d space-time
    Kest = SpaceTime.RipleysKEstimator_spacetime(
        t_max=year, d_max=length, t_min=0, d_min=0, width=W, dt=dt
    )

    # Run tests that output figures when complete
    [D, D_null, D_upper, D_lower] = Kest(
        cutoffs=cutoffs,
        mode="K_st",
        max_search_d=max_search_d,
        max_search_t=max_search_t,
        plotornot=1,
    )
    # normed = (D - D_null) / D_null
    # plt.savefig("results/experiments/"+str(num)+"/"+str(len(cutoffs['time']))+"_cutoffs_Dplot.png", transparent=False, bbox_inches = 'tight')

    # ax = plt.gca()

    # plt.close()
    # print("D = " + str(D))
    # print("null_D = " + str(dt * W * 2))
    # print("clustered_D = " + str(D_upper))
    # print("regular_D = " + str(D_lower))

    # df.to_csv(resultdir)
    df = pd.concat(
        [
            df,
            pd.DataFrame(
                data={
                    "ID": [ID],
                    "D": D[0],
                    "null_D": [dt * W * 2],
                    "D95": D_upper[0],
                    "D5": D_lower[0],
                },
                #    index=ID,
            ),
        ]
    )
    print(str(ID) + " tested")

df.to_csv("results/Ktestresults_v119_1_90th.csv")
