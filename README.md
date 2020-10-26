# Cutoffs

A library for integrating real centerlines with MeanderPy, simulating migration that includes nonlocal effects from cutoffs, and analyzing resulting cutoff distributions. 


## meanderpyalt
This module builds upon MeanderPy, a centerline migration model based on curvature alone [1].   
This is an alternative version of meanderpy which encorporates nonlocal effects from cutoffs into centerline migration.  After a cutoff occurs, a gaussian bump is introduced to the location on the centerline where that cutoff occured. the amplitude of this bump is set by the user, as is its decay rate.  Both of these are constants for now.  The length upstream and downstream that this bump extends is equal to 1.19* the length removed by that cutoff [2].
## cutoffs
functions to..
### initiate meanderpy channel from existing centerline and widths
### retrive cutoffs distributions though space and time from meanderpy simulations
### anaylze cutoff distrubutions as point process or comlete spatial randomness, with Ripley's K function, compare to CSR process with monte carlo simulations

## Dependencies
Numpy, Scipy, Pandas, Astropy, Matplotlib, PIL, skimage, Seaborn 

## References
[1]https://github.com/zsylvester/meanderpy
[2]https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2016GL071670
