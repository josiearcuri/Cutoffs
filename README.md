# Cutoffs

A library for integrating real centerlines with MeanderPy, simulating migration that includes nonlocal effects from cutoffs, and analyzing resulting cutoff distributions. 


## meanderpyalt
This module builds upon MeanderPy[1], a centerline migration model based on local and weigthed upstream curvatures[2].   
This is an alternative version of meanderpy which encorporates nonlocal effects from cutoffs into centerline migration.  After a cutoff occurs, a gaussian bump is introduced to the location on the centerline where that cutoff occured. the amplitude of this bump is set by the user, as is its decay rate.  Both of these are constants for now.  The length upstream and downstream that this bump extends is equal to 1.19* the length removed by that cutoff [3].
## cutoffs
functions to..
* initiate meanderpy channel from existing centerline and widths
* retrive cutoffs distributions though space and time from meanderpy simulations
* analyze cutoff distrubutions as point process for clustering with Ripley's K function, compare to CSR process with monte carlo simulations

## Dependencies
Numpy, Scipy, Pandas, Astropy, Matplotlib, PIL, skimage, Seaborn 

## References
[1]https://github.com/zsylvester/meanderpy
[2]Zoltán Sylvester; Paul Durkin; Jacob A. Covault; High curvatures drive river meandering. Geology (2019) 47 (3): 263–266.
DOI: 10.1130/G45608.1
[3]Jon Schwenk, Efi Foufoula‐Georgiou; Meander cutoffs nonlocally accelerate upstream and downstream migration and channel widening. Geophysical Research Letters (2016) 43 (24): 12,437-12,445. DOI: 10.1002/2016GL071670

