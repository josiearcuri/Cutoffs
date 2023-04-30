# Meander Bend Cutoffs Clustering 
Josie Arcuri,  Doug Edmonds, Scott Robeson
*Indiana University, Department of Earth and Atmospheric Sciences*
> A model for simulating channel centerline migration that accounts for nonlocal effects from cutoffs, based on meanderpy[1]. We also include code to test hypothesis concerning what makes cutoffs cluster along channel centerlines through time. 

# HKplus
This module builds upon a previous implementation of Howard and Knutson's centerline migration model[1], meanderpy, based on local and weigthed upstream curvatures[2]. With its main structure nearly identical to meanderpy, this version considers how cutoffs affect bank migration beyond their impact on curvature. We include a general representation of how cutoffs enhance sediment transport along a cutoff bend for years after occurrence. In calculating a nominal migration rate, each point along a centerline is moved according to a migration rate constant muliplied by local curvature and the sum of any present nonlocal effects from cutoffs at that point. The spatial and temporal extent of each cutoff's nonlocal effects are set according to ours remote sensing analyses and previously measured cutoff nonlocal effects[3]. Following, nominal migration rate modified based on upstream curvature, which has been shown to reasonably predict bank migration rates[4].


# SpaceTime
 A set of functions used to test cutoff distributions for spatiotemporal clustering based on their location and time of occurrence.  This includes a monte carlo sampling method to test cutoff distibutions against randomly-generated point patterns[5].  

# Implementation
For an example of one experiment, Run 

    jupyter notebook cutoffs.ipynb 

To reproduce Figure 3:

    python Figure3.py

# Getting the code

You can download a copy of all the files in this repository by cloning this repository:

    git clone https://github.com/josiearcuri/Cutoffs.git/


## References
[1]Alan Howard;Thomas Knutson; Sufficient Conditions for River Meandering: A Simulation Approach. Water Resources Research (1984) 20: 1659–1667. DOI:10.1029/WR020i011p01659
[2]https://github.com/zsylvester/meanderpy  
[3]Zoltán Sylvester; Paul Durkin; Jacob A. Covault; High curvatures drive river meandering. Geology (2019) 47 (3): 263–266.
DOI: 10.1130/G45608.1  
[4]Jon Schwenk, Efi Foufoula‐Georgiou; Meander cutoffs nonlocally accelerate upstream and downstream migration and channel widening. Geophysical Research Letters (2016) 43 (24): 12,437-12,445. DOI: 10.1002/2016GL071670  
[5]Peter Diggle; Statistical Analysis of Spatial and Spatio-Temporal Point Patterns. CRC Press (2014), Third Edition. ISBN:978-1-4665-6023-9

