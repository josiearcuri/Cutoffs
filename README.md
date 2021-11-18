# Meander Bend Cutoffs Cluster from Self-Induced Migration
Josie Arcuri and Doug Edmonds  
*Indiana University, Department of Earth and Atmospheric Sciences*
> A model for simulating channel centerline migration that accounts for nonlocal effects from cutoffs, based on meanderpy[1]. We also include code to test hypothesis concerning what makes cutoffs cluster along channel centerlines through time. 

# Abstract
Cutoffs are iconic features of meandering rivers which reroute flow and limit the spatial extent of channels. This behavior influences the character of fluvial deposits preserved below the ground, while life on the surface of meandering floodplains relies on cutoff-generated oxbow lakes. Where and how often cutoffs occur along a river informs how floodplain ecosystems evolve and stratigraphy develops. Despite their relevance, what controls cutoff spacing and timing remains unclear. Here we show how nonlocal effects produced by a cutoff temporarily enhances a channel’s lateral migration rate and cause cutoffs to cluster in space and time. We measured nonlocal effects following cutoffs on 9 rivers in the Amazon Basin from Landsat-derived surface water maps.  Meander bends adjacent to cutoffs display lateral migration rates up to 3.5 times greater than elsewhere on the river for up to 12 years. On 5 out of 9 rivers, nonlocal effects appear to reduce time until the next cutoff occurs and cause some to cluster. To understand why only certain rivers exhibit cutoff clustering and others do not we add nonlocal effects to a well-known meandering bank migration model. With 28 experiments that vary migration rate constant, magnitude, and decay time scale of nonlocal effects, cutoff clustering is compared to random events. When nonlocal effects are ignored, cutoffs occur with regular spacing and timing along channel centerlines.  As the magnitude of nonlocal effects rises, cutoffs aggregate, and clusters emerge. Because cutoffs self-organize into clusters, channel deposits and oxbow lakes associated with cutoffs may also inherently aggregate.  

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

# Dependencies

    conda nev create -f cutoffenv.yml

## References
[1]Alan Howard;Thomas Knutson; Sufficient Conditions for River Meandering: A Simulation Approach. Water Resources Research (1984) 20: 1659–1667. DOI:10.1029/WR020i011p01659
[2]https://github.com/zsylvester/meanderpy  
[3]Zoltán Sylvester; Paul Durkin; Jacob A. Covault; High curvatures drive river meandering. Geology (2019) 47 (3): 263–266.
DOI: 10.1130/G45608.1  
[4]Jon Schwenk, Efi Foufoula‐Georgiou; Meander cutoffs nonlocally accelerate upstream and downstream migration and channel widening. Geophysical Research Letters (2016) 43 (24): 12,437-12,445. DOI: 10.1002/2016GL071670  
[5]Peter Diggle; Statistical Analysis of Spatial and Spatio-Temporal Point Patterns. CRC Press (2014), Third Edition. ISBN:978-1-4665-6023-9

