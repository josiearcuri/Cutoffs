import numpy as np
import pandas as pd
import time, sys

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
from matplotlib import cm


from scipy.stats import norm
import scipy.interpolate
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.spatial import distance



def update_progress(progress, start_time):
    """progress bar from https://stackoverflow.com/questions/3160699/python-progress-bar
    update_progress() : Displays or updates a console progress bar
    Accepts a float between 0 and 1. Any int will be converted to a float.
    A value under 0 represents a 'halt'.
    A value at 1 or bigger represents 100%"""
    barLength = 20 # Modify this to change the length of the progress bar
    status = ""
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
        status = "error: progress var must be float\r\n"
    if progress < 0:
        progress = 0
        status = "Halt...\r\n"
    if progress >= 1:
        progress = 1
        status = "Done...\r\n"
    block = int(round(barLength*progress))
    text = "\rPercent: [{0}] {1}% {2}".format( "#"*block + "-"*int(round(barLength-block)), int(progress*100), status)
    sys.stdout.write(text + "-- %s minute(s) --" % int(round((time.time() - start_time)/60)))
    sys.stdout.flush()

class Channel:
    """class for Channel objects"""
    def __init__(self,x,y,W,D):
        """initialize Channel object
        x, y, z  - coordinates of centerline
        W - channel width
        D - channel depth"""
        self.x = x
        self.y = y
        self.W = W
        self.D = D

        

class Cutoff:
    """class for Cutoff objects"""
    def __init__(self,x,y,W,dist, time, cut_len):
        """initialize Cutoff object
        x, y,  - coordinates of centerline
        W - channel width
        """
        self.x = x
        self.y = y
        self.W = W
        self.dist = np.max(dist)
        self.time = time
        self.cut_len = cut_len
        
class ChannelBelt:
    """class for ChannelBelt objects"""
    def __init__(self, channels, cutoffs, cl_times, cutoff_times, cutoff_dists, decay_rate, bump_scale, cut_thresh):
        """initialize ChannelBelt object
        channels - list of Channel objects
        cutoffs - list of Cutoff objects
        cl_times - list of ages of Channel objects
        cutoff_times - list of ages of Cutoff objects
        cutoff_dists - list of cutoff distances downstream
        decay_rate - rate at which nonlocal effects dissipate 
        bump_scale - amplitude of nonlocal effects bump, scalar
        cut_thresh - how many cutoffs to simulate"""
        self.channels = channels
        self.cutoffs = cutoffs
        self.cl_times = cl_times
        self.cutoff_times = cutoff_times
        self.cutoff_dists = cutoff_dists
        self.decay_rate = decay_rate
        self.bump_scale = bump_scale
        self.cut_thresh = cut_thresh

    def migrate_years(self,nit,saved_ts,deltas,pad, crdist,Cf,kl,dt,dens=1000):
        """function for computing migration rates along channel centerlines and moving them, limited by number of iterations
        inputs:
        nit - number of iterations 
        saved_ts - which time steps will be saved
        deltas - distance between nodes on centerline
        pad - padding for upstream bc (number of nodepoints along centerline)
        crdist - threshold distance at which cutoffs occur
        Cf - dimensionless Chezy friction factor
        kl - migration rate constant (m/s)
        dt - time step (s)"""
        start_time = time.time()

        channel = self.channels[-1] # first channel is the same as last channel of input
        x = channel.x; y = channel.y; W = channel.W; D = channel.D; 
        
        k = 1.0 # constant in HK equation
        xc = [] # initialize cutoff coordinates
        yc = []
        cut_dist = []# initialize cutoff distance ds array
        cut_len = []# initialize cutoff length removal array
        # determine age of last channel:
        if len(self.cl_times)>0:
            last_cl_time = self.cl_times[-1]
        else:
            last_cl_time = 0
        dx, dy, ds, s = compute_derivatives(x,y)
        omega = -1.0 # constant in curvature calculation (Howard and Knutson, 1984)
        gamma = 2.5 # from Ikeda et al., 1981 and Howard and Knutson, 1984
        ne = np.zeros_like(x) #array to keep track of nonlocal effects
        ymax = self.bump_scale*kl*2
        itn = 0
        for itn in range(nit): # main loop
            update_progress(itn/nit, start_time)
            ne = update_nonlocal_effects(ne, s, self.decay_rate, self.bump_scale, cut_dist, cut_len) #update array of ne with last itn's cutoff(s) and decay old ne
            klarray = nominal_rate(kl, ne)## compute array of nominal migration rate in m/s with nonlocal effects accounted for
            x, y = migrate_one_step(x,y,W,klarray,dt,k,Cf,D,pad,omega,gamma)
            x,y,xc,yc,cut_dist, cut_len = cut_off_cutoffs(x,y,s,crdist,deltas) # find and execute cutoffs
            x,y,dx,dy,ds,s = resample_centerline(x,y,deltas) # resample centerline
            
            if len(xc)>0: # save cutoff data
                cutoff = Cutoff(xc,yc,W,cut_dist,last_cl_time+(itn)*dt/(365*24*60*60.0), cut_len) # create cutoff object
                #keep track of year cutoff occurs, where it occurs, and save an object. 
                self.cutoff_times.append(last_cl_time+(itn)*dt/(365*24*60*60.0))
                self.cutoff_dists.append(cut_dist)
                self.cutoffs.append(cutoff)
            # saving centerlines:
            if np.mod(itn,saved_ts)==0:
                channel = Channel(x,y,W,D) # create channel object, save year
                self.cl_times.append(last_cl_time+(itn)*dt/(365*24*60*60.0))
                self.channels.append(channel)
    def migrate_cuts(self,saved_ts,deltas,pad, crdist,Cf,kl,dt,dens=1000):
        """function for computing migration rates along channel centerlines and moving them, limited by number of cutoffs the channel experiences
        inputs:
        saved_ts - which time steps will be saved
        deltas - distance between nodes on centerline
        pad - padding for upstream bc (number of nodepoints along centerline)
        crdist - threshold distance at which cutoffs occur
        Cf - dimensionless Chezy friction factor
        kl - migration rate constant (m/s)
        dt - time step (s)"""
        start_time = time.time()

        channel = self.channels[-1] # first channel is the same as last channel of input
        x = channel.x; y = channel.y; W = channel.W; D = channel.D; 
        
        k = 1.0 # constant in HK equation
        xc = [] # initialize cutoff coordinates
        yc = []
        cut_dist = []# initialize cutoff distance ds array
        cut_len = []# initialize cutoff length removal array
        # determine age of last channel:
        if len(self.cl_times)>0:
            last_cl_time = self.cl_times[-1]
        else:
            last_cl_time = 0
        dx, dy, ds, s = compute_derivatives(x,y)
        omega = -1.0 # constant in curvature calculation (Howard and Knutson, 1984)
        gamma = 2.5 # from Ikeda et al., 1981 and Howard and Knutson, 1984
        ne = np.zeros_like(x) #array to keep track of nonlocal effects
        ymax = self.bump_scale*kl*2
        itn = 0
        while len(self.cutoffs)<self.cut_thresh: # main l oop
            itn = itn+1
            update_progress(len(self.cutoffs)/self.cut_thresh, start_time) 
            ne = update_nonlocal_effects(ne, s, self.decay_rate, self.bump_scale, cut_dist, cut_len) #update array of ne with last itn's cutoff(s) and decay old ne
            klarray = nominal_rate(kl, ne)## compute array of nominal migration rate in m/s with nonlocal effects accounted for
            x, y = migrate_one_step(x,y,W,klarray,dt,k,Cf,D,pad,omega,gamma)
            x,y,xc,yc,cut_dist, cut_len = cut_off_cutoffs(x,y,s,crdist,deltas) # find and execute cutoffs
            x,y,dx,dy,ds,s = resample_centerline(x,y,deltas) # resample centerline
            
            if len(xc)>0: # save cutoff data
                cutoff = Cutoff(xc,yc,W,cut_dist,last_cl_time+(itn)*dt/(365*24*60*60.0), cut_len) # create cutoff object
                #keep track of year cutoff occurs, where it occurs, and save an object. 
                self.cutoff_times.append(last_cl_time+(itn)*dt/(365*24*60*60.0))
                self.cutoff_dists.append(cut_dist)
                self.cutoffs.append(cutoff)
            # saving centerlines:
            if np.mod(itn,saved_ts)==0 or len(self.cutoffs)>=self.cut_thresh:
                channel = Channel(x,y,W,D) # create channel object, save year
                self.cl_times.append(last_cl_time+(itn)*dt/(365*24*60*60.0))
                self.channels.append(channel)
    def plot_channels(self):
        cot = np.array(self.cutoff_times)
        sclt = np.array(self.cl_times)
        times = np.unique(np.sort(np.hstack((cot,sclt))))
        
        # set up min and max x and y coordinates of the plot:
        xmin = np.min(self.channels[0].x)
        xmax = 1
        ymax = 1
        for i in range(len(self.channels)):
            ymax = max(ymax, np.max(np.abs(self.channels[i].y)))
            xmax = max(xmax, np.max(np.abs(self.channels[i].x)))

        ymax = ymax+1000# add a bit of space on top and bottom
        ymin = -1*ymax
        # size figure so that its size matches the size of the model:
        fig, ax = plt.subplots(figsize=(10,(ymax-ymin)*10/(xmax-xmin)))
        cmap = cm.get_cmap('gray_r',len(sclt))
        ax.set_xlim([xmin-1000,xmax+1000])
        ax.set_ylim([ymin,ymax])

        plt.axis('equal')
        ax.plot([xmin, xmin+5000],[ymin, ymin], 'k', linewidth=2)
        ax.text(xmin+1500, ymin+200+100, '5 km', fontsize=8)
        order = 0
        for i in range(0,len(times)):
            if times[i] in cot:
                ind = np.where(cot==times[i])[0][0]
                for j in range(0,len(self.cutoffs[ind].x)):
                    x1 = self.cutoffs[ind].x[j]
                    y1 = self.cutoffs[ind].y[j]
                    W = self.channels[-1].W
                    xm, ym = get_channel_banks(x1,y1,W)
                    order += 1
                    if times[i]==cot[-1] and j==(len(self.cutoffs[ind].x)-1):
                        plt.fill(xm,ym,facecolor='r',edgecolor = 'none', alpha = .5,zorder=order, label= 'cutoff')

                    else:
                        plt.fill(xm,ym,facecolor='r',edgecolor = 'none', alpha = .5,zorder=order, label= '_nolegend_')
                    
            if times[i] in sclt:
                ind = np.where(sclt==times[i])[0][0]        
                x1 = self.channels[ind].x
                y1 = self.channels[ind].y
                W = self.channels[ind].W
                xm, ym = get_channel_banks(x1,y1,W)
                order += 1
                plt.fill(xm,ym,facecolor=cmap(ind/len(sclt)),edgecolor = 'k', linewidth=0.1,zorder=order, label= '_nolegend_')
            

        ax.legend(frameon = False, loc = 'lower left',bbox_to_anchor=(7000/(xmax+2000), 0, .1, .1), markerscale = .1)
        ax.axis('off')

        return fig
    
    def cutoff_distributions(self, year, filepath, mode):
        """pull cutoff data from channel belt object and export csv, return dataframe for plotting
        year - last centerline year
        filepath - where cutoff info csv and distrubution plot are saved
        mode - for plotting - "OnlyCurvature or "NonlocalEffects"
        """
        #pull cutoff locations, downstrem distance, time from channel belt as a pandas dataframe
        distances = [i.dist for i in self.cutoffs]
        times = [i.time for i in self.cutoffs]
        x_location = [i.x for i in self.cutoffs]
        y_location = [i.y for i in self.cutoffs]
        cuts = pd.DataFrame({'downstream_distance': distances, 'time': times, 'x': x_location, 'y': y_location})
        
        #save distribution to csv
        newcuts = cuts.to_csv(filepath+mode+str(len(cuts['time']))+"_cutoffs_distribution.csv", index_label = "Cutoff")
        plot_cuts(cuts,self.channels[-1].W, filepath)
        return cuts
def plot_cuts(cuts,W, filepath):
    fig = plt.figure(figsize = (5,5))
    plt.rcParams.update({'font.size': 10})

    plt.scatter(cuts['downstream_distance']/W,cuts['time'], c = 'black', s = 1.5, edgecolor = 'black')
    ncuts = len(cuts['time'])
    plt.ylabel("time (years)")
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    plt.xlabel("distance downstream (ch-w)")
    return fig     
    
def resample_centerline(x,y,deltas):
    dx, dy, ds, s = compute_derivatives(x,y) # compute derivatives
    # resample centerline so that 'deltas' is roughly constant
    # [parametric spline representation of curve]
    tck, u = scipy.interpolate.splprep([x,y],s=0) 
    unew = np.linspace(0,1,1+int(s[-1]/deltas)) # vector for resampling
    out = scipy.interpolate.splev(unew,tck) # resampling
    x, y = out[0], out[1] # assign new coordinate values
    dx, dy, ds, s = compute_derivatives(x,y) # recompute derivatives
    return x,y,dx,dy,ds,s

def nominal_rate(kl, ne):
    """update nominal migration rate with nonlocal effects array"""
    new_kl = kl*(1+ne)
    return new_kl
def migrate_one_step(x,y,W,klarray,dt,k,Cf,D,pad,omega,gamma):
    dx, dy, ds, s = compute_derivatives(x,y)
    curv = W*compute_curvature(x,y)# dimensionless curvature
    R0 = klarray*curv #nominal migration rate with local curvature
    alpha = k*2*Cf/D # exponent for convolution function G
    R1 = compute_migration_rate(pad,len(x),ds,alpha,omega,gamma,R0)
    # calculate new centerline coordinates:
    dy_ds = dy/ds
    dx_ds = dx/ds
    # move x and y coordinates:
    x = x + R1*dy_ds*dt  
    y = y - R1*dx_ds*dt 
    return x,y

def generate_initial_channel(W,D,deltas,pad,n_bends):
    """generate straight Channel object with some noise added that can serve
    as input for initializing a ChannelBelt object
    from MeanderPy
    W - channel width
    D - channel depth
    deltas - distance between nodes on centerline
    pad - padding (number of nodepoints along centerline)
    n_bends - approximate number of bends to be simulated"""
    noisy_len = n_bends*10*W/2.0 # length of noisy part of initial centerline
    pad1 = pad//10
    #padding at upstream end can be shorter than padding on downstream end
    if pad1<5:
        pad1 = 5
    x = np.linspace(0, noisy_len+(pad+pad1)*deltas, int(noisy_len/deltas+pad+pad1)+1) # x coordinate
    y = 10.0 * (2*np.random.random_sample(int(noisy_len/deltas)+1,)-1)
    y = np.hstack((np.zeros((pad1),),y,np.zeros((pad),))) # y coordinate
    return Channel(x,y,W,D)

def load_initial_channel(filepath, W, D, deltas):
    """generate initial channel from centerline csv that can serve
    as input for initializing a ChannelBelt object.  must be fine enough resolution to not warrant smoothing
    filepath - csv with x, y coordinates and no headers
    W - channel width
    D - channel depth
    Sl - channel gradient
    deltas - distance between nodes on centerline"""
    df = pd.read_csv(filepath, sep = ',', header=None).values
    x = df[:,0]
    y = df[:,1]
    return Channel(x,y,W,D)
def generate_channel_from_file(filelist, D_in= 10, smooth_factor=.25, matlab_corr= -1):
    """function for creating a MeanderPy Channel object from an externally-sourced centerline in .csv file format.
        inputs:
        filelist - filelist must be a list of filepaths.  thie first should be a csv containing x and y values for each point on a centerline.  The second should be the widths of each point along the centerline
        slope - average channel slope (default = .01)
        D_in - channel depth (m, default = 10)
        smooth_factor - fraction of centerline points to sample from spline (default = 1/4)
        matlab_corr - 1 if y-axis need not be flipped, -1 for centerlines exported from matlab and need flipping
        
        outputs:
        ch - MeanderPy object of channel centerline
        x - uninterpolated array of x coordinates
        y - uninterpolated array of y coordinates
        cl_len - length of  centerline(m) """
    #use pandas to load x,y coordinates and widths of centerline from csv
    varlist = [pd.read_csv(file, sep = ',', header=None).values for file in filelist]
    x = varlist[0][:,0]*30  ## x-dim array in Sm
    y = varlist[0][:,1]*matlab_corr*30 ##southern hemisphere y-dim array in m
    
    # shift coordinates so all are positive
    if sum(n < 0 for n in y) > 0:
        y = y-min(y)## y-dim earray in m 
    if sum(n < 0 for n in x) > 0:
        x = x-min(x)
        
    #average over widths to get a reach-constant width scalar
    W = np.mean(varlist[1][:,0])*30
    print(W)
    ## water depth scalar#
    D = D_in  
    # Linear length along the line, add a zero for first point:
    points = np.vstack([x, y]).T
    distance = np.cumsum( np.sqrt(np.sum( np.diff(points, axis=0)**2, axis=1 )) )
    distance = np.insert(distance, 0, 0)

    # Build a list of the spline function, one for each dimension:
    splines = [InterpolatedUnivariateSpline(distance, coords) for coords in points.T]

    # Compute the spline for the smoothed(sampled) distances:
    points_fitted = np.vstack([spl(np.linspace(0, distance[-1],round(len(x)*smooth_factor))) for spl in splines])
    
    ## z-dim array, interpolated with constant slope along points of centerline.  assumes centerline points are equidistantly placed along original centerline. 
    #deltas = round(distance[-1]/(len(points_fitted[0])-1)) 
    return Channel(points_fitted[0],points_fitted[1],W,D)


def compute_migration_rate(pad,ns,ds,alpha,omega,gamma,R0):
    """compute migration rate as weighted sum of upstream curvatures
    pad - padding (number of nodepoints along centerline)
    ns - number of points in centerline
    ds - distances between points in centerline
    omega - constant in HK model
    gamma - constant in HK model
    R0 - nominal migration rate (dimensionless curvature * migration rate constant)"""
    R1 = np.zeros(ns) # preallocate adjusted channel migration rate
    pad_up = ns-pad
    #########Periodic Boundary#########################
    for i in range(2,pad):
        si2 = np.hstack((np.array([0]),np.cumsum(np.hstack((ds[i-1::-1], ds[ns-1:pad_up:-1])))))
        G = np.exp(-alpha*si2) # convolution vector for downstream boundary to wrap around 
        
        R1[i] = omega*R0[i] + gamma*np.sum(np.hstack((R0[i::-1], R0[ns-1:pad_up:-1]))*G)/np.sum(G) # main equation, weighted sum of curvatures upstream from downstream boundary - periodic boundary condition
    #####################################################
    for i in range(pad,ns):
        si2 = np.hstack((np.array([0]),np.cumsum(ds[i-1::-1])))  # distance along centerline, backwards from current point 
        G = np.exp(-alpha*si2) # convolution vector
        R1[i] = omega*R0[i] + gamma*np.sum(R0[i::-1]*G)/np.sum(G) # main equation


    return R1

def compute_derivatives(x,y):
    """function for computing first derivatives of a curve (centerline)
    x,y are cartesian coodinates of the curve
    outputs:
    dx - first derivative of x coordinate
    dy - first derivative of y coordinate
    ds - distances between consecutive points along the curve
    s - cumulative distance along the curve"""
    if len(x) < 2:
        dx = [0]; dy = [0]; ds = [0]; s = [0]
    else:
        dx = np.gradient(x) # first derivatives
        dy = np.gradient(y)   
        ds = np.sqrt(dx**2+dy**2)
        s = np.hstack((0,np.cumsum(ds[1:])))
    return dx, dy, ds, s

def compute_curvature(x,y):
    """function for computing first derivatives and curvature of a curve (centerline)
    x,y are cartesian coodinates of the curve
    outputs:
    dx - first derivative of x coordinate
    dy - first derivative of y coordinate
    ds - distances between consecutive points along the curve
    s - cumulative distance along the curve
    curvature - curvature of the curve (in 1/units of x and y)"""
    dx = np.gradient(x) # first derivatives
    dy = np.gradient(y)      
    ddx = np.gradient(dx) # second derivatives 
    ddy = np.gradient(dy) 
    curvature = (dx*ddy-dy*ddx)/((dx**2+dy**2)**1.5)
    return curvature

def kth_diag_indices(a,k):
    """function for finding diagonal indices with k offset
    [from https://stackoverflow.com/questions/10925671/numpy-k-th-diagonal-indices]"""
    rows, cols = np.diag_indices_from(a)
    if k<0:
        return rows[:k], cols[-k:]
    elif k>0:
        return rows[k:], cols[:-k]
    else:
        return rows, cols
    
def find_cutoffs(x,y,crdist,deltas):
    """function for identifying locations of cutoffs along a centerline
    and the indices of the segments that will become part of the oxbows
    from MeanderPy
    x,y - coordinates of centerline
    crdist - critical cutoff distance
    deltas - distance between neighboring points along the centerline"""
    diag_blank_width = int((crdist+20*deltas)/deltas)
    # distance matrix for centerline points:
    dist = distance.cdist(np.array([x,y]).T,np.array([x,y]).T)
    dist[dist>crdist] = np.NaN # set all values that are larger than the cutoff threshold to NaN
    # set matrix to NaN along the diagonal zone:
    for k in range(-diag_blank_width,diag_blank_width+1):
        rows, cols = kth_diag_indices(dist,k)
        dist[rows,cols] = np.NaN
    i1, i2 = np.where(~np.isnan(dist))
    ind1 = i1[np.where(i1<i2)[0]] # get rid of unnecessary indices
    ind2 = i2[np.where(i1<i2)[0]] # get rid of unnecessary indices
    return ind1, ind2 # return indices of cutoff points and cutoff coordinates
def cut_off_cutoffs(x,y,s,crdist,deltas):
    """function for executing cutoffs - removing oxbows from centerline and storing cutoff coordinates
    from MeanderPy
    x,y - coordinates of centerline
    crdist - critical cutoff distance
    deltas - distance between neighboring points along the centerline
    outputs:
    x,y - updated coordinates of centerline
    xc, yc - lists with coordinates of cutoff segments
    cl_dist - distance cutoff occurs down centerline"""
    xc = []
    yc = []
    cl_dist = []
    cut_len = []
    max_curv = []
    ind1, ind2 = find_cutoffs(x,y,crdist,deltas) # initial check for cutoffs
    
    while len(ind1)>0:
        xc.append(x[ind1[0]:ind2[0]+1]) # x coordinates of cutoff
        yc.append(y[ind1[0]:ind2[0]+1]) # y coordinates of cutoff
        dx, dy, ds, s_little = compute_derivatives(x[:ind1[0]+1],y[:ind1[0]+1])#compute derivatives upstream of cutoff
        cl_dist.append(s_little[-1]) #cutoff distance downstream
        #max_curv = np.max(compute_curvature(x[ind1[0]:ind2[0]+1],y[ind1[0]:ind2[0]+1]))  #maximum curvature along cutoff bend
        dx, dy, ds, s_between = compute_derivatives(xc[-1],yc[-1])#compute derivatives along cutoff bend
        cut_len.append(s_between[-1]) #length removed by cutoff
        x = np.hstack((x[:ind1[0]+1],x[ind2[0]:])) # x coordinates after cutoff
        y = np.hstack((y[:ind1[0]+1],y[ind2[0]:])) # y coordinates after cutoff
        
        ind1, ind2 = find_cutoffs(x,y,crdist,deltas) 
    return x,y,xc,yc, cl_dist, cut_len

def get_channel_banks(x,y,W):
    """function for finding coordinates of channel banks, given a centerline and a channel width
    from MeanderPy
    x,y - coordinates of centerline
    W - channel width
    outputs:
    xm, ym - coordinates of channel banks (both left and right banks)"""
    x1 = x.copy()
    y1 = y.copy()
    x2 = x.copy()
    y2 = y.copy()
    ns = len(x)
    dx = np.diff(x); dy = np.diff(y) 
    ds = np.sqrt(dx**2+dy**2)
    x1[:-1] = x[:-1] + 0.5*W*np.diff(y)/ds
    y1[:-1] = y[:-1] - 0.5*W*np.diff(x)/ds
    x2[:-1] = x[:-1] - 0.5*W*np.diff(y)/ds
    y2[:-1] = y[:-1] + 0.5*W*np.diff(x)/ds
    x1[ns-1] = x[ns-1] + 0.5*W*(y[ns-1]-y[ns-2])/ds[ns-2]
    y1[ns-1] = y[ns-1] - 0.5*W*(x[ns-1]-x[ns-2])/ds[ns-2]
    x2[ns-1] = x[ns-1] - 0.5*W*(y[ns-1]-y[ns-2])/ds[ns-2]
    y2[ns-1] = y[ns-1] + 0.5*W*(x[ns-1]-x[ns-2])/ds[ns-2]
    xm = np.hstack((x1,x2[::-1]))
    ym = np.hstack((y1,y2[::-1]))
    return xm, ym
def update_nonlocal_effects(ne, s, decay, scale, cut_dist, cut_len, thresh = .05):
    #reshape array to fit new centerline
    ne_new = np.interp(np.arange(len(s)),np.arange(len(ne)), ne)
   
    ###decay old NE
    ne_new = ne_new*np.exp(-decay)
    ### remove ne that are less than some threshold, default = .05 (1/20 of background rate)
    ne_new[np.where(ne_new<thresh)] = 0

    for k in range(len(cut_dist)): #for each cutoff, add new NE

        #gaussian bump
        mu = cut_dist[k]

        sigma = (cut_len[k]*1.19)/2 # want the whole bump within 1.19*cut_len

        y_bump = norm.pdf(s, mu, sigma)
        ne_new = ne_new + ((scale-1)*y_bump/np.max(y_bump))
  
    return ne_new


       
