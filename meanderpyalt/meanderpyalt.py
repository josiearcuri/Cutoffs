import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import norm
import scipy.interpolate
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.spatial import distance
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
import time, sys
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec


def update_progress(progress):
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
    text = "\rPercent: [{0}] {1}% {2}".format( "#"*block + "-"*(barLength-block), progress*100, status)
    sys.stdout.write(text)
    sys.stdout.flush()

class Channel:
    """class for Channel objects"""
    def __init__(self,x,y,z,W,D):
        """initialize Channel object
        x, y, z  - coordinates of centerline
        W - channel width
        D - channel depth"""
        self.x = x
        self.y = y
        self.z = z
        self.W = W
        self.D = D
        

class Cutoff:
    """class for Cutoff objects"""
    def __init__(self,x,y,z,W,D,dist, time):
        """initialize Cutoff object
        x, y, z  - coordinates of centerline
        W - channel width
        D - channel depth"""
        self.x = x
        self.y = y
        self.z = z
        self.W = W
        self.D = D
        self.dist = np.mean(dist)
        self.time = time
        #print("cutoff at " + str(time))
class ChannelBelt:
    """class for ChannelBelt objects"""
    def __init__(self, channels, cutoffs, cl_times, cutoff_times, cutoff_dists, decay_rate, bump_scale):
        """initialize ChannelBelt object
        channels - list of Channel objects
        cutoffs - list of Cutoff objects
        cl_times - list of ages of Channel objects
        cutoff_times - list of ages of Cutoff objects
        nonlocal_effects - list of magnitudes of additional nominal migration resulting from a cutoff, m by n matrix with m nodes along centerline and n active cutoffs. columns are iteratively added + filled based on cutoff length removal when new cutoffs occur, and removed when max of magnitudes is below threshold.  each dt columns are decayed according to an exponential rate.
        decay_rate - rate at which nonlocal effects dissipate"""
        self.channels = channels
        self.cutoffs = cutoffs
        self.cl_times = cl_times
        self.cutoff_times = cutoff_times
        self.cutoff_dists = cutoff_dists
        self.decay_rate = decay_rate
        self.bump_scale = bump_scale

    def migrate(self,nit,saved_ts,deltas,pad, crdist,Cf,kl,kv,dt,dens,D=[]):
        """function for computing migration rates along channel centerlines and moving the centerlines accordingly
        inputs:
        nit - number of iterations
        saved_ts - which time steps will be saved
        deltas - distance between nodes on centerline
        pad - padding for upstream bc (number of nodepoints along centerline)
        crdist - threshold distance at which cutoffs occur
        Cf - dimensionless Chezy friction factor
        kl - migration rate constant (m/s)
        kv - vertical slope-dependent erosion rate constant (m/s)
        dt - time step (s)
        dens - density of fluid (kg/m3)"""
        channel = self.channels[-1] # first channel is the same as last channel of input
        x = channel.x; y = channel.y; z = channel.z
        W = channel.W; D = channel.D
        
        k = 1.0 # constant in HK equation
        xc = [] # initialize cutoff coordinates
        yc = []
        cut_dist = []
        cut_len = []
        # determine age of last channel:
        if len(self.cl_times)>0:
            last_cl_time = self.cl_times[-1]
        else:
            last_cl_time = 0
        dx, dy, dz, ds, s = compute_derivatives(x,y,z)
        slope = np.gradient(z)/ds
        omega = -1.0 # constant in curvature calculation (Howard and Knutson, 1984)
        gamma = 2.5 # from Ikeda et al., 1981 and Howard and Knutson, 1984
        ne = np.zeros_like(x)
        for itn in range(nit): # main loop
            update_progress(itn/nit) 
            ne = update_nonlocal_effects(ne, s, self.decay_rate, self.bump_scale, cut_dist, cut_len) #update array of ne with last itn's cutoff(s) and decay old ne
            klarray = nominal_rate(kl, ne) ## compute array of nominal migration rate in m/s
            x, y = migrate_one_step(x,y,z,W,klarray,dt,k,Cf,D,pad,omega,gamma)
            x,y,z,xc,yc,zc,cut_dist, cut_len = cut_off_cutoffs(x,y,z,s,crdist,deltas) # find and execute cutoffs
            x,y,z,dx,dy,dz,ds,s = resample_centerline(x,y,z,deltas) # resample centerline
            slope = np.gradient(z)/ds
            
            if len(xc)>0: # save cutoff data
                self.cutoff_times.append(last_cl_time+(itn+1)*dt/(365*24*60*60.0))
                self.cutoff_dists.append(cut_dist)
                cutoff = Cutoff(xc,yc,zc,W,D,cut_dist,last_cl_time+(itn+1)*dt/(365*24*60*60.0)) # create cutoff object
                self.cutoffs.append(cutoff)
            # saving centerlines:
            if np.mod(itn+1,saved_ts)==0:
                self.cl_times.append(last_cl_time+(itn+1)*dt/(365*24*60*60.0))
                channel = Channel(x,y,z,W,D) # create channel object
                self.channels.append(channel)

    def plot(self,plot_type,pb_age,ob_age,*end_time):
        """plot ChannelBelt object
        plot_type - can be either 'strat' (for stratigraphic plot) or 'morph' (for morphologic plot)
        pb_age - age of point bars (in years) at which they get covered by vegetation
        ob_age - age of oxbow lakes (in years) at which they get covered by vegetation
        end_time (optional) - age of last channel to be plotted (in years)"""
        cot = np.array(self.cutoff_times)
        sclt = np.array(self.cl_times)
        if len(end_time)>0:
            cot = cot[cot<=end_time]
            sclt = sclt[sclt<=end_time]

        times = np.sort(np.hstack((cot,sclt)))
        times = np.unique(times)
        order = 0 # variable for ordering objects in plot
        # set up min and max x and y coordinates of the plot:
        xmin = np.min(self.channels[0].x)
        xmax = np.max(self.channels[0].x)
        ymax = 0
        for i in range(len(self.channels)):
            ymax = max(ymax, np.max(np.abs(self.channels[i].y)))
        ymax = ymax+2*self.channels[0].W # add a bit of space on top and bottom
        ymin = -1*ymax
        # size figure so that its size matches the size of the model:
        fig = plt.figure(figsize=(20,(ymax-ymin)*20/(xmax-xmin))) 
        if plot_type == 'morph':
            pb_crit = len(times[times<times[-1]-pb_age])/float(len(times))
            ob_crit = len(times[times<times[-1]-ob_age])/float(len(times))
            green = (106/255.0,159/255.0,67/255.0) # vegetation color
            pb_color = (189/255.0,153/255.0,148/255.0) # point bar color
            ob_color = (15/255.0,58/255.0,65/255.0) # oxbow color
            pb_cmap = make_colormap([green,green,pb_crit,green,pb_color,1.0,pb_color]) # colormap for point bars
            ob_cmap = make_colormap([green,green,ob_crit,green,ob_color,1.0,ob_color]) # colormap for oxbows
            plt.fill([xmin,xmax,xmax,xmin],[ymin,ymin,ymax,ymax],color=(106/255.0,159/255.0,67/255.0))
        for i in range(0,len(times)):
            if times[i] in sclt:
                ind = np.where(sclt==times[i])[0][0]
                x1 = self.channels[ind].x
                y1 = self.channels[ind].y
                W = self.channels[ind].W
                xm, ym = get_channel_banks(x1,y1,W)
                if plot_type == 'morph':
                    if times[i]>times[-1]-pb_age:
                        plt.fill(xm,ym,facecolor=pb_cmap(i/float(len(times)-1)),edgecolor='k',linewidth=0.2)
                    else:
                        plt.fill(xm,ym,facecolor=pb_cmap(i/float(len(times)-1)))
                else:
                    order = order+1
                    plt.fill(xm,ym,sns.xkcd_rgb["light tan"],edgecolor='k',linewidth=0.25,zorder=order)
            if times[i] in cot:
                ind = np.where(cot==times[i])[0][0]
                for j in range(0,len(self.cutoffs[ind].x)):
                    x1 = self.cutoffs[ind].x[j]
                    y1 = self.cutoffs[ind].y[j]
                    xm, ym = get_channel_banks(x1,y1,self.cutoffs[ind].W)
                    if plot_type == 'morph':
                        plt.fill(xm,ym,color=ob_cmap(i/float(len(times)-1)))
                    else:
                        order = order+1
                        plt.fill(xm,ym,sns.xkcd_rgb["ocean blue"],edgecolor='k',linewidth=0.1,zorder=order)
        x1 = self.channels[len(sclt)-1].x
        y1 = self.channels[len(sclt)-1].y
        xm, ym = get_channel_banks(x1,y1,self.channels[len(sclt)-1].W)
        order = order+1
        plt.fill(xm,ym,color=(16/255.0,73/255.0,90/255.0),zorder=order,edgecolor='k')
        plt.axis('equal')

        return fig

def resample_centerline(x,y,z,deltas):
    dx, dy, dz, ds, s = compute_derivatives(x,y,z) # compute derivatives
    # resample centerline so that 'deltas' is roughly constant
    # [parametric spline representation of curve; note that there is *no* smoothing]
    tck, u = scipy.interpolate.splprep([x,y,z],s=0) 
    unew = np.linspace(0,1,1+int(s[-1]/deltas)) # vector for resampling
    out = scipy.interpolate.splev(unew,tck) # resampling
    x, y, z = out[0], out[1], out[2] # assign new coordinate values
    dx, dy, dz, ds, s = compute_derivatives(x,y,z) # recompute derivatives
    return x,y,z,dx,dy,dz,ds,s

def nominal_rate(kl, ne):
    new_kl = kl*(1+ne)
    return new_kl
def migrate_one_step(x,y,z,W,klarray,dt,k,Cf,D,pad,omega,gamma):
    ns=len(x)
    curv = compute_curvature(x,y)
    dx, dy, dz, ds, s = compute_derivatives(x,y,z)
    curv = W*curv # dimensionless curvature
    R0 = klarray*curv
    alpha = k*2*Cf/D # exponent for convolution function G
    R1 = compute_migration_rate(pad,ns,ds,alpha,omega,gamma,R0)
    # calculate new centerline coordinates:
    dy_ds = dy/ds
    dx_ds = dx/ds
    # adjust x and y coordinates (this *is* the migration):
    x = x + R1*dy_ds*dt  
    y = y - R1*dx_ds*dt 
    return x,y

def generate_initial_channel(W,D,Sl,deltas,pad,n_bends):
    """generate straight Channel object with some noise added that can serve
    as input for initializing a ChannelBelt object
    W - channel width
    D - channel depth
    Sl - channel gradient
    deltas - distance between nodes on centerline
    pad - padding (number of nodepoints along centerline)
    n_bends - approximate number of bends to be simulated"""
    noisy_len = n_bends*10*W/2.0 # length of noisy part of initial centerline
    pad1 = int(pad/10.0) # padding at upstream end can be shorter than padding on downstream end
    if pad1<5:
        pad1 = 5
    x = np.linspace(0, noisy_len+(pad+pad1)*deltas, int(noisy_len/deltas+pad+pad1)+1) # x coordinate
    y = 10.0 * (2*np.random.random_sample(int(noisy_len/deltas)+1,)-1)
    y = np.hstack((np.zeros((pad1),),y,np.zeros((pad),))) # y coordinate
    deltaz = Sl * deltas*(len(x)-1)
    z = np.linspace(0,deltaz,len(x))[::-1] # z coordinate
    return Channel(x,y,z,W,D)


def generate_channel_from_file(filelist, slope = .01, D_in= 10, smooth_factor=.25, matlab_corr= -1):
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
        z - array of z coordinates
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
    z = np.interp(np.asarray(range(len(points_fitted[0]))), [1, len(points_fitted[0])], [(slope*distance[-1]), 0]) 
    deltas = round(distance[-1]/(len(points_fitted[0])-1)) 
    return [Channel(points_fitted[0],points_fitted[1],z,W,D), x, y, z, distance[-1], deltas]

def compute_migration_rate(pad,ns,ds,alpha,omega,gamma,R0):
    """compute migration rate as weighted sum of upstream curvatures
    pad - padding (number of nodepoints along centerline)
    ns - number of points in centerline
    ds - distances between points in centerline
    omega - constant in HK model
    gamma - constant in HK model
    R0 - nominal migration rate (dimensionless curvature * migration rate constant)"""
    R1 = np.zeros(ns) # preallocate adjusted channel migration rate

    for i in range(pad,ns):
        si2 = np.hstack((np.array([0]),np.cumsum(ds[i-1::-1])))  # distance along centerline, backwards from current point 
        G = np.exp(-alpha*si2) # convolution vector
        R1[i] = omega*R0[i] + gamma*np.sum(R0[i::-1]*G)/np.sum(G) # main equation

#########Periodic Boundary#########################
    for i in range(0,pad):
        si2 = np.hstack((np.array([0]),np.cumsum(np.hstack((ds[i-1::-1], ds[ns::ns-(10*pad)])))))  # distance along centerline, backwards from corresponding point on downstream boundary
        G = np.exp(-alpha*si2) # convolution vector for downstream boundary to wrap around 
        R1[i] = omega*R0[i] + gamma*np.sum(np.hstack((R0[i::-1], R0[ns::ns-(10*pad)]))*G)/np.sum(G) # main equation, weighted sum of curvatures upstream from downstream boundary - periodic boundary condition
    return R1

def compute_derivatives(x,y,z):
    """function for computing first derivatives of a curve (centerline)
    x,y are cartesian coodinates of the curve
    outputs:
    dx - first derivative of x coordinate
    dy - first derivative of y coordinate
    ds - distances between consecutive points along the curve
    s - cumulative distance along the curve"""
    if len(x) < 2:
        dx = [0]; dy = [0]; dz = [0]; ds = [0]; s = [0]
    else:
        dx = np.gradient(x) # first derivatives
        dy = np.gradient(y)   
        dz = np.gradient(z)   
        ds = np.sqrt(dx**2+dy**2+dz**2)
        s = np.hstack((0,np.cumsum(ds[1:])))
    return dx, dy, dz, ds, s

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

def make_colormap(seq):
    """Return a LinearSegmentedColormap
    seq: a sequence of floats and RGB-tuples. The floats should be increasing
    and in the interval (0,1).
    [from: https://stackoverflow.com/questions/16834861/create-own-colormap-using-matplotlib-and-plot-color-scale]
    """
    seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
    cdict = {'red': [], 'green': [], 'blue': []}
    for i, item in enumerate(seq):
        if isinstance(item, float):
            r1, g1, b1 = seq[i - 1]
            r2, g2, b2 = seq[i + 1]
            cdict['red'].append([item, r1, r2])
            cdict['green'].append([item, g1, g2])
            cdict['blue'].append([item, b1, b2])
    return mcolors.LinearSegmentedColormap('CustomMap', cdict)

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
def cut_off_cutoffs(x,y,z,s,crdist,deltas):
    """function for executing cutoffs - removing oxbows from centerline and storing cutoff coordinates
    x,y - coordinates of centerline
    crdist - critical cutoff distance
    deltas - distance between neighboring points along the centerline
    outputs:
    x,y,z - updated coordinates of centerline
    xc, yc, zc - lists with coordinates of cutoff segments
    cl_dist - distance cutoff occurs down centerline"""
    xc = []
    yc = []
    zc = []
    cl_dist = []
    cut_len = []
    max_curv = []
    ind1, ind2 = find_cutoffs(x,y,crdist,deltas) # initial check for cutoffs
    
    while len(ind1)>0:
        xc.append(x[ind1[0]:ind2[0]+1]) # x coordinates of cutoff
        yc.append(y[ind1[0]:ind2[0]+1]) # y coordinates of cutoff
        zc.append(z[ind1[0]:ind2[0]+1]) # z coordinates of cutoff
        #########JOSIE ADDITIONS###############
        dx, dy, dz, ds, s_little = compute_derivatives(x[:ind1[0]+1],y[:ind1[0]+1],z[:ind1[0]+1])#compute derivatives upstream of cutoff
        cl_dist.append(s_little[-1]/s[-1]) #cutoff distance downstream
        #max_curv = np.max(compute_curvature(x[ind1[0]:ind2[0]+1],y[ind1[0]:ind2[0]+1]))  #maximum curvature along cutoff bend
        dx, dy, dz, ds, s_between = compute_derivatives(xc[-1],yc[-1],zc[-1])#compute derivatives along cutoff bend
        cut_len.append(s_between[-1]) #length removed by cutoff
        ########################
        x = np.hstack((x[:ind1[0]+1],x[ind2[0]:])) # x coordinates after cutoff
        y = np.hstack((y[:ind1[0]+1],y[ind2[0]:])) # y coordinates after cutoff
        z = np.hstack((z[:ind1[0]+1],z[ind2[0]:])) # z coordinates after cutoff
        
        ind1, ind2 = find_cutoffs(x,y,crdist,deltas) 
    return x,y,z,xc,yc,zc, cl_dist, cut_len

def get_channel_banks(x,y,W):
    """function for finding coordinates of channel banks, given a centerline and a channel width
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
def update_nonlocal_effects(ne, s, decay, scale, cut_dist, cut_len, thresh = .0001):
    #reshape array to fit new centerline
    ne_new = np.interp(np.arange(len(s)),np.arange(len(ne)), ne)
   
    ###decay old NE
    ne_new = ne_new*np.exp(-decay)
    ### remove ne that are less than some threshold, default = .001 (1/10th of a percent) 
    ne_new[np.where(ne_new<thresh)] = 0

    for k in range(len(cut_dist)): #for each cutoff, add new NE
        #get distances of upstream and downstream extent to add NE
        US = s[-1]*cut_dist[k] - cut_len[k]*1.19
        DS = s[-1]*cut_dist[k] + cut_len[k]*1.19
        if US<0:
            US = 0
        if DS>s[-1]:
            DS = s[-1]
        #get indices corresponding to these distances
        idx_us = np.where(s<=US)[0][-1]
        idx_ds = np.where(s>=DS)[0][0]

        #gaussian bump
        mu = s[-1]*cut_dist[k]

        sigma = (cut_len[k]*1.19)/4 # want the whole bump within 1.19*cut_len

        y_bump = norm.pdf(s, mu, sigma)
        ne_new = ne_new + (scale*y_bump/np.max(y_bump))
    return ne_new
