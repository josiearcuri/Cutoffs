import numpy as np
import pandas as pd
import time, sys
import numba 
import math
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
from matplotlib import cm
import matplotlib as mpl

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
    def __init__(self,x,y,W,D,MR):
        """initialize Channel object
        x, y, z  - coordinates of centerline
        W - channel width
        D - channel depth
        MR - migration rate from last centerline"""
        self.x = x
        self.y = y
        self.W = W
        self.D = D
        self.MR = MR

        

class Cutoff:
    """class for Cutoff objects"""
    def __init__(self,x,y,W,dist, time, cut_len, rad, cllen):
        """initialize Cutoff object
        x, y,  - coordinates of centerline
        W - channel width
        dist - distance along centerline of upstream most node of cutoff bend. 
        time - year of occurrence
        cut_len - ditance removed from centerline by cutoff
        radius - radius of curvature 
        cllen - centerline length after cutoff
        """
        self.x = x
        self.y = y
        self.W = W
        self.dist = dist
        self.time = time
        self.cut_len = cut_len
        self.radius = rad
        self.cllen = cllen
        
class ChannelBelt:
    """class for ChannelBelt objects"""
    def __init__(self, channels, cutoffs, cl_times, cutoff_times, cutoff_dists, decay_rate, bump_scale, cut_thresh, sinuosity):
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
        self.sinuosity = sinuosity

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
        MR = 0
        DS = 0 
        for itn in range(nit): # main loop
            update_progress(itn/nit, start_time)
            ne = update_nonlocal_effects(ne, s, self.decay_rate, self.bump_scale, cut_dist, cut_len) #update array of ne with last itn's cutoff(s) and decay old ne
            klarray = nominal_rate(kl, ne)## compute array of nominal migration rate in m/s with nonlocal effects accounted for
            curv = compute_curvature(x,y)#get curvature of bends before cutoffs happen
            x, y, R1 = migrate_one_step(x,y,W,klarray,dt,k,Cf,D,pad,omega,gamma)
            x,y,xc,yc,cut_dist, cut_len, ind1, ind2 = cut_off_cutoffs(x,y,s,crdist,deltas) # find and execute cutoffs
            x,y,dx,dy,ds,s = resample_centerline(x,y,deltas) # resample centerline
            
            Sin = get_sinuosity(x,y,s)
            if len(xc)>0: # save cutoff data
                rad = get_radii(curv, ind1, ind2, W)
                cutoff = Cutoff(xc,yc,W,cut_dist[0],last_cl_time+(itn)*dt/(365*24*60*60.0), cut_len[0], rad, s[-1]) # create cutoff object
                #keep track of year cutoff occurs, where it occurs, and save an object. 
                self.cutoff_times.append(last_cl_time+(itn)*dt/(365*24*60*60.0))
                self.cutoff_dists.append(cut_dist)
                self.cutoffs.append(cutoff)
            # saving centerlines:
            if np.mod(itn,saved_ts)==0 or itn == nit-1:
                channel = Channel(x,y,W,D, MR) # create channel object, save year
                self.cl_times.append(last_cl_time+(itn)*dt/(365*24*60*60.0))
                self.channels.append(channel)
                self.sinuosity.append(Sin)

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
        
        while len(self.cutoffs)<self.cut_thresh: # main loop
            itn = itn+1
            update_progress(len(self.cutoffs)/self.cut_thresh, start_time) 
            ne = update_nonlocal_effects(ne, s, self.decay_rate, self.bump_scale, cut_dist, cut_len) #update array of ne with last itn's cutoff(s) and decay old ne
            curv = compute_curvature(x,y)
            klarray = nominal_rate(kl, ne)## compute array of nominal migration rate in m/s with nonlocal effects accounted for
            x, y, MR  = migrate_one_step(x,y,W,klarray,dt,k,Cf,D,pad,omega,gamma)
          
            x,y,xc,yc,cut_dist, cut_len,ind1, ind2 = cut_off_cutoffs(x,y,s,crdist,deltas) # find and execute cutoffs
            x,y,dx,dy,ds,s = resample_centerline(x,y,deltas) # resample centerline
            Sin = get_sinuosity(x,y,s)
            if len(xc)>0: # save cutoff data
                rad = get_radii(curv, ind1, ind2, W)
                cutoff = Cutoff(xc,yc,W,cut_dist,last_cl_time+(itn)*dt/(365*24*60*60.0), cut_len, rad, s[-1]) # create cutoff object
                #keep track of year cutoff occurs, where it occurs, and save an object. 
                self.cutoff_times.append(last_cl_time+(itn)*dt/(365*24*60*60.0))
                self.cutoff_dists.append(cut_dist)
                self.cutoffs.append(cutoff)

            # saving centerlines:
            if np.mod(itn,saved_ts)==0 or len(self.cutoffs)>=self.cut_thresh:
                channel = Channel(x,y,W,D,MR) # create channel object, save year
                self.cl_times.append(last_cl_time+(itn)*dt/(365*24*60*60.0))
                self.channels.append(channel)
                self.sinuosity.append(Sin)
    
    def migrate_bendtracking(self,saved_ts,deltas,pad, crdist,Cf,kl,dt,dens=1000):
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
        
        while len(self.cutoffs)<self.cut_thresh: # main loop
            itn = itn+1
            update_progress(len(self.cutoffs)/self.cut_thresh, start_time) 
            ne = update_nonlocal_effects(ne, s, self.decay_rate, self.bump_scale, cut_dist, cut_len) #update array of ne with last itn's cutoff(s) and decay old ne
            curv = compute_curvature(x,y)
            klarray = nominal_rate(kl, ne)## compute array of nominal migration rate in m/s with nonlocal effects accounted for
            x, y, MR  = migrate_one_step(x,y,W,klarray,dt,k,Cf,D,pad,omega,gamma)
        
            x,y,xc,yc,cut_dist, cut_len,ind1, ind2 = cut_off_cutoffs(x,y,s,crdist,deltas) # find and execute cutoffs
            x,y,dx,dy,ds,s = resample_centerline(x,y,deltas) # resample centerline
            
            if len(xc)>0: # save cutoff data
                rad = get_radii(curv, ind1, ind2, W)
                cutoff = Cutoff(xc,yc,W,cut_dist,last_cl_time+(itn)*dt/(365*24*60*60.0), cut_len, rad, s[-1]) # create cutoff object
                #keep track of year cutoff occurs, where it occurs, and save an object. 
                self.cutoff_times.append(last_cl_time+(itn)*dt/(365*24*60*60.0))
                self.cutoff_dists.append(cut_dist)
                self.cutoffs.append(cutoff)
            # saving centerlines:
            if np.mod(itn,saved_ts)==0 or len(self.cutoffs)>=self.cut_thresh:
        
                channel = Channel(x,y,W,D,MR) # create channel object, save year
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
    
    def plot_channels2(self,ncutoff, directory):
        cot = np.array(self.cutoff_times)
        sclt = np.array(self.cl_times)
      
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

        cmap = cm.get_cmap('gray_r',self.bump_scale)

        order = 0
        figs=[]
    
        cmap = LinearSegmentedColormap.from_list('custom red', ['#ffffff','#de2d26'], N=256)
        ind = np.where(sclt>=cot[ncutoff])[0][0]
        times = np.arange(0, len(sclt))[(ind-1):1+(ind+(4*4))]
        fig, axs = plt.subplots(6,1, figsize = (2.5, 5), constrained_layout = True)

        cx1 = self.cutoffs[ncutoff].x[0]
        cy1 = self.cutoffs[ncutoff].y[0]
        dist = self.cutoffs[ncutoff].cut_len
        loc =  self.cutoffs[ncutoff].dist
     
        bump = self.bump_scale
        decay = self.decay_rate
        steps = np.arange(1, 18, 4)   
        ne = np.zeros_like(self.channels[times[0]].x)
        s = np.zeros_like(self.channels[times[0]].x)
        for t in range(0,len(times)):
                

            
            W = self.channels[-1].W
            x1 = self.channels[times[t]].x
            y1 = self.channels[times[t]].y
            dx, dy, ds, s = compute_derivatives(x1,y1)
            ne = update_nonlocal_effects(ne, s, decay, bump, [], [])    
            xm, ym = get_channel_banks(x1,y1,W)
            order += 1
            central_node = np.where(x1>=cx1[0])[0][0]      
            upstream_node = np.where(s>=s[central_node]-dist)[0][0]
            downstream_node = np.where(s>=s[central_node]+dist)[0][0]
            ne_x = x1[upstream_node-1:downstream_node+1]
            ne_y = y1[upstream_node-1:downstream_node+1]
            ne_mag = ne[upstream_node-1:downstream_node+1] 
            if t==0:
                
         
                ne = update_nonlocal_effects(ne, s, decay, bump,np.array(loc), np.array(dist))
                central_node = np.where(x1>=cx1[0])[0][0]      
                upstream_node = np.where(s>=s[central_node]-dist)[0][0]
                downstream_node = np.where(s>=s[central_node]+dist)[0][0]
                ne_x = x1[upstream_node-1:downstream_node+1]
                ne_y = y1[upstream_node-1:downstream_node+1]
                ne_mag = ne[upstream_node-1:downstream_node+1]  
                axs[0].set_xlim([cx1[0]-dist[0]/1.85,cx1[0]+dist[0]/1.85])
                axs[0].fill(xm,ym,facecolor='k',edgecolor = 'k', alpha = 1, linewidth=0.1,zorder=order)
                #axs[i].fill(xm[idx,y1,facecolor='r',edgecolor = 'none', alpha = .5,zorder=order, label= '_nolegend_')
                #plt.plot(chx1[central], chy1[chx1<=x1[0]], color = 'k', linewidth = 1, zorder = order)
                axs[0].annotate(text ='year '+str(int(t)) ,  xy = (.9,.9), xycoords='axes fraction', fontsize = 10)
                axs[0].axis('off')
                

                
            if t  in steps:

                nextind = np.where(steps==t)[0][0]+1
                
                axs[nextind].fill(xm,ym,facecolor='k',edgecolor = 'k', alpha = 1, linewidth=0.1,zorder=order, label= 'year '+str(2+int((t-1)/4)))
                #axs[i].fill(xm[idx,y1,facecolor='r',edgecolor = 'none', alpha = .5,zorder=order, label= '_nolegend_')
                #plt.plot(chx1[central], chy1[chx1<=x1[0]], color = 'k', linewidth = 1, zorder = order)
                axs[nextind].annotate(text ='year '+str(1+int((t-1)/4)) ,  xy = (.9,.9), xycoords='axes fraction', fontsize = 10)
                axs[nextind].axis('off')
                nemap = axs[nextind].scatter(ne_x[ne_mag>0], ne_y[ne_mag>0], s= 40, c = ne_mag[ne_mag>0], vmin = 0, vmax = bump, cmap = cmap)
                axs[nextind].set_xlim([cx1[0]-dist[0]/1.85,cx1[0]+dist[0]/1.85])
            
                
                
                
  
            #axs[t].set_ylim([ymin,ymax])

            
      
        norm = mpl.colors.Normalize(vmin=0, vmax=bump)
        fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),orientation='horixontal', ax = axs, label = 'nonlocal effect [*]',  shrink = 0.6)
        #axs[-1].set_xlim([cx1[0]-dist[0]/4,cx1[0]+dist[0]/4])
        
    
        return fig
    def cutoff_distributions(self, year, filepath, mode):
        """pull cutoff data from channel belt object and export csv, return dataframe for plotting
        year - last centerline year
        filepath - where cutoff info csv and distrubution plot are saved
        mode - for plotting - "OnlyCurvature or "NonlocalEffects"
        """
        #pull cutoff locations, downstrem distance, time from channel belt as a pandas dataframe
        distances = [i.dist[0] for i in self.cutoffs]
        times = [i.time for i in self.cutoffs]
        radius = [i.radius for i in self.cutoffs]
        cutlen = [i.cut_len[0] for i in self.cutoffs]
        cllen = [i.cllen for i in self.cutoffs]
        cuts = pd.DataFrame({'downstream_distance': distances, 'time': times, 'radius':radius, 'cutlen': cutlen, 'cllen': cllen })
        
        #save distribution to csv
        newcuts = cuts.to_csv(filepath+mode+"_"+str(len(cuts['time']))+"_cutoffs_distribution.csv", index_label = "Cutoff")
        plot_cuts(cuts,self.channels[-1].W, filepath)
        return cuts
    def MR_time(self, filepath):
        MR = [[bend for bend in i.MR] for i in self.channels[1:]]
    
        clt = np.array(self.cl_times[1:])
        
        MRdf= pd.DataFrame(MR).dropna(axis=1, how = 'all').dropna(axis=0, how = 'all')

        MRdf.to_csv(filepath, index=False, header=False)
        MRdf = pd.read_csv(filepath, sep = ',', header=None).to_numpy()
        fig, axs = plt.subplots(1,1, figsize=(6,8))
        plot_segmented_MR(MRdf,np.nanmean(MRdf), fig, axs, 'k', filepath[-16:-18])
        
       

        
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
    
    # calculate bend-by-bend migration rate
    MR = segmented_MR(curv, R1, s)
    # calculate new centerline coordinates:
    dy_ds = dy/ds
    dx_ds = dx/ds
    # move x and y coordinates:
    x = x + R1*dy_ds*dt  
    y = y - R1*dx_ds*dt 
    return x,y, MR

def generate_initial_channel(W,D,deltas,pad, Length):
    """generate straight Channel object with some noise added that can serve
    as input for initializing a ChannelBelt object
    from MeanderPy
    W - channel width
    D - channel depth
    deltas - distance between nodes on centerline
    pad - padding (number of nodepoints along centerline)"""
    cl_length = Length#((50*W)*10)# length of noisy part of initial centerline
    x = np.linspace(0, cl_length+(2*pad)*deltas, int(cl_length/deltas+(2*pad))+1) # x coordinate
    y = 10.0 * (2*np.random.random_sample(int(cl_length/deltas)+1,)-1)
    y = np.hstack((np.zeros((pad),),y,np.zeros((pad),))) # y coordinate
    MR = np.zeros_like(x)
   
    return Channel(x,y,W,D, MR)

def load_initial_channel(filepath, W, D, deltas):
    """generate initial channel from centerline csv that can serve
    as input for initializing a ChannelBelt object.  must be fine enough resolution to not warrant smoothing
    filepath - csv with x, y coordinates and no headers
    W - channel width
    D - channel depth
    deltas - distance between nodes on centerline"""
    df = pd.read_csv(filepath, sep = ',', header=None).values
    x = df[:,0]
    y = df[:,1]
    MR = np.zeros(int(len(x)/30))
    
    return Channel(x,y,W,D,MR)
def generate_channel_from_file(filelist, deltas, matlab_corr= -1):
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

    ## water depth scalar#
    D = np.exp(np.log(W/18.8)/1.41)  
    # Linear length along the line, add a zero for first point:
    points = np.vstack([x, y]).T
    distance = np.cumsum( np.sqrt(np.sum( np.diff(points, axis=0)**2, axis=1 )) )
    distance = np.insert(distance, 0, 0)

    # Build a list of the spline function, one for each dimension:
    splines = [InterpolatedUnivariateSpline(distance, coords) for coords in points.T]

    # Compute the spline for the smoothed(sampled) distances:
    points_fitted = np.vstack([spl(np.linspace(0, distance[-1],int(distance[-1]/(W//2)))) for spl in splines])
    
    ## z-dim array, interpolated with constant slope along points of centerline.  assumes centerline points are equidistantly placed along original centerline. 
    #deltas = round(distance[-1]/(len(points_fitted[0])-1)) 
    MR = np.zeros_like(points_fitted[0])
  
    return Channel(points_fitted[0],points_fitted[1],W,D, MR), D, W

@numba.jit(nopython=True)
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
    from meanderpy
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
    ind1_save = []
    ind2_save = []
    if len(ind1)>0:
        ind1_save = ind1[0]
        ind2_save = ind2[0]
    while len(ind1)>0:
        xc.append(x[ind1[0]:ind2[0]+1]) # x coordinates of cutoff
        yc.append(y[ind1[0]:ind2[0]+1]) # y coordinates of cutoff
        dx, dy, ds, s_little = compute_derivatives(x[:ind1[0]+1],y[:ind1[0]+1])#compute derivatives upstream of cutoff
        cl_dist.append(s_little[-1]) #cutoff distance downstream
        dx, dy, ds, s_between = compute_derivatives(xc[-1],yc[-1])#compute derivatives along cutoff bend
        cut_len.append(s_between[-1]) #length removed by cutoff
    
        x = np.hstack((x[:ind1[0]+1],x[ind2[0]:])) # x coordinates after cutoff
        y = np.hstack((y[:ind1[0]+1],y[ind2[0]:])) # y coordinates after cutoff
        ind1, ind2 = find_cutoffs(x,y,crdist,deltas) 
    return x,y,xc,yc, cl_dist, cut_len, ind1_save, ind2_save

def get_channel_banks(x,y,W):
    """function for finding coordinates of channel banks, given a centerline and a channel width
    from meanderpy
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
    ne = np.interp(np.arange(len(s)),np.arange(len(ne)), ne)
    ###decay old NE
    ne_new = ne*np.exp(-decay)
    ### remove ne that are less than some threshold, default = .05 (1/20 of background rate)
    ne_new[np.where(ne_new<(thresh))[0]] = 0
    for k in range(len(cut_dist)): #for each cutoff, add new NE

        #gaussian bump
        mu = cut_dist[k]
        print(cut_len[k])
        sigma = (cut_len[k]*1.19)/2 # want 88% of the bump within 1.19*cut_len, rest is negligible

        y_bump = norm.pdf(s, mu, sigma)
        ne_new = ne_new + ((scale)*y_bump/np.max(y_bump))
  
    return ne_new

def get_sinuosity(x,y,s):
    dist = np.hypot(abs(x[-1]-x[0]), abs(y[-1]-y[0]))
    Sin = s[-1]/dist
    return Sin

def plot_sinuosity(time, sin):
    fig, ax = plt.subplots(1,1)
    ax.plot(range(0,len(sin)), sin, 'k', label = "mean S = "+ str(np.mean(sin)))
    plt.legend()
    ax.set_title("Sinuosity")
    return fig
def get_radii(c, ind1, ind2, W):
    #unsigned curvature
    radii = abs(c)
    # Width divided by maximum dimensionless curvature of cutoff bend
    max_rad = W/(np.max(radii[ind1:ind2])*W)
               
    return max_rad

def plot_distribution(cuts,W, filepath):
    x = cuts['downstream_distance']/W
    y = cuts['time']- min(cuts['time'])

    # definitions for the axes
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    spacing = 0.005
    
    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom + height + spacing, width, 0.2]
    rect_histy = [left + width + spacing, bottom, 0.2, height]

    # start with a square Figure
    fig, ax = plt.subplots(1,1,figsize=(3, 3))
    plt.rcParams.update({'font.size': 10})
    #ax = fig.add_axes([0,0,1,1])
    #ax_histx = fig.add_axes(rect_histx, sharex=ax)
    #ax_histy = fig.add_axes(rect_histy, sharey=ax)
    ax.scatter(x, y, c = 'black', s = 10, alpha = .8)
    # use the previously defined function
   # scatter_hist(x, y, ax, ax_histx, ax_histy)
    ax.set_ylabel("time (years)")
    #plt.xlim(left=0)
    #plt.ylim(bottom=0)
    ax.set_xlabel("distance downstream (ch-w)")
    ax.set_xlim([0,1500])
    ax.set_ylim([0,250])
    return fig
def scatter_hist(x, y, ax, ax_histx, ax_histy):
    # no labels
    ax_histx.tick_params(axis="x", labelbottom=False)
    ax_histy.tick_params(axis="y", labelleft=False)

    # the scatter plot:
    ax.scatter(x, y, c = 'black', s = 2)

    # now determine nice limits by hand:
    xbinwidth = 25
    ybinwidth = 10
    xmax = np.max(np.abs(x))
    ymax = np.max(np.abs(y))
    xlim = (int(xmax/xbinwidth) + 1) * xbinwidth
    ylim = (int(ymax/ybinwidth) + 1) * ybinwidth

    xbins = np.arange(0, xlim + xbinwidth, xbinwidth)
    ybins = np.arange(0, ylim + ybinwidth, ybinwidth)
    ax_histx.hist(x, bins=xbins)
    ax_histy.hist(y, bins=ybins, orientation='horizontal')
def segmented_MR(curv, R1, s, n=90):
    """
    approximate a bend-by-bend nth percentile lateral migration rate. 
    
    Inputs:
    curve: array of curvature for every node along the centerline.
    R1: array of already computed migration distances for every node along the centerline.
    ds: array of cumulative distance downstream between nodes.
    n: percentile
    
    Output:
    MR: array of nth percentile migration rate for each segment.  
    upstream: distance downstream of each segment start
    downstream: distance downstream of each segment end
    """
    R1 = np.array(R1)*365*24*60*60.0 #m/yr instead of m/sec
    infs1 = np.where(np.logical_and(curv[1:] > 0, curv[:-1] < 0))[0] + 1
    infs2 = np.where(np.logical_and(curv[1:] < 0, curv[:-1] > 0))[0] + 1
    idx = np.sort(np.concatenate((infs1, infs2)))
    #nodes = np.array([(curv[i-2]*curv[i])< 0 and (curv[i+2]*curv[i])> 0 for i in range(2,len(curv)-2)])
    #idx = np.where(nodes==1)[0]
    idx = [idx[i] for i in range(len(idx)-1) if idx[i+1]-idx[i] >5]
    
     
    
    MR = [np.percentile(np.abs(R1[idx[i]:idx[i+1]]), 100) for i in range(len(idx)-1)]

    return MR
def moving_average(matrix, window):
    """
    averages migration rates on each bend over a set time window upwind. 
    
    inputs
    matrix: n years by m bend migration rate numpy array
    window: how many indices to average over, as years in the past
    
    output
    mid: moving average migration rate over window years for each bend
    """

    years = len(matrix)
    bends = np.min([np.count_nonzero(~np.isnan(matrix[i])) for i in range(years)])

    mid = np.zeros(shape = (years, bends))
    mid[0,:bends] = matrix[0][:bends]
    for year in range(1,years):
        if year < window:
            mid[year, :] = np.nanmean(matrix[:year][:bends])
        else:
            mid[year, :] = np.nanmean(matrix[year-window:(year+1)][:bends])
     
    return mid
    
def plot_segmented_MR(MR,obmr, fig, axs, c, name, dt = 1):
   
    mean_all = round(np.nanmean(MR[:, 2:-2]),3)
    max_all = round(np.nanmax(MR[:, 2:-2]), 3)
    
    
    #heatmap = ax1.imshow(MR, cmap = 'gist_heat')
    #heatmap = axs[1].imshow(MR/100, cmap = 'cividis', vmin = 0, vmax =np.nanmax(MR/100), aspect = 'auto', origin = 'lower')
    #axs[1].set_xlabel('distance downstream (bend #)')
  
    #cb = fig.colorbar(heatmap, ax=axs[1])
    #cb.set_label("maximum migration rate (ch-w/yr)")
    
    #mid = moving_average(MR, 2)/100
    #min
    #axs[0].plot(MR/100, range(len(MR)), alpha=.2, c = 'k')
    #mean
   # MR_ddt = np.nanpercentile((MR[1:, :] - MR[:-1, :])/MR[1:, :], 50, axis = 1) 
    #f_spec, Pxx_mid = periodogram(np.nanpercentile(MR[:, 10:], 50, axis = 1),dt, scaling = 'spectrum')
    
    
    
    #f_spec, Pxx_mid = welch(np.nanpercentile(MR[:, :], 50, axis = 1), dt, nperseg=500, scaling = 'spectrum')
    #f, t, Zxx = stft(np.nanpercentile(MR[:, :], 50, axis = 1), nperseg=100)
    #line =axs.pcolormesh(t, f, np.abs(Zxx), vmin=0,vmax=int(np.nanpercentile(Zxx, 90)), shading='gouraud')
    
    #filter
    #Pxx_mid[1:-1]= np.convolve(Pxx_mid, np.ones(3)/3, mode='valid')
    #cumu = np.cumsum(Pxx_mid)
    #midfreq = f_spec[cumu>=(cumu[-1]/2)][0]
    
    #freqscale = (f_spec[1:])/(500/len(MR[:, 0]))
    
    #Pxx_mid = savgol_filter(Pxx_mid,5,3)
   # Pxx_den = medfilt(Pxx_den)
    #upper = np.percentile(Pxx_mid, 99)
    #upper_noise = np.percentile(Pxx_mid[int(len(Pxx_mid)/2):], 95)

    #NP = np.mean(Pxx_mid[-int((len(Pxx_mid)/4)):])
    #peaks, _ = find_peaks(Pxx_mid)
    #[proms, leftbases, rightbases] =  peak_prominences(Pxx_mid, peaks)
    #mids = ((rightbases-leftbases)/2) +leftbases
 
    #newmids = mids[proms>=np.percentile(proms, 90)]
    #print(newmids)
    #input()
    #proms = proms[proms>=np.percentile(proms, 90)]
  
    
    #idx = proms.argsort()[:]
    #proms = proms[idx]
    #newmids = newmids[idx]
    #lower_noise = np.percentile(Pxx_mid[-int(len(Pxx_mid)/5):], 1)
    
    #amp = Pxx_mid.max()
    #stn = amp/NP
    #maxfreq = f_spec[Pxx_mid == Pxx_mid.max()]
    line = axs.hist(MR.flatten(), range=(0,500), bins= range(0, 50),color= c, linewidth = 1.5, histtype='step')
    #line = axs.semilogy(f[1:], NP/(np.power(f, 2))[1:],alpha=.8, linewidth = 1, ls = "--", c = 'brown', label = '_nolegend_')
    #line = axs.semilogy(f_spec[1:], 1/(np.power(f_spec, .8))[1:],alpha=.8, linewidth = 1, ls = "--", c = 'pink', label = '_nolegend_')
    #line = axs.semilogy(f[1:], NP+(f[1:]*0),alpha=.8, linewidth = 1, ls = "--", c = 'purple', label = '_nolegend_')
    #line = axs.semilogy(f_spec[1:]/(500/len(MR[:, 0])), 1/(np.ones(len(f_spec)-1)*(upper_noise/NP)),alpha=.5, linewidth = .5, ls = "--", c = 'grey', label = '_nolegend_')
    #line = axs.semilogy(f_spec[1:]/(500/len(MR[:, 0])), 1/(np.ones(len(f_spec)-1)*(lower_noise/NP)),alpha=.5, linewidth = .5, ls = "--", c = 'grey', label = '_nolegend_')
    #line = axs.loglog((f_spec[1:]), 1/(f_spec[1:]**2),alpha=.8, linewidth = 1, ls = "--", c = 'brown', label = '_nolegend_')
    #line = axs.loglog((f_spec[1:]), 1/(f_spec[1:]**1),alpha=.8, linewidth = 1, ls = "--", c = 'pink', label = '_nolegend_')
    #line = axs.loglog((f_spec[1:]), (np.ones(len(f_spec)-1)),alpha=.8, linewidth = 1, ls = "--", c = 'grey', label = '_nolegend_')
    
    #line = axs.semilogy((f_spec[1:]), Pxx_mid[1:]/NP,alpha=.8, linewidth = 1, c = c)#,  label = name + " SNR = "+ str(round(stn, 3)))
    #axs.scatter(newmids, np.ones(len(newmids)),marker = 'o', color = c,alpha = .8, label = '_nolegend_')
    #axs.set_xlim([0, 100])
    #xtickz = axs.get_xticks()
    #axs.set_xticklabels(1/xtickz)
    #axs.plot(, range(len(MR[:,0])), alpha=1, c = 'r', label = 'max = '+str(max_all))
    axs.legend(loc = 'upper right', fontsize=8, frameon = False)
    print("mean obmr = "+str(round(mean_all,5)))
    return fig, axs, line
       
