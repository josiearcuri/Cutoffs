import numpy as np
import math
import matplotlib.pyplot as plt
import itertools as it
"""
class of functions to statistically compare cutoff 2-d spatiotemporal point processes to simulated complete spatial randomness.  This is an implementation of Ripley's K function for 2-D space time.
"""
class RipleysKEstimator_spacetime:
    def __init__(self,t_max, d_max, t_min, d_min, width):
        """initialize estimator
        t_max - int, last year of model iteration
        d_max - int, longest centerline length
        t_min - int, 0 
        d_min - in, 0
        width - constant channel width in m"""
        self.t_max = t_max # last year
        self.d_max = d_max # max centerline length
        self.t_min = t_min # 0
        self.d_min = d_min # 0
        self.width = width  #150
        

    def __call__(self, cutoffs, mode, max_search_d, max_search_t):
        """
        perform main function
        """
        return self.mc_env(cutoffs = cutoffs, nit=99, mode = mode, max_search_d=max_search_d, max_search_t=max_search_t)
    
    def _pairwise_diffs(self, data):
        """
        compute array of distances between every point in 1D
        data - 1-D array of deltas, space or time
        """
        npts = len(data)
        diff = np.zeros(shape=(npts * (npts - 1) // 2), dtype=np.double) 
        k = 0
        for i in range(npts - 1):
            for j in range(i+1, npts):
                diff[k] = abs(data[i] - data[j])
                k += 1
        return diff
    def _near_neigh(self, data):
        """
        compute array of distance between every point and its nearest neighbor in 1D
        data - 1-D array of deltas, space or time
        """
        npts = len(data)
        diff = np.zeros(shape=npts, dtype=np.double)
        for i in range(npts):
            others= np.hstack((data[i::-1], data[i+1:]))
            mask = np.ones(len(data), dtype=bool)
            mask[i] = False
            others = others[mask]
            diff[i] = np.min(abs(data[i] - others))
        return diff
    
    def evaluate(self, data, dist_space, dist_time, mode):
        """
        
        INPUTS
        data - 2-D array of N by 2 size where N is number of cutoffs in dataset, column [:,0] records distance downstream and           [:,1] is time.
        dist_space - 1-D array of radii to search in for intensity estimates
        dist_time - 1-D array of durations to search in for intensity estimates
        mode - statistical measurement to be made, str, either 'K_st, 'K', 'G', or 'H'
        
        OUTPUTS
        stat_d - 1-D Ripleys K at distances dist_space
        stat_t -  1-D temporal Ripleys K at durations dist_time
        stat_dt - 2-D spatiotemporal Ripleys K of M by N{n} size, where M{t} is dist_space, T is dist_time, and each array member is the intensity of the sampled point process at n search distance and t search time. 
        """
        data = np.asarray(data)
        npts = len(data)
        

        stat_d = np.zeros(len(dist_space)) #1-D spatial statistic
        stat_t = np.zeros(len(dist_time))#1-D temporal statistic
        stat_dt = np.zeros((len(dist_space), len(dist_time))) #2D space-time statistic
        
        null = stat_dt.copy()
        if mode == "H":
            """
            H , the probability of finding neighbors in search dist
            """
            deltaspace = self._pairwise_diffs(data[:,0])
            deltatime = self._pairwise_diffs(data[:,1])
            for i in range(len(dist_space)):
                d_indicator = (deltaspace <=dist_space[i])
                stat_d[i] = (d_indicator).sum()
            for i in range(len(dist_time)):
                t_indicator = (deltatime<=dist_time[i])
                stat_t[i] = (t_indicator).sum()
            stat_t = 2*stat_t/(npts*(npts-1))  
            stat_d = 2*stat_d/(npts*(npts-1)) 
            return (stat_d, stat_t)
        if mode == "G":
            """
            G, probability the nearest neighbor is within search dist
            """
            deltaspace = self._near_neigh(data[:,0])
            deltatime = self._near_neigh(data[:,1])
            for i in range(len(dist_space)):
                d_indicator = (deltaspace <=dist_space[i])
                stat_d[i] = (d_indicator).sum()
            for i in range(len(dist_time)):
                t_indicator = (deltatime<=dist_time[i])
                stat_t[i] = (t_indicator).sum()
           
            stat_t = stat_t/(npts)
            stat_d = stat_d/(npts)
            return (stat_d, stat_t)
        if mode == "K":
            """
            number of additional events near other events on time scales of dist_time and spatial scales of dist_space, 2 1-d plots
            """
            deltaspace = self._pairwise_diffs(data[:,0])
            deltatime = self._pairwise_diffs(data[:,1])
            for i in range(len(dist_space)):
                d_indicator = (deltaspace <=dist_space[i])
                stat_d[i] = (d_indicator).sum()
            for i in range(len(dist_time)):
                t_indicator = (deltatime<=dist_time[i])
                stat_t[i] = (t_indicator).sum()
            stat_t = 2*(self.t_max*stat_t/(npts*(npts-1)))
            stat_d = 2*(self.d_max*stat_d/((npts-1)*npts))
            return (stat_d, stat_t)
        if mode == "K_st":
            """
            number of additional events near other events given sepcific search distances and durations. 2D heatmap
            """
            areas = np.ones_like(stat_dt)
            deltaspace = self._pairwise_diffs(data[:,0])
            deltatime = self._pairwise_diffs(data[:,1])
            for x in range(len(dist_space)):
                for t in range(len(dist_time)):
                    dt_indicator = (deltatime<=dist_time[t])&(deltaspace <=dist_space[x])
                    stat_dt[x,t] = (dt_indicator).sum()
            stat_dt = 2*(self.d_max*self.t_max*stat_dt)/(npts*(npts-1))
           
            return(stat_dt)
         
    def mc_env(self,cutoffs, nit, mode, max_search_d, max_search_t): 
        """
        generate random distibutions in same space + time ranges as data
        """
        rng = np.random.default_rng(seed = 42)
        data = cutoffs[['downstream_distance', 'time']].to_numpy() 
        num_samples = len(cutoffs.time)
        r_time = np.linspace(1,max_search_t, max_search_t)
        r_space = np.linspace(self.width,self.width*max_search_d, max_search_d)
        mc_d = np.zeros((len(r_space), nit))
        mc_t = np.zeros((len(r_time), nit))
        mc_dt = np.zeros((len(r_space), len(r_time), nit))
        z = np.zeros((num_samples, 2))
        for i in range(nit):
            z[:,0] = rng.random(size = num_samples)*self.d_max
            z[:,1] = rng.random(size = num_samples)*self.t_max

            k_dt = self.evaluate(data=z, dist_time=r_time, dist_space=r_space, mode='K_st') 
            mc_dt[:, :, i]= k_dt

            k_d, k_t = self.evaluate(data=z, dist_time=r_time, dist_space=r_space, mode='K') 
            mc_d[:,i] = k_d
            mc_t[:,i] = k_t

    
        if mode == 'K_st':
            ## monte carlo envelope - limits on probable randomness
            upper_dt = np.ma.max(mc_dt, axis = 2)
            lower_dt = np.ma.min(mc_dt, axis = 2)
            middle_dt = np.ma.mean(mc_dt, axis = 2)
            
            upper_d = np.ma.max(mc_d, axis = 1)
            lower_d = np.ma.min(mc_d, axis = 1)
            upper_t = np.ma.max(mc_t, axis = 1)
            lower_t = np.ma.min(mc_t, axis = 1)


            #space-time K
            stat_dt = self.evaluate(data=data, dist_time=r_time, dist_space=r_space, mode=mode)
            dependent_clustered = (stat_dt>np.multiply(upper_d.reshape(len(upper_d),1),upper_t))
            dependent_regular = (stat_dt<np.multiply(lower_d.reshape(len(lower_d),1),lower_t))
            

            #locations of statictically nonrandom, dependent K values
            #significantly more aggregated than upper mc env, and
            clustered = (stat_dt>upper_dt)
            regular =  (stat_dt<lower_dt)
            sig_mask = (clustered+regular)*(dependent_clustered+dependent_regular)
         
            normalized = (stat_dt-middle_dt)/(4*np.multiply(r_space.reshape(len(r_space),1),r_time.reshape(1,len(r_time))))
            D = np.sum(normalized*sig_mask)
            self.plot_st(r_space, r_time, normalized, sig_mask, D)
            
            return [D, np.count_nonzero(normalized*sig_mask)]
        else:
            #monte carlo envelope
            upper_d = np.ma.max(mc_d, axis = 1)/(r_space*2)-1
            upper_t = np.ma.max(mc_t, axis = 1)/(r_time*2)-1
        
            lower_d = np.ma.min(mc_d, axis = 1)/(r_space*2)-1
            lower_t = np.ma.min(mc_t, axis = 1)/(r_time*2)-1
            #Simulated Poisson distrubution
            middle_d = np.ma.mean(mc_d, axis = 1)/(r_space*2)-1
            middle_t = np.ma.mean(mc_t, axis = 1)/(r_time*2)-1
        
            #K values
            stat_d, stat_t = self.evaluate(data=data, dist_time=r_time, dist_space=r_space, mode=mode)
            
            #normalize to what's expected under poisson
            stat_d = (stat_d)/(r_space*2) -1
            stat_t = (stat_t)/(r_time*2) -1
            self.plot(upper_d,upper_t, lower_d, lower_t, middle_d, middle_t, r_space, r_time, stat_d, stat_t, num_samples)
        
    def plot_st(self, r_space, r_time, normalized, sig_mask, D):
        plt.rcParams.update({'font.size': 10})
        fig,ax = plt.subplots(figsize = (8,4))
        
        cmap = plt.get_cmap('PiYG')
        
        #im = ax.imshow(np.ma.masked_values(normalized, 0),origin='lower',vmin = -2, vmax = 2, cmap = cmap)
        im =ax.pcolormesh(np.swapaxes(normalized,0,1), cmap = cmap,vmin = -5, vmax = 5)
        #plt.pcolor(np.ma.masked_values(np.swapaxes(normalized*sig_mask,0,1),0), edgecolors='k', linewidths=4, alpha=0.)
        im2 =ax.pcolormesh(np.ma.masked_values(np.swapaxes(normalized*sig_mask,0,1)/np.swapaxes(normalized*sig_mask,0,1),0), zorder=2, linewidths = .01,facecolor='none', edgecolors='k',cmap='gray')
        plt.title('D ='+str(D), pad = 10)
        cbar = ax.figure.colorbar(im, ax=ax)#, ticks = [-2,-1,0,1,2])
        cbar.ax.set_ylabel("[K-null]/[d*t]", va="bottom", rotation=-90)
        #cbar.ax.set_yticklabels(['<-2', '-1', '0','1','>2']) 
        ax.set_ylim(bottom=0, top=max(r_time))
        ax.set_xlim(left=0, right=max(r_space/self.width))
        ax.set_ylabel('time window (years)')
        ax.set_xlabel('search distance (ch-w)')
        ax.set_xticks(r_space[1::2]/self.width, minor=True)
        ax.set_yticks(r_time[1::2], minor=True)
        ax.set_xticks(r_space[1::2]/self.width-1)
        ax.set_yticks(r_time[1::2] -1)
        ax.set_yticklabels((r_time[1::2]).astype(int))
        ax.set_xticklabels((r_space[1::2]/self.width).astype(int))#, rotation='vertical')
        
        ax.tick_params(axis = 'both', which = 'major', top =False, bottom = False, left = False, right = False)
        #ax.grid(True, which='minor', color='k', linewidth=.1)
        return fig
    def plot(self,upper_d,upper_t, lower_d, lower_t, middle_d, middle_t, r_space, r_time, stat_d, stat_t, num_samples):
        #1-d spatial Ripley's K
        fig = plt.figure()
        #plot CSR envelope
        plt.plot(r_space/self.width, upper_d, color='red', ls=':', label='_nolegend_', linewidth = .5)
        plt.plot(r_space/self.width, lower_d, color='red', ls=':', label='_nolegend_', linewidth = .5)
        plt.plot(r_space/self.width, middle_d, color='red', ls=':', label='CSR', linewidth = .5)
        plt.plot(r_space/self.width, stat_d, color = "black", linewidth = .5,label = str(num_samples)+ ' cutoffs')
        plt.legend(loc = 'lower right')
        plt.xlabel("d along centerline (ch-w)")
        plt.ylabel('K/2-d')
        plt.title("Homegrown 1D space EDF")
        #plt.savefig(resultdir + str(year)+"yrs_Space_Ripley_"+mode+".jpg", dpi = 500)
        plt.show()
    
        
        #1-D Temporal Ripley's K
        fig2 = plt.figure()
        #plot CSR envelope
        plt.plot(r_time, upper_t, color='red', ls=':',linewidth = .5, label='_nolegend_')
        plt.plot(r_time, lower_t, color='red', ls=':',linewidth = .5, label='_nolegend_')
        plt.plot(r_time, middle_t, color='red', ls=':',linewidth = .5, label='CSR')
        plt.plot(r_time, stat_t, color = "black", linewidth = .5, label =str(num_samples)+ ' cutoffs')
        plt.legend(loc = 'lower right')
        plt.xlabel("t in years")
        plt.ylabel('K/2 -t')
        plt.title("Homegrown 1D time EDF")
        plt.show()