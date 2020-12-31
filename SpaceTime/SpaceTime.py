import numpy as np
import math
import matplotlib.pyplot as plt
import itertools as it

class RipleysKEstimator_spacetime:
    def __init__(self,t_max, d_max, t_min, d_min, width):
        self.t_max = t_max # last year
        self.d_max = d_max # max centerline length
        self.t_min = t_min # 0
        self.d_min = d_min # 0
        self.width = width 
        

    def __call__(self, cutoffs, mode):
        return self.mc_env(cutoffs = cutoffs, nit=99, mode = mode)
    
    def _pairwise_diffs(self, data):
        """compute array of distanced between every point in 1D"""
        npts = len(data)
        diff = np.zeros(shape=(npts * (npts - 1) // 2), dtype=np.double) 
        k = 0
        for i in range(npts - 1):
            for j in range(i+1, npts):
                diff[k] = abs(data[i] - data[j])
                k += 1
        return diff
    def _near_neigh(self, data):
        """compute array of distance between every point and its nearest neighbor in 1D"""
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
        data = np.asarray(data)
        npts = len(data)
        

        stat_d = np.zeros(len(dist_space)) #1-D spatial statistic
        stat_t = np.zeros(len(dist_time))#1-D temporal statistic
        stat_dt = np.zeros((len(dist_space), len(dist_time))) #2D space-time statistic
        
        null = stat_dt.copy()
        if mode == "H":
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
            deltaspace = self._pairwise_diffs(data[:,0])
            deltatime = self._pairwise_diffs(data[:,1])
            for i in range(len(dist_space)):
                d_indicator = (deltaspace <=dist_space[i])
                stat_d[i] = (d_indicator).sum()
            for i in range(len(dist_time)):
                t_indicator = (deltatime<=dist_time[i])
                stat_t[i] = (t_indicator).sum()
            stat_t = (self.t_max*stat_t/(npts*(npts-1)))
            stat_d = (self.d_max*stat_d/((npts-1)*npts))
            return (stat_d, stat_t)
        if mode == "K_st":
            areas = np.ones_like(stat_dt)
            deltaspace = self._pairwise_diffs(data[:,0])
            deltatime = self._pairwise_diffs(data[:,1])
            for x in range(len(dist_space)):
                for t in range(len(dist_time)):
                    dt_indicator = (deltatime<=dist_time[t])&(deltaspace <=dist_space[x])
                    stat_dt[x,t] = (dt_indicator).sum()
            stat_dt = (self.d_max*self.t_max*stat_dt)/(npts*(npts-1))
           
            return(stat_dt)
         
    def mc_env(self,cutoffs, nit, mode): 
        #generate random distibutions in same space + time ranges as data
        rng = np.random.default_rng(seed = 42)
        data = cutoffs[['downstream_distance', 'time']].to_numpy() 
        num_samples = len(cutoffs.time)
        r_time = np.linspace(1,50)
        r_space = np.linspace(self.width,self.width*50)
        mc_d = np.zeros((len(r_space), nit))
        mc_t = np.zeros((len(r_time), nit))
        mc_dt = np.zeros((len(r_space), len(r_time), nit))
        z = np.zeros((num_samples, 2))
        for i in range(nit):
            z[:,0] = rng.random(size = num_samples)*self.d_max
            z[:,1] = rng.random(size = num_samples)*self.t_max
            if mode == 'K_st':
                k_dt = self.evaluate(data=z, dist_time=r_time, dist_space=r_space, mode=mode) 
                mc_dt[:, :, i]= k_dt
            else:
                k_d, k_t = self.evaluate(data=z, dist_time=r_time, dist_space=r_space, mode=mode) 
                mc_d[:,i] = k_d
                mc_t[:,i] =k_t

    
        if mode == 'K_st':
      
            upper_dt = np.ma.max(mc_dt, axis = 2)
            lower_dt = np.ma.min(mc_dt, axis = 2)
            middle_dt = np.ma.mean(mc_dt, axis = 2)
            k_d, k_t = self.evaluate(data=data, dist_time=r_time, dist_space=r_space, mode='K') 

            stat_dt = self.evaluate(data=data, dist_time=r_time, dist_space=r_space, mode=mode)
 

            clustered = (stat_dt>upper_dt)*(stat_dt>(k_d*(k_t.reshape(len(r_time),1))))
            regular =  (stat_dt<lower_dt)*(stat_dt<(k_d*(k_t.reshape(len(r_time),1))))
            nonrandom = (clustered+regular)*((stat_dt/(r_space*(r_time.reshape(len(r_time),1))))-2)
            self.plot_st(r_space, r_time, nonrandom)
        else:
            upper_d = np.ma.max(mc_d, axis = 1) -(r_space) 
            upper_t = np.ma.max(mc_t, axis = 1) -(r_time)
        
            lower_d = np.ma.min(mc_d, axis = 1)-(r_space)
            lower_t = np.ma.min(mc_t, axis = 1) -(r_time)
   
            middle_d = np.ma.mean(mc_d, axis = 1) -(r_space)
            middle_t = np.ma.mean(mc_t, axis = 1) -(r_time)
            stat_d, stat_t = self.evaluate(data=data, dist_time=r_time, dist_space=r_space, mode=mode)
            stat_d = (stat_d) -(r_space) 
            stat_t = (stat_t) -(r_time)
            self.plot(upper_d,upper_t, lower_d, lower_t, middle_d, middle_t, r_space, r_time, stat_d, stat_t)
        
    def plot_st(self, r_space, r_time, normalized):
        plt.rcParams.update({'font.size': 8})
        ax = plt.subplot(1,1,1)
        ax.set_ylabel('t (years)')
        ax.set_xlabel('d (ch-w)')
        ax.set_xticks(np.arange(len(r_space))[1::2])
        ax.set_yticks(np.arange(len(r_time))[1::2])
        ax.set_yticklabels((r_time[1::2]).astype(int))
        ax.set_xticklabels((r_space[1::2]/self.width).astype(int))#, rotation='vertical')
        cmap = plt.get_cmap('RdYlGn')
        cmap.set_bad('white')
        #im = ax.imshow(np.ma.masked_values(normalized, 0), origin='lower',vmin = -5, vmax = 5, cmap = cmap)
        im = ax.pcolormesh(np.ma.masked_values(normalized, 0),vmin = -2, vmax = 2, cmap = cmap, edgecolors='k', linewidths=.005, shading='auto')
        plt.title("significant spatiotemporal cutoff intensities")

        cbar = ax.figure.colorbar(im, ax=ax, ticks = [-2,-1,0,1,2])
        cbar.ax.set_ylabel("# additional cutoffs", va="bottom", rotation=-90)
        cbar.ax.set_yticklabels(['<-2', '-1', '0','1','>2']) 
            #
        plt.show()
    def plot(self,upper_d,upper_t, lower_d, lower_t, middle_d, middle_t, r_space, r_time, stat_d, stat_t):
        fig = plt.figure()
     #plot CSR envelope
        plt.plot(r_space/self.width, upper_d, color='red', ls=':', label='_nolegend_', linewidth = .5)
        plt.plot(r_space/self.width, lower_d, color='red', ls=':', label='_nolegend_', linewidth = .5)
        plt.plot(r_space/self.width, middle_d, color='red', ls=':', label='CSR', linewidth = .5)
        plt.plot(r_space/self.width, stat_d, color = "black", linewidth = .5,label = str(num_samples)+ ' cutoffs')
        plt.legend(loc = 'lower right')
        plt.xlabel("d along centerline (ch-w)")
        plt.ylabel(mode)
        plt.title("Homegrown 1D space EDF")
        #plt.savefig(resultdir + str(year)+"yrs_Space_Ripley_"+mode+".jpg", dpi = 500)
        plt.show()
    
        fig2 = plt.figure()
     #plot CSR envelope
        plt.plot(r_time, upper_t, color='red', ls=':',linewidth = .5, label='_nolegend_')
        plt.plot(r_time, lower_t, color='red', ls=':',linewidth = .5, label='_nolegend_')
        plt.plot(r_time, middle_t, color='red', ls=':',linewidth = .5, label='CSR')
        plt.plot(r_time, stat_t, color = "black", linewidth = .5, label =str(num_samples)+ ' cutoffs')
        plt.legend(loc = 'lower right')
        plt.xlabel("t in years")
        plt.ylabel(mode)
        plt.title("Homegrown 1D time EDF")
        plt.show()