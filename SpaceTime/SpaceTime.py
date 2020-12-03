import numpy as np
import math
import matplotlib.pyplot as plt


__all__ = ['RipleysKEstimator_spacetime']

class RipleysKEstimator_spacetime:
    def __init__(self,t_max, d_max, t_min, d_min):
        self.t_max = t_max
        self.d_max = d_max
        self.t_min = t_min
        self.d_min = d_min
        

    def __call__(self, cutoffs, mode):
        return self.mc_env(cutoffs = cutoffs, nit=99, mode = mode)
    
    def _pairwise_diffs(self, data):
        npts = len(data)
        diff = np.zeros(shape=(npts * (npts - 1) // 2), dtype=np.double)
        
        k = 0
        for i in range(npts - 1):
            size = npts - i - 1
            diff[k:k + size] = abs(data[i] - data[i+1:])
            k += size
        
        return diff
    def _near_neigh(self, data):
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

        if not data.shape[1] == 2:
            raise ValueError('data must be an n by 2 array, where n is the '
                             'number of observed points.')

        npts = len(data)
        

        stat_d = np.zeros(len(dist_space))
        stat_t = np.zeros(len(dist_time))
        stat_dt = np.zeros((len(dist_space), len(dist_time)))
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
                    areas[x,t] = dist_space[x]*dist_time[t]
            stat_dt = (self.d_max*self.t_max*stat_dt)/(npts*(npts-1))
           
            return(stat_dt)
         
    def mc_env(self,cutoffs, nit, mode): 
            #generate random distibutions in same space + time ranges as data
        rng = np.random.default_rng(seed = 42)
        data = cutoffs[['downstream_distance', 'time']].to_numpy() 
        num_samples = len(cutoffs.time)
        r_time = np.linspace(1,100, 50)
        r_space = np.linspace(100,40000, 100)
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
            stat_dt = self.evaluate(data=data, dist_time=r_time, dist_space=r_space, mode=mode)
            
            X, Y = np.meshgrid(r_time, r_space)
            fig = plt.figure()
            ax = fig.gca(projection='3d')
            ax.set_xlabel('time (years)')
            ax.set_ylabel('distance (km)')
            ax.set_zlabel('K_st')
            # Plot the surface.

            #urf = ax.plot_surface(X,Y,stat_dt,
             #         linewidth=0)
            surf = ax.plot_surface(X,Y/1000,stat_dt>upper_dt,
                       linewidth=0)
            # Add a color bar which maps values to colors.
            #fig.colorbar(surf)

            plt.show()
        else:
            upper_d = np.ma.max(mc_d, axis = 1)
            upper_t = np.ma.max(mc_t, axis = 1)
        
            lower_d = np.ma.min(mc_d, axis = 1)
            lower_t = np.ma.min(mc_t, axis = 1)
   
            middle_d = np.ma.mean(mc_d, axis = 1)
            middle_t = np.ma.mean(mc_t, axis = 1)
            stat_d, stat_t = self.evaluate(data=data, dist_time=r_time, dist_space=r_space, mode=mode)
            
            fig = plt.figure()
     #plot CSR envelope
            plt.plot(r_space, upper_d, color='red', ls=':', label='_nolegend_', linewidth = .5)
            plt.plot(r_space, lower_d, color='red', ls=':', label='_nolegend_', linewidth = .5)
            plt.plot(r_space, middle_d, color='red', ls=':', label='CSR', linewidth = .5)
            plt.plot(r_space, stat_d, color = "black", linewidth = .5,label = str(num_samples)+ ' cutoffs')
            plt.legend(loc = 'lower right')
            plt.xlabel("d in m along centerline")
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
        
