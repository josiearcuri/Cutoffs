import numpy as np
import math
import matplotlib.pyplot as plt


__all__ = ['RipleysKEstimator_spacetime', 'PreliminaryTesting']


class RipleysKEstimator_spacetime:
    def __init__(self,t_max, d_max, t_min, d_min):
        self.area = (t_max-t_min)*(d_max-d_min)
        self.t_max = t_max
        self.d_max = d_max
        self.t_min = t_min
        self.d_min = d_min
    def __call__(self, data, dist_space, dist_time):
        return self.evaluate(data=data, dist_space=dist_space,dist_time = dist_time)
    
    def _pairwise_diffs(self, data, omega):
        npts = len(data)
        diff = np.zeros(shape=(npts * (npts - 1) // 2), dtype=np.double)
        d_xomega = np.zeros(shape=(npts * (npts - 1) // 2), dtype=np.double)
        d_yomega = np.zeros(shape=(npts * (npts - 1) // 2), dtype=np.double)
        k = 0
        for i in range(npts - 1):
            size = npts - i - 1
            diff[k:k + size] = abs(data[i] - data[i+1:])
            d_xomega[k:k + size] = (omega-data[i])*np.ones((size))
            d_yomega[k:k + size] = omega-data[i+1:]
            k += size
        
        return diff, d_xomega, d_yomega

    def _edge_correction(self, diff, d_xomega, d_yomega):

        k_xy = 1+(diff>d_yomega)
        k_yx = 1+(diff>d_yomega)
        weights = .5*(k_xy +k_yx)

        return weights
    
    def evaluate(self, data, dist_space, dist_time):
        data = np.asarray(data)

        if not data.shape[1] == 2:
            raise ValueError('data must be an n by 2 array, where n is the '
                             'number of observed points.')

        npts = len(data)
        
        intensity_space = npts/((self.d_max - self.d_min))#*(self.d_max - self.d_min))
        intensity_time = npts/(self.t_max - self.t_min)
        #ripley = np.zeros((len(dist_space),len(dist_time)))
        k_d = np.zeros(len(dist_space))
        k_t = np.zeros(len(dist_time))

        deltaspace, d_xomega, d_yomega = self._pairwise_diffs(data[:,0], (self.d_max - self.d_min))
        deltatime, t_xomega, t_yomega = self._pairwise_diffs(data[:,1], (self.t_max - self.t_min))
        d_weights = self._edge_correction(deltaspace, d_xomega, d_yomega) 
        t_weights = self._edge_correction(deltatime, t_xomega, t_yomega)

    
        for d in range(len(dist_space)):
            d_indicator = (deltaspace < dist_space[d])
            k_d[d] = ((d_indicator)).sum()#*d_weights).sum()
            for t in range(len(dist_time)):
                t_indicator = (deltatime<dist_time[t])
                k_t[t] = ((t_indicator)).sum()#*t_weights).sum()
                #dt_indicator = (d_indicator == t_indicator)
                #ripley[d,t]  = (dt_indicator*t_weights*d_weights).sum()

        #ripley = 2*self.area*ripley/(npts*(npts-1)) - (k_t*d_t)
        k_t = (2*k_t/(intensity_time*(npts)))  
        k_d = (2*k_d/(intensity_space*(npts)))      
        return (k_d-(2*dist_space), k_t-(2*dist_time)) 
                
class PreliminaryTesting:
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
        if mode == "K":
            deltaspace = self._pairwise_diffs(data[:,0])
            deltatime = self._pairwise_diffs(data[:,1])
            for i in range(len(dist_space)):
                d_indicator = (deltaspace <=dist_space[i])
                stat_d[i] = (d_indicator).sum()
            for i in range(len(dist_time)):
                t_indicator = (deltatime<=dist_time[i])
                stat_t[i] = (t_indicator).sum()
            stat_t = (self.t_max*stat_t/((npts)*npts*(npts-1)))
            stat_d = (self.d_max*stat_d/((npts)*npts*(npts-1))) 
        return (stat_d, stat_t) 
    def mc_env(self,cutoffs, nit, mode): 
            #generate random distibutions in same space + time ranges as data
        rng = np.random.default_rng(seed = 42)
        data = cutoffs[['downstream_distance', 'time']].to_numpy() 
        num_samples = len(cutoffs.time)
        r_time = np.linspace(1,int(np.sqrt(self.t_max)), 10)
        r_space = np.linspace(1,int(np.sqrt(self.d_max)), 100)
        mc_d = np.zeros((len(r_space), nit))
        mc_t = np.zeros((len(r_time), nit))
        z = np.zeros((num_samples, 2))
        for i in range(nit):
            z[:,0] = rng.random(size = num_samples)*self.d_max
            z[:,1] = rng.random(size = num_samples)*self.t_max
            k_d, k_t = self.evaluate(data=z, dist_time=r_time, dist_space=r_space, mode=mode) 

            mc_d[:,i] = k_d
            mc_t[:,i] =k_t
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
        
