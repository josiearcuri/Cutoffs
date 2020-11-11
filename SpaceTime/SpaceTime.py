import numpy as np
import math
import matplotlib.pyplot as plt


__all__ = ['RipleysKEstimator_spacetime']


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
                
            
