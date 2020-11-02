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
    
    def _pairwise_diffs(self, data):
        npts = len(data)
        diff = np.zeros(shape=(npts * (npts - 1) // 2), dtype=np.double)
        k = 0
        for i in range(npts - 1):
            size = npts - i - 1
            diff[k:k + size] = abs(data[i] - data[i+1:])
            k += size
        
        plt.show()
        return diff

    
    def evaluate(self, data, dist_space, dist_time):
        data = np.asarray(data)

        if not data.shape[1] == 2:
            raise ValueError('data must be an n by 2 array, where n is the '
                             'number of observed points.')

        npts = len(data)
        intensity = npts/self.area
        intensity_space = npts/self.d_max
        intensity_time = npts/self.t_max
        ripley = np.zeros((len(dist_time), len(dist_space)))
        k_d = np.zeros(len(dist_space))
        k_t = np.zeros(len(dist_time))
        deltaspace = self._pairwise_diffs(data[:,0])
        deltatime = self._pairwise_diffs(data[:,1])
 
        for d in range(len(dist_space)):
            idx_space = (deltaspace < dist_space[d])
            k_d[d] = (2*np.count_nonzero(idx_space)/(intensity_space*(npts-1)))
            for t in range(len(dist_time)):
                idx_time = (deltatime < dist_time[t])
                k_t[t] = (2*np.count_nonzero(idx_time)/(intensity_time*(npts-1)))
    
                ripley[t,d]  = (np.count_nonzero(idx_time==idx_space)/(intensity*(npts-1)))
                
        return (ripley, k_d/2 - dist_space, k_t/2-dist_time) 
                
            
