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
    def _edge_correction(self, data,deltaspace,deltatime, npts):
        hor_dist = np.zeros(shape=(npts * (npts - 1)) // 2,
                                dtype=np.double)
        vert_dist = np.zeros(shape=(npts * (npts - 1)) // 2,
                                dtype=np.double)

        for k in range(npts - 1):
            min_hor_dist = min(self.d_max - deltaspace[k],
                                   deltaspace[k] - self.d_min)
            min_ver_dist = min(self.t_max - deltatime[k],
                                   deltatime[k] - self.t_min)
            start = (k * (2 * (npts - 1) - (k - 1))) // 2
            end = ((k + 1) * (2 * (npts - 1) - k)) // 2
            hor_dist[start: end] = min_hor_dist * np.ones(npts - 1 - k)
            vert_dist[start: end] = min_ver_dist * np.ones(npts - 1 - k)

        dist = np.hypot(deltaspace, deltatime)
        dist_ind_dt = dist <= np.hypot(hor_dist, vert_dist)
        dist_ind_d = deltaspace<=hor_dist
        dist_ind_t = deltatime<=vert_dist
        w1_dt = (1 + (np.arccos(np.minimum(vert_dist, dist) / dist) + np.arccos(np.minimum(hor_dist, dist) / dist)) / np.pi)
        w2_dt = (3 / 4 - 0.5 * (np.arccos(vert_dist / dist * ~dist_ind_dt) +np.arccos(hor_dist / dist * ~dist_ind_dt)) / np.pi)
        
        w1_d = (1 -(np.minimum(hor_dist, deltaspace) / deltaspace))
        #w2_d = (3 / 4 - 0.5 * (hor_dist / dist * ~dist_ind_d))
        w1_t = (1 -(np.minimum(vert_dist, deltatime) / deltatime))
        #w2_t = (3 / 4 - 0.5 * (vert_dist / dist * ~dist_ind_d))

        
        d_weights = np.where((dist_ind_d)==0, 1, (dist_ind_d * w1_d) 
        t_weights = dist_ind_t * w1_t 
        print(d_weights)
        return d_weights, t_weights
    
    def evaluate(self, data, dist_space, dist_time):
        data = np.asarray(data)

        if not data.shape[1] == 2:
            raise ValueError('data must be an n by 2 array, where n is the '
                             'number of observed points.')

        npts = len(data)
        intensity_volume = npts/self.area#((self.d_max - self.d_min)*self.area)
        intensity_space = npts/(self.d_max - self.d_min)#*(self.d_max - self.d_min))
        intensity_time = npts/(self.t_max - self.t_min)
        ripley = np.zeros((len(dist_space),len(dist_time)))
        k_d = np.zeros(len(dist_space))
        k_t = np.zeros(len(dist_time))
        deltaspace = self._pairwise_diffs(data[:,0])
        deltatime = self._pairwise_diffs(data[:,1])
        d_weights, t_weights= self._edge_correction(data,deltaspace,deltatime, npts)
        intersec_area_space = ((self.d_max - self.d_min) - deltaspace)
        intersec_area_time =  ((self.t_max - self.t_min) - deltatime)
        intersec_area = (((self.d_max - self.d_min) - deltaspace)*((self.t_max - self.t_min) - deltatime))
    
        for d in range(len(dist_space)):
            d_indicator = (deltaspace < dist_space[d])
            k_d[d] = (d_indicator/d_weights).sum()#((1/intersec_area_space)*d_indicator).sum()
            for t in range(len(dist_time)):
                t_indicator = (deltatime<dist_time[t])
                k_t[t] = (t_indicator/t_weights).sum()#((1/intersec_area_time)*t_indicator).sum()
                dt_indicator = (d_indicator == t_indicator)
                ripley[d,t]  = (dt_indicator/(d_weights*t_weights)).sum()#((1/intersec_area)*dt_indicator).sum()
                
        ripley = ripley/(intensity_volume)
        k_t = k_t/(intensity_time)
        k_d = k_d/(intensity_space)       
        return (ripley - (k_d*k_t), (k_d/2) - dist_space, (k_t/2)-dist_time) 
                
            
