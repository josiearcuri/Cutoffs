import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist


__all__ = ['RipleysKEstimator_spacetime']


class RipleysKEstimator_spacetime:
    def __init__(self,t_max=None, d_max=None, t_min=None, d_min=None):
        self.area = (t_max-t_min)*(d_max-d_min)
        self.t_max = t_max
        self.d_max = d_max
        self.t_min = t_min
        self.d_min = d_min

    @property
    def area(self):
        return self._area

    @area.setter
    def area(self, value):
        if isinstance(value, (float, int)) and value > 0:
            self._area = value
        else:
            raise ValueError('area is expected to be a positive number. '
                             'Got {}.'.format(value))

    @property
    def d_max(self):
        return self._d_max

    @d_max.setter
    def d_max(self, value):
        if value is None or isinstance(value, (float, int)):
            self._y_max = value
        else:
            raise ValueError('d_max is expected to be a real number '
                             'or None. Got {}.'.format(value))

    @property
    def t_max(self):
        return self._t_max

    @t_max.setter
    def t_max(self, value):
        if value is None or isinstance(value, (float, int)):
            self._t_max = value
        else:
            raise ValueError('t_max is expected to be a real number '
                             'or None. Got {}.'.format(value))

    @property
    def d_min(self):
        return self._d_min

    @d_min.setter
    def d_min(self, value):
        if value is None or isinstance(value, (float, int)):
            self._d_min = value
        else:
            raise ValueError('d_min is expected to be a real number. '
                             'Got {}.'.format(value))

    @property
    def t_min(self):
        return self._d_min

    @t_min.setter
    def t_min(self, value):
        if value is None or isinstance(value, (float, int)):
            self._t_min = value
        else:
            raise ValueError('t_min is expected to be a real number. '
                             'Got {}.'.format(value))

    def __call__(self, data, radii_space, radii_time):
        return self.evaluate(data=data, radii_space=radii_space,radii_time = radii_time)
    
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
    
    def evaluate(self, data, radii_space, radii_time):
        data = np.asarray(data)

        if not data.shape[1] == 2:
            raise ValueError('data must be an n by 2 array, where n is the '
                             'number of observed points.')

        npts = len(data)
        ripley = np.zeros((len(radii_time), len(radii_space)))
                
        deltaspace = self._pairwise_diffs(data[:,0])
        deltatime = self._pairwise_diffs(data[:,1])
  
        for d in range(len(radii_space)):
            for t in range(len(radii_time)):
                idx_time = (deltatime < radii_time[t])
                idx_space = (deltaspace < radii_space[d])
                ripley[t,d]  = np.count_nonzero(idx_time==idx_space)
             
        ripley = ripley*self.area/(npts*npts)
        return ripley
                
            
