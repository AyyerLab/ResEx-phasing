import os
import numpy as np
from scipy import interpolate
import mrcfile

# Note: update_support (solvent flattening) not implemented here

class Input():
    def __init__(self, size, vol_object=None):
        self.size = size
        self.num_supp = 0
        self.vol = vol_object
        
        self.obs_mag = None
        self.bragg_calc = None
        self.bragg_mask = None
        self.support = None
        self.supp_loc = None
        self.inverse_cdf = None
        self.local_variation = None
        self._intrad = None

    def parse_intens(self, fname, scale, minsubt=False):
        with mrcfile.open(fname, 'r') as mrc:
            intens = np.fft.ifftshift(mrc.data)
            if intens.dtype != np.dtype('f4'):
                raise TypeError('Intensity file needs to have float32 data')
        self.obs_mag = np.empty_like(intens)
        print("Scale factor = %f" % scale)

        if minsubt:
            print("Positivizing intensities")
            self._calc_intrad()
            radmin = np.ones(self.size, dtype='f4') * np.finfo('f4').max
            sel = (intens != -.1) & (intens > -1.e3)
            np.minimum.at(radmin, self._intrad[sel], intens[sel])
            intens[sel] -= radmin[self._intrad[sel]]

        sel = (intens > 0)
        self.obs_mag[sel] = np.sqrt(intens[sel]) * scale
        self.obs_mag[~sel] = intens[~sel]

    def parse_bragg(self, fname, braggqmax):
        size = self.size
        c = size // 2
        with mrcfile.open(fname, 'r') as mrc:
            self.bragg_calc = np.fft.ifftshift(mrc.data)
            if self.bragg_calc.dtype != np.dtype('c8'):
                raise TypeError('Bragg file needs to have complex64 data')
        
        x, y, z = np.indices((size, size, size))
        self.bragg_calc *= np.exp(-1j * 2. * np.pi * np.fft.ifftshift(x + y + z - 3*c) * c / size)
        self._calc_intrad()
        self.bragg_mask = self._intrad < braggqmax*c
        self.bragg_calc[~self.bragg_mask] = np.finfo('f4').max
        
    def parse_support(self, fname):
        with mrcfile.open(fname, 'r') as mrc:
            self.support = mrc.data
            if self.support.dtype != np.dtype('i1'):
                raise TypeError('Support file needs to have int8 data')

        sx, sy, sz = np.where(self.support > 0)
        self.support_bounds = np.array([sx.min(), sx.max(), sy.min(), sy.max(), sz.min(), sz.max()])
        self.supp_loc = np.where(self.support.ravel() > 0)[0]
        self.num_supp = (self.support>0).sum()
        
        print("num_supp = %ld" % self.num_supp)
        print("Support bounds:", self.support_bounds) 

    def init_iterate(self, model, fname=None, bg_fname=None, do_bg_fitting=False, fixed_seed=False):
        size = self.size
        vol = size*size*size
        do_random_model = False
        do_init_bg = False
        
        if fname is None or not os.path.isfile(fname):
            do_random_model = True
        else:
            with mrcfile.open(fname, 'r') as mrc:
                if mrc.data.dtype != np.dtype('f4'):
                    raise TypeError('Initial model needs to have float32 data')
                model[:vol] = mrc.data
            print('Starting from', fname)
        
        if do_bg_fitting:
            if bg_fname is None or not os.path.isfile(bg_fname):
                do_init_bg = True
                print("...with uniform background")
            else:
                with mrcfile.open(bg_fname, 'r') as mrc:
                    if mrc.data.dtype != np.dtype('f4'):
                        raise TypeError('Initial background model needs to have float32 data')
                    model[vol:] = mrc.data
                print("...with background from", bg_fname)
        
        if fixed_seed:
            print("==== Fixed seed mode ====")
            np.random.seed(0x5EED)
        
        if do_random_model:
            if self.support is None:
                raise ValueError('Need support to generate random model')
            model[:vol] = 0
            model[:vol][self.support > 0] = np.random.random(self.num_supp)

        if do_init_bg:
            model[vol:] = 0
            val = sqrt(vol)
            model[vol:][self.obs_mag > 0] = val

    def parse_histogram(self, fname):
        if self.num_supp == 0:
            raise ValueError('num_supp not defined. Parse support first')
        val, hist = np.loadtxt(fname, skiprows=1, unpack=True)
        cdf = np.insert(np.cumsum(hist, dtype='f8'), 0, 0) * self.num_supp / hist.sum()
        cdf = 0.5 * (cdf[:-1] + cdf[1:])
        ifunc = interpolate.interp1d(cdf, val, copy=False, bounds_error=False, fill_value=(val[0], val[-1]))
        self.inverse_cdf = ifunc(np.arange(self.num_supp))

    def match_histogram(self, in_arr, out_arr):
        '''Histogram matching
           Matches histogram within support volume using inverse_cdf array.
           Can be done in-place
        '''
        supp_val = in_arr[self.supp_loc] 
        sorter = supp_val.argsort()
        out_arr.fill(0)
        out_arr[self.supp_loc[sorter]] = self.inverse_cdf

    def match_bragg(self, fdens, delta=0.):
        '''Match Bragg Fourier components
           delta parameter allows for small distance from input value at each voxel
        '''
        if delta == 0.:
            fdens[self.bragg_mask] = self.bragg_calc[self.bragg_mask]
        else:
            temp = fdens - self.bragg_calc
            mag = np.abs(temp)
            sel = (mag < delta)
            fdens[self.bragg_mask & sel] = self.bragg_calc[self.bragg_mask & sel]
            fdens[self.bragg_mask & ~sel] = self.bragg_calc[self.bragg_mask & ~sel] + delta / mag[self.bragg_mask & ~sel] * temp[self.bragg_mask & ~sel]

    def _calc_intrad(self):
        if self._intrad is not None:
            pass
        elif self.vol is None:
            x, y, z = np.indices(3*(self.size,))
            cen = self.size // 2
            x -= cen
            y -= cen
            z -= cen
            self._intrad = np.sqrt(x*x + y*y + z*z).astype('i4')
            self._intrad = np.fft.ifftshift(self._intrad)
        elif self.vol.intrad is None:
            self.vol.init_radavg()
            self._intrad = self.vol.intrad
        else:
            self._intrad = self.vol.intrad
