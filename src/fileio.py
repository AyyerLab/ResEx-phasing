import os
from scipy import interpolate
import numpy
import mrcfile
try:
    import cupy as np
    CUDA = True
except ImportError:
    import numpy as np
    CUDA = False

# Note: update_support (solvent flattening) not implemented here

class IO():
    def __init__(self, size, vol_object=None):
        self.size = size
        self.vol = vol_object

        self.output_prefix = './'
        self._intrad = None

    def parse_intens(self, proj, fname, scale=1., minsubt=False):
        with mrcfile.open(fname, 'r') as mrc:
            proj.obs_mag = np.fft.ifftshift(mrc.data)
            if proj.obs_mag.dtype != np.dtype('f4'):
                raise TypeError('Intensity file needs to have float32 data')
        print("Scale factor = %f" % scale)

        if minsubt:
            print("Positivizing intensities")
            self._calc_intrad()
            radmin = np.ones(self.size, dtype='f4') * np.finfo('f4').max
            sel = (proj.obs_mag != -1.) & (proj.obs_mag > -1.e3)
            #np.minimum.at(radmin, self._intrad[sel], proj.obs_mag[sel])
            radmin = np.array([float(np.min(proj.obs_mag[sel & (self._intrad==r)])) for r in range(int(self._intrad.max()+1))])
            proj.obs_mag[sel] -= radmin[self._intrad[sel]]

        sel = (proj.obs_mag > 0)
        proj.obs_mag[sel] = np.sqrt(proj.obs_mag[sel]) * scale

    def parse_bragg(self, proj, fname, braggqmax=1.):
        size = self.size
        c = size // 2
        with mrcfile.open(fname, 'r') as mrc:
            proj.bragg_calc = np.fft.ifftshift(mrc.data)
            if proj.bragg_calc.dtype != np.dtype('c8'):
                raise TypeError('Bragg file needs to have complex64 data')
        print('Bragg q_max =', braggqmax)

        x, y, z = numpy.indices((size, size, size), dtype='i4')
        proj.bragg_calc *= np.array(numpy.exp(-1j * 2. * numpy.pi * numpy.fft.ifftshift(x + y + z - 3*c) * c / size))
        self._calc_intrad()
        proj.bragg_mask = (self._intrad < braggqmax*c)
        proj.bragg_calc[~proj.bragg_mask] = np.finfo('f4').max

    def parse_support(self, proj, fname):
        with mrcfile.open(fname, 'r') as mrc:
            proj.support = np.array(mrc.data)
            if proj.support.dtype != np.dtype('i1'):
                raise TypeError('Support file needs to have int8 data')

        suppx, suppy, suppz = np.where(proj.support > 0)
        if proj.do_histogram:
            proj.supp_loc = np.where(proj.support.ravel() > 0)[0]
        proj.num_supp = int((proj.support > 0).sum())

        print("num_supp = %ld" % proj.num_supp)

    def init_iterate(self, proj, model, fname=None, bg_fname=None, do_bg_fitting=False, fixed_seed=False, quiet=False):
        do_random_model = False
        do_init_bg = False

        if fname is None or not os.path.isfile(fname):
            do_random_model = True
            if not quiet: print('Random initial guess')
        else:
            with mrcfile.open(fname, 'r') as mrc:
                if mrc.data.dtype != np.dtype('f4'):
                    raise TypeError('Initial model needs to have float32 data')
                model[0] = np.array(mrc.data)
            if not quiet: print('Starting from', fname)

        if do_bg_fitting:
            if bg_fname is None or not os.path.isfile(bg_fname):
                do_init_bg = True
                if not quiet: print("...with uniform background")
            else:
                with mrcfile.open(bg_fname, 'r') as mrc:
                    if mrc.data.dtype != np.dtype('f4'):
                        raise TypeError('Initial background model needs to have float32 data')
                    model[1] = np.array(mrc.data)
                if not quiet: print("...with background from", bg_fname)

        if fixed_seed:
            if not quiet: print("==== Fixed seed mode ====")
            np.random.seed(0x5EED)

        if do_random_model:
            if proj.support is None:
                raise ValueError('Need support to generate random model')
            model[0] = 0
            model[0][proj.support > 0] = np.random.random(proj.num_supp)

        if do_init_bg:
            model[1] = 0
            val = np.sqrt(model[1].size)
            model[1][proj.obs_mag > 0] = val

    def parse_histogram(self, proj, fname):
        if proj.num_supp == 0:
            raise ValueError('num_supp not defined. Parse support first')
        val, hist = numpy.loadtxt(fname, skiprows=1, unpack=True)
        cdf = numpy.insert(numpy.cumsum(hist, dtype='f8'), 0, 0) * int(proj.num_supp) / hist.sum()
        cdf = 0.5 * (cdf[:-1] + cdf[1:])
        ifunc = interpolate.interp1d(cdf, val, copy=False, bounds_error=False, fill_value=(val[0], val[-1]))
        proj.inverse_cdf = np.array(ifunc(numpy.arange(int(proj.num_supp))))

    def make_recon_folders(self, do_bg_fitting, do_local_variation):
        os.makedirs("%s-slices" % self.output_prefix, 0o755, exist_ok=True)
        os.makedirs("%s-fslices" % self.output_prefix, 0o755, exist_ok=True)
        if do_bg_fitting:
            os.makedirs("%s-radavg" % self.output_prefix, 0o755, exist_ok=True)
        if do_local_variation:
            os.makedirs("%s-support" % self.output_prefix, 0o755, exist_ok=True)

    def _calc_intrad(self):
        if self._intrad is not None:
            pass
        elif self.vol is None:
            x, y, z = numpy.indices(3*(self.size,))
            cen = self.size // 2
            x -= cen
            y -= cen
            z -= cen
            self._intrad = numpy.sqrt(x*x + y*y + z*z).astype('i4')
            self._intrad = np.array(numpy.fft.ifftshift(self._intrad))
        elif self.vol.intrad is None:
            self.vol.init_radavg()
            self._intrad = self.vol.intrad
        else:
            self._intrad = self.vol.intrad

    def dump_slices(self, vol, fname, label='', is_fourier=False, is_support=False):
        '''Save orthogonal central slices to file
        '''
        vsizes = np.ones(3, dtype='f4')
        if len(vol.shape) == 4:
            vol = vol[0]
        if is_fourier:
            temp = np.fft.fftshift(vol)
            vsizes *= -1.
        else:
            temp = vol

        size = vol.shape[-1]
        cen = size // 2
        if is_support:
            slices = np.empty((3, size, size), dtype='i1')
        else:
            slices = np.empty((3, size, size), dtype='f4')
        slices[0] = temp[cen]
        slices[1] = temp[:, cen]
        slices[2] = temp[:, :, cen]

        self.save_as_map(fname, slices, vsizes, label)

    def save_current(self, phas, proj, i, time1, time2, error):
        with open("%s-log.dat" % self.output_prefix, "a") as f:
            f.write("%.4d  %.2e  %.3e\n" % (i, time2 - time1, error))

        self.dump_slices(phas.p1, "%s-slices/%.4d.ccp4" % (self.output_prefix, i), "ResEx-recon p1 %d\n" % i)
        self.dump_slices(proj.exp_mag, "%s-fslices/%.4d.ccp4" % (self.output_prefix, i), "ResEx-recon exp_mag %d\n" % i, is_fourier=True)

        if proj.do_local_variation:
            self.dump_slices(proj.support, "%s-support/%.4d.ccp4" % (self.output_prefix, i), "ResEx-recon support %d\n" % i)
        if proj.do_bg_fitting:
            proj.radavg[self.size//2].tofile("%s-radavg/%.4d.raw" % (self.output_prefix, i))

    def save_output(self, phas, proj, p1, p2):
        rvsizes = np.array([1., 1., 1.], dtype='f4')
        fvsizes = rvsizes * -1.

        self.save_as_map("%s-last.ccp4" % self.output_prefix, phas.iterate, rvsizes, "ResEx-recon Last iteration\n")
        self.save_as_map("%s-pf.ccp4" % self.output_prefix, p1[0], rvsizes, "ResEx-recon average_p1\n")
        self.save_as_map("%s-pd.ccp4" % self.output_prefix, p2[0], rvsizes, "ResEx-recon average_p2\n")

        if proj.do_bg_fitting:
            self.save_as_map("%s-bg.ccp4" % self.output_prefix, phas.p2[1], fvsizes, "ResEx-recon average_p2\n")
            proj.radavg[:self.size//2].tofile("%s-radavg.raw" % self.output_prefix)

        if proj.do_local_variation:
            self.save_as_map("%s-supp.supp" % self.output_prefix, proj.support, rvsizes, "ResEx-recon Refined support\n")

    @staticmethod
    def save_as_map(fname, vol, vsizes, label):
        if CUDA:
            numpy_vol = np.asnumpy(vol)
            numpy_vsizes = np.asnumpy(vsizes)
        else:
            numpy_vol = vol
            numpy_vsizes = vsizes
        mrc = mrcfile.new(fname, overwrite=True, data=numpy_vol)
        mrc.header.cella.x = numpy_vsizes[0]
        mrc.header.cella.y = numpy_vsizes[1]
        mrc.header.cella.z = numpy_vsizes[2]
        num_lines = int(np.ceil(len(label) / 80.))
        for i in range(num_lines):
            mrc.header['label'][i] = label[i*80:(i+1)*80]
        mrc.close()
