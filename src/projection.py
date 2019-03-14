# pylint: disable=unsubscriptable-object
import os
import re
import multiprocessing
try:
    import cupy as np
    import cupy.fft as fft
    from cupyx.scipy import ndimage
    CUDA = True
    PYFFTW = False
except ImportError:
    CUDA = False
    import numpy as np
    from scipy import ndimage
    try:
        import pyfftw
        import pyfftw.interfaces.numpy_fft as fft
        PYFFTW = True
    except ImportError:
        import numpy.fft as fft
        PYFFTW = False

class Projection():
    def __init__(self, size, point_group, num_threads=None):
        self.size = size
        self.point_group = point_group
        self.num_threads = num_threads
        print("Symmetrizing with point group: %s" % self.point_group)

        self.intrad = None
        self.radavg = None
        self.radcount = None
        self.obs_radavg = None

        self.obs_mag = None
        self.bragg_calc = None
        self.bragg_mask = None
        self.support = None
        self.num_supp = 0
        self.supp_loc = None
        self.inverse_cdf = None
        self.local_variation = None
        self.sigma = 0.
        if PYFFTW:
            if self.num_threads is None:
                self.num_threads = multiprocessing.cpu_count()
            self.fft_kwargs = {'planner_effort': 'FFTW_MEASURE', 'threads': self.num_threads}
            pyfftw.interfaces.cache.enable()
            if os.path.exists('data/'):
                wisdom_fname = 'data/wisdom_%d_%d' % (self.size, self.num_threads)
            else:
                wisdom_fname = 'wisdom_%d_%d' % (self.size, self.num_threads)
            if os.path.isfile(wisdom_fname):
                print('Importing PyFFTW wisdom from', wisdom_fname)
                with open(wisdom_fname, 'r') as fptr:
                    wisdom = tuple([(s+')').encode() for s in re.split('\n\)\n', fptr.read())[:-1]])
                pyfftw.import_wisdom(wisdom)
        else:
            self.fft_kwargs = {}

        self.exp_mag = np.empty((size, size, size), dtype='f4')

        self.beta = 1.
        self.do_bg_fitting = False
        self.do_blurring = False
        self.do_histogram = False
        self.do_local_variation = False
        self.do_positivity = False
        self.do_normalize_prtf = False

    def fourier(self, in_arr, out_arr):
        '''Fourier-space projection
           Makes model consistent with data. Modifies Fourier magnitudes and keeps
           phases unchanged. Applies symmetry depending on the specified point group.
           Can be applied in-place (out = in)
        '''
        if self.do_bg_fitting:
            out_arr[1] = in_arr[1]
            if CUDA:
                garr = np.array(in_arr[0])
                fdensity = fft.fftn(garr, **self.fft_kwargs)
            else:
                fdensity = fft.fftn(in_arr[0], **self.fft_kwargs).astype('c8')
            self.symmetrize_incoherent(fdensity, self.exp_mag, out_arr[1])
            self.match_bragg(fdensity, 0.)
            if self.do_blurring:
                self.rotational_blur(self.exp_mag, self.exp_mag)

            sel = (self.obs_mag > 0)
            ratio = self.obs_mag[sel] / self.exp_mag[sel]
            fdensity[sel] *= ratio
            out_arr[1][sel] = in_arr[1][sel] * ratio

            sel = (self.obs_mag == 0.)
            fdensity[sel] = 0
            out_arr[1][sel] = 0

            out_arr[1][self.obs_mag < 0] = 0
        else:
            if CUDA:
                garr = np.array(in_arr[0])
                np.get_default_memory_pool().free_all_blocks()
                fdensity = fft.fftn(garr, axes=(0, 1, 2), **self.fft_kwargs)
                np.get_default_memory_pool().free_all_blocks()
            else:
                fdensity = fft.fftn(in_arr[0], **self.fft_kwargs).astype('c8')
            self.symmetrize_incoherent(fdensity, self.exp_mag)
            self.match_bragg(fdensity, 0.)
            if self.do_blurring:
                self.rotational_blur(self.exp_mag, self.exp_mag)

            sel = (self.obs_mag > 0)
            fdensity[sel] *= self.obs_mag[sel] / self.exp_mag[sel]
            fdensity[self.obs_mag == 0.] = 0

        #out_arr[0] = np.real(fft.ifftn(fdensity, **self.fft_kwargs)) / self.vol
        out_arr[0] = np.real(fft.ifftn(fdensity, **self.fft_kwargs))

    def direct(self, in_arr, out_arr):
        '''Direct-space projection
           Applies the finite support constraint in direct/real-space
           In addition one can apply the following additional constraints by
           setting the following global variables.
               do_bg_fitting - Azimuthally average background part of iterate
               do_positivity - Set negative values inside support to zero
               do_histogram - Project values inside support such that they match a
                              target histogram.
               do_local_variation - Calculate local variation and update support keeping
                                    support size the same
           Can be applied in-place (out = in)
        '''
        if self.do_local_variation:
            self.update_support(in_arr, 2, num_vox=self.num_supp)

        if self.do_histogram:
            self.match_histogram(in_arr, out_arr)
        else:
            out_arr[0] = in_arr[0] * self.support

            if self.do_positivity:
                out_arr[0][out_arr[0] < 0] = 0

        if self.do_bg_fitting:
            self.radial_average(in_arr[1], out_arr[1])

    def symmetrize_incoherent(self, in_arr, out_arr, bg=None, shifted=True):
        '''Symmetrize intensity incoherently according to given point group
           If shifted=True, array is assumed to have q=0 at (0,0,0) instead of in the center of the array
        '''
        in_intens = np.absolute(in_arr)**2
        if shifted:
            in_intens[:] = np.fft.fftshift(in_intens)

        if self.point_group == '222':
            out_arr[:] = 0.25 * (in_intens + in_intens[::-1] + in_intens[:, ::-1] + in_intens[:, :, ::-1])
            if bg is not None:
                bg_intens = np.abs(bg)**2
                bg = 0.25 * (bg_intens + bg_intens[::-1] + bg_intens[:, ::-1] + bg_intens[:, :, ::-1])
                out_arr[:] = np.sqrt(out_arr + bg)
            else:
                out_arr[:] = np.sqrt(out_arr)
        elif self.point_group == "4":
            out_arr[:] = 0.25 * (in_intens +
                                 np.rot90(in_intens, 1, axes=(1, 2)) +
                                 np.rot90(in_intens, 2, axes=(1, 2)) +
                                 np.rot90(in_intens, 3, axes=(1, 2)))
            if bg is not None:
                bg_intens = np.abs(bg)**2
                bg = 0.25 * (bg_intens +
                             np.rot90(bg_intens, 1, axes=(1, 2)) +
                             np.rot90(bg_intens, 2, axes=(1, 2)) +
                             np.rot90(bg_intens, 3, axes=(1, 2)))
                out_arr[:] = np.sqrt(out_arr + bg)
            else:
                out_arr[:] = np.sqrt(out_arr)
        elif self.point_group == "1":
            if bg is not None:
                out_arr[:] = np.sqrt(in_intens + bg**2)
            else:
                out_arr[:] = in_arr
        else:
            raise ValueError("Unrecognized point group: %s\n" % self.point_group)

        if shifted:
            out_arr[:] = np.fft.ifftshift(out_arr)

    def init_radavg(self):
        '''Radial average initialization
           Calculate bin for each voxel and bin occupancy
           Note that q=0 is at (0,0,0) and not in the center of the array
        '''
        size = self.size
        c = size // 2
        self.radavg = np.zeros(size, dtype='f4')
        self.obs_radavg = np.zeros(size, dtype='f4')
        self.radcount = np.zeros(size, dtype='f4')

        self._calc_intrad()
        self.radcount = np.bincount(self.intrad.ravel())
        #self.radcount = np.array([float((self.intrad==r).sum()) for r in range(int(self.intrad.max()+1))])
        #np.add.at(self.radcount, self.intrad, 1)
        self.radcount[self.radcount == 0] = 1

    def radial_average(self, in_arr, out_arr=None, positive=True):
        '''Radial average calculation
           Using previously calculated bins and bin occupancies
           Note that q=0 is at (0,0,0) and not in the center of the array
           If positive=True, radial average forced to be non-negative
        '''
        self.radavg.fill(0)
        self.radavg = np.array([float(in_arr[self.intrad==r].sum()) for r in range(int(self.intrad.max()+1))])
        #np.add.at(self.radavg, self.intrad, in_arr)
        self.radavg /= self.radcount
        if positive:
            self.radavg[self.radavg < 0] = 0

        if out_arr is not None:
            out_arr[:] = self.radavg[self.intrad]

    @staticmethod
    def _rot_from_quat(q):
        q0 = float(q[0])
        q1 = float(q[1])
        q2 = float(q[2])
        q3 = float(q[3])

        rot = np.empty((3,3))
        rot[0] = np.array([1. - 2.*(q2*q2 + q3*q3), 2.*(q1*q2 + q0*q3), 2.*(q1*q3 - q0*q2)])
        rot[1] = np.array([2.*(q1*q2 - q0*q3), 1. - 2.*(q1*q1 + q3*q3), 2.*(q0*q1 + q2*q3)])
        rot[2] = np.array([2.*(q0*q2 + q1*q3), 2.*(q2*q3 - q0*q1), 1. - 2.*(q1*q1 + q2*q2)])

        return rot

    @staticmethod
    def _rot_from_axang(ang, axis):
        x = float(axis[0])
        y = float(axis[1])
        z = float(axis[2])
        c = float(np.cos(ang))
        s = float(np.sin(ang))
        t = 1. - c
        
        rot = np.empty((3,3))
        rot[0] = np.array([t*x*x + c, t*x*y - z*s, t*x*z + y*s])
        rot[1] = np.array([t*x*y + z*s, t*y*y + c, t*y*z - x*s])
        rot[2] = np.array([t*x*z - y*s, t*y*z + x*s, t*z*z + c])

        return rot

    def quat_blur(self, in_arr, out_arr, quat):
        '''Rotate and average intensity distribution using given quaternions and weights
           Note: q=0 is at (c,c,c) and not (0,0,0)
           Cannot set out = in
        '''
        out_arr[:] = 0
        c = int(in_arr.shape[0]) // 2
        for q in quat:
            rot = self._rot_from_quat(q[:4])
            out_arr += ndimage.affine_transform(in_arr*q[4], rot, order=1, offset=np.array([c]*3) - np.dot(rot, np.array([c]*3)))
        if CUDA: np.get_default_memory_pool().free_all_blocks()

    def rotational_blur(self, in_arr, out_arr, num=100):
        if self.sigma == 0.:
            print('WARNING: sigma = 0. Not doing any blurring')
            return
        rot_arr = np.zeros_like(in_arr)
        c = int(in_arr.shape[0]) // 2
        total_w = 0
        for i in range(num):
            angle = np.random.randn(1, dtype='f4') * self.sigma
            w = np.exp(-angle**2 / 2 / self.sigma**2)
            norm = 0
            while norm < 1:
                axis = np.random.rand(3, dtype='f4') * 2 - 1
                norm = np.linalg.norm(axis)
            axis /= norm
            rot = self._rot_from_axang(angle, axis)
            rot_arr += ndimage.affine_transform(in_arr*w, rot, order=1, offset=np.array([c]*3) - np.dot(rot, np.array([c]*3)))
            total_w += w
        if CUDA: np.get_default_memory_pool().free_all_blocks()
        out_arr[:] = rot_arr / total_w

    def match_histogram(self, in_arr, out_arr):
        '''Histogram matching
           Matches histogram within support volume using inverse_cdf array.
           Can be done in-place
        '''
        supp_val = in_arr.ravel()[self.supp_loc]
        sorter = supp_val.argsort()
        out_arr[:] = 0
        out_arr[0].ravel()[self.supp_loc[sorter]] = self.inverse_cdf

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

    def gaussian_blur(self, vol3d, blur):
        if len(vol3d.shape) != 3:
            raise ValueError('gaussian_blur() needs 3D input array')
        if blur == 0:
            return

        fblur = self.size / (2. * np.pi * blur)
        fdensity = fft.fftn(vol3d, axes=(0,1,2))
        if CUDA: np.get_default_memory_pool().free_all_blocks()
        self._calc_intrad()
        if CUDA: np.get_default_memory_pool().free_all_blocks()
        fdensity *= np.exp(-self.radsq / 2. / fblur**2)
        vol3d[:] = np.real(fft.ifftn(fdensity)) / vol3d.size

    def update_support(self, model, blur, threshold=None, num_vox=None):
        if num_vox is None and threshold is None:
            raise ValueError('Need either threshold or num_vox to update support')
        elif num_vox is not None and threshold is not None:
            raise ValueError('Specify only one of threshold or num_vox to update support')

        if len(model.shape) == 4:
            m3d = np.copy(model[0])
        else:
            m3d = np.copy(model)
        if blur > 0.:
            self.gaussian_blur(m3d, blur)

        if threshold is not None:
            self.support[:] = np.where(m3d > threshold, 1, 0)
        elif num_vox is not None:
            self.support[:] = 0
            self.support.ravel()[m3d.argsort(axis=None)[-num_vox:]] = 1

        return int(self.support.sum())

    @staticmethod
    def positive_mode(model):
        '''Get mode of positive values in volume
        '''
        h = np.histogram(model.ravel, bins=np.linspace(0, model.max(), 100))
        mode = h[1][h[0].argmax()]
        mode_err = h[1][1] - h[1][0]
        print("Mode of positive values in volume = %.3e +- %.3e" % (mode, mode_err))
        return mode

    def _calc_intrad(self):
        if self.intrad is not None:
            return
        size = self.size
        cen = size // 2
        ix, iy, iz = np.indices((size, size, size), dtype='f4')
        ix -= cen
        iy -= cen
        iz -= cen
        self.radsq = np.fft.ifftshift(ix*ix + iy*iy + iz*iz)
        self.intrad = np.sqrt(self.radsq).astype('i4')

    def export_wisdom(self):
        if os.path.exists('data/'):
            wisdom_fname = 'data/wisdom_%d_%d' % (self.size, self.num_threads)
        else:
            wisdom_fname = 'wisdom_%d_%d' % (self.size, self.num_threads)
        if not os.path.isfile(wisdom_fname):
            wisdom = pyfftw.export_wisdom()
            with open(wisdom_fname, 'w') as fptr:
                for s in wisdom:
                    fptr.write(s.decode('utf-8'))
            print('Saved pyFFTW wisdom to', wisdom_fname)

'''
    def calc_prtf(self, num_bins):
        model1 = self.average_p1
        model2 = self.average_p2
        c = self.size // 2

        # FFT average_p2 model
        fdensity = fft.fftn(model2[0], **self.fft_kwargs)

        # Calculate exp_mag for average_p2 model
        if self.do_bg_fitting:
            self.volume.symmetrize_incoherent(fdensity, self.exp_mag, model2[1])
        else:
            self.volume.symmetrize_incoherent(fdensity, self.exp_mag, None)

        # Calculate PRTF by comparing with obs_mag
        sel = (self.input.obs_mag > 0)
        prtf = np.ones(self.exp_mag.shape, dtype='f4')*-1
        prtf[sel] = self.exp_mag[sel] / self.input.obs_mag[sel]
        self.p2 = np.fft.fftshift(np.abs(fdensity)**2)

        bin_size = float(c) / num_bins
        x, y, z = numpy.indices(prtf.shape)
        x -= c
        y -= c
        z -= c
        intrad = numpy.around(numpy.sqrt(x*x + y*y + z*z) / bin_size).astype('i4')

        if CUDA:
            numpy_sel = np.asnumpy(sel) & (intrad < num_bins)
            numpy_prtf = np.asnumpy(prtf)
        else:
            numpy_sel = sel & (intrad < num_bins)
            numpy_prtf = prtf
        bin_count = numpy.zeros(num_bins)
        prtf_avg = numpy.zeros(num_bins)
        numpy.add.at(bin_count, intrad[numpy_sel], 1)
        numpy.add.at(prtf_avg, intrad[numpy_sel], numpy_prtf[numpy_sel])
        prtf_avg[bin_count > 0] /= bin_count[bin_count > 0]

        # Save PRTF
        with open("%s-prtf.dat" % self.output_prefix, 'w') as f:
            w = csv.writer(f, delimiter='\t')
            w.writerow(['Fractional q', 'PRTF'])
            w.writerows(zip(np.linspace(0, 1, num_bins), prtf_avg))

        # Save frecon (intensity from model)
        vsizes = np.ones(3, dtype='f4')*-1
        self.volume.save_as_map("%s-frecon.ccp4" % self.output_prefix, self.p2, vsizes, "ResEx-recon Fourier intensities of reconstruction\n")

        # If needed, normalize models by PRTF
        if self.do_normalize_prtf:
            print("Normalizing p2 model by PRTF")
            fdensity[sel] /= prtf[sel]
            rdensity = fft.ifftn(fdensity, **self.fft_kwargs)
            model2[0] = np.real(rdensity) / rdensity.size

            print("Normalizing p1 model by PRTF")
            fdensity = fft.fftn(model1[0], **self.fft_kwargs)
            fdensity[sel] /= prtf[sel]
            model1[0] = np.real(fft.ifftn(fdensity, **self.fft_kwargs)) / fdensity.size
'''
