import numpy as np
import mapio

class Volume():
    def __init__(self, size, point_group):
        self.size = size
        self.point_group = point_group
        print("Symmetrizing with point group: %s" % self.point_group)

        self.intrad = None
        self.radavg = None
        self.radcount = None
        self.obs_radavg = None

    def symmetrize_incoherent(self, in_arr, out_arr, bg=None, shifted=True):
        '''Symmetrize intensity incoherently according to given point group
           If shifted=True, array is assumed to have q=0 at (0,0,0) instead of in the center of the array
        '''
        if shifted:
            in_arr = np.fft.ifftshift(in_arr)
        if self.point_group == '222':
            in_intens = np.abs(in_arr)**2
            out_arr[:] = 0.25 * (in_intens + in_intens[::-1] + in_intens[:,::-1] + in_intens[:,:,::-1])
            if bg is not None:
                bg_intens = np.abs(bg)**2
                bg = 0.25 * (bg_intens + bg_intens[::-1] + bg_intens[:,::-1] + bg_intens[:,:,::-1])
                out_arr[:] = np.sqrt(out_arr + bg)
            else:
                out_arr[:] = np.sqrt(out_arr)
        elif self.point_group == "4":
            in_intens = np.abs(in_arr)**2
            out_arr[:] = 0.25 * (in_intens +
                                 np.rot90(in_intens, 1, axes=(1,2)) +
                                 np.rot90(in_intens, 2, axes=(1,2)) +
                                 np.rot90(in_intens, 3, axes=(1,2)))
            if bg is not None:
                bg_intens = np.abs(bg)**2
                bg = 0.25 * (bg_intens +
                             np.rot90(bg_intens, 1, axes=(1,2)) +
                             np.rot90(bg_intens, 2, axes=(1,2)) +
                             np.rot90(bg_intens, 3, axes=(1,2)))
                out_arr[:] = np.sqrt(out_arr + bg)
            else:
                out_arr[:] = np.sqrt(out_arr)
        elif self.point_group == "1":
            if bg is not None:
                out_arr[:] = np.sqrt(np.abs(in_arr)**2 + bg**2)
            else:
                out_arr[:] = in_arr
        else:
            raise ValueError("Unrecognized point group: %s\n" % self.point_group)

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
        
        ix, iy, iz = np.indices((size, size, size))
        ix -= c
        iy -= c
        iz -= c
        self.intrad = np.fft.ifftshift(np.sqrt(ix*ix + iy*iy + iz*iz).astype('i4'))
        np.add.at(self.radcount, self.intrad, 1)
        self.radcount[self.radcount==0] = 1

    def radial_average(self, in_arr, out_arr=None, positive=True):
        '''Radial average calculation
           Using previously calculated bins and bin occupancies
           Note that q=0 is at (0,0,0) and not in the center of the array
           If positive=True, radial average forced to be non-negative
        '''
        self.radavg.fill(0)
        np.add.at(self.radavg, intrad, in_arr)
        self.radavg /= self.radcount
        if positive:
            self.radavg[self.radavg<0] = 0
        
        if out_arr is not None:
            out_arr[:] = self.radavg[self.intrad]

    def _gen_rot(rot, q):
        q0 = q[0]
        q1 = q[1]
        q2 = q[2]
        q3 = q[3]
        
        q01 = q0*q1
        q02 = q0*q2
        q03 = q0*q3
        q11 = q1*q1
        q12 = q1*q2
        q13 = q1*q3
        q22 = q2*q2
        q23 = q2*q3
        q33 = q3*q3
        
        rot[0][0] = (1. - 2.*(q22 + q33))
        rot[0][1] = 2.*(q12 + q03)
        rot[0][2] = 2.*(q13 - q02)
        rot[1][0] = 2.*(q12 - q03)
        rot[1][1] = (1. - 2.*(q11 + q33))
        rot[1][2] = 2.*(q01 + q23)
        rot[2][0] = 2.*(q02 + q13)
        rot[2][1] = 2.*(q23 - q01)
        rot[2][2] = (1. - 2.*(q11 + q22))

    def rotational_blur(self, in_arr, out_arr, quat):
        '''Rotate and average intensity distribution using given quaternions and weights
           Note: q=0 is at (0,0,0) and not (c,c,c)
           Can set out = in
        '''
        pass

    @staticmethod
    def accumulate(current_arr, sum_arr):
        sum_arr[:] = sum_arr + current_arr

    @staticmethod
    def dump_slices(vol, fname, label='', is_fourier=False, is_support=False):
        '''Save orthogonal central slices to file
        '''
        vsizes = np.ones(3, dtype='f4')
        if is_fourier:
            temp = np.fft.fftshift(vol)
            vsizes *= -1.
        else:
            temp = vol
        
        size = vol.shape[0]
        cen = size // 2
        if is_support:
            slices = np.empty((3, size, size), dtype='i1')
        else:
            slices = np.empty((3, size, size), dtype='f4')
        slices[0] = temp[cen]
        slices[1] = temp[:,cen]
        slices[2] = temp[:,:,cen]
        
        mapio.save_vol_as_map(fname, slices, vsizes, label)

    @staticmethod
    def positive_mode(model):
        '''Get mode of positive values in volume
        '''
        h = np.histogram(model.ravel, bins=np.linspace(0, model.max(), 100))
        mode = h[1][h[0].argmax()]
        mode_err = h[1][1] - h[1][0]
        print("Mode of positive values in volume = %.3e +- %.3e" % (mode, mode_err))
        return mode
