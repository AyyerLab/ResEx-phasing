import os
import multiprocessing
import numpy as np
import pyfftw

class FFT():
    def __init__(self, size, num_threads=None):
        self.size = size
        vol = size*size*size
        self._radsq = None
        
        if num_threads is None:
            self.num_threads = multiprocessing.cpu_count()
        else:
            self.num_threads = num_threads
        self.rdensity = pyfftw.empty_aligned((size, size, size), dtype='c8')
        self.fdensity = pyfftw.empty_aligned((size, size, size), dtype='c8')
        self.wisdom_fname = "data/wisdom_%ld_%d" % (size, self.num_threads)

    def create_plans(self):
        if not os.path.exists(os.path.dirname(self.wisdom_fname)):
            self.wisdom_fname = os.path.basename(self.wisdom_fname)
        if os.path.isfile(self.wisdom_fname):
            with open(self.wisdom_fname, 'r') as f:
                wisdom_string = f.read()
            wlist = wisdom_string.split('\n')
            wisdom = []
            while True:
                try:
                    ind = wlist.index(')')
                    wisdom.append(bytes(''.join(wlist[:ind+1]), 'utf-8'))
                    wlist = wlist[ind+1:]
                except ValueError:
                    break
            wisdom = tuple(wisdom)
            pyfftw.import_wisdom(wisdom)
            self.forward_plan = pyfftw.FFTW(self.rdensity, self.fdensity, axes=(0,1,2), direction='FFTW_FORWARD', flags=['FFTW_MEASURE'], threads=self.num_threads)
            self.inverse_plan = pyfftw.FFTW(self.fdensity, self.rdensity, axes=(0,1,2), direction='FFTW_BACKWARD', flags=['FFTW_MEASURE'], threads=self.num_threads)
        else:
            print('Measuring plans...')
            self.forward_plan = pyfftw.FFTW(self.rdensity, self.fdensity, axes=(0,1,2), direction='FFTW_FORWARD', flags=['FFTW_MEASURE'], threads=self.num_threads)
            self.inverse_plan = pyfftw.FFTW(self.fdensity, self.rdensity, axes=(0,1,2), direction='FFTW_BACKWARD', flags=['FFTW_MEASURE'], threads=self.num_threads)

            wisdom = pyfftw.export_wisdom()
            with open(self.wisdom_fname, "w") as f:
                for s in wisdom:
                    f.write(s.decode('utf-8'))
            print('Created plans. Saved to', self.wisdom_fname)

    def gaussian_blur(self, model, blur):
        size = self.size
        vol = size*size*size
        fblur = size / (2. * np.pi * blur)
        
        self.rdensity[:] = model
        self.forward()
        if self._radsq is None:
            self._calc_radsq()
        self.fdensity *= np.exp(-self._radsq / 2. / fblur**2)
        self.inverse()
        model[:] = np.real(self.rdensity) / vol

    def apply_shrinkwrap(self, model, blur, threshold, support):
        self.gaussian_blur(model, blur)
        sel = model > threshold
        support = np.where(sel, 1, 0)
        return support.sum()

    def forward(self):
        if self.forward_plan is None:
            self.create_plans()
        self.forward_plan()

    def inverse(self):
        if self.forward_plan is None:
            self.create_plans()
        self.inverse_plan()

    def _calc_radsq(self):
        c = self.size // 2
        x, y, z = np.indices((self.size, self.size, self.size))
        x -= c
        y -= c
        z -= c
        self._radsq = np.fft.ifftshift(x*x + y*y + z*z).astype('f4')
