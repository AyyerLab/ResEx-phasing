import os
import configparser
from multiprocessing import cpu_count
try:
    import cupy as np
    import cupy.fft as fft
    CUDA = True
    PYFFTW = False
    print('Using cupy (v%s)' % np.__version__)
except ImportError:
    CUDA = False
    import numpy as np
    try:
        import pyfftw
        import pyfftw.interfaces.numpy_fft as fft
        PYFFTW = True
        print('Using PyFFTW (v%s)' % pyfftw.__version__)
    except ImportError:
        import numpy.fft as fft
        PYFFTW = False
        print('Using numpy FFT')
import projection
import fileio

class Phaser():
    def __init__(self, config_fname, fixed_seed=False):
        self._alg_list = ['ER', 'DM', 'HIO', 'mod-DM', 'RAAR']
        self.iterate = None
        self.p1 = self.p2 = None
        self.r1 = self.r2 = None

        config_dir = os.path.dirname(config_fname)
        config = configparser.ConfigParser()
        config.read(config_fname)

        size = config.getint('parameters', 'size')
        bragg_qmax = config.getfloat('parameters', 'bragg_qmax', fallback=0.)
        scale_factor = config.getfloat('parameters', 'scale_factor', fallback=1.)
        num_threads = config.getfloat('parameters', 'num_threads', fallback=cpu_count())
        point_group = config.get('parameters', 'point_group', fallback='1')

        if point_group not in ('1', '222', '4'):
            raise ValueError("Only '1', '4' and '222' point_group values supported currently")

        self.proj = projection.Projection(size, point_group, num_threads=num_threads)
        self.io = fileio.IO(size)

        intens_fname = os.path.join(config_dir, config.get('files', 'intens_fname', fallback=''))
        bragg_fname = os.path.join(config_dir, config.get('files', 'bragg_fname', fallback=''))
        input_fname = os.path.join(config_dir, config.get('files', 'input_fname', fallback=''))
        inputbg_fname = os.path.join(config_dir, config.get('files', 'inputbg_fname', fallback=''))
        support_fname = os.path.join(config_dir, config.get('files', 'support_fname', fallback=''))
        self.io.output_prefix = os.path.join(config_dir, config.get('files', 'output_prefix', fallback='data/output'))

        algorithm_string = config.get('algorithm', 'algorithm')
        avg_algorithm_string = config.get('algorithm', 'avg_algorithm', fallback=None)
        self.beta = config.getfloat('algorithm', 'beta', fallback=1.)
        self.proj.do_bg_fitting = config.getboolean('algorithm', 'bg_fitting', fallback=False)
        self.proj.do_blurring = config.getboolean('algorithm', 'blurring', fallback=False)
        self.proj.do_histogram = config.getboolean('algorithm', 'histogram', fallback=False)
        self.proj.do_local_variation = config.getboolean('algorithm', 'local_variation', fallback=False)
        self.proj.do_positivity = config.getboolean('algorithm', 'positivity', fallback=False)
        self.proj.do_normalize_prtf = config.getboolean('algorithm', 'normalize_prtf', fallback=False)
        #quat_fname = config.get('algorithm', 'quat_fname', fallback=None)
        #num_div = config.getint('algorithm', 'num_div', fallback=-1)
        hist_fname = config.get('algorithm', 'hist_fname', fallback='')
        #sigma = config.getfloat('algorithm', 'sigma_deg', fallback=0.)

        self.size = size

        self.parse_algorithm_strings(algorithm_string, avg_algorithm_string)
        self.allocate_memory()
        self.io.parse_intens(self.proj, intens_fname, scale_factor, self.proj.do_bg_fitting)
        if CUDA: np.get_default_memory_pool().free_all_blocks()
        self.io.parse_bragg(self.proj, bragg_fname, bragg_qmax)
        if CUDA: np.get_default_memory_pool().free_all_blocks()
        self.io.parse_support(self.proj, support_fname)
        if CUDA: np.get_default_memory_pool().free_all_blocks()
        self.io.init_iterate(self.proj, self.iterate, input_fname, inputbg_fname, self.proj.do_bg_fitting, fixed_seed)
        if CUDA: np.get_default_memory_pool().free_all_blocks()
        if self.proj.do_bg_fitting:
            self.proj.init_radavg()
        if self.proj.do_histogram:
            self.io.parse_histogram(hist_fname)
        if self.proj.do_blurring:
            print('Blurring currently not implemented')
            self.proj.do_blurring = False
        self.io.make_recon_folders(self.proj.do_bg_fitting, self.proj.do_local_variation)

        with open('%s-log.dat' % self.io.output_prefix, 'w') as f:
            f.write("Resolution extension iterative phasing\n")
            f.write("Data: %s %s\n" % (bragg_fname, intens_fname))
            f.write("Support: %s (%ld)\n" % (support_fname, self.proj.num_supp))
            f.write("Algorithm: %s with beta = %.2f\n" % (algorithm_string, self.beta))
            f.write("Averaging algorithm: %s\n" % avg_algorithm_string)
            if self.proj.do_positivity:
                f.write("Assuming electron density is positive\n")
            if self.proj.do_histogram:
                f.write("Applying histogram constraint: %s\n" % hist_fname)
            if self.proj.do_local_variation:
                f.write("Updating support using local variation\n")
            if self.proj.do_bg_fitting:
                f.write("Fitting spherically symmetric background\n")
            #if self.proj.do_blurring:
            #    f.write("Rotationally blurring model with %d orientations\n" % quat.num_rot)
            if self.proj.do_normalize_prtf:
                f.write("Normalizing output by PRTF\n")
            f.write("Output prefix: %s\n" % self.io.output_prefix)
            f.write("-------------------------\n")
            f.write("iter  time (s)  error\n")
            f.write("-------------------------\n")
        if CUDA: np.get_default_memory_pool().free_all_blocks()

    def DM_algorithm(self): # pylint: disable=invalid-name
        r'''Difference Map algorithm
           Update rule (LaTeX syntax)
           x_{n+1} = x_n + \beta \left\{P_D\left[\left(1+\frac{1}{\beta}\right)P_F(x_n) - \frac{x_n}{\beta}\right] - P_F\left[\left(1-\frac{1}{\beta}\right)P_D(x_n) + \frac{x_n}{\beta}\right]\right\}
           Same as all other algorithms for beta = 1
        '''
        self.proj.fourier(self.iterate, self.p1)
        if self.beta != 1.:
            self.proj.direct(self.iterate, self.p2)

        self.r1 = (1. + 1./self.beta) * self.p1 - self.iterate / self.beta
        if self.beta != 1.:
            self.r2 = (1. - 1./self.beta) * self.p2 + self.iterate / self.beta

        self.proj.direct(self.r1, self.p2)
        if self.beta != 1.:
            self.proj.fourier(self.r2, self.p1)

        diff = self.beta * (self.p2 - self.p1)
        self.iterate += diff
        return np.linalg.norm(diff[0]) / np.sqrt(diff[0].size)

    def mod_DM_algorithm(self): # pylint: disable=invalid-name
        r'''Modified Difference Map algorithm
           Update rule (LaTeX syntax)
           x_n' = \beta x_n + (1-\beta) P_F(x_n)
           x_{n+1} = x_n' + P_F\left[2 P_D(x_n') - x_n'\right] - P_D(x_n')
           Same as all other algorithms for beta = 1
        '''
        if self.beta != 1.:
            self.proj.fourier(self.iterate, self.p1)
            self.iterate = self.beta*self.iterate + (1. - self.beta) * self.p1

        self.proj.direct(self.iterate, self.p2)
        self.r1 = 2. * self.p2 - self.iterate
        self.proj.fourier(self.r1, self.p1)

        diff = self.p1 - self.p2
        self.iterate += diff
        return np.linalg.norm(diff[0]) / np.sqrt(diff[0].size)

    def RAAR_algorithm(self): # pylint: disable=invalid-name
        r'''RAAR algorithm
           Update rule (LaTeX syntax)
           x_{n+1} = \beta \left\{x_n + P_D\left[2 P_F(x_n) - x_n\right] - P_F(x_n)\right\} + (1-\beta) P_F(x_n)
         *
           If one does not assume P_D is linear,
           x_{n+1} = \beta \left\{x_n + P_D\left[2 P_F(x_n)\right] + P_D\left[-x_n\right] - P_F(x_n)\right\} + (1-\beta) P_F(x_n)

           Same as all other algorithms for beta = 1
        '''
        self.proj.fourier(self.iterate, self.p1)

        self.r1 = 2. * self.p1
        self.proj.direct(self.r1, self.r2)

        self.r1 = - self.p1
        self.proj.direct(self.r1, self.p2)

        diff = (self.beta - 1.) * self.iterate + self.beta * (self.r2 + self.p2) + (1. - 2. * self.beta) * self.p1
        self.iterate += diff
        return np.linalg.norm(diff[0]) / np.sqrt(diff[0].size)

    def HIO_algorithm(self): # pylint: disable=invalid-name
        r'''HIO algorithm
           Update rule (LaTeX syntax)
           x_{n+1} = x_n + \beta \left\{P_D\left[\left(1+\frac{1}{\beta}\right) P_F(x_n) - \frac{x_n}{\beta}\right] - P_F(x_n)\right\}
           Same as all other algorithms for beta = 1
        '''
        self.proj.fourier(self.iterate, self.p1)
        self.r1 = (1. + 1./self.beta) * self.p1 - self.iterate / self.beta
        self.proj.direct(self.r1, self.p2)

        diff = self.beta * (self.p2 - self.p1)
        self.iterate += diff
        return np.linalg.norm(diff[0]) / np.sqrt(diff[0].size)

    def ER_algorithm(self): # pylint: disable=invalid-name
        r'''Error Reduction algorithm
           Update rule (LaTeX style)
           x_{n+1} = P_D[P_F(x_n)]
           Obviously different from others. Use only in averaging phase.
        '''
        self.proj.fourier(self.iterate, self.p1)
        self.proj.direct(self.p1, self.p2)

        self.iterate = self.p2
        diff = self.p2 - self.p1
        return np.linalg.norm(diff[0]) / np.sqrt(diff[0].size)

    #============================================================

    def allocate_memory(self):
        nmodels = 2 if self.proj.do_bg_fitting else 1
        mshape = (nmodels, self.size, self.size, self.size)

        self.iterate = np.empty(mshape, dtype='f4')
        self.p1 = np.empty(mshape, dtype='f4')
        self.p2 = np.empty(mshape, dtype='f4')
        self.r1 = np.empty(mshape, dtype='f4')
        if self.beta != 1.:
            self.r2 = np.empty(mshape, dtype='f4')

    def run_iteration(self, i):
        if i <= self.num_iter:
            algo = self.algorithms[i-1]
        else:
            algo = self.avg_algorithms[i-self.num_iter-1]

        if algo == "DM":
            error = self.DM_algorithm()
        elif algo == "HIO":
            error = self.HIO_algorithm()
        elif algo == "RAAR":
            error = self.RAAR_algorithm()
        elif algo == "mod-DM":
            error = self.mod_DM_algorithm()
        elif algo == "ER":
            error = self.ER_algorithm()
        else:
            print("Could not understand algorithm name:", algo)
            error = -1.

        if PYFFTW and i == self.num_iter + self.num_avg_iter:
            self.proj.export_wisdom()
        if CUDA: np.get_default_memory_pool().free_all_blocks()
        return error

    def parse_algorithm_strings(self, alg_string, avg_string):
        '''Generates list of algorithms given algorithm and avg_algorithm strings.
        '''
        self.num_iter = 0
        self.num_avg_iter = 0

        tokens = alg_string.split()
        nums = [int(s) for s in tokens[::2]]
        names = tokens[1::2]
        tests = np.array([(n in self._alg_list) for n in names])
        if not tests.all():
            raise ValueError('Could not parse algorithm string')
        self.num_iter = sum(nums)
        print("Total number of normal iterations = %d" % self.num_iter)
        self.algorithms = []
        for num, name in zip(nums, names):
            self.algorithms += [name]*num

        tokens = avg_string.split()
        nums = [int(s) for s in tokens[::2]]
        names = tokens[1::2]
        tests = np.array([(n in self._alg_list) for n in names])
        if not tests.all():
            raise ValueError('Could not parse algorithm string')
        self.num_avg_iter = sum(nums)
        print("Total number of averaging iterations = %d" % self.num_avg_iter)
        self.avg_algorithms = []
        for num, name in zip(nums, names):
            self.avg_algorithms += [name]*num
