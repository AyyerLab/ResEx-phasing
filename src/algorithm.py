import os
import csv
import numpy as np
import input
import volume
import fft

class Algorithm():
    def __init__(self):
        self._alg_list = ['ER', 'DM', 'HIO', 'mod-DM', 'RAAR']
        pass
    
    def proj_fourier(self, in_arr, out_arr):
        '''Fourier-space projection 
           Makes model consistent with data. Modifies Fourier magnitudes and keeps
           phases unchanged. Applies symmetry depending on the specified point group.
           Can be applied in-place (out = in)
        '''
        vol = self.vol
        
        if self.do_bg_fitting:
            self.fft.rdensity[:] = in_arr[0]
            out_arr[1] = in_arr[1]
            
            self.fft.forward()
            self.volume.symmetrize_incoherent(self.fft.fdensity, self.exp_mag, out_arr[1])
            self.input.match_bragg(self.fft.fdensity, 0.)
            #self.volume.rotational_blur(self.exp_mag, self.exp_mag, self.quat)
            
            sel = (self.input.obs_mag > 0)
            ratio = self.input.obs_mag[sel] / self.exp_mag[sel]
            self.fft.fdensity[sel] *= ratio
            out_arr[1][sel] = in_arr[1][sel] * ratio
            
            sel = (self.input.obs_mag == 0.)
            self.fft.fdensity[sel] = 0
            out_arr[1][sel] = 0
            
            out_arr[1][self.input.obs_mag < 0] = 0
        else:
            self.fft.rdensity[:] = in_arr
            self.fft.forward()
            self.volume.symmetrize_incoherent(self.fft.fdensity, self.exp_mag)
            self.input.match_bragg(self.fft.fdensity, 0.)
             
            sel = (self.input.obs_mag > 0)
            self.fft.fdensity[sel] *= self.input.obs_mag[sel] / self.exp_mag[sel]
            self.fft.fdensity[self.input.obs_mag == 0.] = 0
            
        self.fft.inverse()
        out_arr[0] = np.real(self.fft.rdensity) / vol

    def proj_direct(self, in_arr, out_arr):
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
        vol = self.vol
        
        if self.do_local_variation:
            self.input.update_support(in_arr, 2)
        
        if self.do_histogram:
            self.input.match_histogram(in_arr, out_arr)
        else:
            out_arr[0] = in_arr[0] * self.input.support
            
            if self.do_positivity:
                out_arr[0][out_arr[0] < 0] = 0
        
        if self.do_bg_fitting:
            self.volume.radial_average(in_arr[1], out_arr[1])

    def DM_algorithm(self):
        '''Difference Map algorithm
           Update rule (LaTeX syntax)
           x_{n+1} = x_n + \beta \left\{P_D\left[\left(1+\frac{1}{\beta}\right)P_F(x_n) - \frac{x_n}{\beta}\right] - P_F\left[\left(1-\frac{1}{\beta}\right)P_D(x_n) + \frac{x_n}{\beta}\right]\right\}
           Same as all other algorithms for beta = 1
        '''
        self.proj_fourier(self.iterate, self.p1)
        if self.beta != 1.:
            self.proj_direct(self.iterate, self.p2)
        
        self.r1 = (1. + 1./self.beta) * self.p1 - self.iterate / self.beta
        if self.beta != 1.:
            self.r2 = (1. - 1./self.beta) * self.p2 + self.iterate / self.beta
        
        self.proj_direct(self.r1, self.p2)
        if self.beta != 1.:
            self.proj_fourier(self.r2, self.p1)
        
        diff = self.beta * (self.p2 - self.p1)
        self.iterate += diff
        return np.linalg.norm(diff[0])

    def mod_DM_algorithm(self):
        '''Modified Difference Map algorithm
           Update rule (LaTeX syntax)
           x_n' = \beta x_n + (1-\beta) P_F(x_n)
           x_{n+1} = x_n' + P_F\left[2 P_D(x_n') - x_n'\right] - P_D(x_n')
           Same as all other algorithms for beta = 1
        '''
        if self.beta != 1.:
            self.proj_fourier(self.iterate, self.p1)
            self.iterate = self.beta*self.iterate + (1. - self.beta) * self.p1
        
        self.proj_direct(self.iterate, self.p2)
        self.r1 = 2. * self.p2 - self.iterate
        self.proj_fourier(self.r1, self.p1)
        
        diff = self.p1 - self.p2
        self.iterate += diff
        return np.linalg.norm(diff[0])

    def RAAR_algorithm(self):
        '''RAAR algorithm
           Update rule (LaTeX syntax)
           x_{n+1} = \beta \left\{x_n + P_D\left[2 P_F(x_n) - x_n\right] - P_F(x_n)\right\} + (1-\beta) P_F(x_n)
         *
           If one does not assume P_D is linear,
           x_{n+1} = \beta \left\{x_n + P_D\left[2 P_F(x_n)\right] + P_D\left[-x_n\right] - P_F(x_n)\right\} + (1-\beta) P_F(x_n)
           
           Same as all other algorithms for beta = 1
        '''
        self.proj_fourier(self.iterate, self.p1)
        
        self.r1 = 2. * self.p1
        self.proj_direct(self.r1, self.r2)
        
        self.r1 = - self.p1
        self.proj_direct(self.r1, self.p2)
        
        diff = (self.beta - 1.) * self.iterate + self.beta * (self.r2 + self.p2) + (1. - 2. * self.beta) * self.p1
        self.iterate += diff
        return np.linalg.norm(diff[0])

    def HIO_algorithm(self):
        '''HIO algorithm
           Update rule (LaTeX syntax)
           x_{n+1} = x_n + \beta \left\{P_D\left[\left(1+\frac{1}{\beta}\right) P_F(x_n) - \frac{x_n}{\beta}\right] - P_F(x_n)\right\}
           Same as all other algorithms for beta = 1
        '''
        self.proj_fourier(self.iterate, self.p1)
        self.r1 = (1. + 1./self.beta) * self.p1 - self.iterate / self.beta
        self.proj_direct(self.r1, self.p2)
        
        diff = self.beta * (self.p2 - self.p1)
        self.iterate += diff
        return np.linalg.norm(diff[0])

    def ER_algorithm(self):
        '''Error Reduction algorithm
           Update rule (LaTeX style)
           x_{n+1} = P_D[P_F(x_n)]
           Obviously different from others. Use only in averaging phase.
        '''
        self.proj_fourier(self.iterate, self.p1)
        self.proj_direct(self.p1, self.p2)
        
        self.iterate = self.p2
        diff = self.p2 - self.p1
        return np.linalg.norm(diff[0])

    #============================================================

    def make_recon_folders(self):
        os.makedirs("%s-slices" % self.output_prefix, 0o755, exist_ok=True)
        os.makedirs("%s-fslices" % self.output_prefix, 0o755, exist_ok=True)
        if self.do_bg_fitting:
            os.makedirs("%s-radavg" % self.output_prefix, 0o755, exist_ok=True)
        if self.do_local_variation:
            os.makedirs("%s-support" % self.output_prefix, 0o755, exist_ok=True)

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
        
        return error

    def save_current(self, i, t1, t2, error):
        with open("%s-log.dat" % self.output_prefix, "a") as f:
            f.write("%.4d\t%.2f s\t%f\n" % (i, t2 - t1, error))
        
        self.volume.dump_slices(self.p1, "%s-slices/%.4d.ccp4" % (self.output_prefix, i), "ResEx-recon p1 %d\n" % i)
        self.volume.dump_slices(self.exp_mag, "%s-fslices/%.4d.ccp4" % (self.output_prefix, i), "ResEx-recon exp_mag %d\n" % i, is_fourier=True)
        
        if self.do_local_variation:
            self.volume.dump_slices(self.input.support, "%s-support/%.4d.ccp4" % (self.output_prefix, i), "ResEx-recon support %d\n" % i)
        if self.do_bg_fitting:
            self.volume.radavg[self.size//2].tofile("%s-radavg/%.4d.raw" % (self.output_prefix, i))
        
        if i == self.num_iter:
            vsizes = np.array([1., 1., 1.], dtype='f4')
            self.volume.save_as_map("%s-last.ccp4" % self.output_prefix, self.iterate, vsizes, "ResEx-recon Last iteration %d\n" % i)

    def save_output(self):
        rvsizes = np.array([1., 1., 1.], dtype='f4')
        fvsizes = rvsizes * -1.
        
        self.volume.save_as_map("%s-pf.ccp4" % self.output_prefix, self.average_p1, rvsizes, "ResEx-recon average_p1\n")
        self.volume.save_as_map("%s-pd.ccp4" % self.output_prefix, self.average_p2, rvsizes, "ResEx-recon average_p2\n")
        
        if self.do_bg_fitting:
            self.volume.save_as_map("%s-bg.ccp4" % self.output_prefix, self.p2[1], fvsizes, "ResEx-recon average_p2\n")
            self.volume.radavg[:self.size//2].tofile("%s-radavg.raw" % self.output_prefix)
        
        if self.do_local_variation:
            self.volume.save_as_map("%s-supp.supp" % self.output_prefix, self.input.support, rvsizes, "ResEx-recon Refined support\n")

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

    def allocate_memory(self):
        nmodels = 2 if self.do_bg_fitting else 1
        mshape = (nmodels, self.size, self.size, self.size)
        fshape = (self.size, self.size, self.size)
        
        self.iterate = np.empty(mshape, dtype='f4')
        self.exp_mag = np.empty(fshape, dtype='f4')
        
        self.p1 = np.empty(mshape, dtype='f4')
        self.p2 = np.empty(mshape, dtype='f4')
        self.average_p1 = np.zeros(mshape, dtype='f4')
        self.average_p2 = np.zeros(mshape, dtype='f4')
        self.r1 = np.empty(mshape, dtype='f4')
        if self.beta != 1.:
            self.r2 = np.empty(mshape, dtype='f4')

    def calc_prtf(self, num_bins):
        model1 = self.average_p1
        model2 = self.average_p2
        c = self.size // 2
        
        # FFT average_p2 model
        self.fft.rdensity[:] = model2[0]
        self.fft.forward()
        
        # Calculate exp_mag for average_p2 model
        if self.do_bg_fitting:
            self.volume.symmetrize_incoherent(self.fft.fdensity, self.exp_mag, model2[1])
        else:
            self.volume.symmetrize_incoherent(self.fft.fdensity, self.exp_mag, None)
        
        # Calculate PRTF by comparing with obs_mag
        sel = (self.input.obs_mag > 0)
        prtf = np.ones(self.exp_mag.shape, dtype='f4')*-1
        prtf[sel] = self.exp_mag[sel] / self.input.obs_mag[sel]
        self.p2 = np.fft.fftshift(np.abs(self.fft.fdensity)**2)
        
        bin_size = float(c) / num_bins
        x, y, z = np.indices(prtf.shape)
        x -= c
        y -= c
        z -= c
        intrad = np.round(np.sqrt(x*x + y*y + z*z) / bin_size).astype('i4')

        sel = sel & (intrad < num_bins)
        bin_count = np.zeros(num_bins)
        prtf_avg = np.zeros(num_bins)
        np.add.at(bin_count, intrad[sel], 1)
        np.add.at(prtf_avg, intrad[sel], prtf[sel])
        prtf_avg[bin_count>0] /= bin_count[bin_count>0]
        
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
            self.fft.fdensity[sel] /= prtf[sel]
            self.fft.inverse()
            model2[0] = np.real(self.fft.rdensity) / vol
            
            print("Normalizing p1 model by PRTF")
            self.fft.rdensity[:] = model1[0]
            self.fft.forward()
            self.fft.fdensity[sel] /= prtf[sel]
            self.fft.inverse()
            model1[0] = np.real(self.fft.rdensity) / vol
