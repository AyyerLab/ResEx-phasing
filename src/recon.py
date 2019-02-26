import sys
import argparse
import configparser
import time
from multiprocessing import cpu_count
import algorithm
import fft
import input
import volume

def setup(self, config_fname, fixed_seed):
    config = configparser.ConfigParser()
    config.read(config_fname)

    size = config.getint('parameters', 'size')
    bragg_qmax = config.getfloat('parameters', 'bragg_qmax', fallback=0.)
    scale_factor = config.getfloat('parameters', 'scale_factor', fallback=1.)
    num_threads = config.getfloat('parameters', 'num_threads', fallback=cpu_count())
    point_group = config.get('parameters', 'point_group', fallback='1')

    intens_fname = config.get('files', 'intens_fname', fallback=None)
    bragg_fname = config.get('files', 'bragg_fname', fallback=None)
    input_fname = config.get('files', 'input_fname', fallback=None)
    inputbg_fname = config.get('files', 'inputbg_fname', fallback=None)
    support_fname = config.get('files', 'support_fname', fallback=None)
    self.output_prefix = config.get('files', 'output_prefix', fallback='data/output')

    algorithm_string = config.get('algorithm', 'algorithm')
    avg_algorithm_string = config.get('algorithm', 'avg_algorithm', fallback=None)
    self.beta = config.getfloat('algorithm', 'beta', fallback=1.)
    self.do_bg_fitting = config.getboolean('algorithm', 'bg_fitting', fallback=False)
    self.do_blurring = config.getboolean('algorithm', 'blurring', fallback=False)
    self.do_histogram = config.getboolean('algorithm', 'histogram', fallback=False)
    self.do_local_variation = config.getboolean('algorithm', 'local_variation', fallback=False)
    self.do_positivity = config.getboolean('algorithm', 'positivity', fallback=False)
    self.do_normalize_prtf = config.getboolean('algorithm', 'normalize_prtf', fallback=False)
    #quat_fname = config.get('algorithm', 'quat_fname', fallback=None)
    #num_div = config.getint('algorithm', 'num_div', fallback=-1)
    hist_fname = config.get('algorithm', 'hist_fname', fallback=None)
    #sigma = config.getfloat('algorithm', 'sigma_deg', fallback=0.)
    
    if point_group != '1' and point_group != '222' and point_group != '4':
        raise ValueError("Only '1', '4' and '222' point_group values supported currently")
    
    self.size = size
    self.vol = size*size*size
    self.num_vox = 2 * self.vol if self.do_bg_fitting else self.vol
    
    self.fft = fft.FFT(size, num_threads)
    self.fft.create_plans()
    self.input = input.Input(size)
    self.volume = volume.Volume(size, point_group)

    self.parse_algorithm_strings(algorithm_string, avg_algorithm_string)
    self.allocate_memory()
    self.input.parse_intens(intens_fname, scale_factor, self.do_bg_fitting)
    self.input.parse_bragg(bragg_fname, bragg_qmax)
    self.input.parse_support(support_fname)
    if self.do_histogram:
        self.input.parse_histogram(hist_fname)
    if self.do_blurring:
        print('Blurring currently not implemented')
        self.do_blurring = False

    self.input.init_iterate(self.iterate, input_fname, inputbg_fname, self.do_bg_fitting, fixed_seed)
    if self.do_bg_fitting:
        self.volume.init_radavg()
    
    with open('%s-log.dat' % self.output_prefix, 'w') as f:
        f.write("Resolution extension iterative phasing\n")
        f.write("Data: %s %s\n" % (bragg_fname, intens_fname))
        f.write("Support: %s (%ld)\n" % (support_fname, self.input.num_supp))
        f.write("Algorithm: %s with beta = %.2f\n" % (algorithm_string, self.beta))
        f.write("Averaging algorithm: %s\n" % avg_algorithm_string)
        if self.do_positivity:
            f.write("Assuming electron density is positive\n")
        if self.do_histogram:
            f.write("Applying histogram constraint: %s\n" % hist_fname)
        if self.do_local_variation:
            f.write("Updating support using local variation\n")
        if self.do_bg_fitting:
            f.write("Fitting spherically symmetric background\n")
        #if self.do_blurring:
        #    f.write("Rotationally blurring model with %d orientations\n" % quat.num_rot)
        if self.do_normalize_prtf:
            f.write("Normalizing output by PRTF\n")
        f.write("Output prefix: %s\n" % self.output_prefix)
        f.write("-------------------------\n")
        f.write("iter    time    error\n")
        f.write("-------------------------\n")
    
    self.make_recon_folders()

def main():
    parser = argparse.ArgumentParser(description='Resolution extension phasing')
    parser.add_argument('-c', '--config_fname', help='Path to configuration file. Default=config.ini', default='config.ini')
    parser.add_argument('-T', '--testing', help='Flag for whether to run in testing (fixed seed) mode', action='store_true')
    args = parser.parse_args()

    algo = algorithm.Algorithm()
    setup(algo, args.config_fname, args.testing)
    
    for i in range(1, algo.num_iter + algo.num_avg_iter + 1):
        t1 = time.time()
        
        error = algo.run_iteration(i)
        if error < 0:
            sys.exit(1)
        if i > algo.num_iter:
            algo.volume.accumulate(algo.p1, algo.average_p1)
            algo.volume.accumulate(algo.p2, algo.average_p2)
        
        t2 = time.time()
        algo.save_current(i, t1, t2, error)
        sys.stderr.write("\rFinished %d/%d iterations. " % (i, algo.num_iter+algo.num_avg_iter))
        if i > algo.num_iter:
            sys.stderr.write("Now averaging")
    sys.stderr.write("\nCalculating prtf and writing to file.\n")
    
    if algo.num_avg_iter > 0:
        algo.average_p1 /= algo.num_avg_iter
        algo.average_p2 /= algo.num_avg_iter
    else:
        algo.average_p1 = algo.p1
        algo.average_p2 = algo.p2
    
    algo.calc_prtf(100)
    algo.save_output()

if __name__ == '__main__':
    main()
