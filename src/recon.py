import sys
import argparse
import configparser
from multiprocessing import cpu_count

def setup(struct algorithm_data *self, config_fname, fixed_seed):
    self.volume = malloc(sizeof(struct volume_data)) ;
    self.input = malloc(sizeof(struct input_data)) ;
    self.fft = malloc(sizeof(struct fft_data)) ;
    self.quat = malloc(sizeof(struct rotation)) ;
    struct volume_data *volume = self.volume ;
    struct input_data *input = self.input ;
    struct fft_data *fft = self.fft ;
    struct rotation *quat = self.quat ;
    
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
    output_prefix = config.get('files', 'output_prefix', fallback='data/output')

    algorithm_string = config.get('algorithm', 'algorithm')
    avg_algorithm_string = config.get('algorithm', 'avg_algorithm', fallback=None)
    beta = config.getfloat('algorithm', 'beta', fallback=1.)
    do_bg_fitting = config.getboolean('algorithm', 'bg_fitting', fallback=False)
    do_blurring = config.getboolean('algorithm', 'blurring', fallback=False)
    do_histogram = config.getboolean('algorithm', 'histogram', fallback=False)
    do_local_variation = config.getboolean('algorithm', 'local_variation', fallback=False)
    do_positivity = config.getboolean('algorithm', 'positivity', fallback=False)
    do_normalize_prtf = config.getboolean('algorithm', 'normalize_prtf', fallback=False)
    quat_fname = config.get('algorithm', 'quat_fname', fallback=None)
    num_div = config.getint('algorithm', 'num_div', fallback=-1)
    hist_fname = config.get('algorithm', 'hist_fname', fallback=None)
    sigma = config.getfloat('algorithm', 'sigma_deg', fallback=0.)
    
    if point_group != '1' and point_group != '222' and point_group != '4':
        raise ValueError("Only '1', '4' and '222' point_group values supported currently")
    
    self.size = size ;
    self.vol = size*size*size ;
    self.num_vox = self.do_bg_fitting ? self.vol * 2 : self.vol ;
    omp_set_num_threads(num_threads) ;
    
    self.fft = fft.FFT(size, num_threads)
    self.fft.create_plans()
    self.input = input.Input(size)
    self.volume = volume.Volume(size)

    if (parse_algorithm_strings(self, algorithm_string, avg_algorithm_string))
        return 1 ;
    algorithm_allocate_memory(self) ;
    self.input.parse_intens(intens_fname, scale_factor, self.do_bg_fitting)
    self.input.parse_bragg(bragg_fname, bragg_qmax)
    self.input.parse_support(support_fname)
    if do_histogram:
        self.input.parse_histogram(hist_fname)
    if do_blurring:
        print('Blurring currently not implemented')
        
    self.input.init_iterate(self.iterate, input_fname, inputbg_fname, do_bg_fitting, fixed_seed)
    if do_bg_fitting:
        self.volume.init_radavg()
    
    sprintf(line, "%s-log.dat", self.output_prefix) ;
    fp = fopen(line, "w") ;
    fprintf(fp, "Resolution extension iterative phasing\n") ;
    fprintf(fp, "Data: %s %s\n", bragg_fname, intens_fname) ;
    fprintf(fp, "Support: %s (%ld)\n", support_fname, input.num_supp) ;
    fprintf(fp, "Algorithm: %s with beta = %.2f\n", algorithm_string, self.beta) ;
    fprintf(fp, "Averaging algorithm: %s\n", avg_algorithm_string) ;
    if (self.do_positivity)
        fprintf(fp, "Assuming electron density is positive\n") ;
    if (self.do_histogram)
        fprintf(fp, "Applying histogram constraint: %s\n", hist_fname) ;
    if (self.do_local_variation)
        fprintf(fp, "Updating support using local variation\n") ;
    if (self.do_bg_fitting)
        fprintf(fp, "Fitting spherically symmetric background\n") ;
    if (self.do_blurring)
        fprintf(fp, "Rotationally blurring model with %d orientations\n", quat.num_rot) ;
    if (self.do_normalize_prtf)
        fprintf(fp, "Normalizing output by PRTF\n") ;
    fprintf(fp, "Output prefix: %s\n", self.output_prefix) ;
    fprintf(fp, "-------------------------\n") ;
    fprintf(fp, "iter    time    error\n") ;
    fprintf(fp, "-------------------------\n") ;
    fclose(fp) ;
    
    make_recon_folders(self) ;
    
    return 0 ;

def main():
    parser = argparse.ArgumentParser(help='Resolution extension phasing')
    parser.add_argument('-c', '--config_fname', help='Path to configuration file. Default=config.ini', default='config.ini')
    parser.add_argument('-T', '--testing', help='Flag for whether to run in testing (fixed seed) mode', action='store_true')
    args = parser.parse_args()

    long i ;
    int iter, fixed_seed ;
    struct timeval t1, t2 ;
    char config_fname[1024] ;
    
    struct algorithm_data algo ;
    
    if (setup(&algo, args.config_fname, args.testing))
        return 2 ;
    
    for (iter = 1 ; iter <= algo.num_iter+algo.num_avg_iter ; ++iter) {
        gettimeofday(&t1, NULL) ;
        
        float error = run_iteration(&algo, iter) ;
        if (error < 0)
            return 1 ;
        
        if (iter > algo.num_iter) {
            volume_accumulate(algo.p1, algo.average_p1, algo.num_vox) ;
            volume_accumulate(algo.p2, algo.average_p2, algo.num_vox) ;
        }
        
        gettimeofday(&t2, NULL) ;
        save_current(&algo, iter, t1, t2, error) ;
        
        fprintf(stderr, "\rFinished %d/%d iterations. ", iter, algo.num_iter+algo.num_avg_iter) ;
        if (iter > algo.num_iter)
            fprintf(stderr, "Now averaging") ;
    }
    fprintf(stderr, "\nCalculating prtf and writing to file.\n") ;
    
    if (algo.num_avg_iter > 0) {
        for (i = 0 ; i < algo.num_vox ; ++i) {
            algo.average_p1[i] /= algo.num_avg_iter ;
            algo.average_p2[i] /= algo.num_avg_iter ;
        }
    }
    else {
        for (i = 0 ; i < algo.num_vox ; ++i) {
            algo.average_p1[i] = algo.p1[i] ;
            algo.average_p2[i] = algo.p2[i] ;
        }
    }
    
    calc_prtf(&algo, 100) ;
    save_output(&algo) ;
    
    return 0 ;

if __name__ == '__main__':
    main()
