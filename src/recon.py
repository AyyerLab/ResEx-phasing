import sys
import argparse
import time

import algorithm

def main():
    parser = argparse.ArgumentParser(description='Resolution extension phasing')
    parser.add_argument('-c', '--config_fname', help='Path to configuration file. Default=config.ini', default='config.ini')
    parser.add_argument('-T', '--testing', help='Flag for whether to run in testing (fixed seed) mode', action='store_true')
    args = parser.parse_args()

    algo = algorithm.Algorithm(args.config_fname, args.testing)
    
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
