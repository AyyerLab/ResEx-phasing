import sys
import argparse
import time
import numpy
import phaser

def main():
    parser = argparse.ArgumentParser(description='Resolution extension phasing')
    parser.add_argument('-c', '--config_fname',
                        help='Path to configuration file. Default=config.ini',
                        default='config.ini')
    parser.add_argument('-T', '--testing',
                        help='Flag for whether to run in testing (fixed seed) mode',
                        action='store_true')
    args = parser.parse_args()

    phas = phaser.Phaser(args.config_fname, args.testing)
    average_p1 = numpy.zeros(phas.p1.shape, dtype='f4')
    average_p1 = numpy.zeros(phas.p1.shape, dtype='f4')

    for i in range(1, phas.num_iter + phas.num_avg_iter + 1):
        time1 = time.time()

        error = phas.run_iteration(i)
        if error < 0:
            sys.exit(1)
        if i > phas.num_iter:
            average_p1 += phas.p1
            average_p2 += phas.p2

        time2 = time.time()
        phas.io.save_current(phas, phas.proj, i, time1, time2, error)
        sys.stderr.write("\rFinished %d/%d iterations. " % (i, phas.num_iter+phas.num_avg_iter))
        if i > phas.num_iter:
            sys.stderr.write("Now averaging")
    sys.stderr.write("\nCalculating prtf and writing to file.\n")

    if phas.num_avg_iter > 0:
        average_p1 /= phas.num_avg_iter
        average_p2 /= phas.num_avg_iter
    else:
        average_p1 = phas.p1
        average_p2 = phas.p2

    #phas.calc_prtf(100)
    phas.io.save_output(phas, phas.proj)

if __name__ == '__main__':
    main()
