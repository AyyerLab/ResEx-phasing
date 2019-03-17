import sys
import argparse
import time
import numpy
try:
    import cupy as np
    CUDA = True
except ImportError:
    import numpy as np
    CUDA = False
import phaser

def main():
    parser = argparse.ArgumentParser(description='Resolution extension phasing')
    parser.add_argument('-c', '--config_fname',
                        help='Path to configuration file. Default=config.ini',
                        default='config.ini')
    parser.add_argument('-T', '--testing',
                        help='Flag for whether to run in testing (fixed seed) mode',
                        action='store_true')
    parser.add_argument('-D', '--device',
                        help='If CUDA, specify device number to run on. Default: 0',
                        type=int, default=0)
    args = parser.parse_args()

    if CUDA:
        np.cuda.Device(args.device).use()
    phas = phaser.Phaser(args.config_fname, args.testing)
    average_p1 = numpy.zeros(phas.p1.shape, dtype='f4')
    average_p2 = numpy.zeros(phas.p1.shape, dtype='f4')

    for i in range(phas.num_loops * (phas.num_iter + phas.num_avg_iter)):
        time1 = time.time()

        if i > 0 and i % (phas.num_iter+phas.num_avg_iter) == 0:
            phas.allocate_memory()
            phas.io.init_iterate(phas.proj, phas.iterate, do_bg_fitting=phas.proj.do_bg_fitting, quiet=True)
        error = phas.run_iteration(i%(phas.num_iter+phas.num_avg_iter) + 1)
        if error < 0:
            sys.exit(1)
        if (phas.num_loops == 1 and i > phas.num_iter) or (phas.num_loops > 1 and i%phas.num_iter == phas.num_iter - 1):
            if CUDA:
                average_p1 += phas.p1.get()
                average_p2 += phas.p2.get()
            else:
                average_p1 += phas.p1
                average_p2 += phas.p2
            if phas.num_loops > 1:
                numpy.save(phas.io.output_prefix+'-pf-%.4d.npy'%(i/phas.num_iter), phas.p1.get())
                numpy.save(phas.io.output_prefix+'-pd-%.4d.npy'%(i/phas.num_iter), phas.p2.get())

        time2 = time.time()
        phas.io.save_current(phas, phas.proj, i+1, time1, time2, error)
        sys.stderr.write('\r')
        if phas.num_loops > 1:
            sys.stderr.write("Loop %d/%d: " % (i // phas.num_iter + 1, phas.num_loops))
        sys.stderr.write("Finished %d/%d iterations. (%e) " % (i%(phas.num_iter+phas.num_avg_iter) + 1, phas.num_iter+phas.num_avg_iter, error))
        if phas.num_loops == 1 and i > phas.num_iter:
            sys.stderr.write("Now averaging")
    sys.stderr.write("\nCalculating prtf and writing to file.\n")

    if phas.num_avg_iter > 0:
        average_p1 /= phas.num_avg_iter
        average_p2 /= phas.num_avg_iter
    else:
        if CUDA:
            average_p1 = phas.p1.get() / phas.num_loops
            average_p2 = phas.p2.get() / phas.num_loops
        else:
            average_p1 = phas.p1 / phas.num_loops
            average_p2 = phas.p2 / phas.num_loops

    #phas.calc_prtf(100)
    phas.io.save_output(phas, phas.proj, average_p1, average_p2)

if __name__ == '__main__':
    main()
