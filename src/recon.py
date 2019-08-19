import sys
import argparse
import time
import multiprocessing as mp
import numpy
try:
    import cupy as np
    CUDA = True
except ImportError:
    import numpy as np
    CUDA = False
import phaser

def reconstruct_loop(l, phas):
    avg_p1 = avg_p2 = None
    avg_p1_phasor = avg_p2_phasor = None
    for i in range(phas.num_iter + phas.num_avg_iter):
        time1 = time.time()

        # Do iteration
        error = phas.run_iteration(i + 1)
        if error < 0:
            sys.exit(1)
        # Start accumulating from end of num_iter loop
        # [Accumulation also flips and aligns]
        if i >= phas.num_iter - 1:
            avg_p1, avg_p1_phasor = phas.accumulate(phas.p1, avg_p1, avg_p1_phasor)
            avg_p2, avg_p2_phasor = phas.accumulate(phas.p2, avg_p2, avg_p2_phasor)

        time2 = time.time()

        # Logging
        phas.io.save_current(phas, phas.proj, l*(phas.num_iter+phas.num_avg_iter)+i+1, time1, time2, error, slices=False)
        sys.stderr.write('\r')
        sys.stderr.write("Loop %d/%d: " % (l+1, phas.num_loops))
        sys.stderr.write("Finished %d/%d iterations. (%e) " % (i+1, phas.num_iter+phas.num_avg_iter, error))
        if i > phas.num_iter:
            sys.stderr.write("Now averaging")
        else:
            sys.stderr.write(" "*13)
        sys.stderr.flush()
    sys.stderr.write('\n')

    avg_p1 /= phas.num_avg_iter + 1
    avg_p2 /= phas.num_avg_iter + 1

    return avg_p1, avg_p2

def main():
    parser = argparse.ArgumentParser(description='Resolution extension phasing')
    parser.add_argument('-c', '--config_fname',
                        help='Path to configuration file. Default=config.ini',
                        default='config.ini')
    parser.add_argument('-T', '--testing',
                        help='Flag for whether to run in testing (fixed seed) mode',
                        action='store_true')
    parser.add_argument('-D', '--device',
                        help='If CUDA, specify [space-separated] device numbers to run on. Default: 0',
                        type=int, nargs='+', default=[0])
    args = parser.parse_args()

    if CUDA:
        np.cuda.Device(args.device[0]).use()
    phas = phaser.Phaser(args.config_fname, args.testing)
    avg_p1 = avg_p2 = None
    avg_p1_phasor = avg_p2_phasor = None

    for l in range(phas.num_loops):
        if l > 0:
            # Reallocate for new loop after first
            phas.allocate_memory()
            phas.io.init_iterate(phas.proj, phas.iterate, do_bg_fitting=phas.proj.do_bg_fitting, quiet=True)

        recon_p1, recon_p2 = reconstruct_loop(l, phas)
        avg_p1, avg_p1_phasor = phas.accumulate(np.array(recon_p1), avg_p1, avg_p1_phasor)
        avg_p2, avg_p2_phasor = phas.accumulate(np.array(recon_p2), avg_p2, avg_p2_phasor)

    sys.stderr.write("Calculating prtf and writing to file.\n")

    phas.io.save_output(phas, phas.proj, avg_p1, avg_p2)
    if phas.num_avg_iter > 0 or phas.num_loops > 1:
        phas.io.save_prtf(phas, avg_p1_phasor, avg_p2_phasor)

if __name__ == '__main__':
    main()
