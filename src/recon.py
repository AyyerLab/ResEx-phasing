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
    avg_p1 = avg_p2 = None
    avg_p1_phasor = avg_p2_phasor = None

    for i in range(phas.num_loops * (phas.num_iter + phas.num_avg_iter)):
        time1 = time.time()

        if i > 0 and i % (phas.num_iter+phas.num_avg_iter) == 0:
            phas.allocate_memory()
            phas.io.init_iterate(phas.proj, phas.iterate, do_bg_fitting=phas.proj.do_bg_fitting, quiet=True)
        error = phas.run_iteration(i%(phas.num_iter+phas.num_avg_iter) + 1)
        if error < 0:
            sys.exit(1)
        if (phas.num_loops == 1 and i > phas.num_iter) or (phas.num_loops > 1 and i%phas.num_iter == phas.num_iter - 1):
            avg_p1, avg_p1_phasor = phas.accumulate(phas.p1, avg_p1, avg_p1_phasor)
            avg_p2, avg_p2_phasor = phas.accumulate(phas.p2, avg_p2, avg_p2_phasor)

        time2 = time.time()
        phas.io.save_current(phas, phas.proj, i+1, time1, time2, error)
        sys.stderr.write('\r')
        if phas.num_loops > 1:
            sys.stderr.write("Loop %d/%d: " % (i // phas.num_iter + 1, phas.num_loops))
        sys.stderr.write("Finished %d/%d iterations. (%e) " % (i%(phas.num_iter+phas.num_avg_iter) + 1, phas.num_iter+phas.num_avg_iter, error))
        if phas.num_loops == 1 and i > phas.num_iter:
            sys.stderr.write("Now averaging")
    sys.stderr.write("\nCalculating prtf and writing to file.\n")

    phas.io.save_output(phas, phas.proj, avg_p1, avg_p2)
    phas.io.save_prtf(phas, avg_p1_phasor, avg_p2_phasor)

if __name__ == '__main__':
    main()
