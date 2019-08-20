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

def reconstruct_loop(l, args, recon_p1, recon_p2, devno=None):
    if CUDA and devno is not None:
        np.cuda.Device(devno).use()
    phas = phaser.Phaser(args.config_fname, args.testing)
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

    #return avg_p1, avg_p2
    recon_p1.put(avg_p1)
    recon_p2.put(avg_p2)

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

    mp.set_start_method('spawn')
    phas = phaser.Phaser(args.config_fname, args.testing, no_parse=True)
    avg_p1 = avg_p2 = None
    avg_p1_phasor = avg_p2_phasor = None

    recon_p1 = mp.Queue()
    recon_p2 = mp.Queue()

    num_devices = len(args.device)
    if num_devices > 1:
        num_cycles = int(np.ceil(phas.num_loops / num_devices))
        for cycle in range(num_cycles):
            l_min = cycle * num_devices
            l_max = min((cycle+1)*num_devices, phas.num_loops)
            print('%d: %d - %d'%(cycle, l_min, l_max))
            jobs = []
            for l in range(l_min, l_max):
                j = mp.Process(target=reconstruct_loop, args=(l, args, recon_p1, recon_p2, l%num_devices))
                jobs.append(j)
                
            [j.start() for j in jobs]
            for j in jobs:
                avg_p1, avg_p1_phasor = phas.accumulate(np.array(recon_p1.get()), avg_p1, avg_p1_phasor)
                avg_p2, avg_p2_phasor = phas.accumulate(np.array(recon_p2.get()), avg_p2, avg_p2_phasor)
            [j.join() for j in jobs]
    else:
        np.cuda.Device(args.device[0]).use()
        for l in range(phas.num_loops):
            reconstruct_loop(l, args, recon_p1, recon_p2)
            avg_p1, avg_p1_phasor = phas.accumulate(np.array(recon_p1.get()), avg_p1, avg_p1_phasor)
            avg_p2, avg_p2_phasor = phas.accumulate(np.array(recon_p2.get()), avg_p2, avg_p2_phasor)

    sys.stderr.write("Calculating prtf and writing to file.\n")

    phas.io.save_output(phas, phas.proj, avg_p1, avg_p2)
    if phas.num_avg_iter > 0 or phas.num_loops > 1:
        phas.io.save_prtf(phas, avg_p1_phasor, avg_p2_phasor)

if __name__ == '__main__':
    main()
