#!/usr/bin/env python

import algorithm
import mrcfile
import pylab as P
import cupy as np

algo = algorithm.Algorithm('/home/ayyerkar/acads/brcont/config.ini', False)
with mrcfile.open('/home/ayyerkar/acads/brcont/data/convert/4et8_sim-srecon.ccp4', 'r') as f:
    sol = np.array(f.data)
    sol = sol[np.newaxis,:,:,:]
with mrcfile.open('/home/ayyerkar/acads/brcont/data/convert/4et8_sim-cpx.ccp4', 'r') as f:
    brg = np.fft.ifftshift(f.data)
'''
fsol = np.fft.fftn(sol[0])
fsol_c = np.copy(fsol)
fbrg = np.fft.fftshift(np.fft.ifftn(brg))

algo.proj_fourier(sol, algo.p1)
algo.r1 = 2 * algo.p1 - sol
algo.proj_direct(algo.r1, algo.p2)
print(float(np.linalg.norm(algo.p2-algo.p1)))
algo.iterate[:] = sol
print(algo.ER_algorithm())

'''
