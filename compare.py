#!/usr/bin/python

import numpy as np

a = np.fromfile("data/recon.raw", dtype='f4')
b = np.fromfile("data/ps2obj.raw", dtype='f4')

rescale = np.average(a) / np.average(b)
print "rescale =", rescale
print "num_support =", np.count_nonzero(a), np.count_nonzero(b)
print "distance = %.4e" % (np.linalg.norm(b-rescale*a)/np.sum(b))

