#!/usr/bin/env python

from __future__ import print_function
import os
import sys
import glob
import argparse
try:
    from ConfigParser import ConfigParser
except ImportError:
    from configparser import ConfigParser
import numpy as np
import pylab as P
import matplotlib.animation as animation
import mrcfile

parser = argparse.ArgumentParser(description='Plot slices of iterates')
parser.add_argument('-c', '--config', help='Path to config file. Default: config.ini', default='config.ini')
parser.add_argument('-f', '--fourier', help='Show Fourier slices. Default: False', action='store_true', default=False)
parser.add_argument('-s', '--support', help='Show support mask slices. Default: False', action='store_true', default=False)
parser.add_argument('-b', '--begin', help='First iteration to parse (default: 1)', type=int, default=1)
parser.add_argument('-e', '--end', help='Last iteration to parse if different from log file', type=int, default=None)
parser.add_argument('-l', '--loop', help='Loop animation. Default=False', action='store_true', default=False)
parser.add_argument('-j', '--jump', help='Jump. Do only every j iterations. Default=1', type=int, default=1)
args = parser.parse_args()

def update_view(num):
    with mrcfile.open(flist[num], 'r') as f:
        fslices = f.data
        if args.fourier:
            fslices = fslices**0.25
    v[0].set_data(fslices[0,bmin:bmax,bmin:bmax])
    v[1].set_data(fslices[1,bmin:bmax,bmin:bmax])
    v[2].set_data(fslices[2,bmin:bmax,bmin:bmax])
    iter_text.set_text('%d'%(args.begin+num))
    return v[0], v[1], v[2], iter_text

config = ConfigParser()
config.read(args.config)
size = config.getint('parameters', 'size')
prefix = os.path.join(os.path.dirname(args.config), config.get('files', 'output_prefix'))

if args.end is None:
    with open(prefix+'-log.dat', 'r') as f:
        args.end = int(f.readlines()[-1].split()[0])

if args.fourier and args.support:
    print('Cannot show both Fourier and support mask slices at the same time. Pick one.')
    sys.exit(1)
elif args.fourier:
    flist = np.array([prefix+'-fslices/%.4d.ccp4'%i for i in range(args.begin, args.end + 1)])
    bmin = 0
    bmax = None
elif args.support:
    flist = np.array([prefix+'-support/%.4d.ccp4'%i for i in range(args.begin, args.end + 1)])
    bmin = size//3
    bmax = 2*size//3
else:
    flist = np.array([prefix+'-slices/%.4d.ccp4'%i for i in range(args.begin, args.end + 1)])
    bmin = size//3
    bmax = 2*size//3

with mrcfile.open(flist[-1], 'r') as f:
    if args.fourier:
        rangemax = f.data.max()
        rangemin = 0
    elif args.support:
        rangemax = 1
        rangemin = 0
    else:
        #rangemax = f.data.max() / 1.5
        rangemax = f.data.max()
        rangemin=-rangemax*0.25
print('%d slices will be plotted with rangemax %.3f' % (len(flist), rangemax))

fig = P.figure(figsize=(15,5))
fig.subplots_adjust(left=0.0, bottom=0.00, right=0.99, wspace=0.0)

s1 = fig.add_subplot(131)
s2 = fig.add_subplot(132)
s3 = fig.add_subplot(133)

with mrcfile.open(flist[0], 'r') as f:
    fslices = f.data
    if args.fourier:
        fslices = fslices**0.25
v = []
interp = 'gaussian' if not args.support else 'none'
v.append(s1.matshow(fslices[0,bmin:bmax,bmin:bmax], vmax=rangemax, vmin=rangemin, cmap='cubehelix', interpolation=interp))
v.append(s2.matshow(fslices[1,bmin:bmax,bmin:bmax], vmax=rangemax, vmin=rangemin, cmap='cubehelix', interpolation=interp))
v.append(s3.matshow(fslices[2,bmin:bmax,bmin:bmax], vmax=rangemax, vmin=rangemin, cmap='cubehelix', interpolation=interp))

s1.set_title("YZ plane", y = 1.01)
s1.axis('off')
s2.set_title("XZ plane", y = 1.01)
s2.axis('off')
s3.set_title("XY plane", y = 1.01)
s3.axis('off')
iter_text = v[0].axes.text(0.05, 0.95, '1', color='w', transform=v[0].axes.transAxes)

figanim = animation.FuncAnimation(fig, update_view, frames=range(0,len(flist),args.jump), interval=15, blit=True, repeat=args.loop)

#figanim.save('images/%s.gif' % os.path.basename(prefix), writer='imagemagick', fps=30, dpi=120)
P.show()
