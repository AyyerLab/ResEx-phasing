#!/usr/bin/env python

import os
import numpy as np
import glob
import pylab as P
import matplotlib.animation as animation
import argparse
import ConfigParser

parser = argparse.ArgumentParser(description='Plot real-space slices of predicted intensities')
parser.add_argument('-c', '--config', help='Path to config file. Default: config.ini', default='config.ini')
parser.add_argument('-n', '--num', help='Number of iterations to parse if different from log file', type=int, default=None)
parser.add_argument('-l', '--loop', help='Loop animation. Default=False', action='store_true', default=False)
parser.add_argument('-j', '--jump', help='Jump. Do only every j iterations. Default=1', type=int, default=1)
args = parser.parse_args()

def update_view(num):
    fslices = np.fromfile(flist[num], '=f4').reshape(3,size,size)
    v[0].set_data(fslices[0,bmin:bmax,bmin:bmax])
    v[1].set_data(fslices[1,bmin:bmax,bmin:bmax])
    v[2].set_data(fslices[2,bmin:bmax,bmin:bmax])
    iter_text.set_text('%d'%(num+1))
    return v[0], v[1], v[2], iter_text

config = ConfigParser.ConfigParser()
config.read(args.config)
size = config.getint('parameters', 'size')
prefix = config.get('files', 'output_prefix')

with open(prefix+'-log.dat', 'r') as f:
    last_iteration = int(f.readlines()[-1].split()[0])
flist = np.array([prefix+'-slices/%.4d.raw'%i for i in range(1, last_iteration+1)])[:args.num]

rangemax = np.fromfile(flist[-1], '=f4').max()/2.
print '%d slices will be plotted with rangemax %.3f' % (len(flist), rangemax)
bmin = size/3
bmax = 2*size/3

fig = P.figure(figsize=(15,5))
fig.subplots_adjust(left=0.0, bottom=0.00, right=0.99, wspace=0.0)

s1 = fig.add_subplot(131)
s2 = fig.add_subplot(132)
s3 = fig.add_subplot(133)

fslices = np.fromfile(flist[0], '=f4').reshape(3,size,size)
v = []
v.append(s1.matshow(fslices[0,bmin:bmax,bmin:bmax], vmax=rangemax, vmin=0))
v.append(s2.matshow(fslices[1,bmin:bmax,bmin:bmax], vmax=rangemax, vmin=0))
v.append(s3.matshow(fslices[2,bmin:bmax,bmin:bmax], vmax=rangemax, vmin=0))

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
