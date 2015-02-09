#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import sys
import Tkinter as Tk
import os
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

if len(sys.argv) < 2:
	print("Need filename")
	sys.exit()

flag = 0
if len(sys.argv) > 2:
	flag = int(sys.argv[2])

# Start TkInter and set defaults
root = Tk.Tk()
fig = plt.figure(figsize=(15,5))
fig.subplots_adjust(left=0.0, bottom=0.00, right=0.99, wspace=0.0)
mng = plt.get_current_fig_manager()

typestr = 'f4'
typesize = 4
size = 427
rangemax = 1e2

ext = os.path.splitext(os.path.basename(sys.argv[1]))[1]

fname = Tk.StringVar()
sizestr = Tk.StringVar()
rangestr = Tk.StringVar()
imagename = Tk.StringVar()
layernum = Tk.IntVar()
fname.set(sys.argv[1])
imagename.set(os.path.splitext(os.path.basename(fname.get()))[0] + '.png')

image_exists = 0

# Determine extension and resulting data types
def parse_extension(ext_string):
	global typestr, size, typesize, rangemax
	if ext_string == '.raw':
		typestr = 'f4'
		typesize = 4
		size = 427
		rangemax = 1e2
	elif ext_string == '.bin':
		typestr = 'f8'
		typesize = 8
		size = 427
		rangemax = 0.002
	elif ext_string == '.supp':
		print "Support file"
		typestr = 'uint8'
		typesize = 1
		size = 427
		rangemax = 1
	else:
		print "Did not understand data type from extension. Defaulting to float."
		typestr = 'f4'
		typesize = 4
		size = 427
		rangemax = 1

center = int(size/2)

parse_extension(ext)

sizestr.set(str(size))
rangestr.set("%.1e" % rangemax)
layernum.set(center)

old_fname = fname.get()
old_sizestr = sizestr.get()
old_rangestr = rangestr.get()

vol = np.zeros((size,size,size), dtype=typestr)

# Parse file to generate three arrays and plot them
def parse_vol():
	global vol, old_fname, old_sizestr
	size = int(sizestr.get())
	
	s = fname.get()
	
	if os.path.isfile(s):
		f = open(s, "r")
	else:
		print "Unable to open", s
		return
	
	if size == len(vol):
		vol = vol.reshape(size*size*size)
	else:
		vol = np.resize(vol, size*size*size)
	
	vol = np.fromfile(f, dtype=typestr, count=size*size*size)
	vol = vol.reshape((size,size,size))
	
	old_fname = fname.get()
	old_sizestr = sizestr.get()

def plot_vol_slices(layernum):
	global image_exists, flag
	imagename.set(os.path.splitext(os.path.basename(fname.get()))[0] + '.png')
	rangemax = float(rangestr.get())
	
	if flag is 0:
		a = vol[layernum,142:285,142:285]
		b = vol[142:285,layernum,142:285]
		c = vol[142:285,142:285,layernum]
	elif flag is 1:
		a = vol[layernum,:,:]	
		b = vol[:,layernum,:]	
		c = vol[:,:,layernum]
		
	s1 = fig.add_subplot(131)
	s1.matshow(a, vmin=0, vmax=rangemax, cmap='hot')
	plt.title("h = 0, YZ plane", y = 1.01)
	plt.axis('off')
	s2 = fig.add_subplot(132)
	s2.matshow(b, vmin=0, vmax=rangemax, cmap='hot')
	plt.title("k = 0, XZ plane", y = 1.01)
	plt.axis('off')
	s3 = fig.add_subplot(133)
	s3.matshow(c, vmin=0, vmax=rangemax, cmap='hot')
	plt.title("l = 0, XY plane", y = 1.01)
	plt.axis('off')
	
	canvas.show()
	
	image_exists = 1

def parse_and_plot(event=None):
	if image_exists == 0:
		parse_vol()
		plot_vol_slices(layernum.get())
	elif old_fname == fname.get() and old_sizestr == sizestr.get():
		plot_vol_slices(layernum.get())
	else:
		print "Reparsing volume:", fname.get(), sizestr.get()
		slider.configure(from_=0,to=int(sizestr.get()))
		layernum.set(int(int(sizestr.get())/2))
		
		parse_vol()
		plot_vol_slices(layernum.get())

def increment_layer(event=None):
	layernum.set(min(layernum.get()+1, size))
	plot_vol_slices(layernum.get())

def decrement_layer(event=None):
	layernum.set(max(layernum.get()-1, 0))
	plot_vol_slices(layernum.get())

def save_plot(event=None):
	fig.savefig(imagename.get(), bbox_inches='tight')
	print "Saved to", imagename.get()

def quit_(event=None):
	root.destroy()

# Tk GUI
root.bind('<Return>', parse_and_plot)
root.bind('<KP_Enter>', parse_and_plot)
root.bind('<Control-s>', save_plot)
root.bind('<Control-q>', quit_)
root.bind('<Up>', increment_layer)
root.bind('<Down>', decrement_layer)
root.bind('<ButtonRelease-1>', parse_and_plot)

canvas = FigureCanvasTkAgg(fig, root)
canvas.get_tk_widget().grid(rowspan=24)

Tk.Label(root, text="Filename: ").grid(row=0,column=1,sticky=Tk.E)
Tk.Entry(
	root, 
	textvariable = fname,
	width = 35
	).grid(row=0,column=2,columnspan=3,sticky=Tk.W)

Tk.Label(root, text="Size: ").grid(row=1,column=1,sticky=Tk.E)
Tk.Entry(
	root,
	textvariable = sizestr,
	width = 20
	).grid(row=1,column=2,columnspan=2,sticky=Tk.W)

Tk.Label(root, text="Range: ").grid(row=2,column=1,sticky=Tk.E)
Tk.Entry(
	root,
	textvariable = rangestr,
	width = 20
	).grid(row=2,column=2,columnspan=2,sticky=Tk.W)

Tk.Label(root, text="Image name: ").grid(row=3,column=1,sticky=Tk.E)
Tk.Entry(
	root,
	textvariable = imagename,
	width = 30
	).grid(row=3,column=2,columnspan=2,sticky=Tk.W)

Tk.Button(
	root,
	text = "Save",
	command = save_plot
	).grid(row=3,column=4,sticky=Tk.W)

Tk.Button(root,text = "+",command = increment_layer
	).grid(row=4,column=3,sticky=Tk.SW)
slider = Tk.Scale(root,
	from_ = 142, 
	to = 285,
	orient = Tk.HORIZONTAL,
	length = 250,
	width = 20,
	variable = layernum
#	command = parse_and_plot
	)
slider.grid(row=4,column=2,sticky=Tk.W)
Tk.Button(root,text = "-",command = decrement_layer
	).grid(row=4,column=1,sticky=Tk.SE)

Tk.Button(
	root,
	text = "Plot",
	command = parse_and_plot
	).grid(row=5,column=1,sticky=Tk.W)

Tk.Button(
	root,
	text = "Quit",
	command = root.quit
	).grid(row=5,column=2,sticky=Tk.W)

parse_and_plot()

root.mainloop()
