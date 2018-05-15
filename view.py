#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import sys
import Tkinter as Tk
import tkFileDialog
import os
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.colors import LogNorm
import matplotlib.patches as patches

if len(sys.argv) < 2:
    print("Need filename")
    sys.exit()

# Start TkInter and set defaults
root = Tk.Tk()
fig = plt.figure(figsize=(15,5))
fig.subplots_adjust(left=0.0, bottom=0.00, right=0.99, wspace=0.0)
mng = plt.get_current_fig_manager()
s1 = fig.add_subplot(131)
s2 = fig.add_subplot(132)
s3 = fig.add_subplot(133)

typestr = 'f4'
typesize = 4
rangemax = 1e1
rangemin = 0
cmap = 'cubehelix'

fname = Tk.StringVar()
rangeminstr = Tk.StringVar()
rangemaxstr = Tk.StringVar()
imagename = Tk.StringVar()
layernum = Tk.IntVar()
radius_x = Tk.StringVar()
radius_y = Tk.StringVar()
radius_z = Tk.StringVar()
flag = Tk.IntVar()
circleflag = Tk.IntVar()

fname.set(sys.argv[1])
imagename.set('images/' + os.path.splitext(os.path.basename(fname.get()))[0] + '.png')
radius_x.set('200')
radius_y.set('200')
radius_z.set('200')
circleflag.set(0)
flag.set(0)
if len(sys.argv) > 2:
    flag.set(int(sys.argv[2]))

image_exists = 0

# Determine extension and resulting data types
def parse_extension(filename):
    global typestr, typesize, rangemax
    ext_string = os.path.splitext(os.path.basename(filename))[1]
    
    if ext_string == '.raw':
        typestr = 'f4'
        typesize = 4
        rangemax = 1e2
    elif ext_string == '.bin':
        typestr = 'f8'
        typesize = 8
        rangemax = 0.002
    elif ext_string == '.supp':
        print "Support file"
        typestr = 'uint8'
        typesize = 1
        rangemax = 1
    elif ext_string == '.cpx':
        print "Complex file"
        typestr = 'complex64'
        typesize = 8
        rangemax = 1e2
    else:
        print "Did not understand data type from extension. Defaulting to float."
        typestr = 'f4'
        typesize = 4
        rangemax = 1

parse_extension(sys.argv[1])

rangeminstr.set("%.1e" % rangemin)
rangemaxstr.set("%.1e" % rangemax)

old_fname = fname.get()
only_slices = False

# Parse file to generate three arrays and plot them
def parse_vol():
    global vol, old_fname, only_slices
    
    s = fname.get()
    
    if not os.path.isfile(s):
        print "Unable to open", s
        return
    
    vol = np.fromfile(s, dtype=typestr)
    size = int(round(vol.size**(1/3.)))
    slice_size = int(round(np.sqrt(vol.size/3)))
    if size*size*size == len(vol):
        vol = vol.reshape(size,size,size)
    elif slice_size*slice_size*3 == len(vol):
        vol = vol.reshape(3,slice_size,slice_size)
        only_slices = True
        size = slice_size
        print 'Only 3 slices and not full volume'
    else:
        vol = np.resize(vol, size*size*size).reshape(size,size,size)
    
    old_fname = fname.get()
    
    if typestr == 'complex64':
        vol = np.square(np.absolute(vol))
    
    return size

def plot_vol_slices(layernum):
    global image_exists, flag
    imagename.set('images/' + os.path.splitext(os.path.basename(fname.get()))[0] + '.png')
    rangemax = float(rangemaxstr.get())
    rangemin = float(rangeminstr.get())
    size = len(vol[0])
    
    if flag.get() is 0:
        min = int(size/3)
        max = int(2*size/3)
        
        if only_slices:
            a = vol[0,min:max,min:max]
            b = vol[1,min:max,min:max]
            c = vol[2,min:max,min:max]
        else:
            a = vol[layernum,min:max,min:max]
            b = vol[min:max,layernum,min:max]
            c = vol[min:max,min:max,layernum]
    elif flag.get() is 1:
        if only_slices:
            a = vol[0]
            b = vol[1]
            c = vol[2]
        else:
            a = vol[layernum,:,:]    
            b = vol[:,layernum,:]    
            c = vol[:,:,layernum]
    
    s1.matshow(a, vmin=rangemin, vmax=rangemax, cmap=cmap)
    s1.set_title("h = 0, YZ plane", y = 1.01)
    s1.axis('off')
    s2.matshow(b, vmin=rangemin, vmax=rangemax, cmap=cmap)
    s2.set_title("k = 0, XZ plane", y = 1.01)
    s2.axis('off')
    s3.matshow(c, vmin=rangemin, vmax=rangemax, cmap=cmap)
    s3.set_title("l = 0, XY plane", y = 1.01)
    s3.axis('off')
    
    if flag.get() is 1:
        [a.remove() for a in list(set(s1.findobj(patches.Ellipse)))]
        [a.remove() for a in list(set(s2.findobj(patches.Ellipse)))]
        [a.remove() for a in list(set(s3.findobj(patches.Ellipse)))]
    
    if circleflag.get() is 1 and flag.get() is 1:
        rx = 2*float(radius_x.get())
        ry = 2*float(radius_y.get())
        rz = 2*float(radius_z.get())
        s1.add_artist(patches.Ellipse((size/2,size/2), rz, ry, 0, ec='white', fc='none'))
        s2.add_artist(patches.Ellipse((size/2,size/2), rz, rx, 0, ec='white', fc='none'))
        s3.add_artist(patches.Ellipse((size/2,size/2), ry, rx, 0, ec='white', fc='none'))
    
    canvas.show()
    
    image_exists = 1

def parse_and_plot(event=None):
    if image_exists == 0:
        size = parse_vol()
        slider.configure(from_=0,to=size-1)
        layernum.set(size/2)
        plot_vol_slices(layernum.get())
    elif old_fname == fname.get():
        plot_vol_slices(layernum.get())
    else:
        print "Reparsing volume:", fname.get()
        parse_extension(fname.get())
        
        size = parse_vol()
        slider.configure(from_=0,to=size-1)
        layernum.set(size/2)
        plot_vol_slices(layernum.get())

def force_plot(event=None):
    print "Reparsing volume:", fname.get()
    parse_vol()
    plot_vol_slices(layernum.get())

def open_file(event=None):
    filename = tkFileDialog.askopenfilename(initialdir=os.path.dirname(fname.get()), title='Choose file')
    fname.set(filename)
    parse_and_plot()

def increment_layer(event=None):
    layernum.set(min(layernum.get()+1, len(vol)-1))
    plot_vol_slices(layernum.get())

def decrement_layer(event=None):
    layernum.set(max(layernum.get()-1, 0))
    plot_vol_slices(layernum.get())

def flag_changed(event=None):
    plot_vol_slices(layernum.get())

def save_plot(event=None):
    fig.savefig(imagename.get(), bbox_inches='tight', dpi=150)
    print "Saved to", imagename.get()

def quit_(event=None):
    root.quit()

# Tk GUI
root.bind('<Return>', parse_and_plot)
root.bind('<KP_Enter>', parse_and_plot)
root.bind('<Control-s>', save_plot)
root.bind('<Control-q>', quit_)
root.bind('<Up>', increment_layer)
root.bind('<Down>', decrement_layer)
#root.bind('<ButtonRelease-1>', parse_and_plot)
root.rowconfigure(0, weight=1)
root.columnconfigure(0, weight=1)

canvas_frame = Tk.Frame(root)
canvas_frame.grid(row=0, column=0, sticky='news')
canvas_frame.rowconfigure(0, weight=1)
canvas_frame.columnconfigure(0, weight=1)
canvas = FigureCanvasTkAgg(fig, canvas_frame)
canvas.get_tk_widget().pack(fill='both', expand=1)

config_frame = Tk.Frame(root)
config_frame.grid(row=0, column=1,sticky=Tk.N)

Tk.Label(config_frame, text="Filename: ").grid(row=0,column=0,sticky=Tk.E)
Tk.Entry(
    config_frame,
    textvariable = fname,
    width = 35
    ).grid(row=0,column=1,columnspan=2,sticky=Tk.W)
Tk.Button(
    config_frame,
    text='Browse',
    command=open_file
    ).grid(row=0,column=3,sticky=Tk.W)

Tk.Label(config_frame, text="Range: ").grid(row=2,column=0,sticky=Tk.E)
Tk.Entry(
    config_frame,
    textvariable = rangeminstr,
    width = 10
    ).grid(row=2,column=1,columnspan=1,sticky=Tk.W)
Tk.Entry(
    config_frame,
    textvariable = rangemaxstr,
    width = 10
    ).grid(row=2,column=2,columnspan=1,sticky=Tk.W)

Tk.Label(config_frame, text="Image name: ").grid(row=3,column=0,sticky=Tk.E)
Tk.Entry(
    config_frame,
    textvariable = imagename,
    width = 35
    ).grid(row=3,column=1,columnspan=2,sticky=Tk.W)

Tk.Button(
    config_frame,
    text = "Save",
    command = save_plot
    ).grid(row=3,column=3,sticky=Tk.W)

Tk.Button(config_frame,text = "+",command = increment_layer
    ).grid(row=4,column=3,sticky=Tk.SW)
slider = Tk.Scale(config_frame,
    from_ = 0, 
    to = 1,
    orient = Tk.HORIZONTAL,
    length = 350,
    width = 20,
    variable = layernum
#    command = parse_and_plot
    )
slider.grid(row=4,column=1,columnspan=2,sticky=Tk.W)
Tk.Button(config_frame,text = "-",command = decrement_layer
    ).grid(row=4,column=0,sticky=Tk.SE)

Tk.Button(
    config_frame,
    text = "Plot",
    command = parse_and_plot
    ).grid(row=5,column=0,sticky=Tk.W)

Tk.Button(
    config_frame,
    text = "Reparse",
    command = force_plot
    ).grid(row=5,column=1,sticky=Tk.W)

Tk.Button(
    config_frame,
    text = "Quit",
    command = root.quit
    ).grid(row=5,column=2,sticky=Tk.W)

Tk.Checkbutton(
    config_frame,
    text = "Show full volume",
    variable = flag,
    command = flag_changed
    ).grid(row=6,column=0,columnspan=2)

Tk.Checkbutton(
    config_frame,
    text = "Show circles",
    variable = circleflag,
    command = flag_changed
    ).grid(row=6,column=2)

Tk.Entry(
    config_frame,
    textvariable = radius_x,
    width = 10
    ).grid(row=7,column=1,columnspan=1,sticky=Tk.W)
Tk.Entry(
    config_frame,
    textvariable = radius_y,
    width = 10
    ).grid(row=7,column=2,columnspan=1,sticky=Tk.W)
Tk.Entry(
    config_frame,
    textvariable = radius_z,
    width = 10
    ).grid(row=7,column=3,columnspan=1,sticky=Tk.W)

parse_and_plot()

root.mainloop()
