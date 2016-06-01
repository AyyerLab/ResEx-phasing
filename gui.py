#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import sys
import Tkinter as Tk
import ttk
import tkMessageBox
import os
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.patches as patches
import subprocess

class GUI():
    def __init__(self, master, merge_fname, map_fname):
        self.master = master
        master.title('ResEx Phasing GUI')

        self.typestr = 'f4'
        rangemax = 10.
        rangemin = 0.
        self.merge_fname = Tk.StringVar()
        self.map_fname = Tk.StringVar()
        self.rangeminstr = Tk.StringVar()
        self.rangemaxstr = Tk.StringVar()
        self.imagename = Tk.StringVar()
        self.layernum = Tk.IntVar()
        self.radiusmin = Tk.StringVar()
        self.radiusmax = Tk.StringVar()
        self.resedge = Tk.StringVar()
        self.circleflag = Tk.IntVar()
        self.map_image_exists = False
        self.vol_image_exists = False

        self.merge_fname.set(merge_fname)
        self.map_fname.set(map_fname)
        self.fname = merge_fname
        self.imagename.set('images/' + os.path.splitext(os.path.basename(self.fname))[0] + '.png')
        self.radiusmin.set('0')
        self.radiusmax.set('0')
        self.circleflag.set(0)
        self.size = None
        self.old_fname = None
        self.space = None
        self.rangeminstr.set("%.1e" % rangemin)
        self.rangemaxstr.set("%.1e" % rangemax)
        self.layernum.set(0)
        self.zeroed = False
        self.processed_map = False

        self.master.tk.eval('source [file join themes plastik plastik.tcl]')
        self.master.tk.eval('source [file join themes clearlooks clearlooks8.5.tcl]')
        fstyle = ttk.Style()
        #fstyle.theme_use('clearlooks')
        fstyle.theme_use('plastik')
        self.init_UI()

    def init_UI(self):
        self.master.bind('<Return>', self.replot)
        self.master.bind('<KP_Enter>', self.replot)
        self.master.bind('<Control-s>', self.save_plot)
        self.master.bind('<Control-q>', self.quit_)
        self.master.bind('<Up>', self.increment_layer)
        self.master.bind('<Down>', self.decrement_layer)
        self.master.rowconfigure(0, weight=1)
        self.master.columnconfigure(0, weight=1)

        self.fig = plt.figure(figsize=(15,5))
        self.fig.subplots_adjust(left=0.0, bottom=0.00, right=0.99, wspace=0.0)
        canvas_frame = ttk.Frame(self.master)
        canvas_frame.grid(row=0, column=0, sticky='news')
        canvas_frame.rowconfigure(0, weight=1)
        canvas_frame.columnconfigure(0, weight=1)
        self.canvas = FigureCanvasTkAgg(self.fig, canvas_frame)
        self.canvas.show()
        self.canvas.get_tk_widget().pack(fill='both', expand=1)

        config_frame = ttk.Frame(self.master)
        config_frame.grid(row=0, column=1, sticky='news')

        n = ttk.Notebook(config_frame)
        n.pack(fill=Tk.X)
        self.merge_frame = ttk.Frame(n)
        self.map_frame = ttk.Frame(n)
        n.add(self.merge_frame, text='Merge')
        n.add(self.map_frame, text='Map')
        self.gen_merge_tab()
        self.gen_map_tab()

        line = ttk.Frame(config_frame)
        line.pack(fill=Tk.X, expand=1)

        line = ttk.Frame(config_frame)
        line.pack(fill=Tk.X)
        ttk.Label(line,text="Image name: ").pack(side=Tk.LEFT)
        ttk.Entry(line,textvariable=self.imagename,width=30).pack(side=Tk.LEFT,fill=Tk.X,expand=1)
        ttk.Button(line,text="Save",command=self.save_plot).pack(side=Tk.LEFT)

        line = ttk.Frame(config_frame)
        line.pack(fill=Tk.X)
        ttk.Button(line,text="-",command=self.decrement_layer).pack(side=Tk.LEFT)
        self.slider = ttk.Scale(line,from_=0,to=0,orient=Tk.HORIZONTAL,length=300,variable=self.layernum)
        self.slider.pack(side=Tk.LEFT,fill=Tk.X,expand=1)
        ttk.Button(line,text="+",command=self.increment_layer).pack(side=Tk.LEFT)

        line = ttk.Frame(config_frame)
        line.pack(fill=Tk.X)
        ttk.Label(line,text="Range: ").pack(side=Tk.LEFT)
        ttk.Entry(line,textvariable=self.rangeminstr,width=10).pack(side=Tk.LEFT)
        ttk.Entry(line,textvariable=self.rangemaxstr,width=10).pack(side=Tk.LEFT)

        line = ttk.Frame(config_frame)
        line.pack(fill=Tk.X)
        ttk.Button(line,text="Quit",command=self.master.quit).pack(side=Tk.LEFT)

        self.plot_vol()

    def gen_map_tab(self):
        line = ttk.Frame(self.map_frame)
        line.pack(fill=Tk.X)
        ttk.Label(line,text="Map File:").pack(side=Tk.LEFT)
        ttk.Entry(line,textvariable=self.map_fname,width=30).pack(side=Tk.LEFT,fill=Tk.X,expand=1)
        ttk.Button(line,text="Plot Map",command=self.plot_map).pack(side=Tk.LEFT)

        line = ttk.Frame(self.map_frame)
        line.pack(fill=Tk.X)
        ttk.Label(line,text="Res. at edge (A):").pack(side=Tk.LEFT)
        ttk.Entry(line,textvariable=self.resedge,width=5).pack(side=Tk.LEFT)
        ttk.Button(line,text="Process Map",command=self.process_map).pack(side=Tk.LEFT)

    def gen_merge_tab(self):
        line = ttk.Frame(self.merge_frame)
        line.pack(fill=Tk.X)
        ttk.Label(line, text="Merge File: ").pack(side=Tk.LEFT)
        ttk.Entry(line,textvariable=self.merge_fname,width=30).pack(side=Tk.LEFT,fill=Tk.X,expand=1)
        ttk.Button(line,text="Plot Merge",command=self.plot_vol).pack(side=Tk.LEFT)

        line = ttk.Frame(self.merge_frame)
        line.pack(fill=Tk.X)
        ttk.Label(line,text='Cutoff radii:').pack(side=Tk.LEFT)
        ttk.Entry(line,textvariable=self.radiusmin,width=5).pack(side=Tk.LEFT)
        ttk.Entry(line,textvariable=self.radiusmax,width=5).pack(side=Tk.LEFT)
        ttk.Checkbutton(line,text="Show",variable=self.circleflag,command=self.replot).pack(side=Tk.LEFT)

        line = ttk.Frame(self.merge_frame)
        line.pack(fill=Tk.X)
        ttk.Button(line,text='Zero Outer',command=self.zero_outer).pack(side=Tk.LEFT)
        ttk.Button(line,text='Calc. Scale',command=self.calc_scale).pack(side=Tk.LEFT)

    def replot(self, event=None):
        if self.map_image_exists:
            self.plot_map()
        elif self.vol_image_exists:
            self.plot_vol()
        else:
            self.plot_vol() # Default plotting merge

    def parse_extension(self, filename):
        ext_string = os.path.splitext(os.path.basename(filename))[1]
        
        if ext_string == '.raw':
            self.typestr = 'f4'
        elif ext_string == '.bin':
            self.typestr = 'f8'
        elif ext_string == '.supp':
            print "Support file"
            self.typestr = 'uint8'
            self.rangemaxstr.set('%.1e' % 1)
        elif ext_string == '.cpx':
            print "Complex file"
            self.typestr = 'complex64'
        else:
            print "Did not understand data type from extension. Defaulting to float."
            self.typestr = 'f4'

    def parse_vol(self):
        if not os.path.isfile(self.fname):
            print "Unable to open", self.fname 
            return
        self.parse_extension(self.fname)
        self.vol = np.fromfile(self.fname, dtype=self.typestr)
        size = int(round(self.vol.size**(1/3.)))
        self.size = size
        self.vol_size = size
        self.vol = self.vol.reshape(size,size,size)
        
        if self.typestr == 'complex64':
            self.vol = np.square(np.absolute(self.vol))
        self.slider.configure(from_=0,to=self.size-1)
        if not self.vol_image_exists:
            self.rangemaxstr.set('%.1e' % (5*self.vol.mean()))
            self.layernum.set(self.size/2)
            self.radiusmin.set('%d' % (self.size/2/2))
            self.radiusmax.set('%d' % (self.size/2))
        self.old_fname = self.fname

    def parse_map(self):
        f = open(self.fname, 'rb')
        f.seek(28, 0)
        grid = np.fromfile(f, '=i4', count=3)
        nx, ny, nz = tuple(grid)
        f.seek(1024, 0)
        vol = np.fromfile(f, '=f4').reshape(nx, ny, nz)
        vol = np.roll(vol, nx/2, axis=0)
        vol = np.roll(vol, ny/2, axis=1)
        vol = np.roll(vol, nz/2, axis=2)
        s = max(nx, ny, nz)
        self.size = s
        self.vol = np.pad(vol, (((s-nx)/2,s-nx-(s-nx)/2),((s-ny)/2,s-ny-(s-ny)/2),((s-nz)/2,s-nz-(s-nz)/2)), mode='constant', constant_values=0)
        self.slider.configure(from_=0,to=self.size-1)
        if not self.map_image_exists:
            self.layernum.set(self.size/2)
        self.old_fname = self.fname

    def plot_slices(self, layernum, space=None):
        if space is None:
            space = self.space
        self.imagename.set('images/' + os.path.splitext(os.path.basename(self.fname))[0] + '.png')
        rangemax = float(self.rangemaxstr.get())
        rangemin = float(self.rangeminstr.get())
        
        a = self.vol[layernum,:,:]	
        b = self.vol[:,layernum,:]	
        c = self.vol[:,:,layernum]
            
        s1 = self.fig.add_subplot(131)
        s1.matshow(a, vmin=rangemin, vmax=rangemax, cmap='jet')
        if space == 'fourier':
            plt.title("h = 0, YZ plane", y = 1.01)
        elif space == 'real':
            plt.title("YZ plane", y = 1.01)
        plt.axis('off')
        s2 = self.fig.add_subplot(132)
        s2.matshow(b, vmin=rangemin, vmax=rangemax, cmap='jet')
        if space == 'fourier':
            plt.title("k = 0, XZ plane", y = 1.01)
        elif space == 'real':
            plt.title("YZ plane", y = 1.01)
        plt.axis('off')
        s3 = self.fig.add_subplot(133)
        s3.matshow(c, vmin=rangemin, vmax=rangemax, cmap='jet')
        if space == 'fourier':
            plt.title("l = 0, XY plane", y = 1.01)
        elif space == 'real':
            plt.title("XY plane", y = 1.01)
        plt.axis('off')
        
        [a.remove() for a in list(set(s1.findobj(patches.Circle)))]
        [a.remove() for a in list(set(s2.findobj(patches.Circle)))]
        [a.remove() for a in list(set(s3.findobj(patches.Circle)))]
        
        if self.circleflag.get() is 1: 
            rmin = float(self.radiusmin.get())
            rmax = float(self.radiusmax.get())
            s1.add_artist(patches.Circle((self.size/2,self.size/2), rmin, ec='white', fc='none'))
            s1.add_artist(patches.Circle((self.size/2,self.size/2), rmax, ec='white', fc='none'))
            s2.add_artist(patches.Circle((self.size/2,self.size/2), rmin, ec='white', fc='none'))
            s2.add_artist(patches.Circle((self.size/2,self.size/2), rmax, ec='white', fc='none'))
            s3.add_artist(patches.Circle((self.size/2,self.size/2), rmin, ec='white', fc='none'))
            s3.add_artist(patches.Circle((self.size/2,self.size/2), rmax, ec='white', fc='none'))
        
        self.space = space
        self.canvas.show()

    def plot_vol(self, event=None):
        self.fname = self.merge_fname.get()
        if not self.vol_image_exists:
            self.parse_vol()
            #self.rangemaxstr.set('%.1e' % 1.e8)
        elif self.old_fname != self.fname: 
            print "Reparsing volume:", self.fname
            self.parse_vol()
        
        self.plot_slices(self.layernum.get(), space='fourier')
        self.vol_image_exists = True
        self.map_image_exists = False

    def plot_map(self, event=None):
        self.fname = self.map_fname.get()
        if not self.map_image_exists:
            self.parse_map()
            self.rangemaxstr.set('%.1e' % 10)
        elif self.old_fname != self.fname:
            print "Reparsing map:", self.fname
            self.parse_map()
        self.plot_slices(self.layernum.get(), space='real')
        self.map_image_exists = True
        self.vol_image_exists = False

    def zero_outer(self, event=None):
        rmin = int(self.radiusmin.get())
        rmax = int(self.radiusmax.get())
        print '-'*80
        os.system('./utils/zero_outer %s %d %d %d' % (self.merge_fname.get(), self.size, rmin, rmax))
        print '-'*80
        self.fname = os.path.splitext(self.merge_fname.get())[0] + '-zero.raw'
        self.parse_vol()
        self.plot_slices(self.layernum.get(), space='fourier')
        
        if not self.zeroed:
            ttk.Separator(self.merge_frame, orient=Tk.HORIZONTAL).pack(fill=Tk.X, padx=5, pady=5)
            line = ttk.Frame(self.merge_frame)
            line.pack(fill=Tk.X)
            ttk.Label(line,text='Zero-ed volume: %s' % self.fname).pack(side=Tk.LEFT)
        self.zeroed = True

    def process_map(self, event=None):
        mapnoext = os.path.splitext(os.path.basename(self.map_fname.get()))[0]
        if self.processed_map:
            if not tkMessageBox.askyesno('Process Map', 'Already processed. Try again?', default=tkMessageBox.NO, icon=tkMessageBox.QUESTION, parent=self.merge_frame):
                return
        if os.path.isfile('data/'+mapnoext+'.cpx'):
            if not tkMessageBox.askyesno('Process Map', 'Found processed map output. Overwrite?', default=tkMessageBox.NO, icon=tkMessageBox.QUESTION, parent=self.merge_frame):
                self.processed_map = True
                return
        if self.resedge.get() is '':
            print 'Need resolution at edge of volume'
            return
        print '-'*80
        os.system('./process_map.sh %s %d %f' % (self.map_fname.get(), self.vol_size, float(self.resedge.get())))
        print '-'*80
        self.processed_map = True

    def calc_scale(self, event=None):
        rmin = int(self.radiusmin.get())
        rmax = int(self.radiusmax.get())
        mapnoext = os.path.splitext(os.path.basename(self.map_fname.get()))[0]
        sym_model = 'data/'+mapnoext+'-sym.raw'
        cmd = './utils/calc_scale %s %s %d %d %d' % (sym_model, self.merge_fname.get(), self.size, rmin, rmax)
        print '-'*80
        output = subprocess.check_output(cmd.split(), shell=False)
        self.scale_factor = output.split()[4]
        print 'scale_factor =', self.scale_factor
        print '-'*80

    def increment_layer(self, event=None):
        self.layernum.set(min(self.layernum.get()+1, self.size))
        self.plot_slices(self.layernum.get())

    def decrement_layer(self, event=None):
        self.layernum.set(max(self.layernum.get()-1, 0))
        self.plot_slices(self.layernum.get())

    def save_plot(self, event=None):
        self.fig.savefig(self.imagename.get(), bbox_inches='tight', dpi=150)
        print "Saved to", self.imagename.get()

    def quit_(self, event=None):
        self.master.quit()

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Need merge name and map name")
        sys.exit()

    root = Tk.Tk()
    gui = GUI(root, sys.argv[1], sys.argv[2])
    root.mainloop()
