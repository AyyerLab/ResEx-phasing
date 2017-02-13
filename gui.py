#!/usr/bin/env python

import matplotlib
matplotlib.use('TkAgg')
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
import multiprocessing

class GUI():
    def __init__(self, master, merge_fname='', map_fname=''):
        self.master = master
        master.title('ResEx Phasing GUI')

        self.typestr = 'f4'
        rangemax = 10.
        rangemin = 0.
        self.merge_fname = Tk.StringVar()
        self.map_fname = Tk.StringVar()
        self.fname = Tk.StringVar()
        self.rangeminstr = Tk.StringVar()
        self.rangemaxstr = Tk.StringVar()
        self.imagename = Tk.StringVar()
        self.layernum = Tk.IntVar()
        self.radiusmin = Tk.StringVar()
        self.radiusmax = Tk.StringVar()
        self.scaleradmin = Tk.StringVar()
        self.scaleradmax = Tk.StringVar()
        self.resedge = Tk.StringVar()
        self.circleflag = Tk.IntVar()
        self.checkflag = Tk.IntVar()
        self.fslices = Tk.IntVar()
        self.scaleradflag = Tk.IntVar()
        self.suppressflag = Tk.IntVar()
        self.rangelock = Tk.IntVar()
        self.suppradstr = Tk.StringVar()
        self.suppthreshstr = Tk.StringVar()
        self.output_prefix = Tk.StringVar()
        self.point_group = Tk.StringVar()
        self.config_fname = Tk.StringVar()
        self.bgfitting_flag = Tk.IntVar()
        self.variation_flag = Tk.IntVar()
        self.positivity_flag = Tk.IntVar()
        self.histogram_flag = Tk.IntVar()

        self.merge_fname.set(merge_fname)
        self.map_fname.set(map_fname)
        self.fname.set(merge_fname)
        self.imagename.set('images/' + os.path.splitext(os.path.basename(self.fname.get()))[0] + '.png')
        self.radiusmin.set('0')
        self.radiusmax.set('0')
        self.scaleradmin.set('0')
        self.scaleradmax.set('0')
        self.circleflag.set(0)
        self.checkflag.set(0)
        self.fslices.set(0)
        self.scaleradflag.set(0)
        self.suppressflag.set(0)
        self.rangelock.set(0)
        self.suppradstr.set('3.')
        self.suppthreshstr.set('0.1')
        self.output_prefix.set('data/recon/test')
        self.point_group.set('222')
        self.config_fname.set('config.ini')
        self.bgfitting_flag.set(0)
        self.variation_flag.set(0)
        self.positivity_flag.set(0)
        self.histogram_flag.set(0)
        self.size = None
        self.vol = None
        self.rad = None
        self.old_fname = None
        self.space = None
        self.rangeminstr.set("%.1e" % rangemin)
        self.rangemaxstr.set("%.1e" % rangemax)
        self.layernum.set(0)
        self.map_image_exists = False
        self.vol_image_exists = False
        self.zeroed = False
        self.calculated_scale = False
        self.processed_map = False
        self.zoomed = False
        self.added_recon_tab = False

        if sys.platform != 'darwin':
            self.master.tk.eval('source [file join themes plastik plastik.tcl]')
            self.master.tk.eval('source [file join themes clearlooks clearlooks8.5.tcl]')
            fstyle = ttk.Style()
            #fstyle.theme_use('clearlooks')
            fstyle.theme_use('plastik')

        self.init_UI()

    def init_UI(self):
        self.master.bind('<Return>', lambda event: self.replot(zoom='current'))
        self.master.bind('<KP_Enter>', lambda event: self.replot(zoom='current'))
        self.master.bind('<Control-s>', self.save_plot)
        self.master.bind('<Control-q>', self.quit_)
        #self.master.bind('<Right>', self.increment_layer)
        #self.master.bind('<Left>', self.decrement_layer)
        self.master.bind('<Up>', self.increment_layer)
        self.master.bind('<Down>', self.decrement_layer)
        self.master.bind('<Control-m>', lambda event: self.plot_vol(fname=self.merge_fname.get()))
        self.master.bind('<Control-n>', lambda event: self.plot_map())
        self.master.rowconfigure(0, weight=1)
        self.master.columnconfigure(0, weight=1)

        # Canvas frame
        canvas_frame = ttk.Frame(self.master)
        canvas_frame.grid(row=0, column=0, sticky='news')
        canvas_frame.rowconfigure(0, weight=1)
        canvas_frame.columnconfigure(0, weight=1)
        
        line = ttk.Frame(canvas_frame)
        line.pack(fill=Tk.X)
        ttk.Label(line, textvariable=self.fname).pack(side=Tk.LEFT)
        self.fig = plt.figure(figsize=(15,5))
        self.fig.subplots_adjust(left=0.0, bottom=0.00, right=0.99, wspace=0.0)
        self.canvas = FigureCanvasTkAgg(self.fig, canvas_frame)
        self.canvas.show()
        self.canvas.get_tk_widget().pack(fill='both', expand=1)

        # Config frame
        config_frame = ttk.Frame(self.master)
        config_frame.grid(row=0, column=1, sticky='news')

        self.notebook = ttk.Notebook(config_frame)
        self.notebook.pack(fill=Tk.X)
        self.notebook.enable_traversal()
        self.gen_merge_tab()
        self.gen_map_tab()
        self.gen_recon_tab()

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
        self.slider.bind('<ButtonRelease-1>', lambda event: self.replot(zoom='current'))
        ttk.Button(line,text="+",command=self.increment_layer).pack(side=Tk.LEFT)

        line = ttk.Frame(config_frame)
        line.pack(fill=Tk.X)
        ttk.Label(line,text="Range: ").pack(side=Tk.LEFT)
        ttk.Entry(line,textvariable=self.rangeminstr,width=10).pack(side=Tk.LEFT)
        ttk.Entry(line,textvariable=self.rangemaxstr,width=10).pack(side=Tk.LEFT)
        ttk.Checkbutton(line,text='Lock',variable=self.rangelock).pack(side=Tk.LEFT)

        line = ttk.Frame(config_frame)
        line.pack(fill=Tk.X)
        ttk.Button(line,text="Preprocess",command=self.preprocess).pack(side=Tk.LEFT)
        ttk.Button(line,text="Quit",command=self.master.quit).pack(side=Tk.LEFT)

        self.plot_vol()

    def gen_map_tab(self):
        self.map_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.map_frame, text='Map')

        line = ttk.Frame(self.map_frame)
        line.pack(fill=Tk.X)
        ttk.Label(line,text="Map File:").pack(side=Tk.LEFT)
        ttk.Entry(line,textvariable=self.map_fname,width=30).pack(side=Tk.LEFT,fill=Tk.X,expand=1)
        ttk.Button(line,text="Plot Map",command=self.plot_map).pack(side=Tk.LEFT)

        line = ttk.Frame(self.map_frame)
        line.pack(fill=Tk.X)
        ttk.Label(line,text="Res. at edge (A):").pack(side=Tk.LEFT)
        ttk.Entry(line,textvariable=self.resedge,width=5).pack(side=Tk.LEFT)
        ttk.Label(line,text='Point group:').pack(side=Tk.LEFT)
        ttk.Entry(line,textvariable=self.point_group,width=5).pack(side=Tk.LEFT)
        ttk.Button(line,text="Process Map",command=self.process_map).pack(side=Tk.LEFT)
        ttk.Button(line,text='Reset',command=self.reset_map_tab).pack(side=Tk.LEFT)

    def reset_map_tab(self, event=None):
        self.processed_map = True
        self.map_frame.destroy()
        self.gen_map_tab()
        self.notebook.insert(1, self.map_frame)
        self.notebook.select(self.notebook.tabs()[1])

    def gen_merge_tab(self):
        self.merge_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.merge_frame, text='Merge')

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
        ttk.Checkbutton(line,text="Show",variable=self.circleflag,command=lambda: self.replot(zoom=False)).pack(side=Tk.LEFT)

        line = ttk.Frame(self.merge_frame)
        line.pack(fill=Tk.X)
        ttk.Label(line,text='Scale radii:').pack(side=Tk.LEFT)
        ttk.Entry(line,textvariable=self.scaleradmin,width=5).pack(side=Tk.LEFT)
        ttk.Entry(line,textvariable=self.scaleradmax,width=5).pack(side=Tk.LEFT)
        ttk.Checkbutton(line,text="Show",variable=self.scaleradflag,command=lambda: self.replot(zoom=False)).pack(side=Tk.LEFT)
        
        line = ttk.Frame(self.merge_frame)
        line.pack(fill=Tk.X)
        ttk.Button(line,text='Zero Outer',command=self.zero_outer).pack(side=Tk.LEFT)
        ttk.Button(line,text='Calc. Scale',command=self.calc_scale).pack(side=Tk.LEFT)
        ttk.Button(line,text='Reset',command=self.reset_merge_tab).pack(side=Tk.LEFT)

    def reset_merge_tab(self, event=None):
        self.zeroed = False
        self.calculated_scale = False
        self.merge_frame.destroy()
        self.gen_merge_tab()
        self.notebook.insert(0, self.merge_frame)
        self.notebook.select(self.notebook.tabs()[0])

    def parse_extension(self, filename):
        ext_string = os.path.splitext(os.path.basename(filename))[1]
        
        if ext_string == '.raw':
            self.typestr = 'f4'
        elif ext_string == '.bin':
            self.typestr = 'f8'
        elif ext_string == '.supp':
            self.typestr = 'uint8'
            if self.rangelock.get() == 0:
                self.rangemaxstr.set('%.1e' % 1)
        elif ext_string == '.cpx':
            self.typestr = 'complex64'
        else:
            print "Did not understand data type from extension. Defaulting to float."
            self.typestr = 'f4'

    def parse_vol(self):
        if not os.path.isfile(self.fname.get()):
            if self.fname.get() != '':
                print "Unable to open", self.fname.get()
            return
        self.parse_extension(self.fname.get())
        self.vol = np.fromfile(self.fname.get(), dtype=self.typestr)
        size = int(round(self.vol.size**(1/3.)))
        self.size = size
        self.vol_size = size
        self.vol = self.vol.reshape(size,size,size)
        
        if self.typestr == 'complex64':
            self.vol = np.square(np.absolute(self.vol))
        self.slider.configure(from_=0,to=self.size-1)
        if self.rangelock.get() == 0:
            rmax = self.vol[self.vol>0].mean()
            rmax = 5*self.vol[(self.vol>0.01*rmax) & (self.vol<100*rmax)].mean()
            self.rangemaxstr.set('%.1e' % rmax)
        if not self.vol_image_exists:
            self.layernum.set(self.size/2)
            self.radiusmin.set('%d' % (self.size/2/2))
            self.radiusmax.set('%d' % (self.size/2))
            self.scaleradmin.set('%d' % (self.size/2/2*0.9))
            self.scaleradmax.set('%d' % (self.size/2/2*1.1))
        self.old_fname = self.fname.get()

    def parse_map(self):
        with open(self.fname.get(), 'rb') as f:
            f.seek(28, 0)
            grid = np.fromfile(f, '=i4', count=3)
            nx, ny, nz = tuple(grid)
            f.seek(1024, 0)
            vol = np.fromfile(f, '=f4', count=nx*ny*nz).reshape(nx, ny, nz)
        edgesum = (np.abs(vol[:,:,0]).sum() + np.abs(vol[:,:,-1]).sum() + np.abs(vol[:,0]).sum() + np.abs(vol[:,-1]).sum() + np.abs(vol[0]).sum() + np.abs(vol[-1]).sum()) / 6.
        centralsum = (np.abs(vol[:,:,nz/2]).sum() + np.abs(vol[:,ny/2]).sum() + np.abs(vol[nx/2]).sum())/ 3.
        if edgesum > centralsum:
            vol = np.roll(vol, nx/2, axis=0)
            vol = np.roll(vol, ny/2, axis=1)
            vol = np.roll(vol, nz/2, axis=2)
        s = max(nx, ny, nz)
        self.size = s
        self.vol = np.pad(vol, (((s-nx)/2,s-nx-(s-nx)/2),((s-ny)/2,s-ny-(s-ny)/2),((s-nz)/2,s-nz-(s-nz)/2)), mode='constant', constant_values=0)
        self.slider.configure(from_=0,to=self.size-1)
        if not self.map_image_exists:
            self.layernum.set(self.size/2)
        self.old_fname = self.fname.get()

    def plot_slices(self, layernum, space=None, zoom=False, slices=None):
        if space is None:
            space = self.space
        self.imagename.set('images/' + os.path.splitext(os.path.basename(self.fname.get()))[0] + '.png')
        rangemax = float(self.rangemaxstr.get())
        rangemin = float(self.rangeminstr.get())
        if self.vol is None and slices is None:
            return
        
        if zoom == 'current':
            zoom = self.zoomed
        else:
            self.zoomed = zoom
        
        if slices is not None:
            if zoom:
                min = self.size/3
                max = 2*min
                a = slices[0,min:max,min:max]
                b = slices[1,min:max,min:max]
                c = slices[2,min:max,min:max]
            else:
                a = slices[0]
                b = slices[1]
                c = slices[2]
        else:
            if zoom:
                min = self.size/3
                max = 2*min
                a = self.vol[layernum,min:max,min:max]
                b = self.vol[min:max,layernum,min:max]
                c = self.vol[min:max,min:max,layernum]
            else:
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
        
        if self.scaleradflag.get() is 1: 
            rmin = float(self.scaleradmin.get())
            rmax = float(self.scaleradmax.get())
            s1.add_artist(patches.Circle((self.size/2,self.size/2), rmin, ec='white', fc='none', ls='dashed'))
            s1.add_artist(patches.Circle((self.size/2,self.size/2), rmax, ec='white', fc='none', ls='dashed'))
            s2.add_artist(patches.Circle((self.size/2,self.size/2), rmin, ec='white', fc='none', ls='dashed'))
            s2.add_artist(patches.Circle((self.size/2,self.size/2), rmax, ec='white', fc='none', ls='dashed'))
            s3.add_artist(patches.Circle((self.size/2,self.size/2), rmin, ec='white', fc='none', ls='dashed'))
            s3.add_artist(patches.Circle((self.size/2,self.size/2), rmax, ec='white', fc='none', ls='dashed'))
        
        self.space = space
        self.canvas.show()

    def replot(self, event=None, **kwargs):
        if self.map_image_exists:
            self.plot_map(**kwargs)
        elif self.vol_image_exists:
            self.plot_vol(fname=self.fname.get(), **kwargs)
        else:
            self.plot_vol(fname=self.fname.get(), **kwargs) # Default plotting merge

    def plot_vol(self, fname=None, event=None, force=False, sigma=False, **kwargs):
        if fname is None:
            self.fname.set(self.merge_fname.get())
        else:
            self.fname.set(fname)
        if not self.vol_image_exists:
            self.parse_vol()
        elif self.old_fname != self.fname.get() or force == True: 
            print "Reparsing volume:", self.fname.get()
            self.parse_vol()
        
        if sigma and self.suppressflag.get() == 1:
            c = self.vol.shape[0] / 2
            if self.rad is None:
                if os.path.isfile('data/sigma_%d.bin' % self.size):
                    self.rad = np.fromfile('data/rad_%d.bin' % self.size, '=f8').reshape(self.size,self.size,self.size)
                    self.sigma = np.fromfile('data/sigma_%d.bin' % self.size, '=f8').reshape(self.size,self.size,self.size)
                else:
                    x, y, z = np.indices(self.vol.shape)
                    x -= c
                    y -= c
                    z -= c
                    self.rad = np.sqrt(x*x + y*y + z*z)
                    print 'Calculated self.rad'
                    self.sigma = np.exp(-8 * self.rad**2 / (c**2))
                    print 'Calculated self.sigma'
                    self.rad.tofile('data/rad_%d.bin' % self.size)
                    self.sigma.tofile('data/sigma_%d.bin' % self.size)
            self.vol *= (1. - self.sigma)
        self.plot_slices(self.layernum.get(), space='fourier', **kwargs)
        self.vol_image_exists = True
        self.map_image_exists = False

    def plot_map(self, event=None, force=False, sigma=False, **kwargs):
        self.fname.set(self.map_fname.get())
        if not self.map_image_exists:
            self.parse_map()
            if self.rangelock.get() == 0:
                self.rangemaxstr.set('%.1e' % 10)
        elif self.old_fname != self.fname.get() or force == True:
            print "Reparsing map:", self.fname.get()
            self.parse_map()
        self.plot_slices(self.layernum.get(), space='real', **kwargs)
        self.map_image_exists = True
        self.vol_image_exists = False

    def zero_outer(self, event=None):
        rmin = float(self.radiusmin.get())
        rmax = float(self.radiusmax.get())
        print '-'*80
        os.system('./utils/zero_outer %s %d %d %d' % (self.merge_fname.get(), self.size, rmin, rmax))
        print '-'*80
        
        if not self.zeroed:
            #ttk.Separator(self.merge_frame, orient=Tk.HORIZONTAL).pack(fill=Tk.X, padx=5, pady=5)
            zero_fname = os.path.splitext(self.merge_fname.get())[0] + '-zero.raw'
            line = ttk.Frame(self.merge_frame)
            line.pack(fill=Tk.X)
            ttk.Label(line,text='Zero-ed volume: ').pack(side=Tk.LEFT)
            ttk.Button(line,text=zero_fname,command=lambda: self.plot_vol(fname=zero_fname)).pack(side=Tk.LEFT)
        self.zeroed = True
        #if self.calculated_scale and self.processed_map and not self.added_recon_tab:
        #    self.gen_recon_tab()

    def calc_scale(self, event=None):
        rmin = float(self.scaleradmin.get())
        rmax = float(self.scaleradmax.get())
        mapnoext = os.path.splitext(os.path.basename(self.map_fname.get()))[0]
        sym_model = 'data/convert/'+mapnoext+'-sym.raw'
        cmd = './utils/calc_scale %s %s %d %d %d' % (sym_model, self.merge_fname.get(), self.size, rmin, rmax)
        output = subprocess.check_output(cmd.split(), shell=False)
        self.scale_factor = float(output.split()[4])
        if not self.calculated_scale:
            line = ttk.Frame(self.merge_frame)
            line.pack(fill=Tk.X)
            self.scale_label = ttk.Label(line,text='Scale factor = %.6e'%self.scale_factor)
            self.scale_label.pack(side=Tk.LEFT)
        else:
            self.scale_label.config(text='Scale factor = %.6e' % self.scale_factor)
        self.calculated_scale = True
        #if self.zeroed and self.processed_map and not self.added_recon_tab:
        #    self.gen_recon_tab()

    def process_map(self, event=None):
        mapnoext = os.path.splitext(os.path.basename(self.map_fname.get()))[0]
        if os.path.isfile('data/convert/'+mapnoext+'.cpx') and not tkMessageBox.askyesno('Process Map', 'Found processed map output. Overwrite?', default=tkMessageBox.NO, icon=tkMessageBox.QUESTION, parent=self.master):
            with open('results/'+mapnoext+'.log', 'r') as f:
                words = f.read().split()
                warray = np.array(words)
                self.resedge.set(float(words[words.index('./utils/read_map')+2])/(self.vol_size/2))
                self.point_group.set(words[np.where(warray=='data/convert/'+mapnoext+'-srecon.raw')[0][0]+3])
                self.suppradstr.set('%.1f'%float(words[np.where(warray=='./utils/create_support')[0][-1]+3]))
                self.suppthreshstr.set('%.1f'%float(words[np.where(warray=='./utils/create_support')[0][-1]+4]))
        else:
            if self.resedge.get() is '':
                print 'Need resolution at edge of volume'
                return
            print '-'*80
            resedge = float(self.resedge.get())
            supp_rad = float(self.suppradstr.get())
            supp_thresh = float(self.suppthreshstr.get())
            command = './process_map.sh %s %d %f 0 %f %f %s' % (self.map_fname.get(), self.vol_size, resedge, supp_rad, supp_thresh, self.point_group.get())
            print command
            os.system(command)
            print '-'*80
        
        if not self.processed_map:
            self.add_to_map_frame(mapnoext)
        self.processed_map = True
        #if self.zeroed and self.calculated_scale and not self.added_recon_tab:
        #    self.gen_recon_tab()

    def add_to_map_frame(self, mapnoext):
        prefix = 'data/convert/'+mapnoext
        line = ttk.Frame(self.map_frame)
        line.pack(fill=Tk.X)
        ttk.Label(line,text='Complex: ').pack(side=Tk.LEFT)
        ttk.Button(line,text=prefix+'cpx',command=lambda: self.plot_vol(fname=prefix+'.cpx')).pack(side=Tk.LEFT)
        
        line = ttk.Frame(self.map_frame)
        line.pack(fill=Tk.X)
        ttk.Label(line,text='Symmetrized: ').pack(side=Tk.LEFT)
        ttk.Button(line,text=prefix+'-sym.raw',command=lambda: self.plot_vol(fname=prefix+'-sym.raw', sigma=True)).pack(side=Tk.LEFT)
        ttk.Checkbutton(line,text='Suppress low-q',variable=self.suppressflag,command=lambda: self.replot(zoom='current', sigma=True)).pack(side=Tk.LEFT)
        
        line = ttk.Frame(self.map_frame)
        line.pack(fill=Tk.X)
        ttk.Label(line,text='Density: ').pack(side=Tk.LEFT)
        ttk.Button(line,text=prefix+'-srecon.raw',command=lambda: self.plot_vol(fname=prefix+'-srecon.raw', zoom=True)).pack(side=Tk.LEFT)
        
        line = ttk.Frame(self.map_frame)
        line.pack(fill=Tk.X)
        ttk.Label(line,text='Support: ').pack(side=Tk.LEFT)
        ttk.Button(line,text=prefix+'.supp',command=lambda: self.plot_vol(fname=prefix+'.supp', zoom=True)).pack(side=Tk.LEFT)
        
        line = ttk.Frame(self.map_frame)
        line.pack(fill=Tk.X)
        ttk.Label(line,text='Support: Convolution radius = ').pack(side=Tk.LEFT)
        ttk.Entry(line,textvariable=self.suppradstr,width=5).pack(side=Tk.LEFT)
        ttk.Label(line,text='vox. Threshold = ').pack(side=Tk.LEFT)
        ttk.Entry(line,textvariable=self.suppthreshstr,width=5).pack(side=Tk.LEFT)
        
        line = ttk.Frame(self.map_frame)
        line.pack(fill=Tk.X)
        ttk.Button(line,text='Update support',command=self.reprocess_map).pack(side=Tk.LEFT)

    def reprocess_map(self, event=None):
        if self.resedge.get() is '':
            print 'Need resolution at edge of volume'
            return
        mapnoext = os.path.splitext(os.path.basename(self.map_fname.get()))[0]
        resedge = float(self.resedge.get())
        supp_radius = float(self.suppradstr.get())
        supp_thresh = float(self.suppthreshstr.get())
        print '-'*80
        command = './process_map.sh %s %d %f 1 %f %f %s' % (self.map_fname.get(), self.vol_size, resedge, supp_radius, supp_thresh, self.point_group.get())
        print command
        os.system(command)
        print '-'*80
        self.replot(force=True, zoom='current')

    def gen_recon_tab(self, event=None):
        self.recon_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.recon_frame, text='Recon')
        
        line = ttk.Frame(self.recon_frame)
        line.pack(fill=Tk.X)
        ttk.Label(line, text='Output prefix: ').pack(side=Tk.LEFT)
        ttk.Entry(line, textvariable=self.output_prefix).pack(side=Tk.LEFT,fill=Tk.X,expand=1)
        
        line = ttk.Frame(self.recon_frame)
        line.pack(fill=Tk.X)
        ttk.Label(line, text='Config filename: ').pack(side=Tk.LEFT)
        ttk.Entry(line, textvariable=self.config_fname).pack(side=Tk.LEFT,fill=Tk.X,expand=1)
        
        line = ttk.Frame(self.recon_frame)
        line.pack(fill=Tk.X)
        ttk.Button(line, text='Gen. Config', command=self.gen_config).pack(side=Tk.LEFT)
        ttk.Checkbutton(line, text='BG fitting', variable=self.bgfitting_flag).pack(side=Tk.LEFT)
        ttk.Checkbutton(line, text='Variation support', variable=self.variation_flag).pack(side=Tk.LEFT)
        ttk.Checkbutton(line, text='Positivity', variable=self.positivity_flag).pack(side=Tk.LEFT)
        
        line = ttk.Frame(self.recon_frame)
        line.pack(fill=Tk.X)
        ttk.Button(line, text='Launch Recon', command=self.launch_recon).pack(side=Tk.LEFT)
        ttk.Checkbutton(line,text='Keep Checking',variable=self.checkflag,command=self.keep_checking).pack(side=Tk.LEFT)
        ttk.Checkbutton(line,text='Fourier Slices',variable=self.fslices).pack(side=Tk.LEFT)
        
        self.added_recon_tab = True

    def gen_config(self, event=None):
        with open(self.config_fname.get(), 'w') as f:
            f.write('[parameters]\n')
            f.write('size = %d\n' % self.vol_size)
            f.write('bragg_qmax = %f\n' % (float(self.radiusmin.get())/(self.vol_size/2)))
            f.write('scale_factor = %f\n' % self.scale_factor)
            f.write('num_threads = %d\n' % multiprocessing.cpu_count())
            f.write('point_group = %s\n' % self.point_group.get())
            
            mapnoext = os.path.splitext(os.path.basename(self.map_fname.get()))[0]
            f.write('\n[files]\n')
            f.write('intens_fname = %s\n' % (os.path.splitext(self.merge_fname.get())[0].rstrip()+'-zero.raw'))
            f.write('bragg_fname = %s\n' % ('data/convert/'+mapnoext+'.cpx'))
            f.write('support_fname = %s\n' % ('data/convert/'+mapnoext+'.supp'))
            #f.write('input_fname = %s\n')
            f.write('output_prefix = %s\n' % self.output_prefix.get())
            
            f.write('\n[algorithm]\n')
            f.write('# Algorithm choices: DM, HIO, RAAR, mod-DM, ER\n')
            f.write('# With beta = 1, all algorithms except ER are equivalent\n')
            f.write('# By default, the end iterations are averaged. To use ER set avg_algorithm = ER\n')
            f.write('algorithm = DM\n')
            f.write('beta = 1.\n')
            if self.positivity_flag == 1:
                f.write('positivity = 1\n')
            if self.bgfitting_flag == 1:
                f.write('bg_fitting = 1\n')
            if self.variation_flag == 1:
                f.write('local_variation = 1\n')
            if self.histogram_flag == 1:
                f.write('histogram = 1\n')
                f.write('hist_fname = data/3wu2_hist.dat\n')
        print 'Generated %s:' % self.config_fname.get()
        os.system('cat %s' % self.config_fname.get())

    def launch_recon(self, event=None):
        pass

    def keep_checking(self, event=None):
        if self.checkflag.get() == 1:
            if self.fslices.get() == 0:
                self.fname.set(self.output_prefix.get()+'-slices.raw')
            else:
                self.fname.set(self.output_prefix.get()+'-fslices.raw')
            
            if os.path.isfile(self.fname.get()):
                done = False
                while not done:
                    s = np.fromfile(self.fname.get(), '=f4')
                    self.size = int(np.round((s.shape[0]/3)**0.5))
                    try:
                        s = s.reshape(3,self.size,self.size) 
                        done = True
                    except ValueError:
                        pass
                if self.fslices.get() == 0:
                    self.plot_slices(0, slices=s, zoom=True)
                else:
                    self.plot_slices(0, slices=s, zoom=False)
            self.master.after(1000, self.keep_checking)

    def preprocess(self, event=None):
        self.zero_outer()
        self.calc_scale()
        self.process_map()

    def increment_layer(self, event=None):
        self.layernum.set(min(self.layernum.get()+1, self.size))
        self.plot_slices(self.layernum.get(), zoom='current')

    def decrement_layer(self, event=None):
        self.layernum.set(max(self.layernum.get()-1, 0))
        self.plot_slices(self.layernum.get(), zoom='current')

    def save_plot(self, event=None):
        self.fig.savefig(self.imagename.get(), bbox_inches='tight', dpi=150)
        print "Saved to", self.imagename.get()

    def quit_(self, event=None):
        self.master.quit()

if __name__ == '__main__':
    root = Tk.Tk()
    if len(sys.argv) < 3:
        gui = GUI(root)
    else:
        gui = GUI(root, sys.argv[1], sys.argv[2])
    root.mainloop()
