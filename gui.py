#!/usr/bin/env python

from __future__ import print_function
import sys
import os
import subprocess
import multiprocessing
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
# TK Stuff
'''
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import Tkinter as Tk
import ttk
import tkMessageBox
'''
# Qt Stuff
try:
    from PyQt5 import QtCore, QtWidgets, QtGui # pylint: disable=import-error
    import matplotlib
    matplotlib.use('qt5agg')
    from matplotlib.backends.backend_qt5agg import FigureCanvas # pylint: disable=no-name-in-module
    os.environ['QT_API'] = 'pyqt5'
except ImportError:
    import sip
    sip.setapi('QString', 2)
    from PyQt4 import QtCore, QtGui # pylint: disable=import-error
    from PyQt4 import QtGui as QtWidgets # pylint: disable=import-error
    import matplotlib
    matplotlib.use('qt4agg')
    from matplotlib.backends.backend_qt4agg import FigureCanvas # pylint: disable=no-name-in-module
    os.environ['QT_API'] = 'pyqt'

class GUI(QtWidgets.QMainWindow):
    def __init__(self, merge_fname='', map_fname=''):
        super(GUI, self).__init__()
        self.merge_fname = merge_fname
        self.map_fname = map_fname

        self.typestr = 'f4'
        '''
        self.merge_fname = Tk.StringVar()
        self.map_fname = Tk.StringVar()
        self.current_fname = Tk.StringVar()
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
        self.project_flag = Tk.IntVar()
        self.current_angle = Tk.StringVar()

        self.merge_fname.set(merge_fname)
        self.map_fname.set(map_fname)
        self.current_fname.set(merge_fname)
        self.imagename.set('images/' + os.path.splitext(os.path.basename(self.current_fname.get()))[0] + '.png')
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
        self.project_flag.set(0)
        self.current_angle.set('XY')
        self.rangeminstr.set("%.1e" % rangemin)
        self.rangemaxstr.set("%.1e" % rangemax)
        self.layernum.set(0)
        '''
        self.size = None
        self.vol = None
        self.rad = None
        self.old_fname = None
        self.space = None
        self.map_image_exists = False
        self.vol_image_exists = False
        self.zeroed = False
        self.calculated_scale = False
        self.processed_map = False
        self.zoomed = False
        self.added_recon_tab = False
        self.angle_list = ['XY', 'XZ', 'YZ']

        '''
        if sys.platform != 'darwin':
            self.master.tk.eval('source [file join themes plastik plastik.tcl]')
            self.master.tk.eval('source [file join themes clearlooks clearlooks8.5.tcl]')
            fstyle = ttk.Style()
            #fstyle.theme_use('clearlooks')
            fstyle.theme_use('plastik')
        '''

        self.init_UI()

    def init_UI(self):
        with open('style.css', 'r')as f:
            self.setStyleSheet(f.read())
        #self.setWindowFlags(QtCore.Qt.FramelessWindowHint)
        self.setWindowTitle('ResEx Phasing GUI')
        overall = QtWidgets.QWidget()
        self.setCentralWidget(overall)
        layout = QtWidgets.QHBoxLayout(overall)
        layout.setContentsMargins(0, 0, 0, 0)
        splitter = QtWidgets.QSplitter(QtCore.Qt.Horizontal)
        layout.addWidget(splitter)

        # Menu items
        menubar = self.menuBar()
        menubar.setNativeMenuBar(False)

        '''
        # File Menu
        filemenu = menubar.addMenu('&File')
        action = QtWidgets.QAction('&Load Volume', self)
        action.triggered.connect(self._load_volume)
        filemenu.addAction(action)
        action = QtWidgets.QAction('&Save Image', self)
        action.triggered.connect(self._save_plot)
        filemenu.addAction(action)
        action = QtWidgets.QAction('Save Log &Plot', self)
        action.triggered.connect(self._save_log_plot)
        filemenu.addAction(action)
        action = QtWidgets.QAction('&Quit', self)
        action.triggered.connect(self.close)
        filemenu.addAction(action)
        '''

        # Color map picker
        cmapmenu = menubar.addMenu('&Color Map')
        self.color_map = QtWidgets.QActionGroup(self, exclusive=True)
        for i, cmap in enumerate(['cubehelix', 'coolwarm', 'CMRmap', 'gray', 'gray_r', 'jet']):
            action = self.color_map.addAction(QtWidgets.QAction(cmap, self, checkable=True))
            if i == 0:
                action.setChecked(True)
            action.triggered.connect(lambda: self.replot(zoom='current'))
            cmapmenu.addAction(action)

        # Config frame
        config_frame = QtWidgets.QWidget(self)
        config_frame.setMinimumWidth(400)
        config_frame.setObjectName('config')
        vbox = QtWidgets.QVBoxLayout()
        config_frame.setLayout(vbox)
        splitter.addWidget(config_frame)

        hbox = QtWidgets.QHBoxLayout()
        vbox.addLayout(hbox)
        label = QtWidgets.QLabel('ResEx Phasing GUI', self)
        label.setObjectName('heading')
        hbox.addWidget(label)

        self.notebook = QtWidgets.QTabWidget()
        self.merge_tab = QtWidgets.QWidget()
        self.notebook.addTab(self.merge_tab, 'Merge')
        self.map_tab = QtWidgets.QWidget()
        self.notebook.addTab(self.map_tab, 'Map')
        self.recon_tab = QtWidgets.QWidget()
        self.notebook.addTab(self.recon_tab, 'Recon')
        vbox.addWidget(self.notebook)
        '''
        self.notebook = ttk.Notebook(config_frame)
        self.notebook.pack(fill=Tk.X)
        self.notebook.enable_traversal()
        self.gen_merge_tab()
        self.gen_map_tab()
        self.gen_recon_tab()
        '''

        vbox.addStretch(1)

        hbox = QtWidgets.QHBoxLayout()
        vbox.addLayout(hbox)
        label = QtWidgets.QLabel('Image name: ', self)
        hbox.addWidget(label)
        if self.merge_fname != '':
            fname = 'images/' + os.path.splitext(os.path.basename(self.merge_fname))[0] + '.png'
        else:
            fname = ''
        self.image_name = QtWidgets.QLineEdit(fname, self)
        hbox.addWidget(self.image_name, stretch=1)
        button = QtWidgets.QPushButton("Save", self)
        button.clicked.connect(self.save_plot)
        hbox.addWidget(button)

        '''
        line = ttk.Frame(config_frame)
        line.pack(fill=Tk.X)
        ttk.Button(line,text="-",command=self.decrement_layer).pack(side=Tk.LEFT)
        self.slider = ttk.Scale(line,from_=0,to=0,orient=Tk.HORIZONTAL,length=300,variable=self.layernum)
        self.slider.pack(side=Tk.LEFT,fill=Tk.X,expand=1)
        self.slider.bind('<ButtonRelease-1>', lambda event: self.replot(zoom='current'))
        ttk.Button(line,text="+",command=self.increment_layer).pack(side=Tk.LEFT)
        '''

        hbox = QtWidgets.QHBoxLayout()
        vbox.addLayout(hbox)
        label = QtWidgets.QLabel('Range: ', self)
        hbox.addWidget(label)
        self.rangemin = QtWidgets.QLineEdit('0.', self)
        self.rangemin.setFixedWidth(80)
        hbox.addWidget(self.rangemin)
        self.rangemax = QtWidgets.QLineEdit('10.', self)
        self.rangemax.setFixedWidth(80)
        hbox.addWidget(self.rangemax)
        self.rangelock = QtWidgets.QCheckBox('Lock', self)
        hbox.addWidget(self.rangelock)
        hbox.addStretch(1)

        '''
        line = ttk.Frame(config_frame)
        line.pack(fill=Tk.X)
        ttk.Button(line,text="Preprocess",command=self.preprocess).pack(side=Tk.LEFT)
        ttk.Button(line,text="Quit",command=self.master.quit).pack(side=Tk.LEFT)
        '''

        # Canvas frame
        canvas_frame = QtWidgets.QWidget(self)
        vbox = QtWidgets.QVBoxLayout()
        canvas_frame.setLayout(vbox)
        splitter.addWidget(canvas_frame)

        hbox = QtWidgets.QHBoxLayout()
        vbox.addLayout(hbox)
        self.current_fname = QtWidgets.QLabel(self.merge_fname, self)
        hbox.addWidget(self.current_fname)
        hbox.addStretch(1)
        self.project_flag = QtWidgets.QCheckBox('Projection', self)
        self.project_flag.stateChanged.connect(lambda: self.replot(zoom='current'))
        hbox.addWidget(self.project_flag)
        hbox.addStretch(1)
        button = QtWidgets.QPushButton("Prev", self)
        button.clicked.connect(self.next_angle)
        hbox.addWidget(button)
        self.current_angle = QtWidgets.QLabel('XY', self)
        hbox.addWidget(self.current_angle)
        button = QtWidgets.QPushButton("Next", self)
        button.clicked.connect(self.next_angle)
        hbox.addWidget(button)

        self.fig = plt.figure(figsize=(7,7))
        self.fig.subplots_adjust(left=0.0, bottom=0.00, right=0.99, top=0.95, wspace=0.0)
        self.canvas = FigureCanvas(self.fig)
        self.canvas.draw()
        vbox.addWidget(self.canvas)

        self.show()
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
            if self.rangelock.text() == 0:
                self.rangemaxstr.set('%.1e' % 1)
        elif ext_string == '.cpx':
            self.typestr = 'complex64'
        else:
            print("Did not understand data type from extension. Defaulting to float.")
            self.typestr = 'f4'

    def parse_vol(self):
        if not os.path.isfile(self.current_fname.text()):
            if self.current_fname.text() != '':
                print("Unable to open", self.current_fname.text())
            return
        self.parse_extension(self.current_fname.text())
        self.vol = np.fromfile(self.current_fname.text(), dtype=self.typestr)
        size = int(round(self.vol.size**(1/3.)))
        self.size = size
        self.vol_size = size
        self.vol = self.vol.reshape(size,size,size)
        
        if self.typestr == 'complex64':
            self.vol = np.square(np.absolute(self.vol))
        if not self.rangelock.isChecked():
            rmax = self.vol[self.vol>0].mean()
            rmax = 5*self.vol[(self.vol>0.01*rmax) & (self.vol<100*rmax)].mean()
            self.rangemax.setText('%.1e' % rmax)
            '''
        self.slider.configure(from_=0,to=self.size-1)
        if not self.vol_image_exists:
            self.layernum.set(self.size/2)
            self.radiusmin.set('%d' % (self.size/2/2))
            self.radiusmax.set('%d' % (self.size))
            self.scaleradmin.set('%d' % (self.size/2/2*0.9))
            self.scaleradmax.set('%d' % (self.size/2/2*1.1))
        '''
        self.old_fname = self.current_fname.text()

    def parse_map(self):
        with open(self.current_fname.text(), 'rb') as f:
            grid = np.fromfile(f, '=i4', count=3) # Grid size
            f.seek(64, 0)
            ordering = np.fromfile(f, '=i4', count=3) # Axis ordering
            grid = grid[ordering-1]
            nx, ny, nz = tuple(grid)
            f.seek(1024, 0) # End of header
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
        self.old_fname = self.current_fname.text()

    def plot_slices(self, layernum, space=None, zoom=False, slices=None):
        if space is None:
            space = self.space
        self.image_name.setText('images/' + os.path.splitext(os.path.basename(self.current_fname.text()))[0] + '.png')
        rangemax = float(self.rangemax.text())
        rangemin = float(self.rangemin.text())
        project = (self.project_flag.isChecked())
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
                minx = self.size/3
                maxx = 2*minx
                if project:
                    rangemax *= self.size
                    a = self.vol[:,minx:maxx,minx:maxx].sum(0)
                    b = self.vol[minx:maxx,:,minx:maxx].sum(1)
                    c = self.vol[minx:maxx,minx:maxx,:].sum(2)
                else:
                    a = self.vol[layernum,minx:maxx,minx:maxx]
                    b = self.vol[minx:maxx,layernum,minx:maxx]
                    c = self.vol[minx:maxx,minx:maxx,layernum]
            else:
                if project:
                    rangemax *= self.size
                    a = self.vol.sum(0)
                    b = self.vol.sum(1)
                    c = self.vol.sum(2)
                else:
                    a = self.vol[layernum,:,:]	
                    b = self.vol[:,layernum,:]	
                    c = self.vol[:,:,layernum]
        
        self.fig.clear()
        s = self.fig.add_subplot(111)
        if str(self.current_angle.text()) == 'XY':
            view = a
        elif str(self.current_angle.text()) == 'XZ':
            view = b
        else:
            view = c
        cmap = self.color_map.checkedAction().text()
        s.matshow(view, vmin=rangemin, vmax=rangemax, cmap=cmap)
        plt.title(self.current_angle.text(), y=1.01)
        plt.axis('off')
        [a.remove() for a in list(set(s.findobj(patches.Circle)))]
        
        '''
        if self.circleflag.text() is 1: 
            rmin = float(self.radiusmin.text())
            rmax = float(self.radiusmax.text())
            s.add_artist(patches.Circle((self.size/2,self.size/2), rmin, ec='white', fc='none'))
            s.add_artist(patches.Circle((self.size/2,self.size/2), rmax, ec='white', fc='none'))
        
        if self.scaleradflag.text() is 1: 
            rmin = float(self.scaleradmin.text())
            rmax = float(self.scaleradmax.text())
            s.add_artist(patches.Circle((self.size/2,self.size/2), rmin, ec='white', fc='none', ls='dashed'))
            s.add_artist(patches.Circle((self.size/2,self.size/2), rmax, ec='white', fc='none', ls='dashed'))
        '''
        
        self.space = space
        self.canvas.draw()

    def replot(self, event=None, **kwargs):
        if self.map_image_exists:
            self.plot_map(**kwargs)
        elif self.vol_image_exists:
            self.plot_vol(**kwargs)
        else:
            self.plot_vol(**kwargs) # Default plotting merge

    def plot_vol(self, fname=None, event=None, force=False, sigma=False, **kwargs):
        if fname is not None:
            self.current_fname.setText(fname)
        if not self.vol_image_exists:
            self.parse_vol()
        elif self.old_fname != self.current_fname.text() or force == True: 
            print("Reparsing volume:", self.current_fname.text())
            self.parse_vol()
        
        '''
        if sigma and self.suppressflag.text() == 1:
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
                    print('Calculated self.rad')
                    self.sigma = np.exp(-8 * self.rad**2 / (c**2))
                    print('Calculated self.sigma')
                    self.rad.tofile('data/rad_%d.bin' % self.size)
                    self.sigma.tofile('data/sigma_%d.bin' % self.size)
            self.vol *= (1. - self.sigma)
        '''
        #self.plot_slices(self.layernum.text(), space='fourier', **kwargs)
        self.plot_slices(150, space='fourier', **kwargs)
        self.vol_image_exists = True
        self.map_image_exists = False

    def plot_map(self, event=None, force=False, sigma=False, **kwargs):
        self.current_fname.set(self.map_fname.text())
        if not self.map_image_exists:
            self.parse_map()
            if self.rangelock.text() == 0:
                self.rangemaxstr.set('%.1e' % 10)
        elif self.old_fname != self.current_fname.text() or force == True:
            print("Reparsing map:", self.current_fname.text())
            self.parse_map()
        self.plot_slices(self.layernum.text(), space='real', **kwargs)
        self.map_image_exists = True
        self.vol_image_exists = False

    def zero_outer(self, event=None):
        rmin = float(self.radiusmin.text())
        rmax = float(self.radiusmax.text())
        print('-'*80)
        os.system('./utils/zero_outer %s %d %d' % (self.merge_fname.text(), rmin, rmax))
        print('-'*80)
        
        if not self.zeroed:
            #ttk.Separator(self.merge_frame, orient=Tk.HORIZONTAL).pack(fill=Tk.X, padx=5, pady=5)
            zero_fname = os.path.splitext(self.merge_fname.text())[0] + '-zero.raw'
            line = ttk.Frame(self.merge_frame)
            line.pack(fill=Tk.X)
            ttk.Label(line,text='Zero-ed volume: ').pack(side=Tk.LEFT)
            ttk.Button(line,text=zero_fname,command=lambda: self.plot_vol(fname=zero_fname)).pack(side=Tk.LEFT)
        self.zeroed = True

    def calc_scale(self, event=None):
        rmin = float(self.scaleradmin.text())
        rmax = float(self.scaleradmax.text())
        mapnoext = os.path.splitext(os.path.basename(self.map_fname.text()))[0]
        sym_model = 'data/convert/'+mapnoext+'-sym.raw'
        cmd = './utils/calc_scale %s %s %d %d' % (sym_model, self.merge_fname.text(), rmin, rmax)
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

    def process_map(self, event=None):
        mapnoext = os.path.splitext(os.path.basename(self.map_fname.text()))[0]
        if os.path.isfile('data/convert/'+mapnoext+'.cpx') and not tkMessageBox.askyesno('Process Map', 'Found processed map output. Overwrite?', default=tkMessageBox.NO, icon=tkMessageBox.QUESTION, parent=self.master):
            with open('results/'+mapnoext+'.log', 'r') as f:
                words = f.read().split()
                warray = np.array(words)
                self.resedge.set(float(words[words.index('./utils/read_map')+2])/(self.vol_size/2))
                self.point_group.set(words[np.where(warray=='data/convert/'+mapnoext+'-srecon.raw')[0][0]+2])
                self.suppradstr.set('%.1f'%float(words[np.where(warray=='./utils/create_support')[0][-1]+2]))
                self.suppthreshstr.set('%.1f'%float(words[np.where(warray=='./utils/create_support')[0][-1]+3]))
        else:
            if self.resedge.text() is '':
                print('Need resolution at edge of volume')
                return
            print('-'*80)
            resedge = float(self.resedge.text())
            supp_rad = float(self.suppradstr.text())
            supp_thresh = float(self.suppthreshstr.text())
            command = './process_map.sh %s %d %f 0 %f %f %s' % (self.map_fname.text(), self.vol_size, resedge, supp_rad, supp_thresh, self.point_group.text())
            print(command)
            os.system(command)
            print('-'*80)
        
        if not self.processed_map:
            self.add_to_map_frame(mapnoext)
        self.processed_map = True

    def add_to_map_frame(self, mapnoext):
        prefix = 'data/convert/'+mapnoext
        line = ttk.Frame(self.map_frame)
        line.pack(fill=Tk.X)
        ttk.Label(line,text='Complex: ').pack(side=Tk.LEFT)
        ttk.Button(line,text=prefix+'.cpx',command=lambda: self.plot_vol(fname=prefix+'.cpx')).pack(side=Tk.LEFT)
        
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
        if self.resedge.text() is '':
            print('Need resolution at edge of volume')
            return
        mapnoext = os.path.splitext(os.path.basename(self.map_fname.text()))[0]
        resedge = float(self.resedge.text())
        supp_radius = float(self.suppradstr.text())
        supp_thresh = float(self.suppthreshstr.text())
        print('-'*80)
        command = './process_map.sh %s %d %f 1 %f %f %s' % (self.map_fname.text(), self.vol_size, resedge, supp_radius, supp_thresh, self.point_group.text())
        print(command)
        os.system(command)
        print('-'*80)
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
        with open(self.config_fname.text(), 'w') as f:
            f.write('[parameters]\n')
            f.write('size = %d\n' % self.vol_size)
            f.write('bragg_qmax = %f\n' % (float(self.radiusmin.text())/(self.vol_size/2)))
            f.write('scale_factor = %f\n' % self.scale_factor)
            f.write('num_threads = %d\n' % multiprocessing.cpu_count())
            f.write('point_group = %s\n' % self.point_group.text())
            
            mapnoext = os.path.splitext(os.path.basename(self.map_fname.text()))[0]
            f.write('\n[files]\n')
            f.write('intens_fname = %s\n' % (os.path.splitext(self.merge_fname.text())[0].rstrip()+'-zero.raw'))
            f.write('bragg_fname = %s\n' % ('data/convert/'+mapnoext+'.cpx'))
            f.write('support_fname = %s\n' % ('data/convert/'+mapnoext+'.supp'))
            #f.write('input_fname = %s\n')
            f.write('output_prefix = %s\n' % self.output_prefix.text())
            
            f.write('\n[algorithm]\n')
            f.write('# Algorithm choices: DM, HIO, RAAR, mod-DM, ER\n')
            f.write('# With beta = 1, all algorithms except ER are equivalent\n')
            f.write('# algorithm and avg_algorithm are space separated with alternating numbers and names\n')
            f.write('algorithm = 200 DM\n')
            f.write('avg_algorithm = 100 DM\n')
            f.write('beta = 1.\n')
            if self.positivity_flag.text() == 1:
                f.write('positivity = 1\n')
            if self.bgfitting_flag.text() == 1:
                f.write('bg_fitting = 1\n')
            if self.variation_flag.text() == 1:
                f.write('local_variation = 1\n')
            if self.histogram_flag.text() == 1:
                f.write('histogram = 1\n')
                f.write('hist_fname = data/3wu2_hist.dat\n')
        print('Generated %s:' % self.config_fname.text())
        os.system('cat %s' % self.config_fname.text())

    def launch_recon(self, event=None):
        pass

    def keep_checking(self, event=None):
        if self.checkflag.text() == 1:
            if self.fslices.text() == 0:
                self.current_fname.set(self.output_prefix.text()+'-slices.raw')
            else:
                self.current_fname.set(self.output_prefix.text()+'-fslices.raw')
            
            if os.path.isfile(self.current_fname.text()):
                done = False
                while not done:
                    s = np.fromfile(self.current_fname.text(), '=f4')
                    self.size = int(np.round((s.shape[0]/3)**0.5))
                    try:
                        s = s.reshape(3,self.size,self.size) 
                        done = True
                    except ValueError:
                        pass
                if self.fslices.text() == 0:
                    self.plot_slices(0, slices=s, zoom=True)
                else:
                    self.plot_slices(0, slices=s, zoom=False)
            self.master.after(1000, self.keep_checking)

    def preprocess(self, event=None):
        self.zero_outer()
        self.calc_scale()
        self.process_map()

    def increment_layer(self, event=None):
        self.layernum.set(min(self.layernum.text()+1, self.size))
        self.plot_slices(self.layernum.text(), zoom='current')

    def decrement_layer(self, event=None):
        self.layernum.set(max(self.layernum.text()-1, 0))
        self.plot_slices(self.layernum.text(), zoom='current')

    def next_angle(self):
        curr = self.angle_list.index(self.current_angle.text())
        curr = (curr + 1) % 3
        self.current_angle.setText(self.angle_list[curr])
        self.replot(zoom='current')

    def prev_angle(self):
        curr = self.angle_list.index(self.current_angle.text())
        curr = (curr - 1) % 3
        self.current_angle.setText(self.angle_list[curr])
        self.replot(zoom='current')

    def save_plot(self, event=None):
        self.fig.savefig(self.image_name.text(), bbox_inches='tight', dpi=150)
        print("Saved to", self.image_name.text())

    def quit_(self, event=None):
        self.master.quit()

    def keyPressEvent(self, event): # pylint: disable=C0103
        '''Override of default keyPress event handler'''
        key = event.key()
        mod = int(event.modifiers())

        if key == QtCore.Qt.Key_Return or key == QtCore.Qt.Key_Enter:
            self.replot(zoom='current')
        elif key == QtCore.Qt.Key_Up:
            self.increment_layer()
        elif key == QtCore.Qt.Key_Down:
            self.decrement_layer()
        elif QtGui.QKeySequence(mod+key) == QtGui.QKeySequence('Ctrl+Q'):
            self.close()
        elif QtGui.QKeySequence(mod+key) == QtGui.QKeySequence('Ctrl+S'):
            self.save_plot()
        elif QtGui.QKeySequence(mod+key) == QtGui.QKeySequence('Ctrl+M'):
            self.plot_vol(fname=self.merge_fname.text())
        elif QtGui.QKeySequence(mod+key) == QtGui.QKeySequence('Ctrl+N'):
            self.plot_map()
        else:
            event.ignore()

if __name__ == '__main__':
    app = QtWidgets.QApplication([])
    if len(sys.argv) > 2:
        GUI(sys.argv[1], sys.argv[2])
    else:
        GUI()
    sys.exit(app.exec_())
