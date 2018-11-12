#!/usr/bin/env python

from __future__ import print_function
import sys
import os
import time
import subprocess
import multiprocessing
from functools import partial
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
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

class GUIWorker(QtCore.QObject):
    finished = QtCore.pyqtSignal()
    returnval = QtCore.pyqtSignal(str)

    @QtCore.pyqtSlot(str, str, float, float)
    def calc_scale(self, model, merge, rmin, rmax):
        cmd = './utils/calc_scale %s %s %d %d' % (model, merge, rmin, rmax)
        output = subprocess.check_output(cmd.split(), shell=False)
        self.returnval.emit(output.split()[4].decode('utf-8'))
        self.finished.emit()

    @QtCore.pyqtSlot(str, float, float)
    def zero_outer(self, model, rmin, rmax):
        print('-'*80)
        subprocess.call(('./utils/zero_outer %s %d %d' % (model, rmin, rmax)).split())
        print('-'*80)
        self.returnval.emit('')
        self.finished.emit()

    @QtCore.pyqtSlot(str, int, float, bool, float, float, str)
    def process_map(self, map_fname, size, resedge, full_flag, supp_rad, supp_thresh, point_group):
        flag = int(full_flag)
        command = './process_map.sh %s %d %f %d %f %f %s' % (map_fname, size, resedge, flag, supp_rad, supp_thresh, point_group)
        print(command)
        subprocess.call(command.split())
        print('-'*80)
        mapnoext = os.path.splitext(os.path.basename(map_fname))[0]
        self.returnval.emit(mapnoext)
        self.finished.emit()

    @QtCore.pyqtSlot(str)
    def launch_recon(self, fname):
        cmd = './recon -c %s'%fname
        print('-'*80)
        subprocess.call(cmd.split())
        print('-'*80)

class GUI(QtWidgets.QMainWindow):
    def __init__(self, merge_fname='', map_fname=''):
        super(GUI, self).__init__()
        self.input_merge_fname = merge_fname
        self.input_map_fname = map_fname

        self.typestr = 'f4'
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

        self.checker = QtCore.QTimer(self)
        self.checker.timeout.connect(self.keep_checking)
        self.worker = GUIWorker()
        self.thread = QtCore.QThread()
        self.worker.moveToThread(self.thread)
        self.worker.finished.connect(self.thread.quit)

        self.init_UI()

    def init_UI(self):
        QtGui.QFontDatabase.addApplicationFont('Oxygen-Regular.ttf')
        QtGui.QFontDatabase.addApplicationFont('Kalam-Bold.ttf')
        with open('style.css', 'r')as f:
            self.setStyleSheet(f.read())
        #self.setWindowFlags(QtCore.Qt.FramelessWindowHint)
        self.setWindowTitle('ResEx Phasing GUI')
        #self.showMaximized()
        self.resize(1000, 700)
        overall = QtWidgets.QWidget()
        self.setCentralWidget(overall)
        layout = QtWidgets.QHBoxLayout(overall)
        layout.setContentsMargins(0, 0, 0, 0)
        self.splitter = QtWidgets.QSplitter(QtCore.Qt.Horizontal)
        layout.addWidget(self.splitter)

        # Config frame
        self.config_frame = QtWidgets.QWidget(self)
        #self.config_frame.setMinimumWidth(80)
        self.config_frame.setMinimumWidth(300)
        self.config_frame.setObjectName('config')
        vbox = QtWidgets.QVBoxLayout()
        self.config_frame.setLayout(vbox)
        self.splitter.addWidget(self.config_frame)

        hbox = QtWidgets.QHBoxLayout()
        vbox.addLayout(hbox)
        self.title_label = QtWidgets.QLabel('ResEx Phasing GUI', self)
        self.title_label.setObjectName('heading')
        hbox.addWidget(self.title_label)
        hbox.addStretch(1)
        '''
        self.collapse_button = QtWidgets.QPushButton('<', self)
        self.collapse_button.clicked.connect(self.toggle_config)
        hbox.addWidget(self.collapse_button)
        '''

        self.notebook = QtWidgets.QTabWidget()
        vbox.addWidget(self.notebook, stretch=2)

        vbox.addStretch(1)

        self.bottom_frame = QtWidgets.QWidget(self)
        vbox.addWidget(self.bottom_frame)
        
        hbox = QtWidgets.QHBoxLayout()
        vbox.addLayout(hbox)
        hbox.addStretch(1)
        #button = QtWidgets.QPushButton("Preprocess", self)
        #button.clicked.connect(self.preprocess)
        #hbox.addWidget(button)
        button = QtWidgets.QPushButton("Quit", self)
        button.clicked.connect(self.close)
        hbox.addWidget(button)

        # Bottom frame
        vbox = QtWidgets.QVBoxLayout()
        self.bottom_frame.setLayout(vbox)
        
        hbox = QtWidgets.QHBoxLayout()
        vbox.addLayout(hbox)
        label = QtWidgets.QLabel('Image name: ', self)
        hbox.addWidget(label)
        if self.input_merge_fname != '':
            fname = 'images/' + os.path.splitext(os.path.basename(self.input_merge_fname))[0] + '.png'
        else:
            fname = ''
        self.image_name = QtWidgets.QLineEdit(fname, self)
        hbox.addWidget(self.image_name, stretch=1)
        button = QtWidgets.QPushButton("Save", self)
        button.clicked.connect(self.save_plot)
        hbox.addWidget(button)

        hbox = QtWidgets.QHBoxLayout()
        vbox.addLayout(hbox)
        self.layer_slider = QtWidgets.QSlider(QtCore.Qt.Horizontal, self)
        self.layer_slider.setRange(0, 100)
        self.layer_slider.sliderMoved.connect(self.layer_slider_moved)
        self.layer_slider.sliderReleased.connect(lambda: self.replot(zoom='current'))
        hbox.addWidget(self.layer_slider, stretch=1)
        self.layernum = QtWidgets.QSpinBox(self)
        self.layernum.setValue(self.layer_slider.value())
        self.layernum.setMinimum(0)
        self.layernum.setMaximum(100)
        self.layernum.valueChanged.connect(self.layernum_changed)
        self.layernum.editingFinished.connect(self.layernum_changed)
        self.layernum.setFixedWidth(48)
        hbox.addWidget(self.layernum)

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

        # Canvas frame
        canvas_frame = QtWidgets.QWidget(self)
        canvas_frame.setObjectName('canvas')
        vbox = QtWidgets.QVBoxLayout()
        canvas_frame.setLayout(vbox)
        self.splitter.addWidget(canvas_frame)

        hbox = QtWidgets.QHBoxLayout()
        vbox.addLayout(hbox)
        self.current_fname = QtWidgets.QLabel(self.input_merge_fname, self)
        hbox.addWidget(self.current_fname)
        hbox.addStretch(1)
        button = QtWidgets.QPushButton("Prev", self)
        button.clicked.connect(self.next_angle)
        hbox.addWidget(button)
        self.current_angle = QtWidgets.QLabel('XY', self)
        hbox.addWidget(self.current_angle)
        button = QtWidgets.QPushButton("Next", self)
        button.clicked.connect(self.next_angle)
        hbox.addWidget(button)
        hbox.addStretch(1)
        self.project_flag = QtWidgets.QCheckBox('Projection', self)
        self.project_flag.stateChanged.connect(self.toggle_projection)
        hbox.addWidget(self.project_flag)

        self.fig = plt.figure(figsize=(7,7))
        self.fig.subplots_adjust(left=0.0, bottom=0.00, right=0.99, top=0.95, wspace=0.0)
        self.canvas = FigureCanvas(self.fig)
        self.canvas.draw()
        vbox.addWidget(self.canvas)

        self.gen_merge_tab()
        self.gen_map_tab()
        self.gen_recon_tab()
        self.splitter.setSizes([400,796])
        self.show()
        self.plot_vol()

    def gen_merge_tab(self, add=True):
        self.merge_tab = QtWidgets.QWidget()
        if add:
            self.notebook.addTab(self.merge_tab, 'Merge')
        vbox = QtWidgets.QVBoxLayout(self.merge_tab)

        hbox = QtWidgets.QHBoxLayout()
        vbox.addLayout(hbox)
        label = QtWidgets.QLabel('Merge File', self)
        hbox.addWidget(label)
        self.merge_fname = QtWidgets.QLineEdit(self.input_merge_fname, self)
        hbox.addWidget(self.merge_fname, stretch=1)
        button = QtWidgets.QPushButton('Plot Merge', self)
        button.clicked.connect(lambda :self.plot_vol(fname=self.merge_fname.text(), zoom=False))
        hbox.addWidget(button)

        hbox = QtWidgets.QHBoxLayout()
        vbox.addLayout(hbox)
        label = QtWidgets.QLabel('Cutoff radii:', self)
        hbox.addWidget(label)
        self.radiusmin = QtWidgets.QLineEdit('0', self)
        self.radiusmin.setFixedWidth(80)
        hbox.addWidget(self.radiusmin)
        self.radiusmax = QtWidgets.QLineEdit('200', self)
        self.radiusmax.setFixedWidth(80)
        hbox.addWidget(self.radiusmax)
        self.circleflag = QtWidgets.QCheckBox('Show', self)
        self.circleflag.stateChanged.connect(lambda: self.replot(zoom=False))
        hbox.addWidget(self.circleflag)
        hbox.addStretch(1)

        hbox = QtWidgets.QHBoxLayout()
        vbox.addLayout(hbox)
        label = QtWidgets.QLabel('Scale radii:', self)
        hbox.addWidget(label)
        self.scaleradmin = QtWidgets.QLineEdit('60', self)
        self.scaleradmin.setFixedWidth(80)
        hbox.addWidget(self.scaleradmin)
        self.scaleradmax = QtWidgets.QLineEdit('80', self)
        self.scaleradmax.setFixedWidth(80)
        hbox.addWidget(self.scaleradmax)
        self.scaleradflag = QtWidgets.QCheckBox('Show', self)
        self.scaleradflag.stateChanged.connect(lambda: self.replot(zoom=False))
        hbox.addWidget(self.scaleradflag)
        hbox.addStretch(1)
        
        hbox = QtWidgets.QHBoxLayout()
        vbox.addLayout(hbox)
        button = QtWidgets.QPushButton('Zero Outer', self)
        button.clicked.connect(self.zero_outer)
        hbox.addWidget(button)
        button = QtWidgets.QPushButton('Calc. Scale', self)
        button.clicked.connect(self.calc_scale)
        hbox.addWidget(button)
        button = QtWidgets.QPushButton('Reset', self)
        button.clicked.connect(self.reset_merge_tab)
        hbox.addWidget(button)
        
        vbox.addStretch(1)

    def gen_map_tab(self, add=True):
        self.map_tab = QtWidgets.QWidget()
        if add:
            self.notebook.addTab(self.map_tab, 'Map')
        vbox = QtWidgets.QVBoxLayout(self.map_tab)

        hbox = QtWidgets.QHBoxLayout()
        vbox.addLayout(hbox)
        label = QtWidgets.QLabel('Map File', self)
        hbox.addWidget(label)
        self.map_fname = QtWidgets.QLineEdit(self.input_map_fname, self)
        hbox.addWidget(self.map_fname, stretch=1)
        button = QtWidgets.QPushButton('Plot Map', self)
        button.clicked.connect(self.plot_map)
        hbox.addWidget(button)

        hbox = QtWidgets.QHBoxLayout()
        vbox.addLayout(hbox)
        label = QtWidgets.QLabel('Res. at edge (A):', self)
        hbox.addWidget(label)
        self.resedge = QtWidgets.QLineEdit('2.0', self)
        #self.resedge.setFixedWidth(60)
        hbox.addWidget(self.resedge)
        label = QtWidgets.QLabel('Point group:', self)
        hbox.addWidget(label)
        self.point_group = QtWidgets.QLineEdit('222', self)
        #self.point_group.setFixedWidth(60)
        hbox.addWidget(self.point_group)
        hbox.addStretch(1)

        hbox = QtWidgets.QHBoxLayout()
        vbox.addLayout(hbox)
        button = QtWidgets.QPushButton('Process Map', self)
        button.clicked.connect(self.process_map)
        hbox.addWidget(button)
        self.reset_button = QtWidgets.QPushButton('Reset', self)
        self.reset_button.clicked.connect(self.reset_map_tab)
        self.reset_button.setEnabled(False)
        hbox.addWidget(self.reset_button)
        hbox.addStretch(1)

        vbox.addStretch(1)

    def gen_recon_tab(self):
        self.recon_tab = QtWidgets.QWidget()
        self.notebook.addTab(self.recon_tab, 'Recon')
        vbox = QtWidgets.QVBoxLayout(self.recon_tab)

        hbox = QtWidgets.QHBoxLayout()
        vbox.addLayout(hbox)
        label = QtWidgets.QLabel('Output Prefix:', self)
        hbox.addWidget(label)
        self.output_prefix = QtWidgets.QLineEdit('data/recon/output', self)
        hbox.addWidget(self.output_prefix, stretch=1)

        hbox = QtWidgets.QHBoxLayout()
        vbox.addLayout(hbox)
        label = QtWidgets.QLabel('Config file name:', self)
        hbox.addWidget(label)
        self.config_fname = QtWidgets.QLineEdit('config.ini', self)
        hbox.addWidget(self.config_fname, stretch=1)

        hbox = QtWidgets.QHBoxLayout()
        vbox.addLayout(hbox)
        self.bgfitting_flag = QtWidgets.QCheckBox('BG Fitting', self)
        hbox.addWidget(self.bgfitting_flag)
        self.variation_flag = QtWidgets.QCheckBox('Variation support', self)
        hbox.addWidget(self.variation_flag)
        self.positivity_flag = QtWidgets.QCheckBox('Positivity', self)
        hbox.addWidget(self.positivity_flag)
        hbox.addStretch(1)

        hbox = QtWidgets.QHBoxLayout()
        vbox.addLayout(hbox)
        button = QtWidgets.QPushButton('Generate', self)
        button.clicked.connect(self.gen_config)
        hbox.addWidget(button)

        hbox = QtWidgets.QHBoxLayout()
        vbox.addLayout(hbox)
        button = QtWidgets.QPushButton('Launch Recon', self)
        button.clicked.connect(self.launch_recon)
        hbox.addWidget(button)
        self.checkflag = QtWidgets.QCheckBox('Keep Checking', self)
        self.checkflag.stateChanged.connect(self.keep_checking)
        hbox.addWidget(self.checkflag)
        self.fslices= QtWidgets.QCheckBox('Fourier', self)
        hbox.addWidget(self.fslices)
        
        vbox.addStretch(1)
        self.added_recon_tab = True

    def reset_merge_tab(self, event=None):
        self.zeroed = False
        self.calculated_scale = False
        self.notebook.removeTab(0)
        #self.merge_tab.delete()
        self.gen_merge_tab(add=False)
        self.notebook.insertTab(0, self.merge_tab, 'Merge')
        self.notebook.setCurrentIndex(0)
        self.parse_vol(reset=True)
        self.zeroed = False
        self.calculated_scale = False

    def reset_map_tab(self, event=None):
        self.processed_map = True
        self.notebook.removeTab(1)
        #self.map_tab.delete()
        self.map_tab = None
        self.gen_map_tab(add=False)
        self.notebook.insertTab(1, self.map_tab, 'Map')
        self.notebook.setCurrentIndex(1)
        self.processed_map = False
        self.reset_map.setEnabled(False)

    def add_to_map_tab(self, mapnoext):
        if self.processed_map:
            return
        prefix = 'data/convert/'+mapnoext
        vbox = self.map_tab.layout()
        vbox.removeItem(vbox.takeAt(vbox.count()-1))

        hbox = vbox.takeAt(vbox.count()-1).layout()
        hbox.removeItem(hbox.takeAt(hbox.count()-1))
        self.suppressflag = QtWidgets.QCheckBox('Suppress low-q', self)
        self.suppressflag.stateChanged.connect(lambda: self.replot(zoom='current', sigma=True))
        hbox.addWidget(self.suppressflag)
        hbox.addStretch(1)
        vbox.addLayout(hbox)

        grid = QtWidgets.QGridLayout()
        vbox.addLayout(grid)
        label = QtWidgets.QLabel('Complex:', self)
        grid.addWidget(label, 0, 0)
        button = QtWidgets.QPushButton(os.path.basename(prefix + '.cpx'), self)
        button.clicked.connect(lambda: self.plot_vol(fname=prefix + '.cpx', sigma=True))
        grid.addWidget(button, 0, 1)
        label = QtWidgets.QLabel('Symmetrized:', self)
        grid.addWidget(label, 1, 0)
        button = QtWidgets.QPushButton(os.path.basename(prefix + '-sym.raw'), self)
        button.clicked.connect(lambda: self.plot_vol(fname=prefix + '-sym.raw', sigma=True))
        grid.addWidget(button, 1, 1)
        label = QtWidgets.QLabel('Density:', self)
        grid.addWidget(label, 2, 0)
        button = QtWidgets.QPushButton(os.path.basename(prefix + '-srecon.raw'), self)
        button.clicked.connect(lambda: self.plot_vol(fname=prefix + '-srecon.raw', zoom=True))
        grid.addWidget(button, 2, 1)
        label = QtWidgets.QLabel('Support:', self)
        grid.addWidget(label, 3, 0)
        button = QtWidgets.QPushButton(os.path.basename(prefix + '.supp'), self)
        button.clicked.connect(lambda: self.plot_vol(fname=prefix + '.supp', zoom=True, interpolation=None))
        grid.addWidget(button, 3, 1)

        hbox = QtWidgets.QHBoxLayout()
        vbox.addLayout(hbox)
        label = QtWidgets.QLabel('Support: Conv. radius = ', self)
        hbox.addWidget(label)
        self.suppradstr = QtWidgets.QLineEdit('3', self)
        self.suppradstr.setFixedWidth(40)
        hbox.addWidget(self.suppradstr)
        label = QtWidgets.QLabel('vox. Threshold = ', self)
        hbox.addWidget(label)
        self.suppthreshstr = QtWidgets.QLineEdit('1', self)
        self.suppthreshstr.setFixedWidth(40)
        hbox.addWidget(self.suppthreshstr)
        hbox.addStretch(1)

        hbox = QtWidgets.QHBoxLayout()
        vbox.addLayout(hbox)
        button = QtWidgets.QPushButton('Update support', self)
        button.clicked.connect(self.reprocess_map)
        hbox.addWidget(button)
        hbox.addStretch(1)

        vbox.addStretch(1)

        self.cleanup_thread()
        self.reset_button.setEnabled(True)
        self.processed_map = True

    def set_config_size(self, value):
        sizes = self.splitter.sizes()
        w = 80 + 300.*value
        self.splitter.setSizes([w, sum(sizes)-w])

    def toggle_config(self):
        w = self.splitter.sizes()[0]
        self.tl = QtCore.QTimeLine()
        self.tl.setCurveShape(QtCore.QTimeLine.EaseOutCurve)
        self.tl.valueChanged.connect(self.set_config_size)
        self.tl.setDuration(1000)
        if w > 250:
            self.tl.setDirection(QtCore.QTimeLine.Backward)
        self.tl.start()

        if w > 250:
            self.collapse_button.setText('>')
            self.notebook.hide()
            self.bottom_frame.hide()
            self.title_label.setText('RE')
        else:
            self.collapse_button.setText('<')
            self.notebook.show()
            self.bottom_frame.show()
            self.title_label.setText('ResEx Phasing GUI')

    def parse_extension(self, filename):
        ext_string = os.path.splitext(os.path.basename(filename))[1]
        
        if ext_string == '.raw':
            self.typestr = 'f4'
        elif ext_string == '.bin':
            self.typestr = 'f8'
        elif ext_string == '.supp':
            self.typestr = 'uint8'
            if not self.rangelock.isChecked():
                self.rangemax.setText('%.1e' % 1)
        elif ext_string == '.cpx':
            self.typestr = 'complex64'
        else:
            print("Did not understand data type from extension. Defaulting to float.")
            self.typestr = 'f4'

    def parse_vol(self, reset=False):
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
            self.autoset_rangemax(self.vol)
        self.layer_slider.setRange(0, size-1)
        if reset or not self.vol_image_exists or self.layernum.maximum() != self.size-1:
            self.layernum.setMaximum(self.size-1)
            self.layer_slider.setValue(self.size//2)
            self.layer_slider_moved(self.size//2)
            self.radiusmin.setText('%d' % (self.size//2//2))
            self.radiusmax.setText('%d' % (self.size))
            self.scaleradmin.setText('%d' % (self.size//2//2*0.9))
            self.scaleradmax.setText('%d' % (self.size//2//2*1.1))
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
        centralsum = (np.abs(vol[:,:,nz//2]).sum() + np.abs(vol[:,ny//2]).sum() + np.abs(vol[nx//2]).sum())/ 3.
        if edgesum > centralsum:
            vol = np.roll(vol, nx//2, axis=0)
            vol = np.roll(vol, ny//2, axis=1)
            vol = np.roll(vol, nz//2, axis=2)
        s = max(nx, ny, nz)
        self.size = s
        self.vol = np.pad(vol, (((s-nx)//2,s-nx-(s-nx)//2),((s-ny)//2,s-ny-(s-ny)//2),((s-nz)//2,s-nz-(s-nz)//2)), mode='constant', constant_values=0)
        self.layer_slider.setRange(0, self.size-1)
        if not self.map_image_exists:
            self.layernum.setValue(self.size//2)
        self.old_fname = self.current_fname.text()

    def plot_slices(self, layernum, space=None, zoom=False, slices=None, interpolation='gaussian'):
        if space is None:
            space = self.space
        self.image_name.setText('images/' + os.path.splitext(os.path.basename(self.current_fname.text()))[0] + '.png')
        project = (self.project_flag.isChecked())
        if self.vol is None and slices is None:
            print('Nothing to plot')
            return

        if zoom == 'current':
            zoom = self.zoomed
        else:
            self.zoomed = zoom
        
        if zoom:
            minx = self.size//3
            maxx = 2*minx
        else:
            minx = 0
            maxx = None
        if slices is not None:
            a, b, c = tuple(slices[:,minx:maxx,minx:maxx])
        else:
            if project:
                a = self.vol[:,minx:maxx,minx:maxx].mean(0)
                b = self.vol[minx:maxx,:,minx:maxx].mean(1)
                c = self.vol[minx:maxx,minx:maxx,:].mean(2)
            else:
                a = self.vol[layernum,minx:maxx,minx:maxx]
                b = self.vol[minx:maxx,layernum,minx:maxx]
                c = self.vol[minx:maxx,minx:maxx,layernum]

        self.fig.clear()
        s = self.fig.add_subplot(111)
        if str(self.current_angle.text()) == 'XY':
            view = a
        elif str(self.current_angle.text()) == 'XZ':
            view = b
        else:
            view = c

        if project or slices is not None:
            self.autoset_rangemax(view, maxval=True)
        rangemax = float(self.rangemax.text())
        rangemin = float(self.rangemin.text())
        cmap = 'jet'
        s.matshow(view, vmin=rangemin, vmax=rangemax, cmap=cmap, interpolation=interpolation)
        plt.axis('off')

        [a.remove() for a in list(set(s.findobj(patches.Circle)))]
        if self.circleflag.isChecked(): 
            rmin = float(self.radiusmin.text())
            rmax = float(self.radiusmax.text())
            s.add_artist(patches.Circle((self.size//2,self.size//2), rmin, ec='white', fc='none'))
            s.add_artist(patches.Circle((self.size//2,self.size//2), rmax, ec='white', fc='none'))

        if self.scaleradflag.isChecked(): 
            rmin = float(self.scaleradmin.text())
            rmax = float(self.scaleradmax.text())
            s.add_artist(patches.Circle((self.size//2,self.size//2), rmin, ec='white', fc='none', ls='dashed'))
            s.add_artist(patches.Circle((self.size//2,self.size//2), rmax, ec='white', fc='none', ls='dashed'))

        self.space = space
        self.canvas.draw()

    def replot(self, event=None, **kwargs):
        if self.map_image_exists:
            self.plot_map(**kwargs)
        elif self.vol_image_exists:
            self.plot_vol(**kwargs)
        else:
            self.plot_vol(**kwargs) # Default plotting merge
        self.cleanup_thread()

    def plot_vol(self, event=None, fname=None, force=False, sigma=False, **kwargs):
        if fname is not None:
            self.current_fname.setText(fname)
        if not self.vol_image_exists:
            self.parse_vol()
        elif self.old_fname != self.current_fname.text() or force == True: 
            print("Reparsing volume:", self.current_fname.text())
            self.parse_vol()
        
        if sigma and self.suppressflag.isChecked():
            c = self.vol.shape[0] // 2
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
        self.plot_slices(self.layernum.value(), space='fourier', **kwargs)
        self.vol_image_exists = True
        self.map_image_exists = False

    def plot_map(self, event=None, force=False, sigma=False, **kwargs):
        self.current_fname.setText(self.map_fname.text())
        if not self.map_image_exists:
            self.parse_map()
            if not self.rangelock.isChecked():
                self.rangemax.setText('%.1e' % 10)
        elif self.old_fname != self.current_fname.text() or force == True:
            print("Reparsing map:", self.current_fname.text())
            self.parse_map()
        self.plot_slices(self.layernum.value(), space='real', **kwargs)
        self.map_image_exists = True
        self.vol_image_exists = False

    def zero_outer(self, event=None):
        fname = self.merge_fname.text()
        rmin = float(self.radiusmin.text())
        rmax = float(self.radiusmax.text())
        self.thread.started.connect(partial(self.worker.zero_outer, fname, rmin, rmax))
        self.worker.returnval.connect(self.write_zero_line)
        self.thread.start()

    def write_zero_line(self, val):
        if not self.zeroed:
            vbox = self.merge_tab.layout()
            vbox.removeItem(vbox.takeAt(vbox.count()-1))
            hbox = QtWidgets.QHBoxLayout()
            vbox.addLayout(hbox)
            vbox.addStretch(1)
            
            zero_fname = os.path.splitext(self.merge_fname.text())[0] + '-zero.raw'
            label = QtWidgets.QLabel('Zero-ed volume:', self)
            hbox.addWidget(label)
            button = QtWidgets.QPushButton(zero_fname, self)
            button.clicked.connect(lambda: self.plot_vol(fname=zero_fname))
            hbox.addWidget(button)
            hbox.addStretch(1)

        self.cleanup_thread()
        self.zeroed = True

    def calc_scale(self, event=None):
        fname = self.merge_fname.text()
        rmin = float(self.scaleradmin.text())
        rmax = float(self.scaleradmax.text())
        mapnoext = os.path.splitext(os.path.basename(self.map_fname.text()))[0]
        map_fname = 'data/convert/'+mapnoext+'-sym.raw'
        self.thread.started.connect(partial(self.worker.calc_scale, map_fname, fname, rmin, rmax))
        self.worker.returnval.connect(self.write_scale_line)
        self.thread.start()

    def write_scale_line(self, val):
        self.scale_factor = float(val)
        if not self.calculated_scale:
            vbox = self.merge_tab.layout()
            vbox.removeItem(vbox.takeAt(vbox.count()-1))
            hbox = QtWidgets.QHBoxLayout()
            vbox.addLayout(hbox)
            vbox.addStretch(1)
            
            self.scale_label = QtWidgets.QLabel('Scale factor = %.6e'%self.scale_factor, self)
            hbox.addWidget(self.scale_label)
            hbox.addStretch(1)  
        else:
            self.scale_label.setText('Scale factor = %.6e' % self.scale_factor)

        self.cleanup_thread()
        self.calculated_scale = True

    def process_map(self, event=None):
        mapnoext = os.path.splitext(os.path.basename(self.map_fname.text()))[0]
        if os.path.isfile('data/convert/'+mapnoext+'.cpx') and QtWidgets.QMessageBox.question(self, 'Process Map', 'Found processed map output. Overwrite?', QtWidgets.QMessageBox.Yes, QtWidgets.QMessageBox.No) == QtWidgets.QMessageBox.Yes:
            if self.resedge.text() is '':
                print('Need resolution at edge of volume')
                return
            print('-'*80)
            map_fname = self.map_fname.text()
            size = self.vol_size
            resedge = float(self.resedge.text())
            if self.processed_map:
                supp_rad = float(self.suppradstr.text())
                supp_thresh = float(self.suppthreshstr.text())
            else:
                supp_rad = 3.
                supp_thresh = 1.
            point_group = self.point_group.text()
            self.thread.started.connect(partial(self.worker.process_map, map_fname, size, resedge, False, supp_rad, supp_thresh, point_group))
            self.worker.returnval.connect(self.add_to_map_tab)
            self.thread.start()
        else:
            self.add_to_map_tab(mapnoext)
            with open('results/'+mapnoext+'.log', 'r') as f:
                words = f.read().split()
                warray = np.array(words)
                self.resedge.setText(str(float(words[words.index('./utils/read_map')+2])/(self.vol_size//2)))
                self.point_group.setText(words[np.where(warray=='data/convert/'+mapnoext+'-srecon.raw')[0][0]+2])
                self.suppradstr.setText('%.1f'%float(words[np.where(warray=='./utils/create_support')[0][-1]+2]))
                self.suppthreshstr.setText('%.1f'%float(words[np.where(warray=='./utils/create_support')[0][-1]+3]))

    def reprocess_map(self, event=None):
        if self.resedge.text() is '':
            print('Need resolution at edge of volume')
            return
        map_fname = self.map_fname.text()
        size = self.vol_size
        resedge = float(self.resedge.text())
        supp_rad = float(self.suppradstr.text())
        supp_thresh = float(self.suppthreshstr.text())
        point_group = self.point_group.text()
        self.thread.started.connect(partial(self.worker.process_map, map_fname, size, resedge, True, supp_rad, supp_thresh, point_group))
        self.worker.returnval.connect(lambda: self.replot(force=True, zoom='current'))
        self.thread.start()

    def gen_config(self, event=None):
        with open(self.config_fname.text(), 'w') as f:
            f.write('[parameters]\n')
            f.write('size = %d\n' % self.vol_size)
            f.write('bragg_qmax = %f\n' % (float(self.radiusmin.text())/(self.vol_size//2)))
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
            '''
            if self.histogram_flag.text() == 1:
                f.write('histogram = 1\n')
                f.write('hist_fname = data/3wu2_hist.dat\n')
            '''
        print('Generated %s:' % self.config_fname.text())
        with open(self.config_fname.text(), 'r') as f:
            print(f.read())

    def cleanup_thread(self):
        if self.thread.receivers(self.thread.started) > 0:
            self.thread.started.disconnect()
        if self.worker.receivers(self.worker.returnval) > 0:
            self.worker.returnval.disconnect()

    def launch_recon(self, event=None):
        self.thread.started.connect(partial(self.worker.launch_recon, self.config_fname.text()))
        self.worker.returnval.connect(self.cleanup_thread)
        self.thread.start()

    def keep_checking(self, event=None):
        if self.checkflag.isChecked():
            self.update_slices()
            self.checker.start(1000)
        else:
            self.checker.stop()

    def update_slices(self):
        prefix = self.output_prefix.text()
        log_fname = prefix+'-log.dat'
        try:
            with open(log_fname, 'r') as f:
                lines = f.readlines()
                iternum = int(lines[-1].split()[0])
        except (FileNotFoundError, ValueError):
            return

        if self.fslices.isChecked():
            self.current_fname.setText(prefix+'-fslices/%.4d.raw'%iternum)
        else:
            self.current_fname.setText(prefix+'-slices/%.4d.raw'%iternum)

        if os.path.isfile(self.current_fname.text()):
            done = False
            while not done:
                s = np.fromfile(self.current_fname.text(), '=f4')
                self.size = int(np.round((s.shape[0]//3)**0.5))
                try:
                    s = s.reshape(3,self.size,self.size) 
                    done = True
                except ValueError:
                    pass
            if self.fslices.isChecked():
                self.plot_slices(0, slices=s, zoom=False)
            else:
                self.plot_slices(0, slices=s, zoom=True)

    def preprocess(self, event=None):
        self.zero_outer()
        self.calc_scale()
        self.process_map()

    def autoset_rangemax(self, arr, maxval=False):
        if maxval:
            rmax = arr.max()
        else:
            rmax = arr[arr>0].mean()
            rmax = 5*arr[(arr>0.01*rmax) & (arr<100*rmax)].mean()
        self.rangemax.setText('%.1e' % rmax)
        return rmax

    def layer_slider_moved(self, value):
        self.layernum.setValue(value)

    def layernum_changed(self, value=None):
        if value == self.layernum.value():
            self.layer_slider.setValue(value)
        elif value is None:
            self.replot(zoom='current')
    
    def toggle_projection(self):
        if self.project_flag.isChecked():
            self.layer_slider.setEnabled(False)
            self.layernum.setEnabled(False)
        else:
            self.layer_slider.setEnabled(True)
            self.layernum.setEnabled(True)
            self.autoset_rangemax(self.vol)
        self.replot(zoom='current')
    
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
