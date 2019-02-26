from __future__ import print_function
import os
import sys
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
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
import mrcfile

class CanvasPanel(QtWidgets.QWidget):
    def __init__(self, parent, **kwargs):
        super(CanvasPanel, self).__init__(parent, **kwargs)
        self.parent = parent

        self.setObjectName('canvas')
        self.setAttribute(QtCore.Qt.WA_StyledBackground)
        self.typestr = 'f4'
        self.typedict = {'0': 'i1', '2': 'f4', '4': 'c8'}
        self.size = None
        self.vol_size = None
        self.vol = None
        self.rad = None
        self.old_fname = None
        self.space = None
        self.zoomed = False
        self.vol_slices = False
        self.map_image_exists = False
        self.vol_image_exists = False
        self.angle_list = ['XY', 'XZ', 'YZ']

        self.init_UI()

    def init_UI(self):
        vbox = QtWidgets.QVBoxLayout()
        self.setLayout(vbox)

        hbox = QtWidgets.QHBoxLayout()
        vbox.addLayout(hbox)
        self.current_fname = QtWidgets.QLabel('', self)
        self.current_fname.setToolTip('Map file currently being displayed')
        hbox.addWidget(self.current_fname)
        hbox.addStretch(1)
        button = QtWidgets.QPushButton("Prev", self)
        button.clicked.connect(self.next_angle)
        button.setToolTip('Previous orthogonal view')
        hbox.addWidget(button)
        self.current_angle = QtWidgets.QLabel('XY', self)
        self.current_angle.setToolTip('Current view axes (h, v)')
        hbox.addWidget(self.current_angle)
        button = QtWidgets.QPushButton("Next", self)
        button.clicked.connect(self.next_angle)
        button.setToolTip('Next orthogonal view')
        hbox.addWidget(button)
        hbox.addStretch(1)
        self.project_flag = QtWidgets.QCheckBox('Projection', self)
        self.project_flag.stateChanged.connect(self.toggle_projection)
        self.project_flag.setToolTip('When checked, shows projection (sum along third axis) rather than slice')
        hbox.addWidget(self.project_flag)

        self.figure = plt.figure(figsize=(7, 7))
        self.figure.subplots_adjust(left=0.01, bottom=0.01, right=0.99, top=0.99, wspace=0.0)
        self.canvas = FigureCanvas(self.figure)
        self.canvas.draw()
        vbox.addWidget(self.canvas)

        self.show()

    def init_plot_controls(self, input_fname):
        widget = QtWidgets.QWidget()
        vbox = QtWidgets.QVBoxLayout()
        widget.setLayout(vbox)

        hbox = QtWidgets.QHBoxLayout()
        vbox.addLayout(hbox)
        label = QtWidgets.QLabel('Image name: ', self)
        hbox.addWidget(label)
        if input_fname != '':
            fname = 'images/' + os.path.splitext(os.path.basename(input_fname))[0] + '.png'
        else:
            fname = ''
        self.image_name = QtWidgets.QLineEdit(fname, self)
        hbox.addWidget(self.image_name, stretch=1)
        button = QtWidgets.QPushButton("Save", self)
        button.clicked.connect(self.save_plot)
        button.setToolTip('Save current view to file')
        hbox.addWidget(button)

        hbox = QtWidgets.QHBoxLayout()
        vbox.addLayout(hbox)
        self.layer_slider = QtWidgets.QSlider(QtCore.Qt.Horizontal, self)
        self.layer_slider.setRange(0, 100)
        self.layer_slider.sliderMoved.connect(self.layer_slider_moved)
        self.layer_slider.sliderReleased.connect(lambda: self.replot(zoom='current'))
        self.layer_slider.setToolTip('Change slice number in 3D stack')
        hbox.addWidget(self.layer_slider, stretch=1)
        self.layernum = QtWidgets.QSpinBox(self)
        self.layernum.setValue(self.layer_slider.value())
        self.layernum.setMinimum(0)
        self.layernum.setMaximum(100)
        self.layernum.valueChanged.connect(self.layernum_changed)
        self.layernum.editingFinished.connect(self.layernum_changed)
        self.layernum.setFixedWidth(48)
        self.layernum.setToolTip('Change slice number in 3D stack')
        hbox.addWidget(self.layernum)

        hbox = QtWidgets.QHBoxLayout()
        vbox.addLayout(hbox)
        label = QtWidgets.QLabel('Range: ', self)
        hbox.addWidget(label)
        self.rangemin = QtWidgets.QLineEdit('0.', self)
        self.rangemin.setFixedWidth(80)
        self.rangemin.setToolTip('Min value of color scale')
        hbox.addWidget(self.rangemin)
        self.rangemax = QtWidgets.QLineEdit('10.', self)
        self.rangemax.setFixedWidth(80)
        self.rangemax.setToolTip('Max value of color scale')
        hbox.addWidget(self.rangemax)
        self.autorange= QtWidgets.QCheckBox('Auto', self)
        self.autorange.setChecked(True)
        self.autorange.setToolTip('If checked, color range is auto-calculated from data')
        hbox.addWidget(self.autorange)
        hbox.addStretch(1)

        widget.show()
        return widget

    def parse_extension(self, filename):
        ext_string = os.path.splitext(os.path.basename(filename))[1]

        if ext_string == '.raw':
            self.typestr = 'f4'
        elif ext_string == '.bin':
            self.typestr = 'f8'
        elif ext_string == '.supp':
            self.typestr = 'uint8'
        elif ext_string == '.cpx':
            self.typestr = 'complex64'
        else:
            print("Did not understand data type from extension. Defaulting to float.")
            self.typestr = 'f4'

    def parse_vol(self, reset=False, noncube=False):
        if not os.path.isfile(self.current_fname.text()):
            if self.current_fname.text() != '':
                print("Unable to open", self.current_fname.text())
            return False
        with mrcfile.open(self.current_fname.text(), 'r') as f:
            self.vol = f.data
            nx, ny, nz = f.header.nx, f.header.ny, f.header.nz
            self.typestr = self.typedict[str(f.header.mode)]
        s = max(nx, ny, nz)
        self.size = s
        if not noncube:
            self.vol_size = s
        self.vol_slices = (nz <= 3)
        if self.vol_slices:
            print(self.current_fname.text(), 'is image stack')
        if not self.vol_slices:
            edgesum = (np.abs(self.vol[:,:,0]).sum() + np.abs(self.vol[:,0]).sum() + np.abs(self.vol[0]).sum()) / 3.
            centralsum = (np.abs(self.vol[:,:,nz//2]).sum() + np.abs(self.vol[:,ny//2]).sum() + np.abs(self.vol[nx//2]).sum()) / 3.
            if edgesum > centralsum:
                self.vol = np.fft.fftshift(self.vol)
            self.vol = np.pad(self.vol, (((s-nx)//2, s-nx-(s-nx)//2), ((s-ny)//2, s-ny-(s-ny)//2), ((s-nz)//2, s-nz-(s-nz)//2)), mode='constant', constant_values=0)

        if self.typestr == 'c8':
            self.vol = np.square(np.absolute(self.vol))
        if self.autorange.isChecked():
            self.autoset_rangemax(self.vol)
        self.layer_slider.setRange(0, self.size-1)
        if reset or not self.vol_image_exists or self.layernum.maximum() != self.size-1:
            self.layernum.setMaximum(self.size-1)
            self.layer_slider.setValue(self.size//2)
            self.layer_slider_moved(self.size//2)
        self.old_fname = self.current_fname.text()
        return True

    def parse_map(self):
        with open(self.current_fname.text(), 'rb') as f:
            grid = np.fromfile(f, '=i4', count=3) # Grid size
            f.seek(64, 0)
            ordering = np.fromfile(f, '=i4', count=3) # Axis ordering
            grid = grid[ordering-1]
            nx, ny, nz = tuple(grid)
            f.seek(1024, 0) # End of header
            vol = np.fromfile(f, '=f4', count=nx*ny*nz).reshape(nx, ny, nz)
        edgesum = (np.abs(vol[:, :, 0]).sum() + np.abs(vol[:, :, -1]).sum() + np.abs(vol[:, 0]).sum() + np.abs(vol[:, -1]).sum() + np.abs(vol[0]).sum() + np.abs(vol[-1]).sum()) / 6.
        centralsum = (np.abs(vol[:, :, nz//2]).sum() + np.abs(vol[:, ny//2]).sum() + np.abs(vol[nx//2]).sum())/ 3.
        if edgesum > centralsum:
            vol = np.roll(vol, nx//2, axis=0)
            vol = np.roll(vol, ny//2, axis=1)
            vol = np.roll(vol, nz//2, axis=2)
        s = max(nx, ny, nz)
        self.size = s
        self.vol = np.pad(vol, (((s-nx)//2, s-nx-(s-nx)//2), ((s-ny)//2, s-ny-(s-ny)//2), ((s-nz)//2, s-nz-(s-nz)//2)), mode='constant', constant_values=0)
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
            a, b, c = tuple(slices[:, minx:maxx, minx:maxx])
            self.vol_slices = True
        elif not self.vol_slices:
            if project:
                a = self.vol[:, minx:maxx, minx:maxx].mean(0)
                b = self.vol[minx:maxx, :, minx:maxx].mean(1)
                c = self.vol[minx:maxx, minx:maxx, :].mean(2)
            else:
                a = self.vol[layernum, minx:maxx, minx:maxx]
                b = self.vol[minx:maxx, layernum, minx:maxx]
                c = self.vol[minx:maxx, minx:maxx, layernum]
        else:
            return

        self.figure.clear()
        s = self.figure.add_subplot(111)
        if str(self.current_angle.text()) == 'XY':
            view = a
        elif str(self.current_angle.text()) == 'XZ':
            view = b
        else:
            view = c

        if slices is not None:
            self.autoset_rangemax(view, maxval=True)
        rangemax = float(self.rangemax.text())
        rangemin = float(self.rangemin.text())
        cmap = 'jet'
        s.matshow(view, vmin=rangemin, vmax=rangemax, cmap=cmap, interpolation=interpolation)
        plt.axis('off')

        self.space = space
        self.canvas.draw()

    def replot(self, event=None, **kwargs):
        if self.map_image_exists:
            self.plot_map(self.current_fname.text(), **kwargs)
        elif self.vol_image_exists:
            self.plot_vol(**kwargs)
        else:
            self.plot_vol(**kwargs) # Default plotting merge
        #self.cleanup_thread()

    def plot_vol(self, event=None, fname=None, force=False, sigma=False, **kwargs):
        parsed = False
        if fname is not None:
            self.current_fname.setText(fname)
        if not self.vol_image_exists:
            parsed = self.parse_vol()
        elif self.old_fname != self.current_fname.text() or force:
            print("Reparsing volume:", self.current_fname.text())
            parsed = self.parse_vol()

        if sigma:
            c = self.vol.shape[0] // 2
            if self.rad is None:
                if os.path.isfile('data/sigma_%d.bin' % self.size):
                    self.rad = np.fromfile('data/rad_%d.bin' % self.size, '=f8').reshape(self.size, self.size, self.size)
                    self.sigma = np.fromfile('data/sigma_%d.bin' % self.size, '=f8').reshape(self.size, self.size, self.size)
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
        return parsed

    def plot_map(self, fname, event=None, force=False, sigma=False, **kwargs):
        self.current_fname.setText(fname)
        if not self.map_image_exists:
            self.parse_vol(noncube=True)
            if self.autorange.isChecked():
                self.rangemax.setText('%.1e' % 10)
        elif self.old_fname != self.current_fname.text() or force:
            print("Reparsing map:", self.current_fname.text())
            self.parse_vol(noncube=True)
        self.plot_slices(self.layernum.value(), space='real', **kwargs)
        self.map_image_exists = True
        self.vol_image_exists = False

    def update_slices(self, prefix, fslices=False):
        log_fname = prefix+'-log.dat'
        try:
            with open(log_fname, 'r') as f:
                lines = f.readlines()
                iternum = int(lines[-1].split()[0])
        except (FileNotFoundError, ValueError):
            return

        if fslices:
            self.current_fname.setText(prefix+'-fslices/%.4d.ccp4'%iternum)
        else:
            self.current_fname.setText(prefix+'-slices/%.4d.ccp4'%iternum)

        if os.path.isfile(self.current_fname.text()):
            done = False
            while not done:
                s = np.fromfile(self.current_fname.text(), '=f4')
                self.size = int(np.round((s.shape[0]//3)**0.5))
                try:
                    s = s.reshape(3, self.size, self.size)
                    done = True
                except ValueError:
                    pass
            self.plot_slices(0, slices=s, zoom=(not fslices))

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

    def toggle_projection(self):
        if self.project_flag.isChecked():
            self.layer_slider.setEnabled(False)
            self.layernum.setEnabled(False)
            self.autoset_rangemax(self.vol.mean(0))
        else:
            self.layer_slider.setEnabled(True)
            self.layernum.setEnabled(True)
            self.autoset_rangemax(self.vol)
        self.replot(zoom='current')

    def layer_slider_moved(self, value):
        self.layernum.setValue(value)

    def layernum_changed(self, value=None):
        if value == self.layernum.value():
            self.layer_slider.setValue(value)
        elif value is None:
            self.replot(zoom='current')

    def autoset_rangemax(self, arr, maxval=False):
        if maxval:
            rmax = arr.max()
        else:
            rmax = arr[arr>0].mean()
            rmax = 5*arr[(arr>0.01*rmax) & (arr<100*rmax)].mean()
        self.rangemax.setText('%.1e' % rmax)
        return rmax

    def save_plot(self, event=None):
        self.figure.savefig(self.image_name.text(), bbox_inches='tight', dpi=150)
        print("Saved to", self.image_name.text())

