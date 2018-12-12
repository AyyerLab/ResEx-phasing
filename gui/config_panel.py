from __future__ import print_function
import sys
import os
import time
import subprocess
import multiprocessing
import numpy as np
import matplotlib
import matplotlib.patches as patches
try:
    from PyQt5 import QtCore, QtWidgets, QtGui # pylint: disable=import-error
    matplotlib.use('qt5agg')
    from matplotlib.backends.backend_qt5agg import FigureCanvas # pylint: disable=no-name-in-module
    os.environ['QT_API'] = 'pyqt5'
except ImportError:
    import sip
    sip.setapi('QString', 2)
    from PyQt4 import QtCore, QtGui # pylint: disable=import-error
    from PyQt4 import QtGui as QtWidgets # pylint: disable=import-error
    matplotlib.use('qt4agg')
    from matplotlib.backends.backend_qt4agg import FigureCanvas # pylint: disable=no-name-in-module
    os.environ['QT_API'] = 'pyqt'
from . import worker

class IniHighlighter(QtGui.QSyntaxHighlighter):
    def __init__(self, parent):
        super(IniHighlighter, self).__init__(parent)
        self.char_formats = [
            self.make_charformat('red'),
            self.make_charformat('darkGreen', bold=True),
            self.make_charformat('darkBlue'),
            self.make_charformat('gray'),
        ]
        patterns = [r'\[.*\]', r'^(.*)=', r'=([^\n]*)', r'#[^\n]*']
        self.rules = [QtCore.QRegExp(p) for p in patterns]

    def make_charformat(self, color, bold=False, italics=False):
        cf = QtGui.QTextCharFormat()
        cf.setForeground(QtGui.QBrush(QtGui.QColor(color)))
        if bold:
            cf.setFontWeight(QtGui.QFont.Bold)
        cf.setFontItalic(italics)
        return cf

    def highlightBlock(self, text):
        for rule, cf in zip(self.rules, self.char_formats):
            index = rule.indexIn(text, 0)
            while index >= 0:
                length = len(rule.cap(0))
                self.setFormat(index, length, cf)
                index = rule.indexIn(text, index+length)
        self.setCurrentBlockState(0)

class ConfigPanel(QtWidgets.QWidget):
    def __init__(self, parent, **kwargs):
        super(ConfigPanel, self).__init__(parent, **kwargs)
        self.parent = parent
        self.canvas_panel = parent.canvas_panel
        self.input_merge_fname = parent.input_merge_fname
        self.input_map_fname = parent.input_map_fname
        self.setObjectName('config')
        self.setAttribute(QtCore.Qt.WA_StyledBackground)

        self.processed_map = False

        self.checker = QtCore.QTimer(self)
        self.checker.timeout.connect(self.keep_checking)
        self.launcher = worker.Launcher(self, self.parent.launch_cmd)

        self.init_UI()

    def init_UI(self):
        self.setMinimumWidth(300)
        #self.setMinimumWidth(80)
        vbox = QtWidgets.QVBoxLayout()
        self.setLayout(vbox)

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

        # Tabs
        self.notebook = QtWidgets.QTabWidget()
        vbox.addWidget(self.notebook, stretch=2)

        # Bottom panel
        self.plot_controls = self.canvas_panel.init_plot_controls(self.input_merge_fname)
        vbox.addWidget(self.plot_controls)
        
        hbox = QtWidgets.QHBoxLayout()
        vbox.addLayout(hbox)
        hbox.addStretch(1)
        '''
        button = QtWidgets.QPushButton("Preprocess", self)
        button.clicked.connect(self.preprocess)
        hbox.addWidget(button)
        '''
        button = QtWidgets.QPushButton("Quit", self)
        button.clicked.connect(self.parent.close)
        hbox.addWidget(button)

        self.gen_merge_tab()
        self.gen_map_tab()
        self.gen_recon_tab()
        self.show()
        self.plot_vol(fname=self.input_merge_fname)

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
        #self.circleflag.stateChanged.connect(lambda: self.canvas_panel.replot(zoom=False))
        self.circleflag.stateChanged.connect(self.update_rings)
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
        #self.scaleradflag.stateChanged.connect(lambda: self.canvas_panel.replot(zoom=False))
        self.scaleradflag.stateChanged.connect(self.update_rings)
        hbox.addWidget(self.scaleradflag)
        hbox.addStretch(1)
        
        hbox = QtWidgets.QHBoxLayout()
        vbox.addLayout(hbox)
        button = QtWidgets.QPushButton('Zero Outer', self)
        button.clicked.connect(self.launcher.zero_outer)
        hbox.addWidget(button)
        button = QtWidgets.QPushButton('Calc. Scale', self)
        button.clicked.connect(self.launcher.calc_scale)
        hbox.addWidget(button)
        button = QtWidgets.QPushButton('Reset', self)
        button.clicked.connect(self.reset_merge_tab)
        hbox.addWidget(button)
        
        self.zero_outer_line = QtWidgets.QFrame()
        vbox.addWidget(self.zero_outer_line)
        hbox = QtWidgets.QHBoxLayout()
        self.zero_outer_line.setLayout(hbox)
        hbox.setContentsMargins(0, 0, 0, 0)
        zero_fname = os.path.splitext(self.merge_fname.text())[0] + '-zero.raw'
        label = QtWidgets.QLabel('Zero-ed volume:', self)
        hbox.addWidget(label)
        button = QtWidgets.QPushButton(zero_fname, self)
        button.clicked.connect(lambda: self.plot_vol(fname=zero_fname))
        hbox.addWidget(button)
        hbox.addStretch(1)
        self.zero_outer_line.hide()

        self.calc_scale_line = QtWidgets.QFrame()
        vbox.addWidget(self.calc_scale_line)
        hbox = QtWidgets.QHBoxLayout()
        self.calc_scale_line.setLayout(hbox)
        hbox.setContentsMargins(0, 0, 0, 0)
        self.scale_label = QtWidgets.QLabel('Scale factor = %.6e'%0.0, self)
        hbox.addWidget(self.scale_label)
        hbox.addStretch(1)
        self.calc_scale_line.hide()

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
        button.clicked.connect(self.canvas_panel.plot_map)
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
        button = QtWidgets.QPushButton('Show', self)
        button.clicked.connect(self.show_config)
        hbox.addWidget(button)
        button = QtWidgets.QPushButton('Save', self)
        button.clicked.connect(self.save_config)
        hbox.addWidget(button)

        hbox = QtWidgets.QHBoxLayout()
        vbox.addLayout(hbox)
        button = QtWidgets.QPushButton('Launch Recon', self)
        button.clicked.connect(self.launch_recon)
        hbox.addWidget(button, stretch=1)
        self.checkflag = QtWidgets.QCheckBox('Keep Checking', self)
        self.checkflag.stateChanged.connect(self.keep_checking)
        hbox.addWidget(self.checkflag)
        self.fslices = QtWidgets.QCheckBox('Fourier', self)
        hbox.addWidget(self.fslices)
        hbox.addStretch(1)
        
        hbox = QtWidgets.QHBoxLayout()
        vbox.addLayout(hbox, stretch=1)
        self.config_area = QtWidgets.QPlainTextEdit(self)
        self.highlighter = IniHighlighter(self.config_area.document())
        self.config_area.setPlainText('')
        hbox.addWidget(self.config_area)

    def reset_merge_tab(self, event=None):
        self.notebook.removeTab(0)
        #self.merge_tab.delete()
        self.gen_merge_tab(add=False)
        self.notebook.insertTab(0, self.merge_tab, 'Merge')
        self.notebook.setCurrentIndex(0)
        self.parse_vol(reset=True)

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
        self.suppressflag.stateChanged.connect(lambda: self.canvas_panel.replot(zoom='current', sigma=self.suppressflag.ischecked()))
        hbox.addWidget(self.suppressflag)
        hbox.addStretch(1)
        vbox.addLayout(hbox)

        grid = QtWidgets.QGridLayout()
        vbox.addLayout(grid)
        label = QtWidgets.QLabel('Complex:', self)
        grid.addWidget(label, 0, 0)
        button = QtWidgets.QPushButton(os.path.basename(prefix + '.cpx'), self)
        button.clicked.connect(lambda: self.plot_vol(fname=prefix + '.cpx', sigma=self.suppressflag.isChecked()))
        grid.addWidget(button, 0, 1)
        label = QtWidgets.QLabel('Symmetrized:', self)
        grid.addWidget(label, 1, 0)
        button = QtWidgets.QPushButton(os.path.basename(prefix + '-sym.raw'), self)
        button.clicked.connect(lambda: self.plot_vol(fname=prefix + '-sym.raw', sigma=self.suppressflag.isChecked()))
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
        button.clicked.connect(lambda: self.process_map(skip=True))
        hbox.addWidget(button)
        hbox.addStretch(1)

        vbox.addStretch(1)

        self.reset_button.setEnabled(True)
        self.processed_map = True

    def gen_config(self, event=None):
        if self.calc_scale_line.isHidden() or self.zero_outer_line.isHidden():
            print('Need to zero_outer and calc_scale first')
            return
        with open(self.config_fname.text(), 'w') as f:
            f.write('[parameters]\n')
            f.write('size = %d\n' % self.canvas_panel.vol_size)
            f.write('bragg_qmax = %f\n' % (float(self.radiusmin.text())/(self.canvas_panel.vol_size//2)))
            f.write('scale_factor = %f\n' % float(self.scale_label.text().split()[-1]))
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
        print('Generated %s' % self.config_fname.text())
        self.show_config()

    def show_config(self):
        with open(self.config_fname.text(), 'r') as f:
            config_text = f.read()
        self.config_area.setPlainText(config_text)
        self.config_area.show()
        self.highlighter.rehighlight()

    def save_config(self):
        with open(self.config_fname.text(), 'w') as f:
            f.write(self.config_area.toPlainText())

    def plot_vol(self, **kwargs):
        '''Wrapper around canvas_panel.plot_vol'''
        parsed = self.canvas_panel.plot_vol(**kwargs)
        if parsed:
            size = self.canvas_panel.vol_size
            self.radiusmin.setText('%d' % (size//2//2))
            self.radiusmax.setText('%d' % (size))
            self.scaleradmin.setText('%d' % (size//2//2*0.9))
            self.scaleradmax.setText('%d' % (size//2//2*1.1))

    def process_map(self, event=None, skip=False):
        '''Wrapper around launcher.process_map'''
        mapnoext = os.path.splitext(os.path.basename(self.map_fname.text()))[0]
        if skip or (os.path.isfile('data/convert/'+mapnoext+'.cpx') and QtWidgets.QMessageBox.question(self, 'Process Map', 'Found processed map output. Overwrite?', QtWidgets.QMessageBox.Yes, QtWidgets.QMessageBox.No) == QtWidgets.QMessageBox.Yes):
            if self.resedge.text() is '':
                print('Need resolution at edge of volume')
                return
            self.launcher.process_map(skip=skip)
        else:
            self.add_to_map_tab(mapnoext)
            with open('results/'+mapnoext+'.log', 'r') as f:
                words = f.read().split()
                warray = np.array(words)
                self.resedge.setText(str(float(words[words.index('./utils/read_map')+2])/(self.canvas_panel.vol_size//2)))
                self.point_group.setText(words[np.where(warray=='data/convert/'+mapnoext+'-srecon.raw')[0][0]+2])
                self.suppradstr.setText('%.1f'%float(words[np.where(warray=='./utils/create_support')[0][-1]+2]))
                self.suppthreshstr.setText('%.1f'%float(words[np.where(warray=='./utils/create_support')[0][-1]+3]))

    def launch_recon(self, event=None):
        '''Wrapper around launcher.launch_recon'''
        if os.path.isfile(self.output_prefix.text()+'-log.dat') and QtWidgets.QMessageBox.question(self, 'Overwrite Output?', 'Found output with same prefix: %s\nOverwrite?'%self.output_prefix.text(), QtWidgets.QMessageBox.Yes, QtWidgets.QMessageBox.No) == QtWidgets.QMessageBox.Yes:
            self.launcher.launch_recon()

    def keep_checking(self, event=None):
        if self.checkflag.isChecked():
            self.canvas_panel.update_slices(self.output_prefix.text(), fslices=self.fslices.isChecked())
            self.checker.start(1000)
        else:
            self.checker.stop()

    def update_rings(self, event=None):
        s = self.canvas_panel.figure.axes[0]
        size = self.canvas_panel.vol_size

        [a.remove() for a in list(set(s.findobj(patches.Circle)))]
        if self.circleflag.isChecked(): 
            rmin = float(self.radiusmin.text())
            rmax = float(self.radiusmax.text())
            s.add_artist(patches.Circle((size//2,size//2), rmin, ec='white', fc='none'))
            s.add_artist(patches.Circle((size//2,size//2), rmax, ec='white', fc='none'))

        if self.scaleradflag.isChecked(): 
            rmin = float(self.scaleradmin.text())
            rmax = float(self.scaleradmax.text())
            s.add_artist(patches.Circle((size//2,size//2), rmin, ec='white', fc='none', ls='dashed'))
            s.add_artist(patches.Circle((size//2,size//2), rmax, ec='white', fc='none', ls='dashed'))

        self.canvas_panel.canvas.draw()
