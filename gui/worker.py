from __future__ import print_function
import os
import subprocess
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

