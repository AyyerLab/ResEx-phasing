from __future__ import print_function
import os
import subprocess
from functools import partial
try:
    from PyQt5 import QtCore # pylint: disable=import-error
    os.environ['QT_API'] = 'pyqt5'
except ImportError:
    import sip
    sip.setapi('QString', 2)
    from PyQt4 import QtCore # pylint: disable=import-error
    os.environ['QT_API'] = 'pyqt'

class Launcher(object):
    def __init__(self, parent, launch_cmd):
        self.parent = parent
        self.launch_cmd = launch_cmd
        self.worker = GUIWorker()
        self.thread = QtCore.QThread()
        self.worker.moveToThread(self.thread)
        self.worker.finished.connect(self.thread.quit)
        self.make_directories()

    def cleanup_thread(self):
        if self.thread.receivers(self.thread.started) > 0:
            self.thread.started.disconnect()
        if self.worker.receivers(self.worker.returnval) > 0:
            self.worker.returnval.disconnect()

    def make_directories(self):
        os.makedirs('data', exist_ok=True)
        os.makedirs('data/convert', exist_ok=True)
        os.makedirs('data/recon', exist_ok=True)
        os.makedirs('data/logs', exist_ok=True)

    def zero_outer(self, event=None):
        fname = self.parent.merge_fname.text()
        rmin = float(self.parent.radiusmin.text())
        rmax = float(self.parent.radiusmax.text())
        self.thread.started.connect(partial(self.worker.zero_outer, fname, rmin, rmax))
        #self.worker.returnval.connect(self.parent.write_zero_line)
        self.worker.returnval.connect(self.parent.zero_outer_line.show)
        self.worker.returnval.connect(self.cleanup_thread)
        self.thread.start()

    def calc_scale(self, event=None):
        fname = self.parent.merge_fname.text()
        rmin = float(self.parent.scaleradmin.text())
        rmax = float(self.parent.scaleradmax.text())
        mapnoext = os.path.splitext(os.path.basename(self.parent.map_fname.text()))[0]
        map_fname = 'data/convert/'+mapnoext+'-sym.raw'
        self.thread.started.connect(partial(self.worker.calc_scale, map_fname, fname, rmin, rmax))
        self.worker.returnval.connect(self.parent.calc_scale_line.show)
        self.worker.returnval.connect(lambda val: self.parent.scale_label.setText('Scale factor = %e'%float(val)))
        self.worker.returnval.connect(self.cleanup_thread)
        self.thread.start()

    def process_map(self, skip=False):
        map_fname = self.parent.map_fname.text()
        size = self.parent.canvas_panel.vol_size
        resedge = float(self.parent.resedge.text())
        if self.parent.processed_map:
            supp_rad = float(self.parent.suppradstr.text())
            supp_thresh = float(self.parent.suppthreshstr.text())
        else:
            supp_rad = 3.
            supp_thresh = 1.
        point_group = self.parent.point_group.text()
        self.thread.started.connect(partial(self.worker.process_map, map_fname, size, resedge, skip, supp_rad, supp_thresh, point_group))
        self.worker.returnval.connect(self.parent.add_to_map_tab)
        self.worker.returnval.connect(lambda: self.parent.canvas_panel.replot(force=True, zoom='current'))
        self.worker.returnval.connect(self.cleanup_thread)
        self.thread.start()

    def launch_recon(self, event=None):
        self.thread.started.connect(partial(self.worker.launch_recon, self.launch_cmd, self.parent.config_fname.text()))
        self.worker.returnval.connect(self.cleanup_thread)
        self.thread.start()

    def preprocess(self, event=None):
        self.zero_outer()
        self.calc_scale()
        self.process_map()

class GUIWorker(QtCore.QObject):
    finished = QtCore.pyqtSignal()
    returnval = QtCore.pyqtSignal(str)
    rootdir = os.path.dirname(os.path.realpath(__file__))+'/../'

    @QtCore.pyqtSlot(str, str, float, float)
    def calc_scale(self, model, merge, rmin, rmax):
        cmd = os.path.realpath(os.path.join(self.rootdir, 'utils/calc_scale')) + ' %s %s %d %d'%(model, merge, rmin, rmax)
        output = subprocess.check_output(cmd.split(), shell=False)
        self.returnval.emit(output.split()[4].decode('utf-8'))
        self.finished.emit()

    @QtCore.pyqtSlot(str, float, float)
    def zero_outer(self, model, rmin, rmax):
        cmd = os.path.realpath(os.path.join(self.rootdir, 'utils/zero_outer')) + ' %s %d %d'%(model, rmin, rmax)
        print('-'*80)
        subprocess.call(cmd.split())
        print('-'*80)
        self.returnval.emit('')
        self.finished.emit()

    @QtCore.pyqtSlot(str, int, float, bool, float, float, str)
    def process_map(self, map_fname, size, resedge, full_flag, supp_rad, supp_thresh, point_group):
        flag = int(full_flag)
        cmd = os.path.realpath(os.path.join(self.rootdir, 'scripts/process_map.sh')) + ' %s %d %f %d %f %f %s'%(map_fname, size, resedge, flag, supp_rad, supp_thresh, point_group)
        print('-'*80)
        subprocess.call(cmd.split())
        print('-'*80)
        mapnoext = os.path.splitext(os.path.basename(map_fname))[0]
        self.returnval.emit(mapnoext)
        self.finished.emit()

    @QtCore.pyqtSlot(str, str)
    def launch_recon(self, launch_cmd, fname):
        cmd = '%s %s'%(launch_cmd, fname)
        print(cmd)
        print('-'*80)
        subprocess.call(cmd.split())
        print('-'*80)

