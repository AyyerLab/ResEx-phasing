from __future__ import print_function
import sys
import os
try:
    from PyQt5 import QtCore, QtWidgets, QtGui # pylint: disable=import-error
    os.environ['QT_API'] = 'pyqt5'
except ImportError:
    import sip
    sip.setapi('QString', 2)
    from PyQt4 import QtCore, QtGui # pylint: disable=import-error
    from PyQt4 import QtGui as QtWidgets # pylint: disable=import-error
    os.environ['QT_API'] = 'pyqt'
import canvas_panel
import config_panel

class ResExGUI(QtWidgets.QMainWindow):
    def __init__(self, merge_fname='', map_fname=''):
        super(ResExGUI, self).__init__()
        self.input_merge_fname = merge_fname
        self.input_map_fname = map_fname

        self.init_UI()

    def init_UI(self):
        gui_dir = os.path.dirname(os.path.realpath(__file__))
        QtGui.QFontDatabase.addApplicationFont(os.path.join(gui_dir, 'Oxygen-Regular.ttf'))
        QtGui.QFontDatabase.addApplicationFont(os.path.join(gui_dir, 'Kalam-Bold.ttf'))
        with open(os.path.join(gui_dir, 'style.css'), 'r')as f:
            self.setStyleSheet(f.read())
        #self.setWindowFlags(QtCore.Qt.FramelessWindowHint)
        self.setWindowTitle('ResEx Phasing GUI')
        #self.showMaximized()
        self.resize(1200, 700)
        window = QtWidgets.QWidget()
        hbox = QtWidgets.QHBoxLayout()
        hbox.setContentsMargins(0, 0, 0, 0)

        self.splitter = QtWidgets.QSplitter(QtCore.Qt.Horizontal)
        hbox.addWidget(self.splitter)
        self.canvas_panel = canvas_panel.CanvasPanel(self)
        self.config_panel = config_panel.ConfigPanel(self)
        self.splitter.addWidget(self.config_panel)
        self.splitter.addWidget(self.canvas_panel)
        self.splitter.setSizes([400, 696])

        window.setLayout(hbox)
        self.setCentralWidget(window)
        self.show()

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
            self.config_panel.plot_controls.hide()
            self.title_label.setText('RE')
        else:
            self.collapse_button.setText('<')
            self.notebook.show()
            self.config_panel.plot_controls.show()
            self.title_label.setText('ResEx Phasing GUI')

    def keyPressEvent(self, event): # pylint: disable=C0103
        '''Override of default keyPress event handler'''
        key = event.key()
        mod = int(event.modifiers())

        if key in (QtCore.Qt.Key_Return, QtCore.Qt.Key_Enter):
            self.canvas_panel.replot(zoom='current')
        elif QtGui.QKeySequence(mod+key) == QtGui.QKeySequence('Ctrl+Q'):
            self.close()
        elif QtGui.QKeySequence(mod+key) == QtGui.QKeySequence('Ctrl+S'):
            self.canvas_panel.save_plot()
        elif QtGui.QKeySequence(mod+key) == QtGui.QKeySequence('Ctrl+M'):
            self.canvas_panel.plot_vol(fname=self.merge_fname.text())
        elif QtGui.QKeySequence(mod+key) == QtGui.QKeySequence('Ctrl+N'):
            self.canvas_panel.plot_map()
        else:
            event.ignore()

def main():
    app = QtWidgets.QApplication([])
    if len(sys.argv) > 2:
        gui = ResExGUI(sys.argv[1], sys.argv[2])
    else:
        gui = ResExGUI()
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()
