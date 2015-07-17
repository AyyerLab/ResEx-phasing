#!/usr/bin/env python

import numpy as np
import sys
import os
from pyqtgraph.Qt import QtCore, QtGui
import pyqtgraph as pg

if len(sys.argv) < 2:
	print("Need filename")
	sys.exit()

class pgview(QtGui.QWidget):
	def __init__(self):
		super(pgview, self).__init__()
		
		self.flag = 0
		if len(sys.argv) > 2:
			self.flag = int(sys.argv[2])
		
		self.fname = sys.argv[1]
		self.size = 501
		
		self.initUI()
	
	def initUI(self):
		hbox = QtGui.QHBoxLayout()
		vbox = QtGui.QVBoxLayout()
		
		self.iw = pg.ImageView()
		self.iw.setLevels(0,2)
		
		self.file_entry = QtGui.QLineEdit(self.fname, self)
		self.file_entry.editingFinished.connect(self.reparse)
		self.file_entry.resize(self.file_entry.sizeHint())
		
		self.size_entry = QtGui.QLineEdit(str(self.size), self)
		self.size_entry.editingFinished.connect(self.reparse)
		self.size_entry.resize(self.size_entry.sizeHint())
		
		index_slider = QtGui.QSlider(QtCore.Qt.Horizontal, self)
		index_slider.setFocusPolicy(QtCore.Qt.NoFocus)
		if self.flag == 0:
			index_slider.setRange(0, int(self.size/3.))
		else:
			index_slider.setRange(0, self.size)
		index_slider.setTickPosition(QtGui.QSlider.TicksAbove)
		index_slider.setTickInterval(10)
		index_slider.setValue(int(index_slider.maximum()/2))
		index_slider.valueChanged[int].connect(self.index_change)
		index_slider.setFixedSize(300, 20)
		
		hbox.addWidget(self.size_entry)
		hbox.addStretch(1)
		hbox.addWidget(index_slider)
		
		vbox.addWidget(self.iw)
		vbox.addWidget(self.file_entry)
		vbox.addLayout(hbox)
		
		self.setLayout(vbox)
		self.resize(800,800)
		self.show()
	
		self.parse_data()
	
	def parse_data(self):
		if not os.path.isfile(self.fname):
			print "Could not find", self.fname
			return 
		
		typestr, typesize = self.parse_ext(self.fname)
		
		self.data = np.fromfile(self.fname, dtype=typestr, count = self.size**3).reshape(self.size,self.size,self.size)
		if typestr == 'complex64':
			self.data = np.square(np.absolute(self.data))
	
		if self.flag == 0:
			s = int(self.size/3)
			e = int(2*self.size/3)
			self.data = self.data[s:e, s:e, s:e]
			self.iw.setImage(self.data, autoLevels=False)
		else:
			s = self.size
			self.iw.setLevels(0, 1e3)
#			self.iw.setImage(self.data, autoLevels=True)
			self.iw.setImage(self.data, autoLevels=False)
		
		self.iw.jumpFrames(int(s/2))
	
	def parse_ext(self, file):
		ext = os.path.splitext(os.path.basename(file))[1]
		if ext == '.raw':
			typestr = 'f4'
			typesize = 4
			self.size = 501
			self.iw.setLevels(0,10)
		elif ext == '.bin':
			typestr = 'f8'
			typesize = 8
			self.size = 501
			self.iw.setLevels(0,10)
		elif ext == '.supp':
			print "Support file"
			typestr = 'uint8'
			typesize = 1
			self.size = 501
			self.iw.setLevels(0,1)
		elif ext == '.cpx':
			print "Complex file"
			typestr = 'complex64'
			typesize = 8
			self.size = 501
			self.iw.setLevels(0,1e7)
		else:
			print "Did not understand data type from extension. Defaulting to float."
			typestr = 'f4'
			typesize = 4
			self.size = 171
		
		return typestr, typesize
	
	def reparse(self):
		file = str(self.file_entry.text())
		newsize = int(self.size_entry.text())
		if file == self.fname and newsize == self.size:
			return
		else:
			self.fname = file
			print "Reparsing", file
			self.parse_data()
	
	def index_change(self, index):
		self.iw.setCurrentIndex(index)

def main():
	app = QtGui.QApplication(sys.argv)
	sh = pgview()
	sys.exit(app.exec_())

if __name__ == '__main__':
	main()
