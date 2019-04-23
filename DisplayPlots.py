# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 14:46:35 2019

@author: johna
"""
from pyqtgraph.Qt import QtGui
import numpy as np
import pyqtgraph as pg

app = QtGui.QApplication([])

win = pg.GraphicsWindow(title="Basic plotting examples")
win.resize(1000,600)
win.setWindowTitle('Spacecraft Simulation')

layout = pg.LayoutWidget()
topview = layout.addLayout(row=0, col=0, rowspan=-1, colspan=2)
plot11 = layout.addLayout(row=0, col=2)
plot12 = layout.addLayout(row=0, col=3)
plot21 = layout.addLayout(row=1, col=2)
plot22 = layout.addLayout(row=1, col=3)
plot31 = layout.addLayout(row=2, col=2)
plot32 = layout.addLayout(row=2, col=3)
