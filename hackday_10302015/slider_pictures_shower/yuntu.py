#!/usr/bin/env python3
# vim:fileencoding=utf-8

# TODO: 并发下载
# TODO: 下载进度显示
# TODO: 允许加载已经下载但网页上没有的云图
# TODO: 网络作为可选

import os
import sys
import re
import gzip

pic_dir = '.'

try:
  from PySide import QtGui, QtCore
except ImportError:
  from PyQt4 import QtGui, QtCore

def download(pics):
  ret = []
  for p in pics:
    file = os.path.split(p)[1]
    file = os.path.join(pic_dir, file)
    ret.append(file)
  return ret

class YuntuShow(QtGui.QWidget):
  def __init__(self, pics):
    super().__init__()
    self.pics = pics
    self.initUI()

  def initUI(self):
    pic = QtGui.QPixmap(self.pics[-1])
    self.pic = piclabel = QtGui.QLabel(self)
    piclabel.setPixmap(pic)

    slider = QtGui.QSlider(QtCore.Qt.Horizontal, self)
    slider.setTickPosition(QtGui.QSlider.TicksBelow)
    m = len(self.pics) - 1
    slider.setMaximum(m)
    slider.setSliderPosition(m)
    slider.valueChanged[int].connect(self.changePic)

    vbox = QtGui.QVBoxLayout()
    vbox.addWidget(piclabel)
    vbox.addWidget(slider)
    self.setLayout(vbox)

    self.resize(960 + 10, 720 + 50)
    self.setWindowTitle('YuntuShow')
    self.show()

  def keyPressEvent(self, e):
    if e.key() == QtCore.Qt.Key_Q:
      self.close()

  def changePic(self, value):
    pic = QtGui.QPixmap(self.pics[value])
    self.pic.setPixmap(pic)

def main():
  urls = getPics(getPage())
  pics = download(urls)
  app = QtGui.QApplication(sys.argv)
  yt = YuntuShow(pics)
  sys.exit(app.exec_())

if __name__ == '__main__':
  main()
