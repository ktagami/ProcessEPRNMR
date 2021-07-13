from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *
import sys
import time
import threading
import psutil
import os
import numpy as np
import scipy as sp 
import subprocess
from struct import *

def install(name):
    subprocess.call(['pip', 'install', name])

try:
    import nmrglue as ng
except:
    install('nmrglue')
    import nmrglue as ng

try:
    import pyqtgraph as pg
except:
    install('pyqtgrahph')
    import pyqtgraph as pg


class GUI(QWidget):

    def __init__(self, parent=None):
        super(GUI, self).__init__(parent)

        grid = QGridLayout()
        grid.addWidget(self.Load_Data(), 0, 0)
        grid.addWidget(self.Controller_Connect(), 1, 0)
        grid.addWidget(self.plot_trace(), 0, 1)
        grid.addWidget(self.plot_spectra(), 1, 1)
        grid.addWidget(self.plot_data(), 3, 1)
        grid.addWidget(self.process_data(), 3, 0)
        self.setLayout(grid)
        global indx, indy, indz ,blp,lsp,ph,lb,lo,hi,scale,EPR_NMR,ph_method
        indx = 0; indy = 0; indz = 0
        blp = 1
        lsp = 1
        ph = 0
        lb = 0.000001
        lo = int(1)
        hi = int(100)
        scale = 1
        EPR_NMR = 0
        ph_method = 'Manual'
        self.setWindowTitle("Process NMR Data")
        self.resize(900, 700)


    def Load_Data(self):
        self.GB = QGroupBox("Load and Save Data Files")
        layout0 = QGridLayout()

        #Buttons for loading and Saving Data

        self.load_button = QPushButton('Load NMR Data File',self)
        layout0.addWidget(self.load_button, 0, 0)
        self.load_button.clicked.connect(self.on_load_button_clicked)

        self.loadEPR_button = QPushButton('Load EPR Data File',self)
        layout0.addWidget(self.loadEPR_button, 1, 0)
        self.loadEPR_button.clicked.connect(self.on_loadEPR_button_clicked)

        self.save_button = QPushButton('Save Data',self)
        layout0.addWidget(self.save_button, 2, 0)
        self.save_button.clicked.connect(self.on_save_button_clicked)


        self.GB.setLayout(layout0)
        return self.GB


    def on_load_button_clicked(self):
        global data,SW_h,TD,DW,name, dims
        dlg = QFileDialog()
        dlg.setFileMode(QFileDialog.Directory)
        name = dlg.getExistingDirectory(self, 'Open File')
        try:
            dic,data = ng.bruker.read(name)
            re = np.real(data)
            im = np.imag(data)
            dim = np.shape(data)
            data = np.array(data)
            if len(dim) == 1:
                data = data.reshape(1,1,1,dim[0])
                dims = [[2,dim[0],1,1,1]]
            else:
                data = data.reshape(dim[0],1,1,dim[1])
                dims = [[2,dim[1],dim[0],1,1]]
        except:
            print('Could Not Open NMR File')

        try:
            self.plot_trace.clear()
            self.trace_r = pg.PlotDataItem(re[1], pen = pg.mkPen('g',width = 0.9))
            self.trace_i = pg.PlotDataItem(im[1], pen = pg.mkPen('r',width = 0.9))
            self.plot_trace.addItem(self.trace_r)
            self.plot_trace.addItem(self.trace_i)
        except:
            print('Couldnt Plot Trace')

        SW_h = dic['acqus']['SW_h']
        TD = dic['acqus']['TD']
        DW = dic['acqus']['D']


    def on_loadEPR_button_clicked(self):
        global data,re,im, dims
        dlg = QFileDialog()
        name = QFileDialog.getOpenFileName(self, 'Open File')
        try:
            f=open(name[0],"rb")
            data=f.read()
            f.close()
            aq_vars = unpack('I',data[:4])
            form = unpack('I',data[4:8])[0]
            if form == 0:
                dat_form = 'd'; print('data format = double')
                dat_form_len = 8
            if form == 1:
                dat_form = 'f'; print('data format = float') 
                dat_form_len = 4
            aq_vars = aq_vars[0]
            preamble = unpack('%sI' % str(2+6*aq_vars),data[:4*(2+6*aq_vars)])
            dims = preamble[2:]
            dims = [dims[i:(i+6)] for i in 6*np.arange(dims[0])]
            data_len = dims[0][5]
            traces = unpack(str(int(data_len*aq_vars))+dat_form, data[4*(2+6*aq_vars):])
            traces=np.array(traces)
            traces = traces.reshape(dims[0][0],dims[0][2],dims[0][3],dims[0][4],dims[0][1])
            data = traces[0]+traces[1]*1j  
            print('file is .d01')
        except:
            print('Couldnt open .d01 file')
        try:
            data = np.loadtxt(name[0],skiprows=1)
            data = np.transpose(data)
            dim1 = int((len(data)-1)/2)
            trace_length = len(np.unique(data[0]))
            dim2 = int(len(data[0])/trace_length)
            re = []; im = []
            for ind,i in enumerate(data[1:]):
                if ind % 2 == 0:
                    re.append(i)
                else:
                    im.append(i)
            re = np.array(re); im = np.array(im)
            data = re+im*1j
            data = data.reshape(dim1, dim2, 1, trace_length)
            dims = [[2, trace_length, dim1, dim2, 1]]
            print('file is .dat')
        except:
            print('Couldnt open .dat file')

        try:
            self.plot_trace.clear()
            self.trace_r = pg.PlotDataItem(np.real(data[0][0][0]), pen = pg.mkPen('g',width = 0.9))
            self.trace_i = pg.PlotDataItem(np.imag(data[0][0][0]), pen = pg.mkPen('r',width = 0.9))
            self.plot_trace.addItem(self.trace_r)
            self.plot_trace.addItem(self.trace_i)
        except:
            print('Couldnt Plot Trace')

    def on_save_button_clicked(self):
        global savefilename
        savefilename = QFileDialog.getSaveFileName(self,'Save File')
        if savefilename:
            np.savetxt(savefilename[0],int_r)
            np.savetxt(savefilename[0] + '_spectra',re_s)


##############################################################################
    def Controller_Connect(self):
        self.GB = QGroupBox("Process FID and Spectra")
        layout1 = QGridLayout()

        self.settrace_label = QLabel()
        self.settrace_label.setText('Select Trace (x,y,z):')
        layout1.addWidget(self.settrace_label,0,0)

        self.settracex = QDoubleSpinBox(self)
        self.settracex.setDecimals(0)
        self.settracex.setRange(0, 999)
        layout1.addWidget(self.settracex, 0, 1)
        self.settracex.valueChanged.connect(self.on_settracex_clicked)
        self.settracex.valueChanged.connect(self.update_plot)

        self.settracey = QDoubleSpinBox(self)
        self.settracey.setDecimals(0)
        self.settracey.setRange(0, 999)
        layout1.addWidget(self.settracey, 0, 2)
        self.settracey.valueChanged.connect(self.on_settracey_clicked)
        self.settracey.valueChanged.connect(self.update_plot)
        
        self.settracez = QDoubleSpinBox(self)
        self.settracez.setDecimals(0)
        self.settracez.setRange(0, 999)
        layout1.addWidget(self.settracez, 0, 3)
        self.settracez.valueChanged.connect(self.on_settracez_clicked)
        self.settracez.valueChanged.connect(self.update_plot)

        self.phase_method_label = QLabel()
        self.phase_method_label.setText('Phase Type:')
        layout1.addWidget(self.phase_method_label,1,0)

        self.phase_method = QComboBox()
        self.phase_method.addItem('Manual')
        self.phase_method.addItem('Auto')
        self.phase_method.addItem('Magnitude')
        layout1.addWidget(self.phase_method, 1, 1)
        self.phase_method.currentIndexChanged.connect(self.on_phase_method_clicked)
        self.phase_method.currentIndexChanged.connect(self.update_plot)

        self.phase_label = QLabel()
        self.phase_label.setText('Phase:')
        layout1.addWidget(self.phase_label,2,0)

        self.phase = QDoubleSpinBox(self)
        self.phase.setDecimals(1)
        self.phase.setRange(-360, 360)
        layout1.addWidget(self.phase, 2, 1)
        self.phase.valueChanged.connect(self.on_phase_clicked)
        self.phase.valueChanged.connect(self.update_plot)

        self.ls_label = QLabel()
        self.ls_label.setText('Left Shift Points:')
        layout1.addWidget(self.ls_label,3,0)

        self.ls = QDoubleSpinBox(self)
        self.ls.setDecimals(0)
        self.ls.setRange(0, 16000)
        layout1.addWidget(self.ls, 3, 1)
        self.ls.valueChanged.connect(self.on_ls_clicked)
        self.ls.valueChanged.connect(self.update_plot)

        self.em_label = QLabel()
        self.em_label.setText('Line Broadening (Hz):')
        layout1.addWidget(self.em_label,4,0)

        self.em = QDoubleSpinBox(self)
        self.em.setDecimals(6)
        self.em.setRange(0, 16000)
        layout1.addWidget(self.em, 4, 1)
        self.em.valueChanged.connect(self.on_em_clicked)
        self.em.valueChanged.connect(self.update_plot)

        self.bl_label = QLabel()
        self.bl_label.setText('Baseline Correct (%):')
        layout1.addWidget(self.bl_label,5,0)

        self.bl = QDoubleSpinBox(self)
        self.bl.setDecimals(1)
        self.bl.setRange(0,100)
        layout1.addWidget(self.bl, 5, 1)
        self.bl.valueChanged.connect(self.on_bl_clicked)
        self.bl.valueChanged.connect(self.update_plot)

        self.bl_label = QLabel()
        self.bl_label.setText('Integration Method:')
        layout1.addWidget(self.bl_label,6,0)

        self.int_method = QComboBox()
        self.int_method.addItem('FFT')
        self.int_method.addItem('Echo')
        layout1.addWidget(self.int_method, 6, 1)
        self.int_method.currentIndexChanged.connect(self.on_int_method_clicked)
        self.int_method.currentIndexChanged.connect(self.update_plot)

        self.GB.setLayout(layout1)
        return self.GB

    def on_settracex_clicked(self):
        global indx
        indx = self.settracex.value()
        indx = int(indx)
    def on_settracey_clicked(self):
        global indy
        indy = self.settracey.value()
        indy = int(indy)
    def on_settracez_clicked(self):
        global indz
        indz = self.settracez.value()
        indz = int(indz)
        
    def on_phase_method_clicked(self):
        global ph_method
        ph_method = self.phase_method.currentText()

    def on_phase_clicked(self):
        global ph
        ph = self.phase.value()

    def on_ls_clicked(self):
        global lsp
        lsp = self.ls.value()
        lsp = int(lsp)

    def on_em_clicked(self):
        global lb
        lb = self.em.value()
        lb = lb/SW_h

    def on_bl_clicked(self):
        global blp
        blp = self.bl.value()
        blp = int(blp)

    def on_int_method_clicked(self):
        global EPR_NMR
        method = self.int_method.currentText()
        if method == 'FFT':
            EPR_NMR = 0
        else:
            EPR_NMR = 1
##################################################################

##################################################################
    def plot_trace(self):
        self.GB = QGroupBox("Plot FID")
        layout3 = QGridLayout()
        self.plot_trace = pg.PlotWidget()
        self.plot_trace_item = self.plot_trace.plotItem
        self.plot_trace_item.setLabels(left='Intensity')
        self.plot_trace_item.setLabels(bottom='Point Index')
        layout3.addWidget(self.plot_trace, 1, 1)


        self.GB.setLayout(layout3)
        return self.GB

#####################################################################



######################################################################
    def plot_spectra(self):
        self.GB = QGroupBox("Plot Spectra")
        layout3 = QGridLayout()
        self.plot_spectra = pg.PlotWidget()
        self.spectra_item = self.plot_spectra.plotItem
        self.spectra_item.setLabels(left='Intensity')
        self.spectra_item.setLabels(bottom='FFT Index')
        layout3.addWidget(self.plot_spectra, 1, 1)

        self.int_region = pg.LinearRegionItem([200, 300])
        self.plot_spectra.addItem(self.int_region)
        self.int_region.sigRegionChanged.connect(self.update_int)


        self.GB.setLayout(layout3)
        return self.GB
########################################################################
    def update_int(self):
        global lo,hi
        lo,hi = self.int_region.getRegion()
        lo = int(lo); hi = int(hi)

##################################################################
    def plot_data(self):
        self.GB = QGroupBox("Plot Processed Data")
        layout3 = QGridLayout()
        self.plot_data = pg.PlotWidget()
        self.plot_data_item = self.plot_data.plotItem
        self.plot_data_item.setLabels(left='Integral (a.u)')
        self.plot_data_item.setLabels(bottom='Scan Index')
        layout3.addWidget(self.plot_data, 1, 1)


        self.GB.setLayout(layout3)
        return self.GB

#######################################################################

    def process_data(self):
        self.GB = QGroupBox("Process FID and Spectra")
        layout1 = QGridLayout()

        self.settrace_label = QLabel()
        self.settrace_label.setText('Select Trace:')
        layout1.addWidget(self.settrace_label,0,0)

        self.scale_data = QDoubleSpinBox(self)
        self.scale_data.setDecimals(3)
        self.scale_data.setRange(0.001, 100)
        layout1.addWidget(self.scale_data, 0, 1)
        self.scale_data.valueChanged.connect(self.on_scale_data_clicked)
        self.scale_data.valueChanged.connect(self.update_plot)

        self.GB.setLayout(layout1)
        return self.GB

    def on_scale_data_clicked(self):
        global scale
        scale = self.scale_data.value()
#######################################################################
    def update_plot(self):
        global int_r, re_s
        try:
            dat = data.reshape(int(np.prod(np.shape(data))/np.shape(data)[-1]),np.shape(data)[-1])
            dat = ng.proc_bl.cbf(dat,last = blp)
            dat = ng.proc_base.ls(dat,lsp)
            dat = ng.proc_base.em(dat,lb)
            if ph_method == 'Manual':
                dat = ng.proc_base.ps(dat,ph)
            if ph_method == 'Auto':
                dat = ng.proc_autophase.autops(dat,'acme')
            if ph_method == 'Magnitude':
                dat = np.abs(dat)
            else:
                pass
            print(np.shape(dat))
            print(dims)
            dat = dat.reshape(dims[0][2],dims[0][3],dims[0][4],dims[0][1])
            re = np.real(dat); im = np.imag(dat)
            try:
                self.plot_trace.removeItem(self.trace_r)
                self.plot_trace.removeItem(self.trace_i)
                self.plot_trace.removeItem(self.bl_region)
            except:
                pass
            self.bl_region = pg.LinearRegionItem([int((np.shape(re)[1]-lsp)-(np.shape(re)[1]-lsp)*blp/100), np.shape(re)[1]-lsp])
            self.plot_trace.addItem(self.bl_region)
            self.trace_r = pg.PlotDataItem(re[indx][indy][indz], pen = pg.mkPen('g',width = 0.9))
            self.trace_i = pg.PlotDataItem(im[indx][indy][indz], pen = pg.mkPen('r',width = 0.9))
            self.plot_trace.addItem(self.trace_r)
            self.plot_trace.addItem(self.trace_i)
            if EPR_NMR == 0:
                datf = ng.proc_base.fft(dat)
                re_s = np.real(datf); im_s = np.imag(datf)
            else:
                re_s = re; im_s = im
            try:
                self.plot_spectra.removeItem(self.spectra_r)
                self.plot_spectra.removeItem(self.spectra_i)
            except:
                pass
            self.spectra_i = pg.PlotDataItem(im_s[indx][indy][indz], pen = pg.mkPen('r',width = 0.9))
            self.spectra_r = pg.PlotDataItem(re_s[indx][indy][indz], pen = pg.mkPen('g',width = 0.9))
            self.plot_spectra.addItem(self.spectra_r)
            self.plot_spectra.addItem(self.spectra_i)
            ########### Plot Processed Data
            re_si = re_s.reshape(int(np.prod(np.shape(re_s))/np.shape(re_s)[-1]),np.shape(re_s)[-1])
            int_r = [np.sum(i[lo:hi]) for i in re_si]; int_r = int_r/int_r[0]
            #int_i = [np.sum(i[lo:hi]) for i in im_s]; int_i = int_i/int_i[0]
            try:
                self.plot_data.removeItem(self.data_r)
                #self.plot_data.removeItem(self.data_i)
            except:
                pass
            self.data_r = pg.PlotDataItem(int_r,symbol = 'o',symbolSize = 2, pen = pg.mkPen('g',width = 0.9))
            #self.data_i = pg.PlotDataItem(int_i*scale,symbol = 'o', pen = pg.mkPen('r',width = 0.9))
            self.plot_data.addItem(self.data_r)
            #self.plot_data.addItem(self.data_i)


        except:
            print('couldnt plot data')


if __name__ == '__main__':

    app = QApplication(sys.argv)
    G = GUI()
    G.show()
    sys.exit(app.exec_())

def kill_proc_tree(pid, including_parent=True):
    parent = psutil.Process(pid)
    if including_parent:
        parent.kill()

me = os.getpid()
kill_proc_tree(me)
