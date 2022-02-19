import sys
import os
import numpy as np
import pandas as pd
import pyqtgraph as pg
import PyQt5.QtWidgets as QtWidgets
from PyQt5.QtWidgets import QApplication, QWidget, QVBoxLayout, QHBoxLayout, QDoubleSpinBox, QPushButton, QLabel, QComboBox
from PyQt5.QtGui import QIcon
from PyQt5.QtCore import Qt
from pyqtgraph.Qt import QtGui
from scipy import signal
from obspy import read

def butterworthFilter(lowcut, highcut, fs, order):
    fn = 0.5 * fs
    b, a = signal.butter(order, [lowcut/fn, highcut/fn], btype='bandpass')
    return b, a

def bandpassFilter(data, lowcut, highcut, fs, order):
    b, a = butterworthFilter(lowcut, highcut, fs, order=order)
    return signal.lfilter(b, a, data)

def getSNR(data):
    #SNR = round(10 * np.log10(np.mean(data) / np.std(data)), 2)   # as dB
    SNR = round((np.mean(data) / np.std(data)), 2)
    return SNR

def getRMS(data):
    RMS = round((np.sqrt(np.mean( data**2 ))), 2)
    return RMS

def readRecord(path):
        _, file_extension = os.path.splitext(path)
        if file_extension == ".mseed":
            return read(f"records/{path}")[0].data
        elif file_extension == ".gcf":
            return read(f"records/{path}")[0].data
        elif file_extension == ".csv":
            return pd.read_csv(f"records/{path}", names=["Value"]).Value
        else:
            print("Unknown file format")


class Window(QWidget):
    def __init__(self):      # Constructor Method
        super().__init__() 
        self.initUI() 
    def initUI(self):
        self.setWindowTitle("Record Analyzer v0.2")
        self.setWindowIcon(QIcon("logo_low.png"))
        self.setGeometry(100,100,100,100)
        self.mainLayout = QVBoxLayout()
        self.setLayout(self.mainLayout)
        pg.setConfigOptions(antialias=True)              # Apply AntiAliasing

        # ----------------------------------------------------------------
        # Read Records 
        traces = []
        for i in os.listdir("./records"):
            traces.append(i)
        
        self.record1 = readRecord(traces[0])  
        self.record2 = readRecord(traces[1]) 

        # ----------------------------------------------------------------
        #Import Info1 Horizontal Layout
        self.rec1InfoHBoxLayout = QHBoxLayout()
        
        self.rec1Label = QLabel(f" Record 1 : {traces[0]}  ")
        self.rec1Label.setMinimumWidth(100)

        self.Sampling1SpinBox = QDoubleSpinBox()
        self.Sampling1SpinBox.setMinimumWidth(130)
        self.Sampling1SpinBox.setRange(1, 1000)
        self.Sampling1SpinBox.setSingleStep(5)
        self.Sampling1SpinBox.setPrefix("Sampling Rate: ")
        self.Sampling1SpinBox.setValue(100)

        self.Factor1SpinBox = QDoubleSpinBox()
        self.Factor1SpinBox.setMinimumWidth(130)
        self.Factor1SpinBox.setRange(1, 1_000_000_000)
        self.Factor1SpinBox.setPrefix("Inst. Factor: ")
        self.Factor1SpinBox.setValue(1)

        self.Match1SpinBox = QDoubleSpinBox()
        self.Match1SpinBox.setMinimumWidth(130)
        self.Match1SpinBox.setRange(0, 1000)
        self.Match1SpinBox.setPrefix("Shift: ")
        self.Match1SpinBox.setSuffix(" + s")
        self.Match1SpinBox.setValue(0)


        self.spacerItem = QtWidgets.QSpacerItem(0, 0, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        
        self.rec1InfoHBoxLayout.addWidget(self.rec1Label)
        self.rec1InfoHBoxLayout.addWidget(self.Sampling1SpinBox)
        self.rec1InfoHBoxLayout.addWidget(self.Factor1SpinBox)
        self.rec1InfoHBoxLayout.addWidget(self.Match1SpinBox)
        self.rec1InfoHBoxLayout.addItem(self.spacerItem)
        
        self.mainLayout.addLayout(self.rec1InfoHBoxLayout)

        # ----------------------------------------------------------------
        #Import Info2 Horizontal Layout
        self.rec2InfoHBoxLayout = QHBoxLayout()
        
        self.rec2Label = QLabel(f" Record 2 : {traces[1]}  ")
        self.rec2Label.setMinimumWidth(100)

        self.Sampling2SpinBox = QDoubleSpinBox()
        self.Sampling2SpinBox.setMinimumWidth(130)
        self.Sampling2SpinBox.setRange(1, 1000)
        self.Sampling2SpinBox.setSingleStep(5)
        self.Sampling2SpinBox.setPrefix("Sampling Rate: ")
        self.Sampling2SpinBox.setValue(100)

        self.Factor2SpinBox = QDoubleSpinBox()
        self.Factor2SpinBox.setMinimumWidth(130)
        self.Factor2SpinBox.setRange(1, 1_000_000_000)
        self.Factor2SpinBox.setPrefix("Inst. Factor: ")
        self.Factor2SpinBox.setValue(1)
        
        self.Match2SpinBox = QDoubleSpinBox()
        self.Match2SpinBox.setMinimumWidth(130)
        self.Match2SpinBox.setRange(0, 1000)
        self.Match2SpinBox.setPrefix("Shift: ")
        self.Match2SpinBox.setSuffix(" + s")
        self.Match2SpinBox.setValue(0)
        
        self.rec2InfoHBoxLayout.addWidget(self.rec2Label)
        self.rec2InfoHBoxLayout.addWidget(self.Sampling2SpinBox)
        self.rec2InfoHBoxLayout.addWidget(self.Factor2SpinBox)
        self.rec2InfoHBoxLayout.addWidget(self.Match2SpinBox)
        self.rec2InfoHBoxLayout.addItem(self.spacerItem)

        self.mainLayout.addLayout(self.rec2InfoHBoxLayout)

        
        # ----------------------------------------------------------------
        # Bandpass Filter Horizontal Layout
        self.BPFilterHBoxLayout = QHBoxLayout()

        self.bandpassLabel = QLabel("Bandpass Filter ")
        self.bandpassLabel.setMinimumWidth(105)

        self.highpassSpinBox = QDoubleSpinBox()
        self.highpassSpinBox.setMinimumWidth(160)
        self.highpassSpinBox.setRange(0.05, 200)
        self.highpassSpinBox.setSingleStep(0.1)
        self.highpassSpinBox.setPrefix("Highpass: ")
        self.highpassSpinBox.setSuffix("Hz")
        self.highpassSpinBox.setValue(0.2)
        
        self.lowpassSpinBox = QDoubleSpinBox()
        self.lowpassSpinBox.setMinimumWidth(135)
        self.lowpassSpinBox.setRange(1, 200)
        self.lowpassSpinBox.setSingleStep(1)
        self.lowpassSpinBox.setPrefix("Lowpass: ")
        self.lowpassSpinBox.setSuffix("Hz")
        self.lowpassSpinBox.setValue(10)
        
        self.orderSpinBox = QDoubleSpinBox()
        self.orderSpinBox.setMinimumWidth(100)
        self.orderSpinBox.setRange(1, 4)
        self.orderSpinBox.setSingleStep(1)
        self.orderSpinBox.setPrefix("Order: ")
        self.orderSpinBox.setValue(4)

        self.BPApplyPushButton = QPushButton("Apply Filter")
        self.BPApplyPushButton.setMinimumWidth(100)
        self.BPApplyPushButton.setEnabled(False)
        self.BPApplyPushButton.clicked.connect(self.BPApplyFunction)

        self.BPResetPushButton = QPushButton("Reset Filter")
        self.BPResetPushButton.setMinimumWidth(100)
        self.BPResetPushButton.setEnabled(False)
        self.BPResetPushButton.clicked.connect(self.PushButtonTSFunction)

        # ----------------------------------------------------------------
        self.BPFilterHBoxLayout.addWidget(self.bandpassLabel)
        self.BPFilterHBoxLayout.addWidget(self.highpassSpinBox)
        self.BPFilterHBoxLayout.addWidget(self.lowpassSpinBox)
        self.BPFilterHBoxLayout.addWidget(self.orderSpinBox)
        self.BPFilterHBoxLayout.addWidget(self.BPApplyPushButton)
        self.BPFilterHBoxLayout.addWidget(self.BPResetPushButton)
        self.BPFilterHBoxLayout.addItem(self.spacerItem)
        #self.BPFilterHBoxLayout.addWidget(self.referenceInputFunctionSelectionComboBox)

        self.mainLayout.addLayout(self.BPFilterHBoxLayout)

        # ----------------------------------------------------------------
        # Welch's Window Props Horizontal Layout
        self.WWPropHBoxLayout = QHBoxLayout()

        self.WWpropLabel = QLabel("Welch's Method ")
        self.WWpropLabel.setMinimumWidth(105)

        self.windowlengthSpinBox = QDoubleSpinBox()
        self.windowlengthSpinBox.setMinimumWidth(160)
        self.windowlengthSpinBox.setRange(5, 600)
        self.windowlengthSpinBox.setSingleStep(1)
        self.windowlengthSpinBox.setPrefix("Window Size : ")
        self.windowlengthSpinBox.setSuffix("s")
        self.windowlengthSpinBox.setValue(30)
        
        self.overlapSpinBox = QDoubleSpinBox()
        self.overlapSpinBox.setMinimumWidth(135)
        self.overlapSpinBox.setRange(2, 10)
        self.overlapSpinBox.setSingleStep(1)
        self.overlapSpinBox.setPrefix("Overlap : ")
        self.overlapSpinBox.setValue(2)
        
        self.WindowSelectionComboBox = QComboBox()
        self.WindowSelectionComboBox.setMinimumWidth(100)
        #self.referenceInputFunctionSelectionComboBox.currentTextChanged.connect(self.uptadeFunction)
        self.WindowSelectionComboBox.addItem("Hanning")
        self.WindowSelectionComboBox.addItem("Bartlett")
        self.WindowSelectionComboBox.addItem("Blackman")

        self.WPushButton = QPushButton("  Plot Welch's Method Spectra  ")
        self.WPushButton.setMinimumWidth(205)
        self.WPushButton.setEnabled(False)
        self.WPushButton.clicked.connect(self.PushButtonWWFunction)

        self.PushButtonPeaks = QPushButton("Show Peaks")
        self.PushButtonPeaks.setMinimumWidth(100)
        self.PushButtonPeaks.setEnabled(False)
        self.PushButtonPeaks.clicked.connect(self.spectrumShowPeaks)

        self.PushButtonRPeaks = QPushButton("Remove Peaks")
        self.PushButtonRPeaks.setMinimumWidth(100)
        self.PushButtonRPeaks.setEnabled(False)
        self.PushButtonRPeaks.clicked.connect(self.removeSpectrumPeaks)

        # ----------------------------------------------------------------
        self.WWPropHBoxLayout.addWidget(self.WWpropLabel)
        self.WWPropHBoxLayout.addWidget(self.windowlengthSpinBox)
        self.WWPropHBoxLayout.addWidget(self.overlapSpinBox)
        self.WWPropHBoxLayout.addWidget(self.WindowSelectionComboBox)
        self.WWPropHBoxLayout.addWidget(self.WPushButton)
        self.WWPropHBoxLayout.addWidget(self.PushButtonPeaks)
        self.WWPropHBoxLayout.addWidget(self.PushButtonRPeaks)
        self.WWPropHBoxLayout.addItem(self.spacerItem)

        self.mainLayout.addLayout(self.WWPropHBoxLayout)

        # ----------------------------------------------------------------
        # Small Dataset Props Horizontal Layout
        self.SDSPropHBoxLayout = QHBoxLayout()

        self.SDSpropLabel = QLabel("Small Dataset ")
        self.SDSpropLabel.setMinimumWidth(105)

        self.chuncksizeSpinBox = QDoubleSpinBox()
        self.chuncksizeSpinBox.setMinimumWidth(160)
        self.chuncksizeSpinBox.setRange(5, 300)
        self.chuncksizeSpinBox.setSingleStep(5)
        self.chuncksizeSpinBox.setPrefix("Chunk Size : ")
        self.chuncksizeSpinBox.setSuffix("s")
        self.chuncksizeSpinBox.setValue(30)
        
        self.SDSPushButton = QPushButton(" Plot SDS Peaks ")
        self.SDSPushButton.setMinimumWidth(135)
        self.SDSPushButton.setEnabled(False)
        self.SDSPushButton.clicked.connect(self.PushButtonSDSFunction)

        self.SDSbandPushButton = QPushButton(" Show Band  ")
        self.SDSbandPushButton.setMinimumWidth(100)
        self.SDSbandPushButton.setEnabled(False)
        self.SDSbandPushButton.clicked.connect(self.SDSgetMeanStd)

        self.SDSremovebandPushButton = QPushButton(" Remove Band  ")
        self.SDSremovebandPushButton.setMinimumWidth(100)
        self.SDSremovebandPushButton.setEnabled(False)
        self.SDSremovebandPushButton.clicked.connect(self.SDSremoveband)

        # ----------------------------------------------------------------
        self.SDSPropHBoxLayout.addWidget(self.SDSpropLabel)
        self.SDSPropHBoxLayout.addWidget(self.chuncksizeSpinBox)
        self.SDSPropHBoxLayout.addWidget(self.SDSPushButton)
        self.SDSPropHBoxLayout.addWidget(self.SDSbandPushButton)
        self.SDSPropHBoxLayout.addWidget(self.SDSremovebandPushButton)
        self.SDSPropHBoxLayout.addItem(self.spacerItem)
        
        self.mainLayout.addLayout(self.SDSPropHBoxLayout)

        # ----------------------------------------------------------------
        # Add Push Button Items
        self.plotButtonHBoxLayout = QHBoxLayout()

        self.PushButtonTS = QPushButton("Plot Time Series")
        self.PushButtonTS.setMinimumWidth(203)
        self.PushButtonTS.clicked.connect(self.PushButtonTSFunction)

        self.SpectrumSelectionComboBox = QComboBox()
        self.SpectrumSelectionComboBox.setMinimumWidth(203)
        self.SpectrumSelectionComboBox.addItem("Fourier Transform")
        self.SpectrumSelectionComboBox.addItem("Power Spectral Density")
        self.SpectrumSelectionComboBox.addItem("Spectral Density")

        self.PushButtonSP = QPushButton("Plot Spectrum")
        self.PushButtonSP.setMinimumWidth(205)
        self.PushButtonSP.setEnabled(False)
        self.PushButtonSP.clicked.connect(self.PushButtonSPFunction)

        self.plotButtonHBoxLayout.addWidget(self.PushButtonTS)
        self.plotButtonHBoxLayout.addWidget(self.SpectrumSelectionComboBox)
        self.plotButtonHBoxLayout.addWidget(self.PushButtonSP)
        self.plotButtonHBoxLayout.addItem(self.spacerItem)

        self.mainLayout.addLayout(self.plotButtonHBoxLayout)
        
        # ----------------------------------------------------------------
        # Plot 1 Time Series
        self.TSWidget = pg.PlotWidget()
        self.TSWidget.setMinimumSize(960,180) #240
        self.TSWidget.setBackground((220,220,220))

        font=QtGui.QFont()
        font.setPixelSize(13)
        self.TSWidget.getAxis("bottom").tickFont = font
        self.TSWidget.getAxis("bottom").setStyle(tickTextOffset = 8)
        self.TSWidget.setLabel('bottom', "<span style=\"color:#000000;font-size:13px\">"+"Time (second)"+"</span>")
              
        self.TSWidget.getAxis("left").tickFont = font
        self.TSWidget.getAxis("left").setStyle(tickTextOffset = 8)
        self.TSWidget.setLabel('left', "<span style=\"color:#000000;font-size:13px\">"+"Raw Sample Values"+"</span>")
        
        self.TSWidget.showGrid(x=1,y=1)
        self.TSWidget.setTitle("<span style=\"color:#000000;font-size:16px\">"+"Time Series"+"</span>")
        # ---------------------------------------------------------------------
        self.TSWidgetCurve1 = pg.PlotCurveItem(pen=pg.mkPen(color="#001AFF", width=1))
        self.TSWidget.addItem(self.TSWidgetCurve1)
        self.TSWidgetCurve2 = pg.PlotCurveItem(pen=pg.mkPen(color="#Ff1100", width=1))
        self.TSWidget.addItem(self.TSWidgetCurve2)

        # Add Selection Tool
        self.selregion = pg.LinearRegionItem(values=[0, 60])
        self.TSWidget.addItem(self.selregion) 

        # ---------------------------------------------------------------------
        # Plot 2 Spectra
        self.SPWidget = pg.PlotWidget()
        self.SPWidget.setMinimumSize(960,280)
        self.SPWidget.setBackground((220,220,220))
        
        font=QtGui.QFont()
        font.setPixelSize(12)
        self.SPWidget.getAxis("bottom").tickFont = font
        self.SPWidget.getAxis("bottom").setStyle(tickTextOffset = 8)
        self.SPWidget.setLabel('bottom', "<span style=\"color:#000000;font-size:13px\">"+"Frequency"+"</span>")
              
        self.SPWidget.getAxis("left").tickFont = font
        self.SPWidget.getAxis("left").setStyle(tickTextOffset = 8)
        self.SPWidget.setLabel('left', "<span style=\"color:#000000;font-size:13px\">"+"Amplitude"+"</span>")
        
        self.SPWidget.showGrid(x=1,y=1)
        self.SPWidget.setTitle("<span style=\"color:#000000;font-size:16px\">"+"Spectrum"+"</span>")
        # ---------------------------------------------------------------------
        self.SPWidgetCurve1 = pg.PlotCurveItem(pen=pg.mkPen(color="#001AFF", width=1.5))
        self.SPWidget.addItem(self.SPWidgetCurve1)
        self.SPWidgetCurve2 = pg.PlotCurveItem(pen=pg.mkPen(color="#Ff1100", width=1.5))
        self.SPWidget.addItem(self.SPWidgetCurve2)
        #self.SPWidget.plotItem.setLogMode(False, True)

        # ---------------------------------------------------------------------
        # Plot 3 Small Dataset 
        self.SDWidget = pg.PlotWidget()
        self.SDWidget.setMinimumSize(960,280)
        self.SDWidget.setBackground((220,220,220))
        
        font=QtGui.QFont()
        font.setPixelSize(12)
        self.SDWidget.getAxis("bottom").tickFont = font
        self.SDWidget.getAxis("bottom").setStyle(tickTextOffset = 8)
        self.SDWidget.setLabel('bottom', "<span style=\"color:#000000;font-size:13px\">"+"Time (second)"+"</span>")
              
        self.SDWidget.getAxis("left").tickFont = font
        self.SDWidget.getAxis("left").setStyle(tickTextOffset = 8)
        self.SDWidget.setLabel('left', "<span style=\"color:#000000;font-size:13px\">"+"Frequency"+"</span>")
        
        self.SDWidget.showGrid(x=1,y=1)
        self.SDWidget.setTitle("<span style=\"color:#000000;font-size:16px\">"+"Small Dataset"+"</span>")
        
        # ----------------------------------------------------------------
        self.SDWidgetCurve1 = pg.ScatterPlotItem(size=12, symbol="x", pen=pg.mkPen(None), brush=pg.mkBrush(color="#001AFF"))
        self.SDWidget.addItem(self.SDWidgetCurve1)
        self.SDWidgetCurve2 = pg.ScatterPlotItem(size=12, symbol="+", pen=pg.mkPen(None), brush=pg.mkBrush(color="#Ff1100"))
        self.SDWidget.addItem(self.SDWidgetCurve2)
        self.SDWidgetCurve3 = pg.PlotCurveItem(pen=pg.mkPen(color="#001AFF", width=1, style=Qt.DashDotLine))
        self.SDWidget.addItem(self.SDWidgetCurve3)
        self.SDWidgetCurve4 = pg.PlotCurveItem(pen=pg.mkPen(color="#Ff1100", width=1, style=Qt.DashDotLine))
        self.SDWidget.addItem(self.SDWidgetCurve3)
        # ----------------------------------------------------------------
        # Create Plotting Layout
        self.mainLayout.addWidget(self.TSWidget)
        self.mainLayout.addWidget(self.SPWidget)
        self.mainLayout.addWidget(self.SDWidget)

        # ----------------------------------------------------------------
        #Add Legend
        self.legendTS = pg.LegendItem((-20,20), offset=(-20,20))
        self.legendTS.setParentItem(self.TSWidget.graphicsItem())
        self.legendTS.addItem(self.TSWidgetCurve1, f'{traces[0]}')
        self.legendTS.addItem(self.TSWidgetCurve2, f'{traces[1]}')

        self.legendSP = pg.LegendItem((-20,20), offset=(-20,20))
        self.legendSP.setParentItem(self.SPWidget.graphicsItem())
        self.legendSP.addItem(self.SPWidgetCurve1, f'{traces[0]}')
        self.legendSP.addItem(self.SPWidgetCurve2, f'{traces[1]}')

        self.legendSD = pg.LegendItem((-20,20), offset=(-20,20))
        self.legendSD.setParentItem(self.SDWidget.graphicsItem())
        self.legendSD.addItem(self.SDWidgetCurve1, f'{traces[0]}')
        self.legendSD.addItem(self.SDWidgetCurve2, f'{traces[1]}')

        # ----------------------------------------------------------------
        # Show coordinates on Spectra
        self.labelCoor = pg.TextItem(text="Frequency: {} Hz\nAmplitude: {}".format(0, 0))
        self.SPWidget.addItem(self.labelCoor)
        
        self.setMouseTracking(True)
        self.SPWidget.scene().sigMouseMoved.connect(self.onMouseMoved)


    #----------------FUNCTIONS----------------
    def inputs(self): # Inputs
        self.npts1 = int(len(self.record1))                         # number of data points
        self.npts2 = int(len(self.record2))

        self.samplingrate1 = self.Sampling1SpinBox.value()          # sampling rate
        self.samplingrate2 = self.Sampling2SpinBox.value()

        self.fn1 = self.samplingrate1/2                             # Nyquist Frequency
        self.fn2 = self.samplingrate2/2

        self.ds = 10                                                # Downsample divide to


    def onMouseMoved(self, evt):
        if self.SPWidget.plotItem.vb.mapSceneToView(evt):
            point =self.SPWidget.plotItem.vb.mapSceneToView(evt)
            self.labelCoor.setHtml(
                "<p style='color:black'>Frequency : {0} Hz <br> Amplitude : {1}</p>".\
                format(round(point.x(),2), round(point.y(),2)))        

    def BPApplyFunction(self):
        #self.TSWidget.removeItem(self.SNRCoor)
        self.inputs()
        self.y1 =  bandpassFilter(self.y1, self.highpassSpinBox.value(), self.lowpassSpinBox.value(), self.samplingrate1, order = self.orderSpinBox.value())
        self.y2 =  bandpassFilter(self.y2, self.highpassSpinBox.value(), self.lowpassSpinBox.value(), self.samplingrate2, order = self.orderSpinBox.value())
        self.x1_ds = np.array(np.linspace(0 + self.Match1SpinBox.value(), self.npts1/self.samplingrate1 + self.Match1SpinBox.value(), int(self.npts1/self.ds)))   # Downsample
        self.x2_ds = np.array(np.linspace(0 + self.Match2SpinBox.value(), self.npts2/self.samplingrate2 + self.Match2SpinBox.value(), int(self.npts2/self.ds)))
        self.y1_ds = signal.resample(self.y1, int(self.npts1/self.ds))
        self.y2_ds = signal.resample(self.y2, int(self.npts2/self.ds))

        self.TSWidgetCurve1.setData(self.x1_ds, self.y1_ds)
        self.TSWidgetCurve2.setData(self.x2_ds, self.y2_ds)
        self.TSWidget.autoRange()

        # Show SNR and RMS on TS
        #self.SNRCoor = pg.TextItem(text="1 BP SNR, RMS: {} dB, {}\n2 BP SNR, RMS : {} dB, {}".format(getSNR(self.y1), getRMS(self.y1), 
        #                                                                                      getSNR(self.y2), getRMS(self.y2)), color="black", anchor=(1, 1))
        #self.TSWidget.addItem(self.SNRCoor)
        self.BPResetPushButton.setEnabled(True)

    def PushButtonTSFunction(self):
        # Plot TimeSeries
        #self.x = np.array(np.linspace(0, self.npts/self.samplingrate, self.npts))      # Raw
        self.inputs()
        self.y1 = self.record1 / self.Factor1SpinBox.value() #np.array(self.record1["Value"])     read(f"records/{traces[0]}.mseed")
        self.y2 = self.record2 / self.Factor2SpinBox.value()
        self.x1_ds = np.array(np.linspace(0 + self.Match1SpinBox.value(), self.npts1/self.samplingrate1 + self.Match1SpinBox.value(), int(self.npts1/self.ds)))   # Downsample
        self.x2_ds = np.array(np.linspace(0 + self.Match2SpinBox.value(), self.npts2/self.samplingrate2 + self.Match2SpinBox.value(), int(self.npts2/self.ds)))
        self.y1_ds = signal.resample(self.y1, int(self.npts1/self.ds))
        self.y2_ds = signal.resample(self.y2, int(self.npts2/self.ds))

        self.TSWidgetCurve1.setData(self.x1_ds, self.y1_ds)
        self.TSWidgetCurve2.setData(self.x2_ds, self.y2_ds)
        self.TSWidget.autoRange()

        # Show SNR and RMS on TS
        #self.SNRCoor = pg.TextItem(text="1 SNR, RMS: {} dB, {}\n2 SNR, RMS : {} dB, {}".format(getSNR(self.y1), getRMS(self.y1), 
        #                                                                                       getSNR(self.y2), getRMS(self.y2)), color="black", anchor=(1, 0))
        #self.TSWidget.addItem(self.SNRCoor)

        self.PushButtonSP.setEnabled(True)
        self.BPApplyPushButton.setEnabled(True)
        self.WPushButton.setEnabled(True)


    def PushButtonSPFunction(self):
        self.inputs()
        # Trimmed Region Properties
        self.minX, self.maxX = self.selregion.getRegion()
        self.tr_npts = (round(self.maxX)-round(self.minX)) * self.samplingrate2
        self.xmintr, self.xmaxtr = int(round(self.minX)*self.samplingrate2), int(round(self.maxX)*self.samplingrate2)
        
        #Plot Spectra
        self.spectype = self.SpectrumSelectionComboBox.currentText()
        if self.spectype == "Fourier Transform":
            self.f_s1 = np.linspace(0, self.fn1, int((self.tr_npts/2)+1))
            self.f_s2 = np.linspace(0, self.fn2, int((self.tr_npts/2)+1))
            self.a_s1 = 2*abs(np.fft.rfft(self.y1[self.xmintr:self.xmaxtr]))/self.npts1
            self.a_s2 = 2*abs(np.fft.rfft(self.y2[self.xmintr:self.xmaxtr]))/self.npts2
        elif self.spectype == "Power Spectral Density":
            self.f_s1, self.a_s1 = signal.periodogram(self.y1[self.xmintr:self.xmaxtr], self.samplingrate1, scaling="density")
            self.f_s2, self.a_s2 = signal.periodogram(self.y2[self.xmintr:self.xmaxtr], self.samplingrate2, scaling="density")
        elif self.spectype == "Spectral Density":
            self.f_s1, self.a_s1 = signal.periodogram(self.y1[self.xmintr:self.xmaxtr], self.samplingrate1, scaling="density")
            self.f_s2, self.a_s2 = signal.periodogram(self.y2[self.xmintr:self.xmaxtr], self.samplingrate2, scaling="density")
            self.a_s1 = np.sqrt(self.a_s1)
            self.a_s2 = np.sqrt(self.a_s2)
        else:
            pass

        self.SPWidgetCurve1.setData(self.f_s1[1:], self.a_s1[1:])
        self.SPWidgetCurve2.setData(self.f_s2[1:], self.a_s2[1:])
        self.SPWidget.autoRange()

        self.PushButtonPeaks.setEnabled(True)

    def reOrder(self, index, order):
        new_index = []
        for i in order:
            new_index.append(index[i])
        return new_index

    def findPeaks(self, amp, frq, num):
        ind_peaks, _ = signal.find_peaks(amp, prominence = np.mean(amp), distance =20)   #, prominence = np.mean(a_sub), distance = 1
        ind_sort = np.argsort(amp[ind_peaks])[::-1]
        ind_ord = self.reOrder(ind_peaks, ind_sort)
        f_maxspec = []
        for i in range(num):
            f_maxspec.append(frq[ind_ord[i]])
        return f_maxspec

    def spectrumShowPeaks(self):
        peak_num = 1
        r1_peaks = self.findPeaks(self.a_s1, self.f_s1, peak_num)
        r2_peaks = self.findPeaks(self.a_s2, self.f_s2, peak_num)
        for i in range(peak_num):
            self.verticalPeakSpec1 = pg.InfiniteLine( angle=90, pen=pg.mkPen(color="#001AFF", width=1, style=Qt.DashDotLine), movable=False)
            self.verticalPeakSpec1.setValue(r1_peaks[i])
            self.SPWidget.addItem(self.verticalPeakSpec1)
    
            self.verticalPeakSpec2 = pg.InfiniteLine( angle=90, pen=pg.mkPen(color="#Ff1100", width=1, style=Qt.DashDotLine), movable=False)
            self.verticalPeakSpec2.setValue(r2_peaks[i])
            self.SPWidget.addItem(self.verticalPeakSpec2)

        self.PushButtonPeaks.setEnabled(False)
        self.PushButtonRPeaks.setEnabled(True)

    def removeSpectrumPeaks(self):
        self.SPWidget.removeItem(self.verticalPeakSpec1)
        self.SPWidget.removeItem(self.verticalPeakSpec2)
        #self.SPWidget.clear()
        self.PushButtonPeaks.setEnabled(True)
        self.PushButtonRPeaks.setEnabled(False)


    def welchMethod(self, data, window_length, overlap, sampling_rate):
        window_size = int(window_length * sampling_rate)
        if self.windowtype == "Hanning":
            window = signal.windows.hann(window_size)
        elif self.windowtype == "Bartlett":
            window = signal.windows.bartlett(window_size)
        elif self.windowtype == "Blackman":
            window = signal.windows.blackman(window_size)

        w_fft, welch_m = signal.welch(data, fs = sampling_rate, window = window, nperseg = window_size, noverlap = window_size/overlap) 
        return w_fft, welch_m
        
    
    def PushButtonWWFunction(self):
        self.inputs()
        self.windowtype = self.WindowSelectionComboBox.currentText()
        # Trimmed Region Properties
        self.minX, self.maxX = self.selregion.getRegion()
        self.tr_npts = (round(self.maxX)-round(self.minX)) * self.samplingrate1
        self.xmintr, self.xmaxtr = int(round(self.minX)*self.samplingrate1), int(round(self.maxX)*self.samplingrate1)
        #Plot Spectra
        self.f_s1 = self.welchMethod(self.y1[self.xmintr:self.xmaxtr], self.windowlengthSpinBox.value(), self.overlapSpinBox.value(), self.samplingrate1)[0]
        self.f_s2 = self.welchMethod(self.y1[self.xmintr:self.xmaxtr], self.windowlengthSpinBox.value(), self.overlapSpinBox.value(), self.samplingrate2)[0]
        self.a_s1 = self.welchMethod(self.y1[self.xmintr:self.xmaxtr], self.windowlengthSpinBox.value(), self.overlapSpinBox.value(), self.samplingrate1)[1]
        self.a_s2 = self.welchMethod(self.y2[self.xmintr:self.xmaxtr], self.windowlengthSpinBox.value(), self.overlapSpinBox.value(), self.samplingrate2)[1]
        self.SPWidgetCurve1.setData(self.f_s1[1:], self.a_s1[1:])
        self.SPWidgetCurve2.setData(self.f_s2[1:], self.a_s2[1:])
        self.SPWidget.autoRange()

        self.BPResetPushButton.setEnabled(True)
        self.PushButtonPeaks.setEnabled(True)
        self.SDSPushButton.setEnabled(True)

    def frequencyRange(self, npts, samplingrate):
        if npts%2 == 1:
            frq = np.linspace(0, samplingrate/2, round((npts-1)/2)+1)     
        else:
            frq = np.linspace(0, samplingrate/2, round((npts)/2)+1) 
        return frq
    
    def reOrder(self, index, order):
        new_index = []
        for i in order:
            new_index.append(index[i])
        return new_index
    
    def SDS_function(self, data, samplingrate):
        chunk_length = self.chuncksizeSpinBox.value()
        chunks = np.array_split(data, np.floor(len(data)/(chunk_length * samplingrate)))

        chunk_range = np.linspace(chunk_length, len(chunks)*chunk_length, len(chunks))
        f_peaks = []
        for i in chunks:
            f_sub, a_sub = self.welchMethod(i, self.chuncksizeSpinBox.value(), self.overlapSpinBox.value(), samplingrate)

            ind_peaks, _ = signal.find_peaks(a_sub)                 # prominence = np.mean(a_sub), distance = 1
            ind_sort = np.argsort(a_sub[ind_peaks])[::-1]
  
            f_peaks.append(f_sub[self.reOrder(ind_peaks, ind_sort)[0]])         # get only the maxiumum peak

        return chunk_range, f_peaks

    def PushButtonSDSFunction(self):
        self.inputs()
        # Trimmed Region Properties
        self.minX, self.maxX = self.selregion.getRegion()
        self.tr_npts = (round(self.maxX)-round(self.minX)) * self.samplingrate1
        self.xmintr, self.xmaxtr = int(round(self.minX)*self.samplingrate1), int(round(self.maxX)*self.samplingrate1)
        
        #Plot Spectra
        self.max_frqx1 = self.SDS_function(self.y1[self.xmintr:self.xmaxtr], self.samplingrate1)[0]
        self.max_frqx2 = self.SDS_function(self.y2[self.xmintr:self.xmaxtr], self.samplingrate2)[0]
        self.max_frq1 = self.SDS_function(self.y1[self.xmintr:self.xmaxtr], self.samplingrate1)[1]
        self.max_frq2 = self.SDS_function(self.y2[self.xmintr:self.xmaxtr], self.samplingrate2)[1]
        self.SDWidgetCurve1.setData(self.max_frqx1, self.max_frq1)
        self.SDWidgetCurve2.setData(self.max_frqx2, self.max_frq2)

        self.SDSbandPushButton.setEnabled(True)


    def SDSgetMeanStd(self):
        self.horizontalLine1me = pg.InfiniteLine( angle=0, pen=pg.mkPen(color="#001AFF", width=1, style=Qt.DashDotLine), movable=False)
        self.SDSmean1 = np.mean(self.max_frq1)
        self.horizontalLine1me.setValue(self.SDSmean1)
        self.SDWidget.addItem(self.horizontalLine1me)

        self.horizontalLine2me = pg.InfiniteLine( angle=0, pen=pg.mkPen(color="#Ff1100", width=1, style=Qt.DashDotLine), movable=False)
        self.SDSmean2 = np.mean(self.max_frq2)
        self.horizontalLine2me.setValue(np.mean(self.SDSmean2))
        self.SDWidget.addItem(self.horizontalLine2me)


        self.SDSstdmin1 = self.SDSmean1 - np.std(self.max_frq1)
        self.SDSstdmax1 = self.SDSmean1 + np.std(self.max_frq1)
        
        self.horizontalLine3st = pg.PlotCurveItem([0, len(self.max_frq1)*self.chuncksizeSpinBox.value()], [self.SDSstdmin1, self.SDSstdmin1], pen=pg.mkPen(color="#001AFF", width=2, style=Qt.DotLine))
        self.horizontalLine4st = pg.PlotCurveItem([0, len(self.max_frq1)*self.chuncksizeSpinBox.value()], [self.SDSstdmax1, self.SDSstdmax1], pen=pg.mkPen(color="#001AFF", width=2, style=Qt.DotLine))
        self.horizontalfill1st = pg.FillBetweenItem(self.horizontalLine3st, self.horizontalLine4st, brush=pg.mkBrush(color="#001AFF"))
        self.horizontalfill1st.setOpacity(0.2)
        self.SDWidget.addItem(self.horizontalLine3st)
        self.SDWidget.addItem(self.horizontalLine4st)
        self.SDWidget.addItem(self.horizontalfill1st)

        
        self.SDSstdmin2 = self.SDSmean2 - np.std(self.max_frq2)
        self.SDSstdmax2 = self.SDSmean2 + np.std(self.max_frq2)

        self.horizontalLine5st = pg.PlotCurveItem([0, len(self.max_frq2)*self.chuncksizeSpinBox.value()], [self.SDSstdmin2, self.SDSstdmin2], pen=pg.mkPen(color="#Ff1100", width=2, style=Qt.DotLine))
        self.horizontalLine6st = pg.PlotCurveItem([0, len(self.max_frq2)*self.chuncksizeSpinBox.value()], [self.SDSstdmax2, self.SDSstdmax2], pen=pg.mkPen(color="#Ff1100", width=2, style=Qt.DotLine))
        self.horizontalfill2st = pg.FillBetweenItem(self.horizontalLine5st, self.horizontalLine6st, brush=pg.mkBrush(color="#Ff1100"))
        self.horizontalfill2st.setOpacity(0.2)
        self.SDWidget.addItem(self.horizontalLine5st)
        self.SDWidget.addItem(self.horizontalLine6st)
        self.SDWidget.addItem(self.horizontalfill2st)

        self.SDSremovebandPushButton.setEnabled(True)
        self.SDSbandPushButton.setEnabled(False)
        

    def SDSremoveband(self):
        self.SDWidget.removeItem(self.horizontalLine1me)
        self.SDWidget.removeItem(self.horizontalLine2me)
        self.SDWidget.removeItem(self.horizontalLine3st)
        self.SDWidget.removeItem(self.horizontalLine4st)
        self.SDWidget.removeItem(self.horizontalLine5st)
        self.SDWidget.removeItem(self.horizontalLine6st)
        self.SDWidget.removeItem(self.horizontalfill1st)
        self.SDWidget.removeItem(self.horizontalfill2st)

        self.SDSbandPushButton.setEnabled(True)
        self.SDSremovebandPushButton.setEnabled(False)



if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = Window()
    window.show()
    sys.exit(app.exec_())


