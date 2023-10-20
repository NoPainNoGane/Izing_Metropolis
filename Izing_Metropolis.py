import sys
import time
import random
from numpy.random import rand
from PyQt5 import QtWidgets, QtGui, QtCore
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
#from matplotlib.figure import Figure
from PyQt5.QtWidgets import QMessageBox, QProgressBar, QDialog, QLabel

defoult = {
    'nt'      : 32,          #  number of temperature points
    'N'       : 10,          #  size of the lattice, N x N
    'eqSteps' : 2**8,        #  number of MC sweeps for equilibration
    'mcSteps' : 2**10,       #  number of MC sweeps for calculation
    'T1'      : 1.53,
    'T2'      : 3.28
   
    }

defoult2 = {
    'T'       : np.linspace(1.53, 3.28, defoult['nt']), 
    'E' : np.zeros(defoult['nt']),
    'M' : np.zeros(defoult['nt']),
    'C' : np.zeros(defoult['nt']),
    'X' : np.zeros(defoult['nt']),
    'n1' : 1.0/(defoult['mcSteps']*defoult['N']*defoult['N']),
    'n2' : 1.0/(defoult['mcSteps']*defoult['N']*defoult['N'])
    }

defoult.update(defoult2)

choice = {'Simulation', 'Configurations'}


class IsingMetropolis():
    
    def __init__(self, nt, N, eqSteps, mcSteps, T1, T2, name = 'MainCalc'):
        self.nt = nt
        self.N = N
        self.eqSteps = eqSteps
        self.mcSteps = mcSteps
        self.T1 = T1
        self.T2 = T2
        self.T = np.linspace(self.T1, self.T2, self.nt)
        self.E = np.zeros(self.nt)
        self.M = np.zeros(self.nt)
        self.C = np.zeros(self.nt)
        self.X = np.zeros(self.nt)
        self.n1 = 1.0/(self.mcSteps*self.N*self.N)
        self.n2 = 1.0/(self.mcSteps*self.mcSteps*self.N*self.N)
        #self.config = np.empty([self.N, self.N])
        #self.data_conf = []
        
        self.data_conf = []
        self.t = np.empty(6)
        self.temp = .4
        if name == 'Simulation':
            self.name = 'MainCalc'
        else:
            self.name = 'simulate'
        self.checkMethods()
        #self.MainCalc()
    
    
    #----------------------------------------------------------------------
    ##  BLOCK OF FUNCTIONS USED IN THE MAIN CODE
    #----------------------------------------------------------------------
    def checkMethods(self):
        f = getattr(self, self.name)
        f()
    
    def simulate(self):
        pass
    
    def initialstate(self):   
        ''' 
        Generates a random spin configuration for initial condition
        '''
        state = 2*np.random.randint(2, size=(self.N,self.N))-1
        #self.config = state
        return state
        
    
    
    
    def mcmove(self, beta, config):
        '''
        Monte Carlo move using Metropolis algorithm 
        '''
        
        for i in range(self.N):
            for j in range(self.N):
                    a = np.random.randint(0, self.N)
                    b = np.random.randint(0, self.N)
                    s = config[a, b]
                    nb = config[(a + 1) % self.N, b] + config[a, (b + 1) % 
                            self.N] + config[(a - 1) % self.N, b] + config[a, (b - 1) % self.N]
                    cost = 2 * s * nb
                    
                    if cost < 0:
                        s *= -1
                    elif rand() < np.exp(-cost * beta):
                        s *= -1
                    config[a, b] = s
                    
        return config
    
    
    
    def calcEnergy(self, config):
        '''
        Energy of a given configuration
        '''
        energy = 0 
        
        for i in range(len(config)):
            for j in range(len(config)):
                S = config[i,j]
                nb = config[(i + 1) % self.N, j] + config[i, (j + 1) % self.N] + config[(i - 1) %
                        self.N, j] + config[i, (j - 1) % self.N]
                energy += -nb*S
        return energy/2.  # to compensate for over-counting
    
    
    
    def calcMag(self, config):
        '''
        Magnetization of a given configuration
        '''
        mag = np.sum(config)
        return mag
    

        
        
        #----------------------------------------------------------------------
        #  MAIN PART OF THE CODE
        #----------------------------------------------------------------------
    def MainCalc(self):
    
        for tt in range(self.nt):
            config = self.initialstate()         # initialise
    
            E1 = M1 = E2 = M2 = 0
            iT=1.0/self.T[tt]
            iT2=iT*iT
            
            for i in range(self.eqSteps):         # equilibrate
                self.mcmove(iT, config)           # Monte Carlo moves
    
            for i in range(self.mcSteps):
                self.mcmove(iT, config)           
                Ene = self.calcEnergy(config)     # calculate the energy
                Mag = self.calcMag(config)        # calculate the magnetisation
    
                E1 = E1 + Ene
                M1 = M1 + Mag
                M2 = M2 + Mag*Mag 
                E2 = E2 + Ene*Ene
    
    
            # divide by number of sites and iteractions to obtain intensive values    
            self.E[tt] = self.n1 * E1
            self.M[tt] = self.n1 * M1
            self.C[tt] = (self.n1 * E2 - self.n2 * E1 * E1) * iT2
            self.X[tt] = (self.n1 * M2 - self.n2 * M1 * M1) * iT
            
            
            
        #----------------------------------------------------------------------
        #  plot the calculated values   
        #----------------------------------------------------------------------
        
    def draw(self):
        #self.MainCalc()
        f = plt.figure(figsize=(18, 10)); #  
    
    
        sp =  f.add_subplot(2, 2, 1 );
        plt.scatter(self.T, self.E, s=50, marker='o', color='IndianRed')
        plt.xlabel("Temperature (T)", fontsize=20);
        plt.ylabel("Energy ", fontsize=20);         plt.axis('tight');
    
    
        sp =  f.add_subplot(2, 2, 2 );
        plt.scatter(self.T, abs(self.M), s=50, marker='o', color='RoyalBlue')
        plt.xlabel("Temperature (T)", fontsize=20); 
        plt.ylabel("Magnetization ", fontsize=20);   plt.axis('tight');
    
    
        sp =  f.add_subplot(2, 2, 3 );
        plt.scatter(self.T, self.C, s=50, marker='o', color='IndianRed')
        plt.xlabel("Temperature (T)", fontsize=20);  
        plt.ylabel("Specific Heat ", fontsize=20);   plt.axis('tight');   
    
    
        sp =  f.add_subplot(2, 2, 4 );
        plt.scatter(self.T, self.X, s=50, marker='o', color='RoyalBlue')
        plt.xlabel("Temperature (T)", fontsize=20); 
        plt.ylabel("Susceptibility", fontsize=20);   plt.axis('tight');   
        
        
        plt.savefig("fig1.png")


        
   




class Ising():
    ''' Simulating the Ising model ''' 
    
    
    ## monte carlo moves
    def mcmove(self, config, N, beta):
        ''' This is to execute the monte carlo moves using 
        Metropolis algorithm such that detailed
        balance condition is satisified'''
        for i in range(N):
            for j in range(N):            
                    a = np.random.randint(0, N)
                    b = np.random.randint(0, N)
                    s =  config[a, b]
                    nb = config[(a+1)%N,b] + config[a,(b+1)%N] + config[(a-1)%N,b] + config[a,(b-1)%N]
                    cost = 2*s*nb
                    if cost < 0:	
                        s *= -1
                    elif rand() < np.exp(-cost*beta):
                        s *= -1
                    config[a, b] = s
        return config
    

    
    def simulate(self, f):   
        ''' This module simulates the Ising model'''
        N, temp     = 64, .4        # Initialse the lattice
        config = 2*np.random.randint(2, size=(N,N))-1
        #f = plt.figure(figsize=(15, 15), dpi=80);    
        self.configPlot(f, config, 0, N, 1);
        
        msrmnt = 1001
        for i in range(msrmnt):
            self.mcmove(config, N, 1.0/temp)
            if i == 1:       self.configPlot(f, config, i, N, 2);
            if i == 4:       self.configPlot(f, config, i, N, 3);
            if i == 32:      self.configPlot(f, config, i, N, 4);
            if i == 100:     self.configPlot(f, config, i, N, 5);
            if i == 1000:    self.configPlot(f, config, i, N, 6);
                 
                    
    def configPlot(self, f, config, i, N, n_):
        ''' This modules plts the configuration once passed to it along with time etc '''
        X, Y = np.meshgrid(range(N), range(N))
        sp =  f.add_subplot(3, 3, n_ )  
        plt.setp(sp.get_yticklabels(), visible=False)
        plt.setp(sp.get_xticklabels(), visible=False)      
        plt.pcolormesh(X, Y, config, cmap=plt.cm.RdBu);
        plt.title('Time=%d'%i); plt.axis('tight')    
    plt.show()







# =============================================================================
# class Progress(QDialog):
#     def __init__(self):
#         super().__init__()
#         self.initUI()
#         
#  
#     def initUI(self):
#         n = 500
#         layout = QtWidgets.QVBoxLayout()
#         self.progressBar = QProgressBar()
#         self.progressBar.setMinimum(0)
#         self.progressBar.setMaximum(n)
#         
#         layout.addWidget(self.progressBar)
# 
#         self.setLayout(layout)
#         self.show()
#         
#     def run(self):
#         self.show()
#         n = 500
#         for i in range(n):
#             self.show()
#             time.sleep(0.001)
#             self.progressBar.setValue(i+1)
#             self.show()
#         self.show()
# =============================================================================

class Ui_Form(object):
    def setupUi(self, Form):
        Form.setObjectName("Form")
        Form.resize(1075, 84)
        self.progressBar = QtWidgets.QProgressBar(Form)
        self.progressBar.setGeometry(QtCore.QRect(30, 30, 1000, 35))
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.progressBar.sizePolicy().hasHeightForWidth())
        self.progressBar.setSizePolicy(sizePolicy)
        self.progressBar.setMinimumSize(QtCore.QSize(1000, 35))
        self.progressBar.setMaximumSize(QtCore.QSize(1000, 35))
        self.progressBar.setProperty("value", 0)
        self.progressBar.setObjectName("progressBar")
    
        self.retranslateUi(Form)
        QtCore.QMetaObject.connectSlotsByName(Form)

    def retranslateUi(self, Form):
        _translate = QtCore.QCoreApplication.translate
        Form.setWindowTitle(_translate("Form", "Progress bar"))
    
class ProgressBar(QtWidgets.QDialog, Ui_Form):
    def __init__(self, desc = None, parent=None):
        super(ProgressBar, self).__init__(parent)
        self.setupUi(self)
        self.show()

        if desc != None:
            self.setDescription(desc)

    def setValue(self, val): # Sets value
        self.progressBar.setProperty("value", val)

    def setDescription(self, desc): # Sets Pbar window title
        self.setWindowTitle(desc)
    
    


class HelpWindow(QDialog):
    
     def __init__(self):
         
         super().__init__()
         self.initUI()
 
     def initUI(self):
        
         lbl = QLabel('''The Ising model plays a central role in the theory of phase transitions. 
 It is a mathematical model of ferromagnetism (for example, iron can be magnetized 
in a magnetic field, but if heated, it loses magnetization beyond Curie temperature).
This blog contains Python code and a detailed algorithm for the Monte Carlo simulation
of the Ising model. A corresponding numerical simulation of a continuum description
of the Ising model is also provided below.
The Ising model is named after Ernst Ising, Ph.D. in Physics (1924) from the University 
of Hamburg under the supervision of Wilhelm Lenz. Ising solved the one-dimensional (1D) 
Ising model exactly to find no phase transition. He also provided arguments on why there
would not be a phase transition in higher dimensions either. In 1936, Peierls argued that 
both 2D and 3D Ising models admit phase transitions. ''', self)
         lbl.move(35, 40)
        
         qbtn = QtWidgets.QPushButton('OK', self)
         qbtn.resize(qbtn.sizeHint())
         qbtn.move(520, 170)
         qbtn.clicked.connect(self.close)
        
         self.setWindowTitle ('Help')
         self.setFixedSize(600, 200)





class MyWindow(QtWidgets.QMainWindow):
    def __init__(self, parent=None):
        super(MyWindow, self).__init__()
        self.central = None
        self.initUI()

    def initUI(self):
        self.setWindowTitle("Ising Model")
        self.resize(1920, 1010)

        self.create_menu()
        self.create_status_bar()
        
        self.central = CentralWidget(self)
        self.setCentralWidget(self.central)
        self.create_heading()
        # self.status_text.setText(self.central.time_text)
        # self.update()
        
    


    def create_menu(self):
        self.file_menu = self.menuBar().addMenu("&File")
        save_plot_file = self.create_action("&Save plot", shortcut='Ctrl+S', slot=self.save_plot,
                                            tip='Save the plot')

        exit_file = self.create_action("&Quit", shortcut='Ctrl+x', slot=self.close, tip='Exit')
        self.add_action(self.file_menu, (save_plot_file, None))

        self.help_menu = self.menuBar().addMenu("&Help")
       
        about = self.create_action("&About programm", shortcut='Ctrl+6', slot=self.helpWindow,
                                            tip='About')
        self.add_action(self.help_menu, (about,))

    def create_action(self, text, slot=None, shortcut=None,
                      icon=None, tip=None, checkable=False):

        newAction = QtWidgets.QAction(text, self)

        if shortcut is not None:
            newAction.setShortcut(shortcut)
        if tip is not None:
            newAction.setToolTip(tip)
            newAction.setStatusTip(tip)
        if slot is not None:
            newAction.triggered.connect(slot)
        if checkable:
            newAction.setCheckable(True)
        return newAction

    @staticmethod
    def add_action(target, actions):
        for action in actions:
            if action is None:
                target.addSeparator()
            else:
                target.addAction(action)

    def create_status_bar(self):
        self.coord_text = QtWidgets.QLabel("X: % 3f, Y: % 3f" % (0, 0))
        self.statusBar().addWidget(self.coord_text, -1)
        
        self.status_text = QtWidgets.QLabel("Waiting")
        self.statusBar().addWidget(self.status_text, 1)

    def save_plot(self):
        file_choices = "PNG (*.png)|*.png"
        path = QtWidgets.QFileDialog.getSaveFileName(self,
                                           'Save file', '',
                                           file_choices)
        if path:
            print(path[0])
            self.central.canvas.print_figure(path[0])
            self.statusBar().showMessage('Saved to %s' % str(path), 1)

        
    def closeEvent(self, event):
     
        reply = QtWidgets.QMessageBox.question(self, 'Уведомление',
            "Вы уверены, что хотите закрыть окно?", QtWidgets.QMessageBox.Yes |
            QtWidgets.QMessageBox.No, QtWidgets.QMessageBox.No)
    
        if reply == QtWidgets.QMessageBox.Yes:
            
            event.accept()
            
        else:
            
            event.ignore()
            
    def create_heading(self):
        #  Заголовок ----------------------------------------------------------
        self.label = QtWidgets.QLabel(self.central)
        self.label.setGeometry(QtCore.QRect(0, 0, 801, 71))
        font = QtGui.QFont()
        font.setPointSize(20)
        font.setBold(True)
        font.setItalic(False)
        font.setUnderline(False)
        font.setWeight(75)
        font.setStrikeOut(False)
        font.setKerning(True)
        self.label.setFont(font)
        self.label.setLayoutDirection(QtCore.Qt.LeftToRight)
        self.label.setAutoFillBackground(False)
        self.label.setStyleSheet("background-color: rgb(106, 106, 106);\n"
"color: rgb(255, 255, 255);")
        self.label.setScaledContents(True)
        self.label.setWordWrap(False)
        self.label.setObjectName("label")
        _translate = QtCore.QCoreApplication.translate
        self.label.setText(_translate("MainWindow", "Модель Изинга и алгоритм Метрополиса"))
        #  --------------------------------------------------------------------
        
    def helpWindow(self):
        
        self.chile_Win = HelpWindow()
        self.chile_Win.show()
        self.chile_Win.exec_()
        


class CentralWidget(QtWidgets.QWidget):

    def __init__(self, parent=None):
        QtWidgets.QWidget.__init__(self, parent)
        self.parent = parent

        self.initUI()
        self.reset()
        
        

    def initUI(self):
        self.fig = plt.figure(figsize=(18, 10), dpi=80)
        self.canvas = FigureCanvasQTAgg(self.fig)

        self.canvas.mpl_connect("motion_notify_event", self.statusbar_coord)


        self.textbox_num_temp_points = QtWidgets.QLineEdit()
        self.textbox_num_temp_points_text = QtWidgets.QLabel("Number of tempreture points: ")
        self.textbox_num_temp_points.setMaximumWidth(100)
        self.textbox_num_temp_points_text.setAlignment(QtCore.Qt.AlignRight)
        

        self.textbox_size_lattice = QtWidgets.QLineEdit()
        self.textbox_size_lattice_text = QtWidgets.QLabel("Size of the lattice: ")
        self.textbox_size_lattice.setMaximumWidth(100)
        self.textbox_size_lattice_text.setAlignment(QtCore.Qt.AlignRight)


        self.textbox_eqSteps = QtWidgets.QLineEdit()
        self.textbox_eqSteps_text = QtWidgets.QLabel("Number of MC sweeps for equilibration (2^8): ")
        self.textbox_eqSteps.setMaximumWidth(100)
        self.textbox_eqSteps_text.setAlignment(QtCore.Qt.AlignRight)

        self.textbox_mcSteps = QtWidgets.QLineEdit()
        self.textbox_mcSteps_text = QtWidgets.QLabel("Number of MC sweeps for calculation (2^10): ")
        self.textbox_mcSteps.setMaximumWidth(100)
        self.textbox_mcSteps_text.setAlignment(QtCore.Qt.AlignRight)


        self.combo_box_methods = QtWidgets.QComboBox()
        self.combo_box_methods.setMaximumWidth(100)
        self.combo_box_methods.addItems(choice)
        self.combo_box_methods_text = QtWidgets.QLabel("Simulation/Configurations")
        self.combo_box_methods_text.setAlignment(QtCore.Qt.AlignRight)


        self.spin_T1_text = QtWidgets.QLabel("T1: ")
        self.spin_T1 = QtWidgets.QDoubleSpinBox()
        self.spin_T1_text.setAlignment(QtCore.Qt.AlignRight)

        self.spin_T2_text = QtWidgets.QLabel("T2: ")
        self.spin_T2 = QtWidgets.QDoubleSpinBox()
        self.spin_T2_text.setAlignment(QtCore.Qt.AlignRight)


        self.reset_button = QtWidgets.QPushButton("Reset")
        self.reset_button.clicked.connect(self.reset)
        self.reset_button.setCursor(QtGui.QCursor(QtCore.Qt.PointingHandCursor))
        self.reset_button.setStyleSheet("background-color: rgb(51, 255, 0);")

        self.draw_button = QtWidgets.QPushButton("Calculate && Draw")
        self.draw_button.clicked.connect(self.on_draw)
        self.draw_button.clicked.connect(self.progress_click)
        self.draw_button.setCursor(QtGui.QCursor(QtCore.Qt.PointingHandCursor))
        self.draw_button.setStyleSheet("background-color: rgb(51, 255, 0);")
        
        self.qbtn = QtWidgets.QPushButton('Clear', self)
        self.qbtn.resize(self.qbtn.sizeHint())
        self.qbtn.setToolTip('<b>Alt+1<b>')
        self.qbtn.clicked.connect(self.close)
        self.qbtn.setStatusTip('Exiting the application')
        self.qbtn.setStyleSheet("background-color: rgb(255, 163, 71);")
        
        grid = QtWidgets.QGridLayout()
        grid.setColumnStretch(1, 1)
        #grid.setRowStretch(1, 1)
        
        hbox = QtWidgets.QVBoxLayout()
        hbox.addWidget(self.canvas)
        hbox.addLayout(grid, 1)
        
        

        grid.addWidget(self.textbox_num_temp_points, 0, 1)
        grid.addWidget(self.textbox_num_temp_points_text, 0, 0)

        grid.addWidget(self.textbox_size_lattice, 1, 1)
        grid.addWidget(self.textbox_size_lattice_text, 1, 0)

        grid.addWidget(self.textbox_eqSteps, 2, 1)
        grid.addWidget(self.textbox_eqSteps_text, 2, 0)

        grid.addWidget(self.textbox_mcSteps, 3, 1)
        grid.addWidget(self.textbox_mcSteps_text, 3, 0)

        grid.addWidget(self.combo_box_methods, 4, 1)
        grid.addWidget(self.combo_box_methods_text, 4, 0)

        grid.addWidget(self.spin_T1_text, 5, 0)
        grid.addWidget(self.spin_T1, 5, 1)

        grid.addWidget(self.spin_T2_text, 6, 0)
        grid.addWidget(self.spin_T2, 6, 1)
        

        grid.addWidget(self.reset_button, 7, 1, 1, 1)
        grid.addWidget(self.draw_button, 8, 1, 1, 1)

        grid.addWidget(QtWidgets.QLabel(""), 9, 1, 2, 2)
        
        grid.addWidget(self.qbtn, 10, 1, 2, 1)

        self.setLayout(hbox)
        
    
    def closeEvent(self, event):
     
        reply = QtWidgets.QMessageBox.question(self, 'Уведомление',
            "Вы уверены, что хотите закрыть окно?", QtWidgets.QMessageBox.Yes |
            QtWidgets.QMessageBox.No, QtWidgets.QMessageBox.No)
    
        if reply == QtWidgets.QMessageBox.Yes:
            
            event.accept()
            
        else:
            
            event.ignore()
            
# =============================================================================
#     def progressWindow(self):
#         
#         self.progressBar = Progress()
#         #self.progressBar.show()
#         #self.progressBar.initUI()
#         self.progressBar.run()
#         #self.progressBar.show()
#         self.progressBar.exec_()
# =============================================================================
    def progress_click(self):
        pb = ProgressBar()

        for i in range(0, 100):
            time.sleep(0.01)
            pb.setValue(((i + 1) / 100) * 100)
            QtWidgets.QApplication.processEvents()    

    def reset(self):
        self.textbox_num_temp_points.setText(str(defoult['nt']))
        self.textbox_size_lattice.setText(str(defoult['N']))
        self.textbox_eqSteps.setText(str(defoult['eqSteps']))
        self.textbox_mcSteps.setText(str(defoult['mcSteps']))

        self.spin_T1.setRange(-20, 20)
        self.spin_T2.setRange(-20, 20)
        self.spin_T1.setValue(defoult['T1'])
        self.spin_T1.setSingleStep(0.01)
        self.spin_T2.setValue(defoult['T2'])
        self.spin_T2.setSingleStep(0.01)
        self.on_draw()

    def on_draw(self):
        self.fig.clf()
        nt = int(self.textbox_num_temp_points.text())
        N = int(self.textbox_size_lattice.text())
        eqSteps = int(self.textbox_eqSteps.text())
        mcSteps = int(self.textbox_mcSteps.text())
        T1 = float(self.spin_T1.text().replace(',', '.'))
        T2 = float(self.spin_T2.text().replace(',', '.'))
        
        self.IM = IsingMetropolis(nt, N, eqSteps, mcSteps, T1, T2, self.combo_box_methods.currentText())
    
        if self.combo_box_methods.currentText() == 'Simulation':
            self.fig.add_subplot(2, 2, 1 );
            plt.scatter(self.IM.T, self.IM.E, s=50, marker='o', color='IndianRed')
            plt.xlabel("Temperature (T)", fontsize=20);
            plt.ylabel("Energy ", fontsize=20);         plt.axis('tight');
        
        
            self.fig.add_subplot(2, 2, 2 );
            plt.scatter(self.IM.T, abs(self.IM.M), s=50, marker='o', color='RoyalBlue')
            plt.xlabel("Temperature (T)", fontsize=20); 
            plt.ylabel("Magnetization ", fontsize=20);   plt.axis('tight');
        
        
            self.fig.add_subplot(2, 2, 3 );
            plt.scatter(self.IM.T, self.IM.C, s=50, marker='o', color='IndianRed')
            plt.xlabel("Temperature (T)", fontsize=20);  
            plt.ylabel("Specific Heat ", fontsize=20);   plt.axis('tight');   
        
        
            self.fig.add_subplot(2, 2, 4 );
            plt.scatter(self.IM.T, self.IM.X, s=50, marker='o', color='RoyalBlue')
            plt.xlabel("Temperature (T)", fontsize=20); 
            plt.ylabel("Susceptibility", fontsize=20);   plt.axis('tight');
            
            self.canvas.draw_idle()
        else:
            IS = Ising()
            IS.simulate(self.fig)
            self.canvas.draw_idle()

            
            #self.canvas.draw_idle()
            
      
        
        
        #self.canvas.draw_idle()

    def statusbar_coord(self, event):
        # Show coordinates time in statusbar
        if event.inaxes is not None:
            text = "X: % .3f, Y % .3f" % (event.xdata, event.ydata)
            self.parent.coord_text.setText(text)
    
        
    

if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    window = MyWindow()
    window.show()
    sys.exit(app.exec_())




