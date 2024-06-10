import os
import platform
import sys
from PyQt6.QtCore import QDate, QRegularExpression, QTime
from PyQt6.QtGui import QRegularExpressionValidator
from PyQt6.QtWidgets import QApplication, QCheckBox, QDialog, QFileDialog, QGridLayout, QLineEdit, QMessageBox, QRadioButton,  QWidget, QPushButton, QLabel, QVBoxLayout, QComboBox
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import imageio
import webbrowser
app = QApplication([])
window = QWidget() 
window.setWindowTitle('CubicSplinePyCalc')
window.setMaximumSize(1300,550)
window.setMinimumSize(1250,500)
layout = QGridLayout() 

global a0,b0,c0,d0,figura
figura=plt.figure(figsize=(14, 4))
a0=1
b0=1
c0=1
d0=1
global a,b,c,d, flagcustomizare, x, fx
flagcustomizare=0


def natcubspline(x,flagcustomizare):
    global fx
    if flagcustomizare==1:
        fx=np.array([f(x[i]) for i in range(len(x))])

    a=fx
    alfa=np.zeros(len(x))
    h=np.diff(x) 
    for i in range(1,len(x)-1):
        alfa[i]=3*(a[i+1]-a[i])/h[i]-3*(a[i]-a[i-1])/h[i-1] 
    l=np.ones(len(x))
    u=np.zeros(len(x)) 
    z=np.zeros(len(x))
    for i in range(1,len(x)-1):
        l[i]=2*(x[i+1]-x[i-1])-h[i-1]*u[i-1]
        u[i]=h[i]/l[i]
        z[i]=(alfa[i]-h[i-1]*z[i-1])/l[i]
    c=np.zeros(len(x))   
    b=np.zeros(len(x))    
    d=np.zeros(len(x)) 
    for j in range(len(x)-2,-1,-1):
        c[j]=z[j]-u[j]*c[j+1]
        b[j]=(a[j+1]-a[j])/h[j]-h[j]*(c[j+1]+2*c[j])/3
        d[j]=(c[j+1]-c[j])/(3*h[j])
    for j in range(len(x)):
        a[j]=round(a[j],2)
        b[j]=round(b[j],2)
        c[j]=round(c[j],2)
        d[j]=round(d[j],2)
    return a,b,c,d

def f(x):
    return a0*x**3+b0*x**2+c0*x+d0

def fderiv(x):
    return 3*a0*x**2+2*b0*x+c0

def clampcubspline(x, flagcustomizare):
    global fx
    if flagcustomizare==1:
        fx=np.array([f(x[i]) for i in range(len(x))])
        FPO=fderiv(x[0])
        FPN=fderiv(x[-1])
    else:
        FPO=float(boxFPO.text())
        FPN=float(boxFPN.text())
    a=fx
    alfa=np.zeros(len(x))
    h=np.diff(x) #x[i+1]-x[i]
    alfa[0]=3*(a[1]-a[0])/h[0]-3*FPO
    alfa[-1]=3*FPN-3*(a[-1]-a[-2])/h[-1]
    
    for i in range(1,len(x)-1):
        alfa[i]=3*(a[i+1]-a[i])/h[i]-3*(a[i]-a[i-1])/h[i-1]
    l=np.ones(len(x))
    u=np.zeros(len(x)) 
    z=np.zeros(len(x))
    l[0]=2*h[0]
    u[0]=0.5
    z[0]=alfa[0]/l[0]
    for i in range(1,len(x)-1):
        l[i]=2*(x[i+1]-x[i-1])-h[i-1]*u[i-1]
        u[i]=h[i]/l[i]
        z[i]=(alfa[i]-h[i-1]*z[i-1])/l[i]
        
        
    c=np.zeros(len(x))   
    b=np.zeros(len(x))    
    d=np.zeros(len(x)) 
    l[-1]=h[-2]*(2-u[-2])
    z[-1]=(alfa[-1]-h[-2]*z[-2])/l[-1]
    c[-1]=z[-1]
    #pas7
    for j in range(len(x)-2,-1,-1):
        c[j]=z[j]-u[j]*c[j+1]
        b[j]=(a[j+1]-a[j])/h[j]-h[j]*(c[j+1]+2*c[j])/3
        d[j]=(c[j+1]-c[j])/(3*h[j])
    for j in range(len(x)):
        for j in range(len(x)):
            a[j]=round(a[j],2)
            b[j]=round(b[j],2)
            c[j]=round(c[j],2)
            d[j]=round(d[j],2)
    return a,b,c,d

def valorispline(x, x_values, a, b, c, d): 
    spline_values=np.zeros(len(x))
    for i in range(len(x)-1): 
        x_segment = x_values[(x_values>=x[i]) & (x_values<=x[i+1])] 
        spline_segment=a[i]+b[i]*(x_segment-x[i])+c[i]*(x_segment-x[i])**2+d[i]*(x_segment-x[i])**3
        spline_values[(x_values>=x[i]) & (x_values<=x[i+1])]=spline_segment
    return spline_values

j=0
label = QLabel('Pick spline type and its data')
layout.addWidget(label, 0, 0)


selectieTip=QComboBox()
selectieTip.addItems(['Natural Cubic Spline', 'Clamped Cubic Spline'])
layout.addWidget(selectieTip, 1, 0,1,2)


dateCustom=QRadioButton()
dateCustom.setText('Custom data')
layout.addWidget(dateCustom, 1, 4)

dateFisier=QRadioButton()
dateFisier.setText('Data from file')
layout.addWidget(dateFisier, 1, 5)

dialogSelectie=QFileDialog()
dialogSelectie.setFileMode(QFileDialog.FileMode.ExistingFile)
dialogSelectie.setNameFilter("NumPy Archive (*.npz)")
butdiag=QPushButton('Open file') 
layout.addWidget(butdiag, 5,5)

boxFPO=QLineEdit(placeholderText='Enter FPO')
boxFPO.setValidator(QRegularExpressionValidator(QRegularExpression("^-?\d+(\.\d+)?$"),boxFPO))
boxFPN=QLineEdit(placeholderText='Enter FPN') 
boxFPN.setValidator(QRegularExpressionValidator(QRegularExpression("^-?\d+(\.\d+)?$"),boxFPN))
boxFPO.setText("0.5")
layout.addWidget(boxFPO,3,5)
boxFPN.setText("1.5")
layout.addWidget(boxFPN,4,5)
labelFP=QLabel('Enter FPO and FPN')
layout.addWidget(labelFP,2,5)
  
label=QLabel('Input coefficients for the custom function: f(x)=a0*x^3+b0*x^2+c0*x+d0')
layout.addWidget(label,2,4)
boxA=QLineEdit(placeholderText='Enter a0')
boxA.setValidator(QRegularExpressionValidator(QRegularExpression("^-?\d+(\.\d+)?$"),boxA))
layout.addWidget(boxA,3,4)
boxB=QLineEdit(placeholderText='Enter b0')
boxB.setValidator(QRegularExpressionValidator(QRegularExpression("^-?\d+(\.\d+)?$"),boxB))
layout.addWidget(boxB,4,4)
boxC=QLineEdit(placeholderText='Enter c0')
boxC.setValidator(QRegularExpressionValidator(QRegularExpression("^-?\d+(\.\d+)?$"),boxC))
layout.addWidget(boxC,5,4)
boxD=QLineEdit(placeholderText='Enter d0')
boxD.setValidator(QRegularExpressionValidator(QRegularExpression("^-?\d+(\.\d+)?$"),boxD))
layout.addWidget(boxD ,6,4)
label2=QLabel('Input the values of x, separated by spaces')
layout.addWidget(label2,7,4)
boxX=QLineEdit(placeholderText='Enter val. x')
boxX.setText("1 1.5 2 2.25 3 4.15")
boxX.setValidator(QRegularExpressionValidator(QRegularExpression("^(-?\d+(\.\d+)?)(\s-?\d+(\.\d+)?)*$"),boxX))
layout.addWidget(boxX,8,4)

saver = QCheckBox('Save data to file')
layout.addWidget(saver,9,4)

saverdate=QCheckBox('Add current date to file name')
layout.addWidget(saverdate,9,5)   

savername=QLineEdit(placeholderText='Enter file name')
savername.setText("solution")
layout.addWidget(savername,10,4,1,2)

npzcreator=QPushButton('NPZ Creator')
layout.addWidget(npzcreator,6,5)

button = QPushButton('Draw cubic spline')
layout.addWidget(button,11,4,2,2) 

butanim=QPushButton('Export cubic spline animation')
layout.addWidget(butanim,12,4)

loopImagine=QCheckBox("Loop animation")
layout.addWidget(loopImagine,12,5)

biesire=QPushButton("Leave")
layout.addWidget(biesire,12,6,2,1)

#design NPZ
npz=QDialog()
npz.setWindowTitle('NPZ Creator')
layoutnpz=QVBoxLayout()


dateCustomNPZ=QRadioButton()
dateCustomNPZ.setText('Custom data')
layoutnpz.addWidget(dateCustomNPZ)

dateFisierNPZ=QRadioButton()
dateFisierNPZ.setText('Data from text file')
layoutnpz.addWidget(dateFisierNPZ)

labelx=QLabel('x, values separated by spaces')
layoutnpz.addWidget(labelx)
liniex=QLineEdit(placeholderText='Input val. x')
liniex.setValidator(QRegularExpressionValidator(QRegularExpression("^(-?\d+(\.\d+)?)(\s-?\d+(\.\d+)?)*$"),liniex))
layoutnpz.addWidget(liniex)

labelfx=QLabel('fx, values separated by spaces')
layoutnpz.addWidget(labelfx)
liniefx=QLineEdit(placeholderText='Input val. fx')
liniefx.setValidator(QRegularExpressionValidator(QRegularExpression("^(-?\d+(\.\d+)?)(\s-?\d+(\.\d+)?)*$"),liniefx))
layoutnpz.addWidget(liniefx)

labelfis=QLabel('Enter file name')
layoutnpz.addWidget(labelfis)

linienumefis=QLineEdit(placeholderText='File name')
layoutnpz.addWidget(linienumefis)

saverdateNPZ=QCheckBox('Add current date to file name')
layoutnpz.addWidget(saverdateNPZ)  

salvarenpz = QPushButton('Save data')
layoutnpz.addWidget(salvarenpz)   

exportnpz = QPushButton('Export data')
layoutnpz.addWidget(exportnpz)   

dialogSelectie2=QFileDialog()
dialogSelectie2.setFileMode(QFileDialog.FileMode.ExistingFile)
dialogSelectie2.setNameFilter("Text file (*.txt);;All files (*)")

 



def btn_pushed():
    try:
        global a0,b0,c0,d0, x, fx, figura
        if dateCustom.isChecked():
            flagcustomizare=1
            a0=float(boxA.text())
            b0=float(boxB.text())
            c0=float(boxC.text())
            d0=float(boxD.text())
            x=np.array([float(nr) for nr in boxX.text().split(' ')])
        else:
            flagcustomizare=0
            
        if selectieTip.currentText()=='Natural Cubic Spline':
            a,b,c,d=natcubspline(x,flagcustomizare)
        else:
            a,b,c,d=clampcubspline(x,flagcustomizare)
        if saver.isChecked():
            suprascriere=QMessageBox()
            suprascriere.setIcon(QMessageBox.Icon.Warning)
            suprascriere.setWindowTitle('Warning')
            suprascriere.setText('Save data to file? This will overwrite any existing files with the same name.')
            suprascriere.setStandardButtons(QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No)
            suprascriere.setDefaultButton(QMessageBox.StandardButton.No)
            rezultat=suprascriere.exec()
            if rezultat==QMessageBox.StandardButton.Yes:
                if saverdate.isChecked():
                    numefis=savername.text()+'_'+QDate.currentDate().toString('dd-MM-yyyy')+'_'+QTime.currentTime().toString('hh-mm-ss')
                else:
                    numefis=savername.text()
                np.savez_compressed(numefis,x=x,fx=fx)
                QMessageBox.information(window, 'Success!', 'Data saved successfully')
            
        x_values=np.linspace(min(x), max(x), len(x))
        y_values=valorispline(x, x_values, a, b, c, d) 
        figmpl.figure.clear() #cleanup
        ax=figura.add_subplot(111)
        ax.scatter(x, fx, label='Data to interpolate')
        ax.plot(x_values, y_values, label='Obtained spline', color='orange')
        ax.set_xlabel('x')
        ax.set_ylabel('f(x)')
        if selectieTip.currentText()=='Natural Cubic Spline':
            ax.set_title('Natural Cubic Spline')
        else:
            ax.set_title('Clamped Cubic Spline')
        ax.legend()
        figmpl.figure= figura
        figmpl.draw()

    except Exception as e:
        QMessageBox.warning(window, 'Invalid data!', f'Error: {e}')
    
def dinFisier():
    boxA.setEnabled(False)
    boxB.setEnabled(False)
    boxC.setEnabled(False)
    boxD.setEnabled(False)
    boxX.setEnabled(False)
    butdiag.setEnabled(True)
    saver.setEnabled(False)
    saver.setChecked(False)
    if selectieTip.currentText()=='Clamped Cubic Spline':
        boxFPO.setEnabled(True)
        boxFPN.setEnabled(True)
    else:
        boxFPO.setEnabled(False)
        boxFPN.setEnabled(False)


def dinProgram():
    boxA.setEnabled(True)
    boxB.setEnabled(True)
    boxC.setEnabled(True)
    boxD.setEnabled(True)
    boxX.setEnabled(True)
    butdiag.setEnabled(False)
    boxFPO.setEnabled(False)
    boxFPN.setEnabled(False)
    saver.setEnabled(True)
    


def selectFisier():
    if dialogSelectie.exec():
        try:
            global x, fx
            date=np.load(dialogSelectie.selectedFiles()[0])
            x=date['x']
            fx=date['fx']  
        except Exception as e:
            QMessageBox.warning(dialogSelectie, 'Warning', f'Error: {e}')
            

def ResetCandSel():
    dateCustom.setChecked(True)
    
    
#FCT NPZ    
def execNPZ():
    npz.exec()
def dinFisierNPZ():
    exportnpz.setEnabled(True)
    salvarenpz.setEnabled(False)
    liniefx.setEnabled(False)
    liniex.setEnabled(False)

def dinProgramNPZ():
    exportnpz.setEnabled(False)
    salvarenpz.setEnabled(True)
    liniefx.setEnabled(True)
    liniex.setEnabled(True)
    
    
def saveNPZ(): 
    suprascriere=QMessageBox()
    suprascriere.setIcon(QMessageBox.Icon.Warning)
    suprascriere.setWindowTitle('Warning')
    suprascriere.setText('Save data to file? This will overwrite any existing files with the same name.')
    suprascriere.setStandardButtons(QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No)
    suprascriere.setDefaultButton(QMessageBox.StandardButton.No)
    rezultat=suprascriere.exec()
    if rezultat==QMessageBox.StandardButton.Yes:
        try:
            x=np.fromstring(liniex.text(), dtype=float, sep=' ')
            fx=np.fromstring(liniefx.text(), dtype=float, sep=' ')
            if saverdateNPZ.isChecked():
                numefis=linienumefis.text()+'_'+QDate.currentDate().toString('dd-MM-yyyy')+'_'+QTime.currentTime().toString('hh-mm-ss')
            else:
                numefis=linienumefis.text()
                np.savez_compressed(numefis,x=x,fx=fx)
            QMessageBox.information(window, 'Success!', 'Data saved successfully')
        except Exception as e:
                QMessageBox.warning(npz, 'Error', f'File could not be saved... Error: {e}')
    
def exportNPZ(): 
    suprascriere=QMessageBox()
    suprascriere.setIcon(QMessageBox.Icon.Warning)
    suprascriere.setWindowTitle('Warning')
    suprascriere.setText('Save data to file? This will overwrite any existing files with the same name.')
    suprascriere.setStandardButtons(QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No)
    suprascriere.setDefaultButton(QMessageBox.StandardButton.No)
    rezultat=suprascriere.exec()
    if rezultat==QMessageBox.StandardButton.Yes:
        if dialogSelectie2.exec():
            try:
                filePath=dialogSelectie2.selectedFiles()[0]
                date=np.loadtxt(filePath, delimiter=' ', unpack=False)
                x=date[0,:]
                fx=date[1,:]
                if saverdateNPZ.isChecked():
                    numefis=linienumefis.text()+'_'+QDate.currentDate().toString('dd-MM-yyyy')+'_'+QTime.currentTime().toString('hh-mm-ss')
                else:
                    numefis=linienumefis.text()
                np.savez_compressed(numefis,x=x,fx=fx)
                QMessageBox.information(window, 'Success!', 'Data saved successfully')
            except Exception as e:
                QMessageBox.warning(npz, 'Error', f'File could not be saved... Error: {e}')
                
                              
def optiuniSalvare():
    if saver.isChecked():
        savername.setEnabled(True)
        saverdate.setEnabled(True)
    else:
        savername.setEnabled(False)
        saverdate.setEnabled(False)
        

def animatie():    
        try:
            global a0,b0,c0,d0, x, fx
            if dateCustom.isChecked():
                flagcustomizare=1
                a0=float(boxA.text())
                b0=float(boxB.text())
                c0=float(boxC.text())
                d0=float(boxD.text())
                x=np.array([float(nr) for nr in boxX.text().split(' ')])
            else:
                flagcustomizare=0
                
            if selectieTip.currentText()=='Natural Cubic Spline':
                a,b,c,d=natcubspline(x,flagcustomizare)
            else:
                a,b,c,d=clampcubspline(x,flagcustomizare)
                
                
            for i in range(1,len(x)+1): #salveaza cate o imagine pt fiecare portiune a spline-ului
                x_values=np.linspace(min(x[:i]), max(x[:i]), len(x[:i]))
                y_values=valorispline(x[:i], x_values, a, b, c, d) 
                figura=plt.figure(figsize=(14, 4))
                ax=figura.add_subplot(111)
                ax.scatter(x[:i], fx[:i], label='Data to interpolate')
                ax.plot(x_values, y_values, label='Obtained spline', color='orange')
                ax.set_xlabel('x')
                ax.set_ylabel('f(x)')
                if selectieTip.currentText()=='Natural Cubic Spline':
                    ax.set_title('Natural Cubic Spline')
                else:
                    ax.set_title('Clamped Cubic Spline')
                ax.legend()
                ax.grid()
                figura.savefig(f"plot_{i}.jpg")
                fisiere=[f"plot_{i}.jpg" for i in range(1,len(x)+1)]
                nume='anim'+'_'+QDate.currentDate().toString('dd-MM-yyyy')+'_'+QTime.currentTime().toString('hh-mm-ss')+'.gif'
            if loopImagine.isChecked(): 
                with imageio.get_writer(nume, mode='I', duration=300, loop=0) as writer:
                    for fisier in fisiere:
                        imagine=imageio.imread(fisier)
                        writer.append_data(imagine)
            else:
                with imageio.get_writer(nume, mode='I', duration=400) as writer:
                    for fisier in fisiere:
                        imagine=imageio.imread(fisier)
                        writer.append_data(imagine)            
            for fisier in fisiere: #cleanup
                os.remove(fisier)
            webbrowser.open('file://' + os.path.realpath(nume))     
                
            

        except Exception as e:
            QMessageBox.warning(dialogSelectie2, 'Eroare', f'Nu s-a putut salva fisierul... Error: {e}')


  
def iesire():
    print("Platform data for debugging: ")
    print(sys.version_info)
    print(platform.system())
    print(platform.version())
    print(platform.architecture())
    sys.exit(0)   
    
    
#setup NPZ
npz.setLayout(layoutnpz)
npz.setWindowOpacity(0.95)
npz.setStyleSheet("color: brown")
npz.setFixedSize(460,305)
linienumefis.setText("file_name")
numefisier=str(linienumefis.text())
dateCustomNPZ.setChecked(True)
dateCustomNPZ.toggled.connect(dinProgramNPZ)
dateFisierNPZ.toggled.connect(dinFisierNPZ)
salvarenpz.clicked.connect(saveNPZ)
exportnpz.clicked.connect(exportNPZ)
exportnpz.setEnabled(False)
saverdateNPZ.setChecked(True)

#fig setup
figmpl=FigureCanvas(figura)
layout.addWidget(figmpl,2,0,10,4)
navigare=NavigationToolbar(figmpl,window)
layout.addWidget(navigare,12,0,1,4)  

         
#window setup
button.clicked.connect(btn_pushed) 
dateFisier.toggled.connect(dinFisier)
dateCustom.toggled.connect(dinProgram)
saver.toggled.connect(optiuniSalvare)
selectieTip.activated.connect(ResetCandSel)    
butanim.clicked.connect(animatie)
dateCustom.setChecked(True)
butdiag.clicked.connect(selectFisier)
npzcreator.clicked.connect(execNPZ)
biesire.clicked.connect(iesire)
savername.setEnabled(False)
saverdate.setEnabled(False)
saverdate.setChecked(True)
window.setStyleSheet("color:darkblue")
window.setLayout(layout)  
window.show()
app.exec() 

