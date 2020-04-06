
#import cora.core
#import cora.core.calibration
import glob
import json
import numpy as np
import os
import re
import time
import matplotlib.path as mpath
import datetime
import sys
from CoreControllGui import Ui_MainWindow
from assembling_orca import Orca#, MeasurementMode
from cora.core.interpolation import get_evenspaced_spectrum, match_x_axes
from cora.core.calibration import relative_intensity_correction, y_axis
from PyQt5.QtCore import QCoreApplication
from PyQt5 import QtCore, uic, QtWidgets
from PyQt5.QtWidgets import QMainWindow, QApplication,QSlider, QWidget, QPushButton, QAction, QLineEdit, QMessageBox, QMenu, QVBoxLayout, QSizePolicy
from PyQt5.QtGui import QIcon
from USBswitch import set_channel
from scipy.optimize import curve_fit
from datetime import datetime
    ######################################################
def numericalSort(value):
        numbers = re.compile(r'(\d+)')
        parts = numbers.split(value)
        parts[1::2] = map(int, parts[1::2])
        return parts 
    ######################################################
def status(form,bar):  
        barname = "progressBar_"+str(bar)   
        zaehler = 0
        while zaehler < 100:
            zaehler += 0.0001
            getattr(form, barname).setValue(zaehler)                  
        zaehler = 100
        while zaehler > 0:
            zaehler -= 0.0001
            getattr(form, barname).setValue(zaehler)         
    ######################################################
def statushalf(form,bar):  
        barname = "progressBar_"+str(bar)   
        zaehler = 0
        while zaehler < 100:
            zaehler += 0.0001
            getattr(form, barname).setValue(zaehler)                        
def exp_fit(x, offset, x0, xe):
        x = np.asarray(x)
        return offset + np.exp((x - x0) / xe)
    ######################################################
def line_fit(x, y0, slope):
        x = np.asarray(x)
        return x * slope + y0
    ######################################################
def save(form,name,canvas,bar):
        #Name=zB. "Spektrallamp_{}.jpg"
        #Canvas=self.canvasNr
        counter = 1
        filename = name  
        while os.path.isfile(filename.format(counter)):
            counter += 1
        filename = filename.format(counter)
        canvas.print_figure(filename, dpi=150)                
        status(form,bar)    
    ######################################################
def Gauss(npx, a, x0, sigma):  
    return a*np.exp(-(npx-x0)**2/(2*sigma**2))
    ######################################################
def lorentzian(npx, a, x0, sigma):  
    return a*sigma**2/(sigma**2+(npx-x0)**2)   
    ######################################################
def create785ax(form,y):   
    #print(y)      
    peak_height_threshold = 1000
    old_laser_wl = 784.5839375
    old_x_coeff = [748.261575786785, 0.142727638865543, -6.09516632365329E-06, -8.92332882686976E-10]
    loadref1(form)
    nominal_peak_data = form.refdata
    calc_laser_wl=False
    method="auto"
    coreresults = cora.core.calibration.calibrate_x(y, peak_height_threshold, calc_laser_wl, nominal_peak_data, method, old_laser_wl, old_x_coeff)
    xcalcoeff = coreresults['coefficients']     # coefficients from coracore
    lwl = coreresults['new_laser_wl']                         ##############
    z = np.arange(len(y)).astype(float)    
    achse785 = xcalcoeff[0]+xcalcoeff[1]*z+xcalcoeff[2]*z**2+xcalcoeff[3]*z**3
    form.wellenzahl = 10**7/lwl-10**7/achse785
    print(form.wellenzahl)
    ######################################################
def create1064ax(form,y):   
    #print(y)      
    peak_height_threshold = 400
    old_laser_wl = 1063.442875
    old_x_coeff = [1078.69466551906, 1.42956620256726, -0.000403740288097827, 1.06210060979639E-06]
    loadref1(form)
    nominal_peak_data = form.refdata
    calc_laser_wl=False
    method="auto"
    coreresults = cora.core.calibration.calibrate_x(y, peak_height_threshold, calc_laser_wl, nominal_peak_data, method, old_laser_wl, old_x_coeff)
    xcalcoeff = coreresults['coefficients']     # coefficients from coracore
    lwl = coreresults['new_laser_wl']                         ##############
    z = np.arange(len(y)).astype(float)    
    achse1064 = xcalcoeff[0]+xcalcoeff[1]*z+xcalcoeff[2]*z**2+xcalcoeff[3]*z**3
    form.wellenzahl = 10**7/lwl-10**7/achse1064
    print(form.wellenzahl)
    ######################################################
def loadref1(form):                        
    status(form,15)    
    refpath = "C:\\Users\\martin.hell\\Cora5x\\csv\\reference\\"       
    os.chdir(refpath)
    fileList=glob.glob("*.json")                                                                                            
    for filename in sorted(fileList, key=numericalSort):     
        with open(filename) as reffile:
            refdatar = json.load(reffile)
           # form.refdata = next(subst for subst in refdatar if subst['name'] == 'cyclohexane')  
            if form.checkBox_7.isChecked(): 
                form.refdata = next(subst for subst in refdatar if subst['name'] == 'cyclohexane')
            if form.checkBox_8.isChecked(): 
                form.refdata = next(subst for subst in refdatar if subst['name'] == 'benzonitrile')
            #print(form.refdata)
    os.chdir(form.loadallpath)  
    ######################################################
def sliderset(form, hs1, hs2, vs1,vs2,hl,vl,fkt):
    #CoraClassHell.sliderset(self,14,15,13,14,2000,2100,self.plotcydarkstd)
    hslider1 = "horizontalSlider_"+str(hs1)
    hslider2 = "horizontalSlider_"+str(hs2)
    vslider1 = "verticalSlider_"+str(vs1)
    vslider2 = "verticalSlider_"+str(vs2)
    QApplication.processEvents()   
    getattr(form,hslider1).setMaximum(hl)
    getattr(form,hslider2).setMaximum(hl)
    getattr(form,vslider1).setMaximum(vl)
    getattr(form,vslider2).setMaximum(vl)
    getattr(form,hslider2).setValue(hl)
    getattr(form,vslider2).setValue(vl)
    getattr(form,hslider1).valueChanged.connect(fkt) 
    getattr(form,hslider2).valueChanged.connect(fkt)
    getattr(form,vslider1).valueChanged.connect(fkt)
    getattr(form,vslider2).valueChanged.connect(fkt)
def slidercall(form, hs1, hs2, vs1,vs2,axis):
    #CoraClassHell.slidercall(self,14,15,13,14,self.ax17) 
    hslider1 = "horizontalSlider_"+str(hs1)
    hslider2 = "horizontalSlider_"+str(hs2)
    vslider1 = "verticalSlider_"+str(vs1)
    vslider2 = "verticalSlider_"+str(vs2)
    form.xap = getattr(form,hslider1).value()
    form.xep = getattr(form,hslider2).value()
    form.yap = getattr(form,vslider1).value()
    form.yep = getattr(form,vslider2).value()
    axis.set_xlim(form.xap,form.xep)
    axis.set_ylim(form.yap,form.yep)           
def spa_files(form, y):   
    form.path = "C:\\Users\\martin.hell\\Cora5x\\csv"
    os.chdir(form.path) 
    form.fileList=glob.glob("*.spa") 
    if form.fileList:                
        y=[]                   
        for filename in sorted(form.fileList, key=numericalSort):
            dat = spafile.load(filename)                                                                    #  import spa file
            raw_counts = np.array( [ capture['Spectrum']['Raw_Counts'] for capture in dat['Captures']])     #######                
            intime = np.array( [ capture['Spectrum']['Set_Integration_Time'] for capture in dat['Captures']])
            i=0                 
            form.spectra = []                
            while i<len(raw_counts):                                   ##############
                form.spectra.append(raw_counts[i]-raw_counts[i+1])     # spectra minus dark      
                i+=2                  
            deleter = []                                   ################
            for i in range (0,len(intime)):                #
                if (i+1)%2 == 0:                           #
                    deleter.append(i)                      #  delete dark entries from intime
            form.intimes = np.delete(intime, deleter)      #
            pos_dic = {}   
            form.elementlist = []                               ######      
            for (index,element) in enumerate(form.intimes):           #
                if element in pos_dic:                          #
                    pos_dic[element].append(index)              #  search integration time duplicates
                else:                                           #
                    pos_dic.update({element:[index]})           #
                    form.elementlist.append(element)                 #############   liste der integrationszeiten                                                                  
            pos_list = list(pos_dic.values())                   #  liste mit arrays der  locations der duplicate
            form.adiz = len(form.elementlist)                   # anzahl der integrationszeiten
            form.spectrasum = []                                                                #
            for i in range (0, len(pos_list)):    
                meanspec0 =0                                      #
                for j in range (0,len(pos_list[i])):                                    #  aufsummieren der duplicate und durchschnitt erstellen
                    val = pos_list[i][j]                                                                       #  -> neue liste spectrasum
                    meanspec0 = np.add(form.spectra[val],meanspec0) 
                meanspec=meanspec0/len(pos_list[i])   
                form.spectrasum.append(meanspec)
            spektren.append(form.spectrasum)   #Liste aller Spektren: spectrasum
def flipspec(form,spektren):   
    for i in range (0,len(spektrum)):
        platzhalter3 = []
        platzhalter4 = spektrum[i]
        for j in range (0,len(platzhalter4)):
            platzhalter1 = platzhalter4[j]
            ende = len(platzhalter1)-1
            platzhalter2 = []                                               
            for k in range (0,len(platzhalter1)):                                                    
                platzhalter2.append(platzhalter1[ende-k])                                
            platzhalter3.append(platzhalter2)        
        spektren.append(platzhalter3)
def spadarkextract(form,inputtime,inputspectra): #takes a spa file and makes 1 spectra for each int.time
    pos_dic = {}   
    form.elementlist = []   
                                ######      
    for (index,element) in enumerate(inputtime):           #
        if element in pos_dic:                          #
            pos_dic[element].append(index)              #  search integration time duplicates
        else:                                           #
            pos_dic.update({element:[index]})           #
            form.elementlist.append(element)                 #############   liste der integrationszeiten                                                                  
    pos_list = list(pos_dic.values())       
    print(pos_list)            #  liste mit arrays der  locations der duplicate
    form.adiz = len(form.elementlist)                   # anzahl der integrationszeiten
    outputlist = []    
                                                            #
    for i in range (0, len(pos_list)):    
        meanspec0 =0                                      #
        for j in range (0,len(pos_list[i])):                                    #  aufsummieren der duplicate und durchschnitt erstellen
            val = pos_list[i][j]                                                                       #  -> neue liste spectrasum
            meanspec0 = np.add((inputspectra[val]),meanspec0) 
                    
        meanspec=meanspec0/len(pos_list[i])                                       
        outputlist.append(meanspec)  
        
    print(form.elementlist)      
    result = {'spectrum': outputlist,'zeiten': form.elementlist}
    return(result) 

def CoreSet(form,adr):
    IP=adr
    #IP="10.131.51.106"
    form.cora = Orca.remote_direct(IP, ignore_attenuation=False)
    n_profiles = form.cora.firmware.ProfileCountGet()['profile_count']      
def Setwavelength(form,x): 
    form.cora.set_profile_by_wavelength(x)   
def CoreAquireSpectra(form,timer):
    try:
        print("-- capturing spectra --")        
        form.spectra=form.cora.acquire_raw1d_spectrum(timer) 
       # SetlaserOff(form) 
        print("Spectra finished")  
        return(form.spectra)  
    except Exception as err:
            tb = sys.exc_info()[2]           
            print("something went in line:",tb.tb_lineno)  
            print("Error:",err)
def CoreAquireDark(form,timer):
    try:
        print("-- capturing dark--")
        form.dark=form.cora.acquire_raw1d_spectrum(timer) 
        print("Dark finished")   
        return(form.dark)         
    except Exception as err:
            tb = sys.exc_info()[2]           
            print("something went in line:",tb.tb_lineno)  
            print("Error:",err)  
def Set532laserOn(form):
    try:
        form.cora.firmware.LaserArmSet(True)
        form.cora.firmware.LaserElectricalPowerPercentSet(form.lineEdit_20.text())
        form.cora.firmware.LaserOpticalPowerSet(form.lineEdit_25.text())
        form.cora.firmware.LaserEnableSet(True)  
        time.sleep(2)  
    except Exception as err:
            tb = sys.exc_info()[2]           
            print("something went in line:",tb.tb_lineno)  
            print("Error:",err) 
def Set785laserOn(form):
    try:
        form.cora.firmware.LaserArmSet(True)
        form.cora.firmware.LaserElectricalPowerPercentSet(form.lineEdit_28.text())
        form.cora.firmware.LaserOpticalPowerSet(form.lineEdit_33.text())
        form.cora.firmware.LaserEnableSet(True)  
        time.sleep(2)  
    except Exception as err:
            tb = sys.exc_info()[2]           
            print("something went in line:",tb.tb_lineno)  
            print("Error:",err)
def Set1064laserOn(form):
    try:
        form.cora.firmware.LaserArmSet(True)
        form.cora.firmware.LaserElectricalPowerPercentSet(form.lineEdit_36.text())
        form.cora.firmware.LaserOpticalPowerSet(form.lineEdit_37.text())
        form.cora.firmware.LaserEnableSet(True)  
        time.sleep(2)  
    except Exception as err:
            tb = sys.exc_info()[2]           
            print("something went in line:",tb.tb_lineno)  
            print("Error:",err)       
def SetlaserOff(form):
    try:
        form.cora.firmware.LaserArmSet(False)
        form.cora.firmware.LaserEnableSet(False)  
        time.sleep(2)
    except Exception as err:
            tb = sys.exc_info()[2]           
            print("something went in line:",tb.tb_lineno)  
            print("Error:",err) 
def singlederivation(form,spektrum):  
    form.std=np.std(spektrum)
    return(form.std)   
def diffderivation(form,spektrum,a):             
    diff=spektrum[a]-spektrum[a-1]
    form.sstd=np.std(diff)   
    return(form.sstd)               
def savefile(form,file,name):  
    np.savetxt(name,list(file), fmt="%s")#,fmt='%1.3f')
def starmarker(form):
    form.star = mpath.Path.unit_regular_star(6)                                  ########## 
    form.circle = mpath.Path.unit_circle()                                       # 
    form.verts = np.concatenate([form.circle.vertices, form.star.vertices[::-1, ...]])     #
    form.codes = np.concatenate([form.circle.codes, form.star.codes])                      # define a star for plot
    form.cut_star = mpath.Path(form.verts, form.codes)    
def peakfit_old(form,spektrum):    
    rangediff=form.range2-form.range1-1
    form.npx=np.arange(form.range1,form.range2).astype(float)
    form.m = spektrum                       
    npy1 = form.m[form.range1:form.range2]              
    form.x = [form.npx[0],form.npx[rangediff]]
    form.y = [npy1[0],npy1[rangediff]]
    coefficients = np.polyfit(form.x,form.y,1)
    polynom = np.poly1d(coefficients) 
    form.npy = npy1 - polynom(form.npx)                                                                                           
    mean = sum(form.npx * form.npy) / sum(form.npy)                        
    sigma = np.sqrt(abs(sum(form.npy * (form.npx - mean)**2) / sum(form.npy)))             ##für gauss/lorentz                   
    form.popt,form.pcov = curve_fit(lorentzian,form.npx,form.npy,p0=[max(form.npy),mean,sigma])  
    form.peakmax=max(lorentzian(form.npx,*form.popt))
    form.maxvalue=max(npy1)
    form.maxvaluewb=max(form.npy)
def peakfit(form,spektrum):   
    bg(form,spektrum) 
    rangediff=form.range2-form.range1-1
    form.npx=np.arange(form.range1,form.range2).astype(float)
    form.m = spektrum                       
    npy1 = form.m[form.range1:form.range2]  
    form.npy = npy1 - form.polynom(form.npx)                                                                                           
    mean = sum(form.npx * form.npy) / sum(form.npy)                          
    sigma = np.sqrt(abs(sum(form.npy * (form.npx - mean)**2) / sum(form.npy)))             ##für gauss/lorentz                   
    form.popt,form.pcov = curve_fit(lorentzian,form.npx,form.npy,p0=[max(form.npy),mean,sigma])  
    form.peakmax=max(lorentzian(form.npx,*form.popt))
    form.maxvalue=max(npy1)
    form.maxvaluewb=max(form.npy)    
def bg(form,spektrum):  
    range1l=form.range1-50 
    range2h=form.range2+50 
    rangediff1=form.range1-range1l-1
    rangediff2=range2h-form.range2-1
    partrange1=np.arange(range1l,form.range1).astype(float)
    partrange2=np.arange(form.range2,range2h).astype(float)
    form.fullXrange=partrange1.tolist() + partrange2.tolist()
    form.m = spektrum           
    form.Xrange=np.arange(range1l,range2h).astype(float)            
    npy1 = form.m[range1l:form.range1] 
    npy2 = form.m[form.range2:range2h] 
    form.fullYrange=npy1.tolist() + npy2.tolist()          
    coefficients = np.polyfit(form.fullXrange,form.fullYrange,int(form.lineEdit_38.text()))
    form.polynom = np.poly1d(coefficients) 
    form.bgnpy = form.polynom(form.Xrange)                                                                                             
def fitrange532c(form):
    form.range1=410
    form.range2=450   
def fitrange785c(form):
    form.range1=625
    form.range2=700
def fitrange1064c(form):
    form.range1=65
    form.range2=75 
def fitrange532p(form):
    form.range1=495
    form.range2=508   
def fitrange785p(form):
    form.range1=760
    form.range2=772
def fitrange1064p(form):
    form.range1=65
    form.range2=75 
def fitrangeSRM(form):
    form.range1=1
    form.range2=2000
def fitrangeSRM1064(form):
    form.range1=1
    form.range2=250
def createfolder(form,type):
    datum=datetime.now().strftime("%d-%m-%Y_%I-%M-%S_%p")
    form.newpath="C:\\Users\APDE1.RD\\Desktop\\Cora5x\\Recorded\\"+type+datum  
    try:
        os.mkdir(form.newpath)
    except OSError:
        print ("Creation of the directory failed" )
    else:
        print ("Successfully created the directory" )        
    os.chdir(form.newpath)   
def readprofiles(form,pro):
    #CoreSet(form) 
    form.cora.firmware.ProfileSelectSet(pro)
    profile=form.cora.firmware.ProfileWavelengthNominalGet(pro)
    form.wellenlaenge = profile['wavelength_nm'] 
def grap(form,fenster,ort):
        pix=fenster.grab()
        counter = 1  
        while os.path.isfile(ort.format(counter)):
            counter += 1
        ort = ort.format(counter)
        pix.save(ort) 