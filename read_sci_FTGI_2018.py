# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 11:15:44 2017

@author: rjh
"""
from tkinter import Tk
from tkinter.filedialog import askopenfilename
import numpy as np

Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
filename = askopenfilename(filetypes=[("binary scientific","*.sci2")]) # show an "Open" dialog box and return the path to the selected file
print(filename)

data=np.empty(shape=(124,60*4000),dtype=float)
dat=np.empty(4000,dtype='i8')
with open(filename, "rb") as f:
    for k in range(60):
        print(k)
        for i in range(124):
            campo1= np.fromfile(f,count=1, dtype='>a2')
            campo2= np.fromfile(f,count=1, dtype='>u2')            
            campo3= np.fromfile(f,count=1, dtype='>u2')
            campo4= np.fromfile(f,count=1, dtype='>u1')
            campo5= np.fromfile(f,count=1, dtype='>u1')
            campo6= np.fromfile(f,count=1, dtype='>u1')
            campo7= np.fromfile(f,count=1, dtype='>u1')
            Tstamp= np.fromfile(f,count=1, dtype='>u4')
            Time=   np.fromfile(f,count=1, dtype='>u4')
            SID=    np.fromfile(f,count=1, dtype='>u2')
            AqErr=  np.fromfile(f,count=1, dtype='>u1')
            channel_ID=  np.fromfile(f,count=1, dtype='>u1')
            cal_flag=  np.fromfile(f,count=1, dtype='>u1')
            cal_sw=  np.fromfile(f,count=1, dtype='>u1')
            ph_sw=  np.fromfile(f,count=1, dtype='>u1') 
            repuesto1=  np.fromfile(f,count=1, dtype='>u1') 
            Nphstates=  np.fromfile(f,count=1, dtype='>u4')
            PHseq= np.fromfile(f,count=16, dtype='>u1')
            Processtype= np.fromfile(f,count=1, dtype='>u1')
            repuesto2= np.fromfile(f,count=1, dtype='>u1')
            Pack_count=  np.fromfile(f,count=1, dtype='>u2')
            Sc_factor=np.fromfile(f,count=1, dtype='>f8')
            Samprate=  np.fromfile(f,count=1, dtype='>u4')
            NSample=  np.fromfile(f,count=1, dtype='u4')
            dat= np.fromfile(f,count=4000, dtype='>i4')
            dat=dat*Sc_factor
            data[i,(k*4000):(k*4000)+4000]=dat
           
            print(k,i)      
            print('campo1 = ' + str(campo1))
            print('campo2 = ' + str(campo2))
            print('campo3 = ' + str(campo3))
            print('campo4 = ' + str(campo4))
            print('campo5 = ' + str(campo5))
            print('campo6 = ' + str(campo6))
            print('campo7 = ' + str(campo7))
            print('Tstamp = ' + str(Tstamp))
            print('Time = ' + str(Time))
            print('SID = ' + str(SID))
            print('AqErr = ' + str(AqErr)) 
            print('channel_ID = ' + str(channel_ID))
            print('cal_flag = ' + str(cal_flag))
            print('cal__sw = ' + str(cal_sw))
            print('ph__sw = ' + str(ph_sw) + ' (1- 16KHz, 2- 8KHz)')
            print('repuesto1 = ' + str(repuesto1))
            print('Nphstates = ' + str(Nphstates))
            print('Phase sequence = ' + str(PHseq))
            print('Process type = ' + str(Processtype))
            print('repuesto2 = ' + str(repuesto2))
            print('packet counter = ' + str(Pack_count))
            print('Scale factor = ' + str(Sc_factor))
            print('Sampling rate = ' + str(Samprate))
            print('No of samples = ' + str(NSample))
            print(dat[0:10])
            
