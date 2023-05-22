# -*- coding: utf-8 -*-
"""
Created on Mon May 22 15:04:07 2023

@author: user
"""

from scipy import fftpack
import math as m
from numpy import sqrt, mean, absolute, real, conj
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
def std(vecstd):
  #Funcion para la desviacion estandar
  return np.std(vecstd)
def mean(meanvec):
  #Funcion para el promedio 
  return np.mean(meanvec)

def counth(vech):
  #Funcion para hacer la cuenta de los valores encima de 1
  vecn1=abs(np.array(vech))
  count1=0
  for i in vecn1:
    if abs(i)>1:
      count1=count1+1
  return count1
def countlow(vecl):
  #Funcion para contar las veces que el valor es menor a 0.1
  veclow=abs(np.array(vecl))
  
  clow=0
  for i in veclow:
    if abs(i)<=0.1:
      clow=clow+1
  return clow
def countinter(vint):
  #Funcion para contar las veces que se encuentra en el intervalo 0.1 y 1
  vinter=abs(np.array(vint))
  cinter=0
  for i in vinter:
    if abs(i)>0.1 and abs(i)<=1:
          cinter=cinter+1
  return cinter

###Le haremos la FFT para sacar una caracteristica llamada THD
def RMS(Vrms1):
  Vrms= np.abs(Vrms1)
  return np.sqrt(sum(np.array(Vrms)**2))
def THD(Vthd):
  Vthd= np.abs(Vthd)
  
  return np.sqrt(((RMS(Vthd)/list(Vthd)[1])**2)-1)*100


def THDIA(vTHDIA1):
  "Esta funcion me saca el THD promedio de todos los ciclos, el menor y el mayor"
  #Dividimos el vector en su numero de ciclos en una lista 
  vTHDIA=list(vTHDIA1)
  interval=[]
  ffi=[]
  THDr=[]
  
  for i in range(1,m.floor(len(vTHDIA)/53.3)-1):
      num1=m.floor(53*(i-1))
      num2=m.floor(53*(i))
      
      interval=vTHDIA[num1:num2]
      #Hacemos la transformada de fourier y lo añadimos a una lista
      ffi=fftpack.fft(interval)
      
      #Le aplicamos el THD
      THDr.append(THD(list(ffi)[0:7]))
      #Le aplicamos el LHD que es THD pero considerando solo los primeros 7 frecuencias
  
  #Vamos a retornar una lista de 6 numeros 
  mp=[np.mean(THDr),min(THDr),max(THDr)]
  return mp
def preproceso(dfs):
  x=[]
  for i in range(len(dfs.T)-1):
    ldfs=list(np.transpose(np.array(dfs[str(i)])))
    carac=[] 
    carac.append(std(ldfs))
    carac.append(mean(ldfs))
    carac.append(counth(ldfs))
    carac.append(countlow(ldfs))
    carac.append(countinter(ldfs))
    TH=THDIA(ldfs)
    carac.append(TH[0])
    carac.append(TH[1])
    carac.append(TH[2])
    x.append(carac)
  return x
def data(dataS,dataN):
  xS=preproceso(dataS)
  xN=preproceso(dataN)
  xData=[]
  yData=[]
  
  for i in range (len(xS)):
    xData.append(xS[i])
    xData.append(xN[i])
    yData.append(1)
    yData.append(0)
  
  return xData, yData
"""
Ntest=pd.read_csv("Normal1.csv")
Stest=pd.read_csv("sag2.csv")
xtest,ytest=data(Stest,Ntest)
print(xtest) 
print("ok")
"""