# -*- coding: utf-8 -*-
"""
Created on Sat Feb 12 15:32:42 2022

@author: mkzay
"""
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

L = 5  # m
H = 3  # m
T0 = 400  # K
Tx = 45  # K
Ty = 35  # K
Txy = 27.5  # K
ax = 1/3.
ay = 1/4.
axy = 1/2.

# Récupérer les MMS et les dériver
x,y = sp.symbols('x y')
T_MMS=T0 + Tx*sp.cos((ax*np.pi*x)/L)+Ty*sp.sin((ay*np.pi*y)/H)+Txy*sp.sin((axy*np.pi*x*y)/(L*H))
f_T_MMS = sp.lambdify([x,y], T_MMS, "numpy")
source = sp.diff(sp.diff(T_MMS,x),x)+sp.diff(sp.diff(T_MMS,y),y)
f_source = sp.lambdify([x,y], source, "numpy")

# Test de la MMS
GRIDx = np.linspace(0,L,100)
GRIDy = np.linspace(0,H,100)
T=np.zeros((len(GRIDx),len(GRIDy)))

for i in range(len(GRIDx)):
        for j in range(len(GRIDy)):
            T[i,j]=f_T_MMS(GRIDx[i],GRIDy[j])
  
#Comparaison des solutions
c = plt.pcolor(GRIDx, GRIDy, T)
plt.colorbar(c)

