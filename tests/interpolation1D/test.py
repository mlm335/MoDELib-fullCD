# /opt/local/bin/python3.13 test.py
import sys, string, os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
plt.rcParams['text.usetex'] = True
sys.path.append("../../python")
from modlibUtils import *
sys.path.append("../../build/tools/pyMoDELib")
import pyMoDELib


#XiYi=np.array([[0.0,1.0],
#               [30.0,10.0],
#               [54.0,0.5],
#               [90.0,45.0],
#               [180-54.0,0.5],
#               [180-30.0,17.0]
#               ])
               
XiYi=np.random.rand(6, 2)
dX=np.max(XiYi[:,0])-np.min(XiYi[:,0])

extrapolation=pyMoDELib.ExtrapolationMethod(2,1.1*dX) # 0=default, 1=linear, 2=periodic
lagrangeInterp=pyMoDELib.LagrangeInterpolant(XiYi,extrapolation)
linearInterp=pyMoDELib.LinearInterpolant(XiYi,extrapolation)

X=np.linspace(np.min(XiYi[:,0])-2, np.max(XiYi[:,0])+2, num=10000, endpoint=True)
Y1=np.zeros(len(X))
Y2=np.zeros(len(X))
for i in range(len(X)):
    Y1[i]=lagrangeInterp.f(X[i])
    Y2[i]=linearInterp.f(X[i])
#    Y1[i]=lagrangeInterp.atPeriodic(X[i],np.min(XiYi[:,0])-0.2,np.max(XiYi[:,0])+0.2)
#    Y2[i]=linearInterp.atPeriodic(X[i],np.min(XiYi[:,0])-0.2,np.max(XiYi[:,0])+0.2)


fig1 = plt.figure()
ax1=plt.subplot(1,1,1)
ax1.plot(XiYi[:,0], XiYi[:,1],'o')
#ax1.plot(X, Y1,label='lagrange')
ax1.plot(X, Y2,label='linear',linestyle='dashed')
plt.legend()
ax1.grid()
plt.show()
