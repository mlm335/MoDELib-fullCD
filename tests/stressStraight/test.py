# /opt/local/bin/python3.13 test.py
import sys, string, os, time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
plt.rcParams['text.usetex'] = True
sys.path.append("../../python")
from modlibUtils import *
sys.path.append("../../build/tools/pyMoDELib")
import pyMoDELib
    
################################
# main
################################
materialFile="../../Library/Materials/W.txt"
mat=pyMoDELib.PolycrystallineMaterialBase(materialFile,300) # material object

P0=np.array([0.0,0.0,0.0])
P1=np.array([1.0,0.0,0.0])
b=np.array([1.0,0.0,0.0])
x=np.array([1.0,1.0,1.0])
ss=pyMoDELib.StressStraight(mat,P0,P1,b,0.0)
start_time = time.perf_counter() # Use perf_counter() for high resolution
for i in range(0,10000000):
    sigma=ss.stress(x)
end_time = time.perf_counter()
execution_time = end_time - start_time
print(f"Execution time: {execution_time} seconds")
print(sigma)
