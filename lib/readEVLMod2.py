import sys, string, os
import numpy as np
# evl file

class EVL:
    nodes=np.empty([0,0])
    
def readEVLtxt(filename):
    evlFile = open(filename+'.txt', "r")
    numNetNodes=int(evlFile.readline().rstrip())
    numLoops=int(evlFile.readline().rstrip())
    numLinks=int(evlFile.readline().rstrip())
    numLoopNodes=int(evlFile.readline().rstrip())
    numSpInc=int(evlFile.readline().rstrip())
    numPolyInc=int(evlFile.readline().rstrip())
    numPolyIncNodes=int(evlFile.readline().rstrip())
    numPolyIncEdges=int(evlFile.readline().rstrip())
    numEDrow=int(evlFile.readline().rstrip())
    numCDrow=int(evlFile.readline().rstrip())
    evl=EVL();
    evl.nodes=np.empty([numNetNodes, 3])
    for k in np.arange(numNetNodes):
        data=np.fromstring(evlFile.readline().rstrip(), sep=' ')
        evl.nodes[k,:]=data[1:4]
    return evl
