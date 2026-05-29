import sys, string, os
import numpy as np
# evl file

class EVLM1:
    nodes=np.empty([0,0])
    loops=np.empty([0])

def readEVLM1txt(filename):
    evlFile = open(filename+'.txt', "r")
    numNodes=int(evlFile.readline().rstrip())
    numLoops=int(evlFile.readline().rstrip())
    numLinks=int(evlFile.readline().rstrip())
    evl=EVLM1();
    evl.nodes=np.empty([numNodes, 3])
    evl.loops=np.empty([numLoops, 1])
    for k in np.arange(numNodes):
        data=np.fromstring(evlFile.readline().rstrip(), sep=' ')
        evl.nodes[k,:]=data[1:4]
    for k in np.arange(numLoops):
        data=np.fromstring(evlFile.readline().rstrip(), sep=' ')
        evl.loops[k]=data[-1]
    return evl



