import sys, string, os
import numpy as np
# evl file

class EVL:
    nodes=np.empty([0,0])
    loopsArea=np.empty([0])
    burgers=np.empty([0,0])
    normal=np.empty([0,0])
    position=np.empty([0,0])
    
def readEVLtxt(filename):
    evlFile = open(filename+'.txt', "r")
    numNodes=int(evlFile.readline().rstrip())
    numLoops=int(evlFile.readline().rstrip())
    numLinks=int(evlFile.readline().rstrip())
    evl=EVL();
    evl.nodes=np.empty([numNodes, 3])
    evl.loopsArea=np.empty([numLoops, 1])
    evl.burgers=np.empty([numLoops, 3])
    evl.normal=np.empty([numLoops, 3])
    evl.position=np.empty([numLoops, 3])
    for k in np.arange(numNodes):
        data=np.fromstring(evlFile.readline().rstrip(), sep=' ')
        evl.nodes[k,:]=data[1:4]
    for k in np.arange(numLoops):
        data=np.fromstring(evlFile.readline().rstrip(), sep=' ')
        evl.loopsArea[k]=data[-1]
        evl.burgers[k,0:3]=data[1:4]
        evl.normal[k,0:3]=data[4:7]
        evl.position[k,0:3]=data[7:10]
    return evl
