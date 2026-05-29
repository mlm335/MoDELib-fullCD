import sys, string, os
import numpy as np
# evl file

class NODES:
    nodesPos=np.empty([0,0])
    nodesV=np.empty([0,0])
    nodesVOld=np.empty([0,0])
    nodesClimbVScalar=np.empty([0,0])
    loopsArea=np.empty([0])

def readNODEStxt(filename):

    evlFile = open(filename+'.txt', "r")
    numNodes=int(evlFile.readline().rstrip())
    numLoops=int(evlFile.readline().rstrip())
    numLinks=int(evlFile.readline().rstrip())
    
    nd=NODES();
    
    nd.nodesPos=np.empty([numNodes, 3])
    nd.nodesV=np.empty([numNodes, 3])
    nd.nodesVOld=np.empty([numNodes, 3])
    nd.nodesClimbVScalar=np.empty([numNodes, 2])
    nd.loopsArea=np.empty([numLoops, 1])

    for k in np.arange(numNodes):
        data=np.fromstring(evlFile.readline().rstrip(), sep=' ')
        nd.nodesPos[k,:]=data[1:4]
        nd.nodesV[k,:]=data[4:7]
        nd.nodesVOld[k,:]=data[7:10]
        nd.nodesClimbVScalar[k,:]=data[11:13]

    for k in np.arange(numLoops):
        data=np.fromstring(evlFile.readline().rstrip(), sep=' ')
        nd.loopsArea[k]=data[-1]
#        nd.burgers[k,0:3]=data[1:4]
#        nd.normal[k,0:3]=data[4:7]
    return nd
