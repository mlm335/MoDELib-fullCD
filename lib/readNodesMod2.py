import sys, string, os
import numpy as np

class NODES:
    nodesPos=np.empty([0,0])
    nodesV=np.empty([0,0])
    nodesVOld=np.empty([0,0])
    nodesClimbVScalar=np.empty([0,0])

    loopsNumber=np.empty(0)
    loopsArea=np.empty(0)
    loopsBurger=np.empty([0,0])
    loopsNormal=np.empty([0,0])
    loopsPos=np.empty([0,0])

    loopNumber=np.empty(0)
    loopNodeNumber=np.empty(0)
    loopNodePos=np.empty([0,0])

def readNODEStxt(filename):
    try: 
        with open(filename+'.txt', "r") as evlFile:
            numNodes=int(evlFile.readline().rstrip())
            numLoops=int(evlFile.readline().rstrip())
            numLinks=int(evlFile.readline().rstrip())
            numLoopNodes=int(evlFile.readline().rstrip())
            numSpInc=int(evlFile.readline().rstrip())
            numPolyInc=int(evlFile.readline().rstrip())
            numPolyIncNodes=int(evlFile.readline().rstrip())
            numPolyIncEdges=int(evlFile.readline().rstrip())
            numEDrow=int(evlFile.readline().rstrip())
            numCDrow=int(evlFile.readline().rstrip())
        
            nd=NODES();
            
            nd.nodesPos=np.empty([numNodes, 3])
            nd.nodesV=np.empty([numNodes, 3])
            nd.nodesVOld=np.empty([numNodes, 3])
            nd.nodesClimbVScalar=np.empty([numNodes, 2])

            nd.loopsNumber=np.empty(numLoops)
            nd.loopsArea=np.empty(numLoops)
            nd.loopsBurger=np.empty([numLoops, 3])
            nd.loopsNormal=np.empty([numLoops, 3])
            nd.loopsPos=np.empty([numLoops, 3])

            nd.loopNumber=np.empty(numLoopNodes)
            nd.loopNodeNumber=np.empty(numLoopNodes)
            nd.loopNodePos=np.empty([numLoopNodes,3])

            # node data
            for k in range(numNodes):
                data = np.fromstring(evlFile.readline().strip(), sep=' ')
                nd.nodesPos[k, :] = data[1:4]
                nd.nodesV[k, :] = data[4:7]
                nd.nodesClimbVScalar[k, :] = data[7:9]
                
            # loop data
            for k in range(numLoops):
                data = np.fromstring(evlFile.readline().strip(), sep=' ')
                nd.loopsNumber[k] = data[0]
                nd.loopsArea[k] = data[-1]
                nd.loopsBurger[k, 0:3] = data[1:4]
                nd.loopsNormal[k, 0:3] = data[4:7]
                nd.loopsPos[k, 0:3] = data[7:10]

            # loop links
            for k in range(numLinks):
                data = np.fromstring(evlFile.readline().strip(), sep=' ')

            # loop nodes
            for k in range(numLoopNodes):
                data = np.fromstring(evlFile.readline().strip(), sep=' ')
                nd.loopNumber[k]=data[0]
                nd.loopNodeNumber[k]=data[1]
                nd.loopNodePos[k, 0:3]=data[2:5]

        return nd

    except Exception as e:
        print(f"Error reading file {filename}: {e}")
    return None
