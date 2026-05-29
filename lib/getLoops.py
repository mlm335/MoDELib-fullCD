import os,glob,sys
import numpy as np


class node:
    def __init__(self,p):
        self.p = np.array([s for s in p])
        
class loop:
    def __init__(self,b,n):
        self.b = np.array([s for s in b])
        self.n = np.array([s for s in n])
        self._c = np.array([0., 0., 0.])
        self._r = 0.0
        self._links = []
        
    def addLinks(self,source,sink):
        self._links.append( [source,sink] )
        
    def updateCenter(self,nodelist):
        for link in self._links:
            for n in link:
                self._c += nodelist[n].p/len(self._links)/len(link)
    
    def getRadius(self,nodelist):
        for link in self._links:
            for n in link:
                self._r += np.sqrt(np.sum(np.square(nodelist[n].p-self._c)))/len(self._links)/len(link)


def getLoops(filename):
#----------------------------------------------------------------------------#
#			        Get input file
#----------------------------------------------------------------------------#
    ff = filename+'.txt'
#----------------------------------------------------------------------------#
#				  Get three lists
#----------------------------------------------------------------------------#

    nlist = {} # Point list
    llist = {} # Loop list

    fcon = open(ff)
    ntotal = int(fcon.readline().strip())
    ltotal = int(fcon.readline().strip())
    stotal = int(fcon.readline().strip())


    for line_number in range(ntotal):
        line = fcon.readline().strip().split()
        nlist.update({ int(line[0]): node(np.array([float(s) for s in line[1:4]])) })

    for line_number in range(ltotal):
        line = fcon.readline().strip().split()
        llist.update({int(line[0]): loop(np.array([float(s) for s in line[1:4]]),np.array([float(s) for s in line[4:7]]))})
        
    for line_number in range(stotal):
        line = fcon.readline().strip().split()
        llist[int(line[0])].addLinks(int(line[1]),int(line[2]))
        
#----------------------------------------------------------------------------#
#				  Update loop center and get radius
#----------------------------------------------------------------------------#
    for l in llist:
        llist[l].updateCenter(nlist)
        llist[l].getRadius(nlist)
        
    return llist
