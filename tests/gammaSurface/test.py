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

def plotGammaSurface(mat,gammaSurface,lab,pts,nc,nR,C):
    energyDensityScaling=mat.mu_SI*mat.b_SI
    fig, axs = plt.subplots(2,2)
    axs[0, 0].set_aspect("equal")
    axs[0, 1].set_aspect("equal")
    a0=gammaSurface.A[:,0]
    a1=gammaSurface.A[:,1]
    r0=gammaSurface.B[:,0]
    r1=gammaSurface.B[:,1]
    X=np.zeros((pts+1,pts+1))
    Y=np.zeros((pts+1,pts+1))
    GSFE=np.zeros((pts+1,pts+1))
    for i in range(pts+1):
        for j in range(pts+1):
            sL=i/pts*a0+j/pts*a1 # slip vector
            #sL=planeValue.localSlipVector(s) # slip vector in the plane local coordinates
            X[j,i]=sL[0]
            Y[j,i]=sL[1]
            GSFE[j,i]=gammaSurface.misfitEnergy(sL)*energyDensityScaling
    
    pmsh=axs[0, 0].pcolormesh (X, Y, GSFE,cmap=cm.coolwarm,shading='gouraud')
    axs[0, 0].arrow(0, 0, a0[0], a0[1])
    axs[0, 0].arrow(0, 0, a1[0], a1[1])
    axs[0, 1].arrow(0, 0, r0[0], r0[1])
    axs[0, 1].arrow(0, 0, r1[0], r1[1])
    for i in range(-nR,nR):
        for j in range(-nR,nR):
            r=i*r0+j*r1
            axs[0, 1].scatter(r[0],r[1],0.1,color='k')
            rcirc=plt.Circle((0, 0), np.linalg.norm(r), color='0.5',fill=False,linewidth=0.25)
            axs[0, 1].add_patch(rcirc)
    axs[0, 1].scatter(gammaSurface.waveVectors[:,0],gammaSurface.waveVectors[:,1],color='m',s=1)
    axs[0, 0].scatter(gammaSurface.points[:,0],gammaSurface.points[:,1],color='k',s=1)
    
    e=np.zeros((nc))
    c=gammaSurface.A@C
    for j in range(0,np.shape(c)[1]):
        cj=c[:,j]
        for i in range(0,nc):
            e[i]=gammaSurface.misfitEnergy(cj/nc*i)*energyDensityScaling
        axs[1, 0].plot(np.linspace(0, np.linalg.norm(cj), num=nc, endpoint=False),e)
    fig.savefig(mat.materialName+lab+".pdf", bbox_inches='tight', format='pdf')

##########################################################################
# Main
materialFile="../../Library/Materials/W.txt"
absoluteTemperature=0.0
mat=pyMoDELib.PolycrystallineMaterialBase(materialFile,absoluteTemperature)
energyDensityScaling=mat.mu_SI*mat.b_SI
C2G=np.array([[1,0,0],[0,1,0],[0,0,1]])
singleCrystal=pyMoDELib.SingleCrystalBase(mat,C2G)

# initialize
nR=4
pts=50 # GSFE resolution in each direction
nc=100 # GSFE cuts resolution

C=np.array([[1,1],
            [0,1]])

for planeKey, planeValue in singleCrystal.planeNormals().items(): # loop over lattice planes
    #b0=planeValue.primitiveVectors[0].cartesian() # first primitive vector of the matrix lattice plane
    #b1=planeValue.primitiveVectors[1].cartesian() # second primitive vector of the matrix lattice plane
    # Matrix
    if planeValue.gammaSurface is not None:
        plotGammaSurface(mat,planeValue.gammaSurface,"_plane_"+str(planeKey)+"_matrix",pts,nc,nR,C)
    # Second Phases
    for phaseKey, phaseValue in singleCrystal.secondPhases().items():
        if phaseValue.gammaSurface(planeKey) is not None:
            plotGammaSurface(mat,phaseValue.gammaSurface(planeKey),"_plane_"+str(planeKey)+"_"+phaseValue.name,pts,nc,nR,C)
