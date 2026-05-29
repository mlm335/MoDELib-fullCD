# /opt/local/bin/python3.14 test.py
import sys, string, os
import numpy as np
import pathlib as pt
import math
import matplotlib.pyplot as plt
from matplotlib import cm
plt.rcParams['text.usetex'] = True
sys.path.append("../../python")
from modlibUtils import *
sys.path.append("../../build/tools/pyMoDELib")
import pyMoDELib

def plotMobility(m,mat,singleCrystal,slipSystem,ax1,ax2):
    s0=0.01 # stress amplitude
    S=np.outer(m,m)*s0 # stress tensor

    theta=np.linspace(0, 360, num=360, endpoint=True)*np.pi/180
    vel=np.zeros(len(theta))

    for k in range(len(theta)):
        xi=angleAxis(theta[k],slipSystem.unitNormal)@slipSystem.unitSlip # rotate b about n to define the line tangent
        vel[k]=slipSystem.velocity(S,xi,mat.T)
    ax1.plot(theta*180/np.pi, vel,label=str(mat.T)+'K')
    
    stress=np.linspace(0, 0.02, num=100)
    vels=np.zeros(len(stress))
    vele=np.zeros(len(stress))
    for k in range(len(stress)):
        S=np.outer(m,m)*stress[k] # stress tensor
        xi=angleAxis(0.0,slipSystem.unitNormal)@slipSystem.unitSlip # rotate b about n to define the line tangent
        vels[k]=slipSystem.velocity(S,xi,mat.T)
        xi=angleAxis(0.5*np.pi,slipSystem.unitNormal)@slipSystem.unitSlip # rotate b about n to define the line tangent
        vele[k]=slipSystem.velocity(S,xi,mat.T)

    ax2.plot(stress, vels,label=str(mat.T)+'K')
    ax2.plot(stress, vele,'--',label=str(mat.T)+'K')

def selectMobility(materialFile,Temps,slipSystemID,theta_deg,phi_deg):
    fig1 = plt.figure()
    ax1=plt.subplot(1,2,1)
    ax2=plt.subplot(1,2,2)
    for T in Temps:
        mat=pyMoDELib.PolycrystallineMaterialBase(materialFile,T) # material object
        singleCrystal=pyMoDELib.SingleCrystalBase(mat,np.identity(3))
        for ssKey, slipSystem in singleCrystal.slipSystems().items(): # loop over lattice planes
            if ssKey==slipSystemID:
                m=slipSystem.unitNormal*np.cos(theta_deg*math.pi/180.0)+slipSystem.unitSlip*np.sin(theta_deg*math.pi/180.0)*np.cos(phi_deg*math.pi/180.0)+np.cross(slipSystem.unitNormal, slipSystem.unitSlip)*np.sin(theta_deg*math.pi/180.0)*np.sin(phi_deg*math.pi/180.0)
                plotMobility(m,mat,singleCrystal,slipSystem,ax1,ax2)
    ax1.grid()
    ax1.legend()
    plt.xlabel('character angle $[deg]$')
    plt.ylabel(' v / c_s [-]')
    ax2.grid()
    ax2.legend()
    plt.xlabel('stress/mu $[-]$')
    plt.ylabel(' v / c_s [-]')
    fig1.savefig(materialName+"_"+str(slipSystemID)+".pdf", bbox_inches='tight')
#    ssCounter=ssCounter+1


    
    
#    if mat.crystalStructure=="FCC" or mat.crystalStructure=="CubicFluorite":
#        b=np.array([0,-1,1]) # burgers vector
#        b=b/np.linalg.norm(b)
#        n=np.array([1,1,1]) # plane normal
#        n=n/np.linalg.norm(n)
#        mob=mobilitySelector.getMobility(mat,mobilityType)
#        if mob is not None:
#            plotMobility(m,b,n,mat,mob,ax1,ax2)
#    elif mat.crystalStructure=="BCC":
#        b=np.array([1,1,1]) # burgers vector
#        b=b/np.linalg.norm(b)
#        n=np.array([-1,0,1]) # plane normal
#        n=n/np.linalg.norm(n)
#        mob=mobilitySelector.getMobility(mat,mobilityType)
#        if mob is not None:
#            plotMobility(m,b,n,mat,mob,ax1,ax2)
#    elif mat.crystalStructure=="HEX":
#        b=np.array([1,1,1]) # burgers vector
#        b=b/np.linalg.norm(b)
#        n=np.array([-1,0,1]) # plane normal
#        n=n/np.linalg.norm(n)
#        mob=mobilitySelector.getMobility(mat,mobilityType)
#        if mob is not None:
#            plotMobility(m,b,n,mat,mob,ax1,ax2)

################################
# main
################################
slipSystemID=0
Temps=np.array([100,300,1000]) # Kelvin
phi_deg=0.0
theta_deg=45 # m = n*cos(theta)+s*sin(theta)*cos(phi)+(nxs)*sin(theta)*sin(phi)
#m=np.array([0,0,1]) # stress axis
#m=m/np.linalg.norm(m); # normalized stress axis

for file in pt.Path('../../Library/Materials').glob('*.txt'):
    if file.is_file():
        materialFile=str(file)
        print(materialFile)
#        materialFile="../../Library/Materials/UO2.txt"
        Tm=getValueInFile(materialFile,"Tm")
        materialName=getStringInFile(materialFile,"materialName")
        selectMobility(materialFile,Temps,slipSystemID,theta_deg,phi_deg)


#        enabledSlipSystems=getStringInFile(materialFile,"enabledSlipSystems")
#        mobilitySelector=pyMoDELib.DislocationMobilitySelector()

#        for ssKey, ssValue in singleCrystal.planeNormals().items(): # loop over lattice planes
#            #b0=planeValue.primitiveVectors[0].cartesian() # first primitive vector of the matrix lattice plane
#            #b1=planeValue.primitiveVectors[1].cartesian() # second primitive vector of the matrix lattice plane
#            # Matrix
#            if planeValue.gammaSurface is not None:
#                plotGammaSurface(mat,planeValue.gammaSurface,"_plane_"+str(planeKey)+"_matrix",pts,nc,nR,C)
#            # Second Phases
#            for phaseKey, phaseValue in singleCrystal.secondPhases().items():
#                if phaseValue.gammaSurface(planeKey) is not None:
#                    plotGammaSurface(mat,phaseValue.gammaSurface(planeKey),"_plane_"+str(planeKey)+"_"+phaseValue.name,pts,nc,nR,C)



#        ssCounter=0
#        for ss in enabledSlipSystems.split():
#            fig1 = plt.figure()
#            ax1=plt.subplot(1,2,1)
#            ax2=plt.subplot(1,2,2)
#            mobilityType=getStringInFile(materialFile,"mobility_"+ss)
#            print("mobilityType="+mobilityType)
#            for T in Temps:
#                selectMobility(materialFile,mobilityType,T,m,Temps,ax1,ax2)
#            ax1.grid()
#            ax1.legend()
#            plt.xlabel('character angle $[deg]$')
#            plt.ylabel(' v / c_s [-]')
#            ax2.grid()
#            ax2.legend()
#            plt.xlabel('stress/mu $[-]$')
#            plt.ylabel(' v / c_s [-]')
#            fig1.savefig(materialName+"_"+str(ssCounter)+".pdf", bbox_inches='tight')
#            ssCounter=ssCounter+1
