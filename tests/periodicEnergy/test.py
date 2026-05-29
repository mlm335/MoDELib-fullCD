# /opt/local/bin/python3.13 test.py
import sys, string, os
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
sys.path.insert(0, '../../python')
from modlibUtils import *

compute=True
periodicImageMatrix=np.array([[0,0,0],[1,1,1],[2,2,2],[4,4,4]])
EwaldLengthDivisions=2
#periodicImageMatrix=np.array([[0,0,0],[1,1,0],[2,2,0],[4,4,0],[4,4,4]])
#EwaldLengthFactors=np.array([0.0,10.0])


if compute:
    curDir=os.getcwd()
    microstrucureFile=getStringInFile('inputFiles/initialMicrostructure.txt','microstructureFile')
    microstrucureFileName=curDir+'/inputFiles/'+microstrucureFile
    ddFileName=curDir+'/inputFiles/DD.txt'
    
    for k in range(0,np.size(periodicImageMatrix,0)):
        periodicImageSize=periodicImageMatrix[k]
        #EwaldLengthFactors=np.array([0,1,periodicImageSize[0]])
        #EwaldLengthFactors=np.array([periodicImageSize[0]])
#        EwaldLengthFactors=np.linspace(0, periodicImageSize[0], num=EwaldLengthDivisions)
        EwaldLengthFactors=np.linspace(0.0, 1.0, num=EwaldLengthDivisions)
        print("periodicImageSize="+str(periodicImageSize))
        print("EwaldLengthFactors="+str(EwaldLengthFactors))


        for c in range(0,np.size(EwaldLengthFactors)):
            #periodicImageSize=np.array([1,1,1])
            EwaldLengthFactor=EwaldLengthFactors[c];
            setInputVariable(ddFileName,'EwaldLengthFactor',str(EwaldLengthFactor))
            setInputVector(ddFileName,'periodicImageSize',periodicImageSize,'number of periodic images along each period vector')

            if microstrucureFile=='periodicDipoleIndividual.txt':
                setInputVariable(microstrucureFileName,'periodicDipoleExitFaceIDs',str(1)) # 1=edge, 0=screw
                #X=np.linspace(0, 200, num=40)
                X=np.linspace(0, 110, num=20)
#                X=np.linspace(0, 110, num=10)

                E=X*0.0;
                S11=X*0.0;
                S12=X*0.0;
                S13=X*0.0;
                S22=X*0.0;
                S23=X*0.0;
                S33=X*0.0;
                for k in range(0,len(X)):
                    x=X[k]
                    os.system('rm -r F/*.*')
                    os.system('rm -r evl/*.*')
                    setInputVariable(microstrucureFileName,'glideSteps',str(x))
                    os.system('../../build/tools/MicrostructureGenerator/microstructureGenerator '+ curDir)
                    os.system('../../build/tools/DDomp/DDomp '+ curDir)
                    F,Flabels=readFfile('./F')
                    E[k]=getFarray(F,Flabels,'dislocation elastic energy [mu b^3]')
                    aux =readAUXtxt('evl/ddAux_0')
                    S11[k]=np.mean(aux.gaussPoints[:,10])
                    S12[k]=np.mean(aux.gaussPoints[:,11])
                    S13[k]=np.mean(aux.gaussPoints[:,12])
                    S22[k]=np.mean(aux.gaussPoints[:,14])
                    S23[k]=np.mean(aux.gaussPoints[:,15])
                    S33[k]=np.mean(aux.gaussPoints[:,18])
            else:
                print('unsupported option')

            data=np.empty([np.size(X),8])
            data[:,0]=X
            data[:,1]=E
            data[:,2]=S11
            data[:,3]=S22
            data[:,4]=S33
            data[:,5]=S23
            data[:,6]=S13
            data[:,7]=S12
            #data=np.array([X.transpose(),E.transpose(),S11.transpose(),S22.transpose(),S12.transpose()])
            #data=np.array([X,E,S11,S22,S12]).reshape(np.size(X),5)
#            if periodicImageCentered:
#                outFilePrefix='centered'
#            else:
#                outFilePrefix='uncentered'

            outFileName='data'+'_'+str(periodicImageSize[0])+'_'+str(periodicImageSize[1])+'_'+str(periodicImageSize[2])+'_'+str(EwaldLengthFactor)+'.txt'
            np.savetxt(outFileName, data)

dataDict=dict()
#uncentered=dict()
for k in range(0,np.size(periodicImageMatrix,0)):
    periodicImageSize=periodicImageMatrix[k]
    #EwaldLengthFactors=np.array([0,1,periodicImageSize[0]])
#    EwaldLengthFactors=np.array([periodicImageSize[0]])
#    EwaldLengthFactors=np.linspace(0, periodicImageSize[1], num=EwaldLengthDivisions)
    EwaldLengthFactors=np.linspace(0.0, 1.0, num=EwaldLengthDivisions)


    for c in range(0,np.size(EwaldLengthFactors)):
        EwaldLengthFactor=EwaldLengthFactors[c];
        key=str(periodicImageSize[0])+'_'+str(periodicImageSize[1])+'_'+str(periodicImageSize[2])+'_'+str(EwaldLengthFactor)
        print(key)
        FileName='data_'+key+'.txt'
#    uFileName='uncentered_'+str(periodicImageSize[0])+'_'+str(periodicImageSize[1])+'_'+str(periodicImageSize[2])+'.txt'
        dataDict[key]=np.loadtxt(FileName);
#    uncentered[key]=np.loadtxt(uFileName);
    
    
fields=['energy','sigma11','sigma22','sigma33','sigma23','sigma13','sigma12']
fieldLabels=['elastic energy $[\mu b^3]$','$\sigma_{11}$ $[\mu]$','$\sigma_{22}$ $[\mu]$','$\sigma_{33}$ $[\mu]$','$\sigma_{23}$ $[\mu]$','$\sigma_{13}$ $[\mu]$','$\sigma_{12}$ $[\mu]$']
for k in range(0,len(fields)):
    # plots with centered images
    fig1 = plt.figure()
    ax1=plt.subplot(1,1,1)
    for key in dataDict:
        data=dataDict[key]
        ax1.plot(data[:,0], data[:,k+1],label=key)
    ax1.grid()
    ax1.legend()
    plt.xlabel('dipole distance $[b]$')
    plt.ylabel(fieldLabels[k])
    fig1.savefig(fields[k], bbox_inches='tight')
    
    # plots with uncentered images
#    fig1 = plt.figure()
#    ax1=plt.subplot(1,1,1)
#    for key in uncentered:
#        data=uncentered[key]
#        ax1.plot(data[:,0], data[:,k+1],label='  u'+key)
#    ax1.grid()
#    ax1.legend()
#    plt.xlabel('dipole distance $[b]$')
#    plt.ylabel(fieldLabels[k])
#    fig1.savefig(fields[k]+'_uncentered', bbox_inches='tight')
