import sys
sys.path.append("../../python/")
from modlibUtils import *

# Create folder structure
folders=['evl','F','inputFiles']
for x in folders:
    if not os.path.exists(x):
        os.makedirs(x)

# Make a local copy of DD parameters file and modify that copy if necessary
DDfile='DD.txt'
DDfileTemplate='../../Library/DislocationDynamics/'+DDfile
print("\033[1;32mCreating  DDfile\033[0m")
shutil.copy2(DDfileTemplate,'inputFiles/'+DDfile)
setInputVariable('inputFiles/'+DDfile,'Lmax','20')  # max segment length (in Burgers vector units)
setInputVariable('inputFiles/'+DDfile,'outputFrequency','1')  # output frequency
setInputVariable('inputFiles/'+DDfile,'outputQuadraturePoints','1')  # output quadrature data
setInputVariable('inputFiles/'+DDfile,'computeElasticEnergyPerLength','1')  # output quadrature data
setInputVariable('inputFiles/'+DDfile,'Nsteps','1')  # number of simulation steps
#setInputVector('inputFiles/'+DDfile,'nodalVelocityConstraints',np.array([0.0,1.0,0.0]),'Velocty constrained in the x direction')
setInputVariable('inputFiles/'+DDfile,'quadPerLength','1.0')
setInputVector('inputFiles/'+DDfile,'periodicImageSize',np.array([0,2,0]),'n of images in each direction')
setInputVariable('inputFiles/'+DDfile,'EwaldLengthFactor','1.0')
setInputVariable('inputFiles/'+DDfile,'maxJunctionIterations','0')
setInputVariable('inputFiles/'+DDfile,'remeshFrequency','0')


materialFile='W.txt';
materialFileTemplate='../../Library/Materials/'+materialFile;
shutil.copy2(materialFileTemplate,'inputFiles/'+materialFile)
b_SI=getValueInFile('inputFiles/'+materialFile,'b_SI')
h_SI=b_SI/np.sqrt(2.0/3.0) # plane heigh in BCC


# Make a local copy of ElasticDeformation file, and modify that copy if necessary
elasticDeformatinoFile='ElasticDeformation.txt';
elasticDeformatinoFileTemplate='../../Library/ElasticDeformation/'+elasticDeformatinoFile;
print("\033[1;32mCreating  elasticDeformatinoFile\033[0m")
shutil.copy2(elasticDeformatinoFileTemplate,'inputFiles/'+elasticDeformatinoFile)
#setInputVector('inputFiles/'+elasticDeformatinoFile,'ExternalStress0',np.array([0.0,0.0,0.0,0.0,0.0,0.0]),'applied stress')
setInputVector('inputFiles/'+elasticDeformatinoFile,'ExternalStressRate',np.array([0.0,0.0,0.0,0.0,0.0,0.0]),'applied strain rate')
setInputVector('inputFiles/'+elasticDeformatinoFile,'stiffnessRatio',np.array([0.0,0.0,0.0,0.0,0.0,0.0]),'stiffness ratio')

meshFile='unitCube24.msh';
meshFileTemplate='../../Library/Meshes/'+meshFile;
print("\033[1;32mCreating  polycrystalFile\033[0m")
shutil.copy2(meshFileTemplate,'inputFiles/'+meshFile)
pf=PolyCrystalFile(materialFile);
pf.absoluteTemperature=300;
pf.meshFile=meshFile
pf.alignToSlipSystem0=1
pf.boxScaling=np.array([300e-10,180e-10,400e-10])/pf.boxEdgesLatticeLengths/b_SI # must be a vector of integers
pf.X0=np.array([0,0,0]) # Centering unitCube mesh. Mesh nodes X are mapped to x=F*(X-X0)
pf.periodicFaceIDs=np.array([-1])
pf.write('inputFiles')


# make a local copy of microstructure file, and modify that copy if necessary
microstructureFile1='periodicDipoleIndividual.txt';
microstructureFileTemplate='../../Library/Microstructures/'+microstructureFile1;
print("\033[1;32mCreating  microstructureFile\033[0m")
shutil.copy2(microstructureFileTemplate,'inputFiles/'+microstructureFile1) # target filename is /dst/dir/file.ext
setInputVector('inputFiles/'+microstructureFile1,'slipSystemIDs',np.array([0]),'slip system IDs for each dipole')
setInputVector('inputFiles/'+microstructureFile1,'exitFaceIDs',np.array([1]),'4 is for edge, 2 for screw')
setInputVector('inputFiles/'+microstructureFile1,'dipoleCenters',np.array([0.0,0.0,0.0]),'dipole center')
setInputVector('inputFiles/'+microstructureFile1,'nodesPerLine',np.array([0]),'number of extra nodes on each dipole')
setInputVector('inputFiles/'+microstructureFile1,'dipoleHeights',np.array([np.round(200e-10/h_SI)]),'height of each dipole, in number of planes')
setInputVector('inputFiles/'+microstructureFile1,'glideSteps',np.array([1.0]),'step of each dipole in the glide plane')

## make a local copy of microstructure file, and modify that copy if necessary
#microstructureTemplate='periodicDipole.txt';
#shutil.copy2('../../../tutorials/DislocationDynamics/MicrostructureLibrary/'+microstructureTemplate, '.') # target filename is /dst/dir/file.ext
#setInputVector(microstructureTemplate,'periodicDipoleSlipSystemIDs',np.array([0]),'slip system IDs for each dipole')
#setInputVector(microstructureTemplate,'periodicDipoleExitFaceIDs',np.array([1]),'1 is for edge, 0 for screw')
#setInputVector(microstructureTemplate,'periodicDipolePoints',np.array([0.0,0.0,0.0]),'Center of the dipole')
#setInputVector(microstructureTemplate,'periodicDipoleNodes',np.array([0]),'number of extra nodes on each dipole')
#setInputVector(microstructureTemplate,'periodicDipoleHeights',np.array([np.round(200e-10/h_SI)]),'height of each dipole, in number of planes')
#setInputVector(microstructureTemplate,'periodicDipoleGlideSteps',np.array([1.0]),'step of each dipole in the glide plane')


print("\033[1;32mCreating  initialMicrostructureFile\033[0m")
with open('inputFiles/initialMicrostructure.txt', "w") as initialMicrostructureFile:
    initialMicrostructureFile.write('microstructureFile='+microstructureFile1+';\n')
#    initialMicrostructureFile.write('microstructureFile='+microstructureFile2+';\n')
