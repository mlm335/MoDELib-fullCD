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
setInputVariable('inputFiles/'+DDfile,'useFEM','0')
setInputVariable('inputFiles/'+DDfile,'useDislocations','1')
setInputVariable('inputFiles/'+DDfile,'useInclusions','0')
setInputVariable('inputFiles/'+DDfile,'useElasticDeformation','1')
setInputVariable('inputFiles/'+DDfile,'useClusterDynamics','0')
setInputVariable('inputFiles/'+DDfile,'timeSteppingMethod','adaptive') # adaptive or fixed
setInputVariable('inputFiles/'+DDfile,'dtMax','1e25')
setInputVariable('inputFiles/'+DDfile,'dxMax','1') # max nodal displacement for when timeSteppingMethod=adaptive
setInputVariable('inputFiles/'+DDfile,'use_velocityFilter','0') # don't filter velocity if noise is enabled
setInputVariable('inputFiles/'+DDfile,'use_stochasticForce','0') # Langevin thermal noise enabled
setInputVariable('inputFiles/'+DDfile,'alphaLineTension','1')
setInputVariable('inputFiles/'+DDfile,'Lmin','5')  # min segment length (in Burgers vector units)
setInputVariable('inputFiles/'+DDfile,'Lmax','20')  # max segment length (in Burgers vector units)
setInputVariable('inputFiles/'+DDfile,'outputFrequency','100')  # output frequency
setInputVariable('inputFiles/'+DDfile,'outputQuadraturePoints','1')  # output quadrature data
setInputVariable('inputFiles/'+DDfile,'computeElasticEnergyPerLength','1')  # output quadrature data
setInputVariable('inputFiles/'+DDfile,'glideSolverType','Galerkin')  # type of glide solver, or none
setInputVariable('inputFiles/'+DDfile,'climbSolverType','none')  # type of clim solver, or none
setInputVariable('inputFiles/'+DDfile,'Nsteps','10000')  # number of simulation steps
setInputVector('inputFiles/'+DDfile,'nodalVelocityConstraints',np.array([0.0,1.0,0.0]),'Velocty constrained in the x direction')
setInputVariable('inputFiles/'+DDfile,'remeshFrequency','0')
setInputVariable('inputFiles/'+DDfile,'maxJunctionIterations','0')
setInputVariable('inputFiles/'+DDfile,'quadPerLength','1.0')
setInputVector('inputFiles/'+DDfile,'periodicImageSize',np.array([2,2,2]),'n of images in each direction')
setInputVariable('inputFiles/'+DDfile,'EwaldLengthFactor','1.0')

gridSize=1024;
gridSpacing=1e-10;
boxSize=gridSize*gridSpacing;

# Make a local copy of noise file, and modify that copy if necessary
noiseFile='AnalyticalSolidSolutionNoise.txt'
noiseFileTemplate='../../Library/GlidePlaneNoise/'+noiseFile;
print("\033[1;32mCreating  noiseFile\033[0m")
shutil.copy2(noiseFileTemplate,'inputFiles/'+noiseFile) # target filename is /dst/dir/file.ext
#setInputVariable(noiseFile,'seed','1')
setInputVariable('inputFiles/'+noiseFile,'MSSS_SI','0.45e18')
setInputVariable('inputFiles/'+noiseFile,'seed','6')
setInputVector('inputFiles/'+noiseFile,'gridSize',np.array([gridSize,gridSize,64]),'noise grid size')


#noiseFile='MDSolidSolution.txt';
#noiseFileTemplate='../../Library/GlidePlaneNoise/'+noiseFile;
#print("\033[1;32mCreating  noiseFile\033[0m")
#shutil.copy2(noiseFileTemplate,'inputFiles/'+noiseFile) # target filename is /dst/dir/file.ext
#setInputVariable('inputFiles/'+noiseFile,'type','MDSolidSolutionNoise')
#setInputVariable('inputFiles/'+noiseFile,'tag','0')
#setInputVariable('inputFiles/'+noiseFile,'seed','0')
#setInputVariable('inputFiles/'+noiseFile,'correlationFile_L',os.path.abspath('../../Library/GlidePlaneNoise/MDSolidSolutionCorrelations_L_MoNbTi.vtk'))
#setInputVariable('inputFiles/'+noiseFile,'correlationFile_T',os.path.abspath('../../Library/GlidePlaneNoise/MDSolidSolutionCorrelations_T_MoNbTi.vtk'))
#setInputVector('inputFiles/'+noiseFile,'gridSize',np.array([100,100,1]),'number of grid points in each direction')
#setInputVector('inputFiles/'+noiseFile,'gridSpacing_SI',np.array([1.12e-10,1.06e-10,1e-10]),'grid spacing in each direction')
#setInputVariable('inputFiles/'+noiseFile,'a_cai_SI','0')

# Make a local copy of material file, and modify that copy if necessary
materialFile='AlMg15.txt';
materialFileTemplate='../../Library/Materials/'+materialFile;
print("\033[1;32mCreating  materialFile\033[0m")
shutil.copy2(materialFileTemplate,'inputFiles/'+materialFile)
#setInputVariable('inputFiles/'+materialFile,'enabledSlipSystems','Shockley')
setInputVariable('inputFiles/'+materialFile,'glidePlaneNoise',noiseFile)
b_SI=getValueInFile('inputFiles/'+materialFile,'b_SI')

# Make a local copy of ElasticDeformation file, and modify that copy if necessary
elasticDeformatinoFile='ElasticDeformation.txt';
elasticDeformatinoFileTemplate='../../Library/ElasticDeformation/'+elasticDeformatinoFile;
print("\033[1;32mCreating  elasticDeformatinoFile\033[0m")
shutil.copy2(elasticDeformatinoFileTemplate,'inputFiles/'+elasticDeformatinoFile)
#setInputVector('inputFiles/'+elasticDeformatinoFile,'ExternalStress0',np.array([0.0,0.0,0.0,0.0,0.0,0.0]),'applied stress')
setInputVector('inputFiles/'+elasticDeformatinoFile,'ExternalStrainRate',np.array([0.0,0.0,0.0,0.0,0.0,1e-8]),'applied strain rate')
setInputVector('inputFiles/'+elasticDeformatinoFile,'stiffnessRatio',np.array([0.0,0.0,0.0,0.0,0.0,1e20]),'stiffness ratio')

# Create polycrystal.txt using local material file
meshFile='unitCube24.msh';
meshFileTemplate='../../Library/Meshes/'+meshFile;
print("\033[1;32mCreating  polycrystalFile\033[0m")
shutil.copy2(meshFileTemplate,'inputFiles/'+meshFile)
pf=PolyCrystalFile(materialFile);
pf.absoluteTemperature=300;
pf.meshFile=meshFile
pf.alignToSlipSystem0=1
pf.boxScaling=np.array([boxSize,boxSize,boxSize])/b_SI # length of box edges in Burgers vector units
pf.periodicFaceIDs=np.array([-1])
pf.write('inputFiles')

# make a local copy of microstructure file, and modify that copy if necessary
microstructureFile1='periodicDipoleIndividual.txt';
microstructureFileTemplate='../../Library/Microstructures/'+microstructureFile1;
print("\033[1;32mCreating  microstructureFile\033[0m")
shutil.copy2(microstructureFileTemplate,'inputFiles/'+microstructureFile1) # target filename is /dst/dir/file.ext
setInputVector('inputFiles/'+microstructureFile1,'slipSystemIDs',np.array([0,-1]),'slip system IDs for each dipole')
setInputVector('inputFiles/'+microstructureFile1,'exitFaceIDs',np.array([4,4]),'4 is for edge, 2 for screw')
setInputMatrix('inputFiles/'+microstructureFile1,'dipoleCenters',np.array([[0.0,0.0,0.0],[0.0,0.0,0.0]]))
setInputVector('inputFiles/'+microstructureFile1,'nodesPerLine',np.array([0,10]),'number of extra nodes on each dipole')
setInputVector('inputFiles/'+microstructureFile1,'dipoleHeights',np.array([boxSize/b_SI/4/(np.sqrt(6)/3),boxSize/b_SI/4/(np.sqrt(6)/3)]),'height of each dipole, in number of planes')
setInputVector('inputFiles/'+microstructureFile1,'glideSteps',np.array([100.0,30.0]),'step of each dipole in the glide plane')

print("\033[1;32mCreating  initialMicrostructureFile\033[0m")
with open('inputFiles/initialMicrostructure.txt', "w") as initialMicrostructureFile:
    initialMicrostructureFile.write('microstructureFile='+microstructureFile1+';\n')
#    initialMicrostructureFile.write('microstructureFile='+microstructureFile2+';\n')
