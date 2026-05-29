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
setInputVariable('inputFiles/'+DDfile,'useInclusions','1')
setInputVariable('inputFiles/'+DDfile,'useElasticDeformation','1')
setInputVariable('inputFiles/'+DDfile,'useClusterDynamics','0')
setInputVariable('inputFiles/'+DDfile,'timeSteppingMethod','adaptive') # adaptive or fixed
setInputVariable('inputFiles/'+DDfile,'dtMax','1e25')
setInputVariable('inputFiles/'+DDfile,'dxMax','1') # max nodal displacement for when timeSteppingMethod=adaptive
setInputVariable('inputFiles/'+DDfile,'use_velocityFilter','0') # don't filter velocity if noise is enabled
setInputVariable('inputFiles/'+DDfile,'use_stochasticForce','0') # Langevin thermal noise enabled
setInputVariable('inputFiles/'+DDfile,'alphaLineTension','0.1') # dimensionless scale factor in for line tension forces
setInputVariable('inputFiles/'+DDfile,'Lmin','5')  # min segment length (in Burgers vector units)
setInputVariable('inputFiles/'+DDfile,'Lmax','20')  # max segment length (in Burgers vector units)
setInputVariable('inputFiles/'+DDfile,'outputFrequency','10')  # output frequency
setInputVariable('inputFiles/'+DDfile,'outputQuadraturePoints','1')  # output quadrature data
setInputVariable('inputFiles/'+DDfile,'computeElasticEnergyPerLength','1')  # output quadrature data
setInputVariable('inputFiles/'+DDfile,'glideSolverType','Galerkin')  # type of glide solver, or none
setInputVariable('inputFiles/'+DDfile,'climbSolverType','none')  # type of clim solver, or none
setInputVariable('inputFiles/'+DDfile,'Nsteps','1000')  # number of simulation steps



# Make a local copy of material file, and modify that copy if necessary
materialFile='Ni.txt';
materialFileTemplate='../../Library/Materials/'+materialFile;
print("\033[1;32mCreating  materialFile\033[0m")
shutil.copy2(materialFileTemplate,'inputFiles/'+materialFile)
setInputVariable('inputFiles/'+materialFile,'enabledSlipSystems','Shockley')
b_SI=getValueInFile('inputFiles/'+materialFile,'b_SI')

# Make a local copy of ElasticDeformation file, and modify that copy if necessary
elasticDeformatinoFile='ElasticDeformation.txt';
elasticDeformatinoFileTemplate='../../Library/ElasticDeformation/'+elasticDeformatinoFile;
print("\033[1;32mCreating  elasticDeformatinoFile\033[0m")
shutil.copy2(elasticDeformatinoFileTemplate,'inputFiles/'+elasticDeformatinoFile)
setInputVector('inputFiles/'+elasticDeformatinoFile,'ExternalStress0',np.array([0.0,0.0,0.0,0.0,0.0,0.01]),'applied stress')

# Create polycrystal.txt using local material file
meshFile='unitCube24.msh';
meshFileTemplate='../../Library/Meshes/'+meshFile;
print("\033[1;32mCreating  polycrystalFile\033[0m")
shutil.copy2(meshFileTemplate,'inputFiles/'+meshFile)
pf=PolyCrystalFile(materialFile);
pf.absoluteTemperature=300;
pf.meshFile=meshFile
pf.grain1globalX1=np.array([0,1,1])     # global x1 axis. Overwritten if alignToSlipSystem0=true
pf.grain1globalX3=np.array([-1,1,-1])    # global x3 axis. Overwritten if alignToSlipSystem0=true
pf.boxEdges=np.array([[0,1,1],[2,1,-1],[-1,1,-1]]) # i-throw is the direction of i-th box edge
#pf.boxScaling=np.array([200,200,100]) # length of box edges in Burgers vector units
pf.boxScaling=np.array([50e-9,50e-9,50e-9])/b_SI # length of box edges in Burgers vector units
pf.X0=np.array([0,0,0]) # Centering unitCube mesh. Mesh nodes X are mapped to x=F*(X-X0)
pf.periodicFaceIDs=np.array([-1])
pf.write('inputFiles')

# make a local copy of microstructure files, and modify those copies as necessary
microstructureFile1='periodicDipoleIndividual.txt';
microstructureFileTemplate='../../Library/Microstructures/'+microstructureFile1;
print("\033[1;32mCreating  microstructureFile\033[0m")
shutil.copy2(microstructureFileTemplate,'inputFiles/'+microstructureFile1) # target filename is /dst/dir/file.ext
setInputVector('inputFiles/'+microstructureFile1,'slipSystemIDs',np.array([0,1]),'slip system IDs for each dipole')
setInputVector('inputFiles/'+microstructureFile1,'exitFaceIDs',np.array([4,4]),'4 is for edge, 2 for screw')
setInputMatrix('inputFiles/'+microstructureFile1,'dipoleCenters',np.array([[0.0,0.0,0.0],[0.0,0.0,0.0]]))
setInputVector('inputFiles/'+microstructureFile1,'nodesPerLine',np.array([10,10]),'number of extra nodes on each dipole')
setInputVector('inputFiles/'+microstructureFile1,'dipoleHeights',np.array([200,200]),'height of each dipole, in number of planes')
setInputVector('inputFiles/'+microstructureFile1,'glideSteps',np.array([10.0,30.0]),'step of each dipole in the glide plane')

microstructureFile2='polyhedronInclusionIndividual.txt';
microstructureFileTemplate='../../Library/Microstructures/'+microstructureFile2;
print("\033[1;32mCreating  microstructureFile\033[0m")
shutil.copy2(microstructureFileTemplate,'inputFiles/'+microstructureFile2) # target filename is /dst/dir/file.ext
setInputMatrix('inputFiles/'+microstructureFile2,'F',np.array([[50.0,0.0,0.0],[0.0,50.0,0.0],[0.0,0.0,30.0]]))
setInputVariable('inputFiles/'+microstructureFile2,'phaseID','0') # ID of the phase of the precipitate
setInputVector('inputFiles/'+microstructureFile2,'eigenDistortion',np.array([0.03,0,0,0,0.03,0,0,0,0.03]),'eigenDistortion')
setInputVariable('inputFiles/'+microstructureFile2,'mshFile','../../../Library/Meshes/cylinder2D.msh') # ID of the phase of the precipitate
#setInputVariable('inputFiles/'+microstructureFile2,'mshFile','../../../Library/Meshes/unitCube8_2D.msh.msh') # ID of the phase of the precipitate

print("\033[1;32mCreating  initialMicrostructureFile\033[0m")
with open('inputFiles/initialMicrostructure.txt', "w") as initialMicrostructureFile:
    initialMicrostructureFile.write('microstructureFile='+microstructureFile1+';\n')
    initialMicrostructureFile.write('microstructureFile='+microstructureFile2+';\n')
