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
setInputVariable('inputFiles/'+DDfile,'useFEM','1')
setInputVariable('inputFiles/'+DDfile,'useDislocations','1')
setInputVariable('inputFiles/'+DDfile,'useInclusions','0')
setInputVariable('inputFiles/'+DDfile,'useElasticDeformation','1')
setInputVariable('inputFiles/'+DDfile,'useClusterDynamics','1')
setInputVariable('inputFiles/'+DDfile,'useClusterDynamicsFEM','1')
setInputVariable('inputFiles/'+DDfile,'use_directSolver_FEM','0')
setInputVariable('inputFiles/'+DDfile,'EwaldLengthFactor','1')
setInputVariable('inputFiles/'+DDfile,'timeSteppingMethod','adaptive') # adaptive or fixed
setInputVariable('inputFiles/'+DDfile,'dtMax','1e25')
setInputVariable('inputFiles/'+DDfile,'dxMax','1') # max nodal displacement for when timeSteppingMethod=adaptive
setInputVariable('inputFiles/'+DDfile,'use_velocityFilter','0') # don't filter velocity if noise is enabled
setInputVariable('inputFiles/'+DDfile,'use_stochasticForce','0') # Langevin thermal noise enabled
setInputVariable('inputFiles/'+DDfile,'quadPerLength','0.1')
setInputVariable('inputFiles/'+DDfile,'alphaLineTension','1.0') # dimensionless scale factor in for line tension forces
setInputVariable('inputFiles/'+DDfile,'Lmin','100')  # min segment length (in Burgers vector units)
setInputVariable('inputFiles/'+DDfile,'Lmax','200')  # max segment length (in Burgers vector units)
setInputVariable('inputFiles/'+DDfile,'outputFrequency','1')  # output frequency
setInputVariable('inputFiles/'+DDfile,'outputQuadraturePoints','1')  # output quadrature data
setInputVariable('inputFiles/'+DDfile,'glideSolverType','none')  # type of glide solver, or none
setInputVariable('inputFiles/'+DDfile,'climbSolverType','Galerkin')  # type of clim solver, or none
setInputVariable('inputFiles/'+DDfile,'Nsteps','500')  # number of simulation steps
setInputVariable('inputFiles/'+DDfile,'outputLoopLength','1')  # crossSlipModel
setInputVariable('inputFiles/'+DDfile,'outputDislocationDensityPerSlipSystem','0')
setInputVariable('inputFiles/'+DDfile,'outputPlasticDistortionPerSlipSystem','0')
setInputVariable('inputFiles/'+DDfile,'remeshFrequency','1')
setInputVariable('inputFiles/'+DDfile,'absoluteAreaThreshold','100.0')


# Make a local copy of material file, and modify that copy if necessary
materialFile='Al.txt';
materialFileTemplate='../../Library/Materials/'+materialFile;
print("\033[1;32mCreating  materialFile\033[0m")
shutil.copy2(materialFileTemplate,'inputFiles/'+materialFile)
setInputVariable('inputFiles/'+materialFile,'enabledSlipSystems','Shockley')
setInputVariable('inputFiles/'+materialFile,'mobileSpeciesSurvivingEfficiency','0.01')
setInputVariable('inputFiles/'+materialFile,'discreteDislocationBias','1.0 1.0')
setInputVariable('inputFiles/'+materialFile,'ISF_SI','0.0e-3')
b_SI=getValueInFile('inputFiles/'+materialFile,'b_SI')



# Make a local copy of ElasticDeformation file, and modify that copy if necessary
elasticDeformatinoFile='ElasticDeformation.txt';
elasticDeformatinoFileTemplate='../../Library/ElasticDeformation/'+elasticDeformatinoFile;
print("\033[1;32mCreating  elasticDeformatinoFile\033[0m")
shutil.copy2(elasticDeformatinoFileTemplate,'inputFiles/'+elasticDeformatinoFile)
setInputVector('inputFiles/'+elasticDeformatinoFile,'ExternalStress0',np.array([0.0,0.0,0.0,0.0,0.0,0.0]),'applied stress of 10 MPa')
setInputVector('inputFiles/'+elasticDeformatinoFile,'stiffnessRatio',np.array([0.0,0.0,0.0,0.0,0.0,0.0]),'applied stiffness ratio')

# Create polycrystal.txt using local material file
meshFile='unitCube8188.msh';
meshFileTemplate='../../Library/Meshes/'+meshFile;
print("\033[1;32mCreating  polycrystalFile\033[0m")
shutil.copy2(meshFileTemplate,'inputFiles/'+meshFile)
pf=PolyCrystalFile(materialFile);
pf.absoluteTemperature=470;
pf.meshFile=meshFile
pf.grain1globalX1=np.array([1,0,0])     # global x1 axis. Overwritten if alignToSlipSystem0=true
pf.grain1globalX3=np.array([0,0,1])    # global x3 axis. Overwritten if alignToSlipSystem0=true
pf.boxEdges=np.array([[1,0,0],[0,1,0],[0,0,1]]) # i-throw is the direction of i-th box edge
#pf.boxScaling=np.array([200,200,100]) # length of box edges in Burgers vector units
pf.boxScaling=np.array([1.5e-6,1.5e-6,1.5e-6])/b_SI # length of box edges in Burgers vector units
pf.X0=np.array([0,0,0]) # Centering unitCube mesh. Mesh nodes X are mapped to x=F*(X-X0)
pf.periodicFaceIDs=np.array([])
pf.write('inputFiles')



# make a local copy of microstructure file, and modify that copy if necessary
#microstructureFile1='shearLoopsDensity.txt';
#microstructureFileTemplate='../../Library/Microstructures/'+microstructureFile1;
#print("\033[1;32mCreating  microstructureFile\033[0m")
#shutil.copy2(microstructureFileTemplate,'inputFiles/'+microstructureFile1) # target filename is /dst/dir/file.ext
#setInputVariable('inputFiles/'+microstructureFile1,'targetDensity_SI','1e11')
#setInputVariable('inputFiles/'+microstructureFile1,'radiusDistributionMean_SI','1e-05')
#setInputVariable('inputFiles/'+microstructureFile1,'radiusDistributionStd_SI','0.0')
#setInputVariable('inputFiles/'+microstructureFile1,'numberOfSides','4')
#setInputVector('inputFiles/'+microstructureFile1,'allowedGrainIDs',np.array([-1]),'set of grain IDs where loops are allowed. Use -1 for all grains')
#setInputVector('inputFiles/'+microstructureFile1,'allowedSlipSystemIDs',np.array([-1]),'set of slip system IDs whose Burgers vector are allowed to be the prism axis. Use -1 for all slip systems')

microstructureFile='frankLoopsIndividual.txt';
microstructureFileTemplate='../../Library/Microstructures/'+microstructureFile;

#microstructureFile1='prismaticLoopsIndividual.txt';
#print("\033[1;32mCreating  microstructureFile\033[0m")
#shutil.copy2(microstructureFileTemplate,'inputFiles/'+microstructureFile1) # target filename is /dst/dir/file.ext
#setInputVariable('inputFiles/'+microstructureFile1,'radiusDistributionMean_SI','1.0e-6')
#setInputVariable('inputFiles/'+microstructureFile1,'targetDensity_SI','1.0e13')

microstructureFile1='frankLoopsIndividual.txt';
print("\033[1;32mCreating  microstructureFile\033[0m")
shutil.copy2(microstructureFileTemplate,'inputFiles/'+microstructureFile1) # target filename is /dst/dir/file.ext
setInputVariable('inputFiles/'+microstructureFile1,'planeIDs','0')
setInputVariable('inputFiles/'+microstructureFile1,'loopRadii_SI','330e-10')
setInputVariable('inputFiles/'+microstructureFile1,'loopSides','12')
setInputVector('inputFiles/'+microstructureFile1,'loopCenters',np.array([0, 0, 0]), 'Loop Centers')
setInputVariable('inputFiles/'+microstructureFile1,'isVacancyLoop','1')

print("\033[1;32mCreating  initialMicrostructureFile\033[0m")
with open('inputFiles/initialMicrostructure.txt', "w") as initialMicrostructureFile:
    initialMicrostructureFile.write('microstructureFile='+microstructureFile1+';\n')
