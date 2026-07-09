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
setInputVariable('inputFiles/'+DDfile,'useElasticDeformationFEM','1')
setInputVariable('inputFiles/'+DDfile,'useDislocations','1')
setInputVariable('inputFiles/'+DDfile,'useInclusions','0')
setInputVariable('inputFiles/'+DDfile,'useClusterDynamics','1')
setInputVariable('inputFiles/'+DDfile,'useClusterDynamicsFEM','1')
setInputVariable('inputFiles/'+DDfile,'Nsteps','120')  # number of simulation steps
setInputVariable('inputFiles/'+DDfile,'timeSteppingMethod','fixed') # adaptive or fixed
setInputVariable('inputFiles/'+DDfile,'dtMax','3.4e19')
setInputVariable('inputFiles/'+DDfile,'dxMax','5') # max nodal displacement for when timeSteppingMethod=adaptive
setInputVariable('inputFiles/'+DDfile,'use_velocityFilter','0') # don't filter velocity if noise is enabled
setInputVariable('inputFiles/'+DDfile,'use_stochasticForce','0') # Langevin thermal noise enabled
setInputVariable('inputFiles/'+DDfile,'useSubCycling','0') # Langevin thermal noise enabled
setInputVariable('inputFiles/'+DDfile,'glideSolverType','none')  # type of glide solver, or none
setInputVariable('inputFiles/'+DDfile,'climbSolverType','Galerkin')  # type of clim solver, or none
setInputVariable('inputFiles/'+DDfile,'quadPerLength','0.0001')  # quadrature points per unit length (in Burgers vector)
setInputVariable('inputFiles/'+DDfile,'alphaLineTension','1.0') # dimensionless scale factor in for line tension forces
setInputVariable('inputFiles/'+DDfile,'remeshFrequency','1')  # min segment length (in Burgers vector units)
setInputVariable('inputFiles/'+DDfile,'Lmin','450')  # min segment length (in Burgers vector units)
setInputVariable('inputFiles/'+DDfile,'Lmax','500')  # max segment length (in Burgers vector units)
setInputVariable('inputFiles/'+DDfile,'outputFrequency','1')  # output frequency
setInputVariable('inputFiles/'+DDfile,'outputQuadraturePoints','1')  # output quadrature data



# Make a local copy of material file, and modify that copy if necessary
materialFile='Zr4_Fitted.txt';
materialFileTemplate='../../Library/Materials/'+materialFile;
print("\033[1;32mCreating  materialFile\033[0m")
shutil.copy2(materialFileTemplate,'inputFiles/'+materialFile)
setInputVariable('inputFiles/'+materialFile,'enabledSlipSystems','fullBasal fullPrismatic')
setInputMatrix('inputFiles/'+materialFile,'reactionPrefactorMap',np.array([[-1, -1, 0.0],
                                                                           [-1, 1, 2.856],
                                                                           [-1, 2, 1.0],
                                                                           [-1, 3, 1.0],
                                                                           [1, 1, 1.0],
                                                                           [1, 2, 1.0],
                                                                           [1, 3, 0.0],
                                                                           [2, 2, 0.0],
                                                                           [2, 3, 0.0],
                                                                           [3, 3, 0.0],])) # Reaction map between species (species,species,reaction prefactor)
setInputVector('inputFiles/'+materialFile,'initLoopSinks_SI',np.array([4e21,4e21,4e21,4e21,2e-5,12e-5,12e-5,12e-5]),'') # [m^-3] and [-] for initial loop sinks
setInputVector('inputFiles/'+materialFile,'loopClusterFraction',np.array([0,0,0,0]),'') # [-] route lost i+3i/2i+2i content into loop families
setInputVariable('inputFiles/'+materialFile,'loopClusterNucDefects','4') # [-] loop-loop Avrami capture factor
setInputVector('inputFiles/'+materialFile,'loopCascadeFraction',np.array([0,0,0,0]),'')  # [-] direct cascade loop seed per [c,a1,a2,a3]
setInputVector('inputFiles/'+materialFile,'loopNucDefects',np.array([8,4,4,4]),'') # [-] defects per nucleated cascade loop per [c,a1,a2,a3]
setInputVector('inputFiles/'+materialFile,'loopContentTau_SI',np.array([1e30,1e30,1e30,1e30]),'') # [s] loop-content removal lifetime per [c,a1,a2,a3]; large = off
setInputVector('inputFiles/'+materialFile,'loopAnnealTau0_SI',np.array([1e30,1e30,1e30,1e30]),'')  # [s] Arrhenius prefactor for loop annealing; large = off
setInputVector('inputFiles/'+materialFile,'loopAnnealEa_eV',np.array([0,0,0,0]),'') # [eV] activation energy for loop annealing per [c,a1,a2,a3]
setInputVector('inputFiles/'+materialFile,'loopCoalLL',np.array([0,0,0,0]),'') # [-] loop-loop coalescence rate scale
setInputVector('inputFiles/'+materialFile,'loopCoalLN',np.array([0,0,0,0]),'') # [-] loop-network coalescence rate scale
setInputVariable('inputFiles/'+materialFile,'loopCoalKappaLL','1') # [-] loop-loop Avrami capture factor
setInputVariable('inputFiles/'+materialFile,'loopCoalKappaLN','1') # [-] loop-network Avrami capture factor
setInputVariable('inputFiles/'+materialFile,'loopCoalNetwork_SI','0.0') # [m^-2] network density used by loop-network coalescence 



# Make a local copy of ElasticDeformation file, and modify that copy if necessary
elasticDeformatinoFile='ElasticDeformation.txt';
elasticDeformatinoFileTemplate='../../Library/ElasticDeformation/'+elasticDeformatinoFile;
print("\033[1;32mCreating  elasticDeformatinoFile\033[0m")
shutil.copy2(elasticDeformatinoFileTemplate,'inputFiles/'+elasticDeformatinoFile)
setInputVector('inputFiles/'+elasticDeformatinoFile,'ExternalStress0',np.array([0.0,0.0,0.0,0.0,0.0,0.0]),'applied stress')



# Create polycrystal.txt using local material file
meshFile='poly30_100K.msh';
meshFileTemplate='../../Library/Meshes/'+meshFile;
print("\033[1;32mCreating  polycrystalFile\033[0m")
shutil.copy2(meshFileTemplate,'inputFiles/'+meshFile)
pf=PolyCrystalFile(materialFile);
pf.absoluteTemperature=553;
pf.meshFile=meshFile
pf.singleCrystal = False
pf.numberGrains = 30
pf.f_param = [0.6, 0.2] # index 0 => % of grains: [001]-axis || global x, index 1 => % of grains: [001]-axis || global y, diff(0,1) =>% of grains: [001] || global z
pf.grain1globalX1=np.array([1,0,0])     # global x1 axis. Overwritten if alignToSlipSystem0=true (only for single crystal)
pf.grain1globalX3=np.array([0,0,1])    # global x3 axis. Overwritten if alignToSlipSystem0=true (only for single crystal)
pf.boxEdges=np.array([[1,0,0],[0,1,0],[0,0,1]]) # i-throw is the direction of i-th box edge (only for single crystal)
pf.boxScaling=np.array([3093,3093,3093]) # must be a vector of integers
pf.X0=np.array([0,0,0]) # Centering unitCube mesh. Mesh nodes X are mapped to x=F*(X-X0)
pf.periodicFaceIDs=np.array([])
pf.write('inputFiles')



# make a local copy of microstructure file, and modify that copy if necessary
microstructureFile1='frankLoopsDensity.txt';
microstructureFileTemplate1='../../Library/Microstructures/'+microstructureFile1;
print("\033[1;32mCreating  microstructureFile\033[0m")
shutil.copy2(microstructureFileTemplate1,'inputFiles/'+microstructureFile1) # target filename is /dst/dir/file.ext
setInputVariable('inputFiles/'+microstructureFile1,'targetDensity','0e11')
setInputVector('inputFiles/'+microstructureFile1,'planeIDs',np.array([0]),'')
setInputVariable('inputFiles/'+microstructureFile1,'radiusDistributionMean','30e-9') # [m] mean of loop radii
setInputVariable('inputFiles/'+microstructureFile1,'radiusDistributionStd','0e-8') # [m] std of loop radii
setInputVariable('inputFiles/'+microstructureFile1,'numberOfSides','20') # [-] number of sides in each loop polygon
setInputVariable('inputFiles/'+microstructureFile1,'burgersFactor','0.5') # [-] scaling of burgers vector
setInputVariable('inputFiles/'+microstructureFile1,'areVacancyLoops','1') # 1=vacancy-type loops, 0=interstiatial-type loops



print("\033[1;32mCreating  initialMicrostructureFile\033[0m")
with open('inputFiles/initialMicrostructure.txt', "w") as initialMicrostructureFile:
    # initialMicrostructureFile.write('microstructureFile='+microstructureFile1+';\n')
    initialMicrostructureFile.write('')