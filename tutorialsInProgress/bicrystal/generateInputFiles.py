# /opt/local/bin/python3.12 test.py
import sys, string, os
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
sys.path.append("../../python")
from modlibUtils import *
sys.path.append("../../build/tools/pyMoDELib")
import pyMoDELib

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
setInputVariable('inputFiles/'+DDfile,'alphaLineTension','0.1') # dimensionless scale factor in for line tension forces
setInputVariable('inputFiles/'+DDfile,'Lmin','25')  # min segment length (in Burgers vector units)
setInputVariable('inputFiles/'+DDfile,'Lmax','150')  # max segment length (in Burgers vector units)
setInputVariable('inputFiles/'+DDfile,'outputFrequency','10')  # output frequency
setInputVariable('inputFiles/'+DDfile,'outputQuadraturePoints','0')  # output quadrature data
setInputVariable('inputFiles/'+DDfile,'glideSolverType','Galerkin')  # type of glide solver, or none
setInputVariable('inputFiles/'+DDfile,'climbSolverType','none')  # type of clim solver, or none
setInputVariable('inputFiles/'+DDfile,'Nsteps','10000')  # number of simulation steps
setInputVariable('inputFiles/'+DDfile,'crossSlipModel','0')  # crossSlipModel


# Make a local copy of material file, and modify that copy if necessary
materialFile='Al.txt';
materialFileTemplate='../../Library/Materials/'+materialFile;
print("\033[1;32mCreating  materialFile\033[0m")
shutil.copy2(materialFileTemplate,'inputFiles/'+materialFile)
setInputVariable('inputFiles/'+materialFile,'enabledSlipSystems','full')
b_SI=getValueInFile('inputFiles/'+materialFile,'b_SI')

# Make a local copy of ElasticDeformation file, and modify that copy if necessary
elasticDeformatinoFile='ElasticDeformation.txt';
elasticDeformatinoFileTemplate='../../Library/ElasticDeformation/'+elasticDeformatinoFile;
print("\033[1;32mCreating  elasticDeformatinoFile\033[0m")
shutil.copy2(elasticDeformatinoFileTemplate,'inputFiles/'+elasticDeformatinoFile)
setInputVector('inputFiles/'+elasticDeformatinoFile,'ExternalStress0',np.array([0.0,0.0,0.0,0.0,0.0,0.0]),'applied stress')

# Create polycrystal.txt using local material file
#meshFile='polycrystalCube_10_periodicX.msh';
#meshFile='polycrystalCube_50_periodicX.msh';
meshFile='bicrystal_48.msh';
meshFileTemplate='../../Library/Meshes/'+meshFile;
print("\033[1;32mCreating  polycrystalFile\033[0m")
shutil.copy2(meshFileTemplate,'inputFiles/'+meshFile)
pf=PolyCrystalFile(materialFile);
pf.absoluteTemperature=300;
pf.meshFile=meshFile;
pf.boxScaling=np.array([1e-7,1e-7,1e-7])/b_SI # length of box edges in Burgers vector units
#pf.boxScaling=np.array([2000,2000,2000]) # length of box edges in Burgers vector units
pf.X0=np.array([0,0,0]) # Centering unitCube mesh. Mesh nodes X are mapped to x=F*(X-X0)
pf.periodicFaceIDs=np.array([]) # no periodicity
pf.write('inputFiles')


N = 50
C2G = []
for i in range(N):
#    R = random_output(i)
    C2G.append(np.identity(3))
with open('inputFiles/polycrystal.txt', 'a') as fID:
    for i in range(N):
        fID.write(f'C2G{i + 1} =')
        fID.write(' '.join([f'{val:.15f}' for val in C2G[i][0]]) + '\n')
        fID.write(' '.join([f'{val:.15f}' for val in C2G[i][1]]) + '\n')
        fID.write(' '.join([f'{val:.15f}' for val in C2G[i][2]]) + ';\n\n')

# make a local copy of microstructure file, and modify that copy if necessary
simulationDir=os.path.abspath(".")
ddBase=pyMoDELib.DislocationDynamicsBase(simulationDir)
microstructureGenerator=pyMoDELib.MicrostructureGenerator(ddBase)

spec=pyMoDELib.ShearLoopIndividualSpecification()
spec.slipSystemIDs=np.array([0])
spec.loopRadii=np.array([1e-8])
spec.loopCenters=np.array([100.0,0.0,0.0])
spec.loopSides=np.array([10])
#microstructureGenerator.addShearLoopIndividual(spec)

microstructureGenerator.writeConfigFiles(0) # write evel_0.txt (optional)


defectiveCrystal=pyMoDELib.DefectiveCrystal(ddBase)
defectiveCrystal.initializeConfiguration(microstructureGenerator.configIO)
#defectiveCrystal.runSteps()

# Exctract displacement and stress fields on the xy-plane through the center of the mesh
mesh=ddBase.mesh
xMax=mesh.xMax(); # vector [max(x1) max(x2) max(x3)] in the mesh
xMin=mesh.xMin(); # vector [min(x1) min(x2) min(x3)] in the mesh

n=200
x=np.linspace(xMin[0], xMax[0], num=n) # grid x-range
y=np.linspace(xMin[1], xMax[1], num=n) # grid x-range
z=0.5*(xMin[2]+xMax[2])
points = np.zeros((0, 3))

for i in range(0,x.size):
    for j in range(0,y.size):
        pnt=np.array([x[i],y[j],z])
        points=np.vstack((points, pnt))

u=defectiveCrystal.displacement(points)*b_SI # displacement field at points

# Place u1 on a grid
u1=np.empty([n, n])
k=0
for i in range(0,x.size):
    for j in range(0,y.size):
        u1[i,j]=u[k,0]
        k=k+1

# Plot u1
fig=plt.figure()
plt.imshow(u1,origin='lower',cmap='jet')
plt.colorbar()
plt.show()

