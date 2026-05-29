import os
import sys
#import vtk
import glob
import string
import numpy as np
import matplotlib.pyplot as plt

#plt.rcParams["text.usetex"] = True
import matplotlib as mpl
import matplotlib.transforms as mtransforms

sys.path.append("../../build/tools/pyMoDELib/")
import pyMoDELib
from scipy.signal import correlate2d
from scipy.stats import norm
from matplotlib import cm
from matplotlib.colors import Normalize

sys.path.append("../../python")
from modlibUtils import *
from pathlib import Path

# Define all file templates (using pathlib for cross-platform paths)
file_templates = {
    "dd_file": Path("../../Library/DislocationDynamics/DD.txt"),
    "noise_file": Path("../../Library/GlidePlaneNoise/MDStackingFault.txt"),
    "material_file": Path("../../Library/Materials/AlMg5.txt"),
    "elastic_deformation_file": Path(
        "../../Library/ElasticDeformation/ElasticDeformation.txt"
    ),
    "mesh": Path("../../Library/Meshes/unitCube24.msh"),
    "microstructure": Path(
        "../../Library/Microstructures/periodicDipoleIndividual.txt"
    ),
}

MDStackingFault_parameters = {
    "variables": {
        # Scalar parameters (use setInputVariable)
        "type": "MDStackingFaultNoise",
        "tag": 1,
        "seed": 2,
        "outputNoise": 1,
        "noiseFile": "AlMg5_ISF.vtk",
        # File paths (convert to absolute paths)
        "correlationFile": Path(
            "../../Library/GlidePlaneNoise/AlMg5_Cx_R100_ISF.vtk"
        ).resolve(),
    },
    "vectors": {
        # Vector parameters with comments (use setInputVector)
        "gridSize": {
            "value": np.array([60, 60, 1]),
            "comment": "number of grid points in each direction",
        },
        "gridSpacing_SI": {
            "value": np.array([1.0e-10, 1.0e-10, 1e-10]),
            "comment": "grid spacing in each direction",
        },
    },
    "copy_to": "inputFiles/MDStackingFault.txt",
}

Material_parameters = {
    "variables": {
        "enabledSlipSystems": "Shockley",
        #'glidePlaneNoise': ['MDSolidSolution.txt', 'MDStackingFault.txt'],
        "glidePlaneNoise": "MDStackingFault.txt",
        "atomsPerUnitCell": "1",
        "dislocationMobilityType": "default",
    },
    "copy_to": "inputFiles/AlMg5.txt",
}

Polycrystal_parameters = {
    "parameters": {
        "absoluteTemperature": 1,
        "grain1globalX1": np.array([0, 1, 1]),
        "grain1globalX3": np.array([-1, 1, -1]),
        "boxEdges": np.array([[0, 1, 1], [-2, -1, 1], [-1, 1, -1]]),
        "boxScaling": np.array([200, 200, 200]),
        #'X0': np.array([0.5, 0.5, 0.5]),
        "X0": np.array([0.0, 0.0, 0.0]),
        "periodicFaceIDs": np.array([0, 1, 2, 3, 4, 5]),
        "gridSpacing_SI": np.array([1.0e-10, 1.0e-10]),
    }
}


def processInitialConfigurations() -> None:
    # Copy all template files
    print("\033[1;32mCopying template files...\033[0m")
    for key, src_path in file_templates.items():
        dest = f"inputFiles/{src_path.name}"
        shutil.copy2(src_path.resolve(), dest)
        # shutil.copy2(src_path, dest)
        print(f"Created {dest}")

    # print("\033[1;32mCreating ddFile\033[0m")
    ## Apply all DD.txt parameters from the dictionary
    # for param, value in DD_parameters['variables'].items():
    #    setInputVariable(DD_parameters['copy_to'], param, str(value))

    print("\033[1;32mCreating  noiseFile\033[0m")
    # copy noise file into inputFiles
    shutil.copy2(
        MDStackingFault_parameters["variables"]["correlationFile"], "inputFiles/"
    )
    # copy declared variables to the configuration txt file
    for param, value in MDStackingFault_parameters["variables"].items():
        if "correlationFile" in param:
            relativePath = f"{str(value).split('/')[-1]}"
            setInputVariable(MDStackingFault_parameters["copy_to"], param, relativePath)
        else:
            setInputVariable(MDStackingFault_parameters["copy_to"], param, str(value))
    for param, data in MDStackingFault_parameters["vectors"].items():
        setInputVector(
            MDStackingFault_parameters["copy_to"], param, data["value"], data["comment"]
        )

    # Process material file
    print("\033[1;32mCreating materialFile...\033[0m")
    # print(Material_parameters['variables'].items())
    matParams = Material_parameters["variables"]
    setInputVariable(
        Material_parameters["copy_to"],
        "enabledSlipSystems",
        matParams["enabledSlipSystems"],
    )

    # for gPlaneNoise in matParams['glidePlaneNoise']:
    #    setInputVariable(Material_parameters['copy_to'],'glidePlaneNoise',gPlaneNoise)
    setInputVariable(
        Material_parameters["copy_to"], "glidePlaneNoise", matParams["glidePlaneNoise"]
    )

    setInputVariable(
        Material_parameters["copy_to"],
        "atomsPerUnitCell",
        matParams["atomsPerUnitCell"],
    )
    setInputVariable(
        Material_parameters["copy_to"],
        "dislocationMobilityType",
        matParams["dislocationMobilityType"],
    )

    # Process polycrystal
    # pf = PolyCrystalFile('FeCrAl_Fe.txt')
    pf = PolyCrystalFile(Material_parameters["copy_to"].split("/")[-1])
    # pf.meshFile='unitCube24.msh'
    pf.meshFile = f'{str(file_templates["mesh"]).split("/")[-1]}'
    for param, value in Polycrystal_parameters["parameters"].items():
        setattr(pf, param, value)
    pf.write("inputFiles")


def constructDotGrid(NX, NY, unitVec1, unitVec2):
    # Construct non-orthogonal dot grid based on the non-minimized atom structure
    firstNearNeighborDist = 2.8
    dotGrid = np.zeros((NX * NY, 2))  # each row contains x and y position of the dot
    vec1 = unitVec1 * firstNearNeighborDist
    vec2 = unitVec2 * firstNearNeighborDist
    for yIdx in range(NY):
        for xIdx in range(NX):
            dotGrid[yIdx * NX + xIdx] = vec1 * xIdx + vec2 * yIdx
    return dotGrid


def initUnitVectors(originalCorrVTK: str):
    latVector1 = None
    latVector2 = None
    copyFlag1 = False
    copyFlag2 = False
    with open(originalCorrVTK) as o:
        for line in o:
            if copyFlag1:
                latVector1 = np.array(line.strip().split(), dtype=float)
                copyFlag1 = False
            elif copyFlag2:
                latVector2 = np.array(line.strip().split(), dtype=float)
                copyFlag2 = False
            elif "VECTORS lattice_basis1 double" in line:
                copyFlag1 = True
            elif "VECTORS lattice_basis2 double" in line:
                copyFlag2 = True

    unitVec1 = latVector1 / np.linalg.norm(latVector1)
    unitVec2 = latVector2 / np.linalg.norm(latVector2)
    unitVec1 = unitVec1[:2]  # drop the z axis value
    unitVec2 = unitVec2[:2]
    return unitVec1, unitVec2


def circularShift(corrArray: np.ndarray, gridSize: np.ndarray) -> np.ndarray:
    NX, NY = np.squeeze(gridSize[0]), np.squeeze(gridSize[1])
    shiftY = NY / 2
    shiftX = NX / 2
    tempArr = np.zeros(NX * NY)
    for y in np.arange(NY):
        newY = int((y + shiftY) % NY)
        for x in np.arange(NX):
            newX = int((x + shiftX) % NX)
            tempArr[newY * NX + newX] += corrArray[y * NX + x]
    return tempArr


def main() -> int:
    # Preparing input files
    folders = ["inputFiles"]
    for x in folders:
        # remove existing data
        if os.path.exists(x):
            shutil.rmtree(x)
        # create necessary folder structure for the simulation
        os.makedirs(x)

    # set simulation parameters in inputFiles
    processInitialConfigurations()

    # setInputVariable(MDStackingFault_parameters['copy_to'], 'testNoiseSampling', str(1))
    simulationDir = os.path.abspath(".")
    ddBase = pyMoDELib.DislocationDynamicsBase(simulationDir)

    matFile = "inputFiles/AlMg5.txt"
    absoluteTemp = 1
    mat = pyMoDELib.PolycrystallineMaterialBase(matFile, absoluteTemp)
    b_SI = getValueInFile(f"{matFile}", "b_SI")
    mu0_SI = getValueInFile(f"{matFile}", "mu0_SI")
    rho_SI = getValueInFile(f"{matFile}", "rho_SI")

    tag = "0"
    originalCorrFile = "inputFiles/AlMg5_Cx_R100_ISF.vtk"
    seed = 1
    gridSize = np.array([100, 100, 1])
    gridSpacing = np.array([1e-10, 1e-10, 1e-10]) / b_SI
    latticeBasis = np.array([[1, 1], [1, 1]]).reshape(2, 2)

    sfNoise = pyMoDELib.MDStackingFaultNoise(
        mat,
        tag,
        originalCorrFile,
        seed,
        gridSize.reshape(3, 1),
        gridSpacing.reshape(3, 1),
        latticeBasis,
    )

    realizationNum = 100
    sampledNoises = np.array(sfNoise.sampleAverageNoise(realizationNum))
    std = np.std(sampledNoises)
    print(f"std = {std} J/m2")

    fig, axs = plt.subplots(1, 1, figsize=(8,6), dpi=200)
    # (returns a veiw if possible), flattens the array for histogram
    sampledNoises = sampledNoises.ravel()
    axs.hist(sampledNoises, density=True, alpha=0.3, bins='auto', label=f"R{realizationNum}")
    axs.grid(True)
    axs.set_title(f"sampled SF noise")
    axs.set_xlabel(f"")
    axs.set_ylabel(f"")
    axs.legend()
    fig.savefig(f'sampledNoiseDistribution_R{realizationNum}.png')

    return 0


if __name__ == "__main__":
    main()
