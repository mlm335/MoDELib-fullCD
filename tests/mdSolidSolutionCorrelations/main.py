import os
import sys
import vtk
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
    "noise_file": Path("../../Library/GlidePlaneNoise/MDSolidSolution.txt"),
    "material_file": Path("../../Library/Materials/FeCrAl_Fe.txt"),
    "elastic_deformation_file": Path(
        "../../Library/ElasticDeformation/ElasticDeformation.txt"
    ),
    "mesh": Path("../../Library/Meshes/unitCube24.msh"),
    "microstructure": Path(
        "../../Library/Microstructures/periodicDipoleIndividual.txt"
    ),
}

MDSolidSolution_parameters = {
    "variables": {
        # Scalar parameters (use setInputVariable)
        "type": "MDSolidSolutionNoise",
        "tag": 1,
        "seed": 2,
        "outputNoise": 1,
        "a_cai_SI": 1e-16,
        # File paths (convert to absolute paths)
        "correlationFile_xz": Path(
            "../../Library/GlidePlaneNoise/MoDELCompatible_FeCr8_xz.vtk"
        ).resolve(),
        "correlationFile_yz": Path(
            "../../Library/GlidePlaneNoise/MoDELCompatible_FeCr8_yz.vtk"
        ).resolve(),
    },
    "vectors": {
        # Vector parameters with comments (use setInputVector)
        "gridSize": {
            "value": np.array([100, 100, 1]),
            "comment": "number of grid points in each direction",
        },
        "gridSpacing_SI": {
            "value": np.array([1.0e-10, 1.0e-10, 1e-10]),
            "comment": "grid spacing in each direction",
        },
    },
    "copy_to": "inputFiles/MDSolidSolution.txt",
}

Material_parameters = {
    "variables": {
        "enabledSlipSystems": "Shockley",
        #'glidePlaneNoise': ['MDSolidSolution.txt', 'MDStackingFault.txt'],
        "glidePlaneNoise": "MDSolidSolution.txt",
        "atomsPerUnitCell": "1",
        "dislocationMobilityType": "default",
    },
    "copy_to": "inputFiles/FeCrAl_Fe.txt",
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


def readVTKnoise(fname: str):
    reader = vtk.vtkGenericDataObjectReader()
    reader.SetFileName(fname)  # declare the vtk filename
    reader.ReadAllVectorsOn()  # read all vector data
    reader.ReadAllScalarsOn()  # read all scalar data
    reader.Update()  # update to new file
    # read data dimensions
    dims = [0, 0, 0]
    reader.GetOutput().GetDimensions(dims)
    NX, NY, NZ = dims
    # NX, NY, NZ = reader.GetOutput().GetDimensions()
    # creates dynamic point data object
    pointData = reader.GetOutput().GetPointData()

    # Initialize a list to store scalar arrays
    scalarArrays = []
    # Iterate over the number of arrays in pointData
    for i in range(pointData.GetNumberOfArrays()):
        array = pointData.GetArray(i)
        # Check if it's a scalar array
        if array.GetNumberOfComponents() == 1:
            scalarArrays.append(np.array(array))
    noiseScalars = scalarArrays[0]
    noiseScalars = np.reshape(noiseScalars, shape=(NX, NY))
    return noiseScalars, dims


def processInitialConfigurations() -> None:
    # Copy all template files
    print("\033[1;32mCopying template files...\033[0m")
    for key, src_path in file_templates.items():
        dest = f"inputFiles/{src_path.name}"
        shutil.copy2(src_path.resolve(), dest)
        # shutil.copy2(src_path, dest)
        print(f"Created {dest}")

    print("\033[1;32mCreating  noiseFile\033[0m")
    # copy noise file into inputFiles
    shutil.copy2(
        MDSolidSolution_parameters["variables"]["correlationFile_xz"], "inputFiles/"
    )
    shutil.copy2(
        MDSolidSolution_parameters["variables"]["correlationFile_yz"], "inputFiles/"
    )
    # copy declared variables to the configuration txt file
    for param, value in MDSolidSolution_parameters["variables"].items():
        if "correlationFile" in param:
            relativePath = f"{str(value).split('/')[-1]}"
            setInputVariable(MDSolidSolution_parameters["copy_to"], param, relativePath)
        else:
            setInputVariable(MDSolidSolution_parameters["copy_to"], param, str(value))
    for param, data in MDSolidSolution_parameters["vectors"].items():
        setInputVector(
            MDSolidSolution_parameters["copy_to"], param, data["value"], data["comment"]
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

    simulationDir = os.path.abspath(".")
    ddBase = pyMoDELib.DislocationDynamicsBase(simulationDir)

    matFile = "inputFiles/FeCrAl_Fe.txt"
    absoluteTemp = 1
    mat = pyMoDELib.PolycrystallineMaterialBase(matFile, absoluteTemp)
    b_SI = getValueInFile(f"{matFile}", "b_SI")
    mu0_SI = getValueInFile(f"{matFile}", "mu0_SI")
    rho_SI = getValueInFile(f"{matFile}", "rho_SI")

    tag = "1"
    originalCorrFile_xz = "./inputFiles/MoDELCompatible_FeCr8_xz.vtk"
    originalCorrFile_yz = "./inputFiles/MoDELCompatible_FeCr8_yz.vtk"
    seed = 1
    gridSize = np.array([100, 100, 1])
    gridSpacing = np.array([1e-10, 1e-10, 1e-10]) / b_SI
    latticeBasis = np.array([[1, 0], [0, 1]]).reshape(2, 2)
    a_cai_SI = 1e-10 / b_SI

    ssNoise = pyMoDELib.MDSolidSolutionNoise(
        mat,
        tag,
        originalCorrFile_xz,
        originalCorrFile_yz,
        seed,
        gridSize.reshape(3, 1),
        gridSpacing.reshape(3, 1),
        latticeBasis,
        a_cai_SI,
    )
    realizationNum = 100
    corr_xz, corr_yz = np.array(ssNoise.averageNoiseCorrelation(realizationNum)) * mu0_SI**2

    corr_xz = circularShift(corr_xz, gridSize)
    corr_yz = circularShift(corr_yz, gridSize)
    corr_xz = corr_xz.reshape(100, 100)
    corr_yz = corr_yz.reshape(100, 100)

    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 6), dpi=200)
    # First plot with colorbar
    axes[0].set_title("corr_xz")
    im0 = axes[0].imshow(corr_xz, cmap="coolwarm", origin="lower")
    fig.colorbar(im0, ax=axes[0], label="Pa2")
    axes[0].set_xlabel("X")
    axes[0].set_ylabel("Y")

    # Second plot with colorbar
    axes[1].set_title("corr_yz")
    im1 = axes[1].imshow(corr_yz, cmap="coolwarm", origin="lower")
    fig.colorbar(im1, ax=axes[1], label="Pa2")
    axes[1].set_xlabel("X")
    axes[1].set_ylabel("Y")
    plt.tight_layout()
    fig.savefig(f"ensemCorrSS_R{realizationNum}.png", transparent=True)

    return 0


if __name__ == "__main__":
    main()
