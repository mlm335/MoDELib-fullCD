import os
import sys
import glob
import string
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams["text.usetex"] = True
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
    "noise_file": Path("../../Library/GlidePlaneNoise/AnalyticalSolidSolutionNoise.txt"),
    "material_file": Path("../../Library/Materials/AlMg5.txt"),
    "elastic_deformation_file": Path(
        "../../Library/ElasticDeformation/ElasticDeformation.txt"
    ),
    "mesh": Path("../../Library/Meshes/unitCube24.msh"),
    "microstructure": Path(
        "../../Library/Microstructures/periodicDipoleIndividual.txt"
    ),
}

AnalyticalSolidSolution_parameters = {
    "variables": {
        # Scalar parameters (use setInputVariable)
        "type": "AnalyticalSolidSolutionNoise",
        "tag": 1,
        "seed": 2,
        "outputNoise": 1,
        "a_cai_SI": 1,
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
    "copy_to": "inputFiles/AnalyticalSolidSolutionNoise.txt",
}

Material_parameters = {
    "variables": {
        "enabledSlipSystems": "Shockley",
        #'glidePlaneNoise': ['MDSolidSolution.txt', 'MDStackingFault.txt'],
        "glidePlaneNoise": "AnalyticalSolidSolutionNoise.txt",
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
    # copy declared variables to the configuration txt file
    for param, value in AnalyticalSolidSolution_parameters["variables"].items():
        if "correlationFile" in param:
            relativePath = f"{str(value).split('/')[-1]}"
            setInputVariable(AnalyticalSolidSolution_parameters["copy_to"], param, relativePath)
        else:
            setInputVariable(AnalyticalSolidSolution_parameters["copy_to"], param, str(value))
    for param, data in AnalyticalSolidSolution_parameters["vectors"].items():
        setInputVector(
            AnalyticalSolidSolution_parameters["copy_to"], param, data["value"], data["comment"]
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


def plotComparisonHistogram(realizationNumbers: list):
    # Plotting on the first subplot
    fig, axs = plt.subplots(1, 1, figsize=(8, 6), dpi=200)  # Adjusted for two subplots

    matFile = str(file_templates["material_file"]).split("/")[-1]
    b_SI = getValueInFile(f"inputFiles/{matFile}", "b_SI")
    mu0_SI = getValueInFile(f"inputFiles/{matFile}", "mu0_SI")

    inputPath = "./inputFiles/"
    for i, realizationNum in enumerate(realizationNumbers):
        data = np.loadtxt(f"{inputPath}/noiseDistributionR{realizationNum}.txt")
        data *= mu0_SI * b_SI
        axs.hist(data, density=True, alpha=0.3, bins="auto", label=f"R{realizationNum}")

        std = np.std(data)
        # Annotate in the top right, offsetting each annotation
        axs.text(
            0.98,
            0.95 - 0.05 * i,
            f"R{realizationNum} = STD {std:.3g} J/m2",
            transform=axs.transAxes,
            ha="right",
            va="top",
        )

    axs.grid(True)
    axs.set_title(f"")
    axs.set_xlabel(f"")
    axs.set_ylabel(f"")
    axs.legend()
    fig.savefig(f"sampledNoiseDistribution.png")
    plt.close()


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

    matFile = "inputFiles/AlMg5.txt"
    absoluteTemp = 1
    mat = pyMoDELib.PolycrystallineMaterialBase(matFile, absoluteTemp)
    b_SI = getValueInFile(f"{matFile}", "b_SI")
    mu0_SI = getValueInFile(f"{matFile}", "mu0_SI")
    rho_SI = getValueInFile(f"{matFile}", "rho_SI")

    tag = "1"
    seed = 0
    gridSize = np.array([256, 256, 64])
    gridSpacing = np.array([1e-10, 1e-10, 1e-10])/b_SI
    latticeBasis = np.array([[1, 0], [0, 1]]).reshape(2, 2)

    msss_SI=0.9e18 /mu0_SI**2; #[Pa^2] Mean Square Shear Stress 
    a = 1e-10 /b_SI
    a_cai_SI = 1e-09 /b_SI # in meter
    #a_cai_SI = 
    anssNoise = pyMoDELib.AnalyticalSolidSolutionNoise(
        tag,
        seed,
        gridSize.reshape(3, 1),
        gridSpacing.reshape(3, 1),
        latticeBasis,
        a,
        a_cai_SI,
        msss_SI
    )

    realizationNum = 100
    corr_xz, corr_yz = np.array(anssNoise.averageNoiseCorrelation(realizationNum)) * mu0_SI**2

    corr_xz = circularShift(corr_xz, gridSize)
    corr_yz = circularShift(corr_yz, gridSize)
    corr_xz = corr_xz.reshape(gridSize[0], gridSize[1])
    corr_yz = corr_yz.reshape(gridSize[0], gridSize[1])

    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 6), dpi=200)
    # First plot with colorbar
    axes[0].set_title("corr_xz")
    im0 = axes[0].imshow(corr_xz, cmap="coolwarm", origin="lower")
    fig.colorbar(im0, ax=axes[0], label="Noise Amplitude")
    axes[0].set_xlabel("X (grid units)")
    axes[0].set_ylabel("Y (grid units)")
    
    # Second plot with colorbar
    axes[1].set_title("corr_yz")
    im1 = axes[1].imshow(corr_yz, cmap="coolwarm", origin="lower")
    fig.colorbar(im1, ax=axes[1], label="Noise Amplitude")
    axes[1].set_xlabel("X (grid units)")
    axes[1].set_ylabel("Y (grid units)")
    plt.tight_layout()
    fig.savefig(f"ensemCorr_R{realizationNum}.png", transparent=True)

    return 0


if __name__ == "__main__":
    main()
