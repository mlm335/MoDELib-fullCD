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
        "a_cai_SI": 1,
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
    seed = 0
    gridSize = np.array([100, 100, 1])
    gridSpacing = np.array([1e-10, 1e-10, 1e-10]) / b_SI
    latticeBasis = np.array([[1, 0], [0, 1]]).reshape(2, 2)

    a_cai_SI = 1e-10 / b_SI  # in meter
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
    sampledNoises = np.array(ssNoise.sampleAverageNoise(realizationNum)) * mu0_SI**2
    # Get all even indices (0, 2, 4,...) for xz
    xzNoises = sampledNoises[::2]
    # Get all odd indices (1, 3, 5,...) for yz
    yzNoises = sampledNoises[1::2]

    fig, axs = plt.subplots(1, 2, figsize=(16,6), dpi=200)
    # (returns a veiw if possible), flattens the array for histogram
    xzNoises = xzNoises.ravel()
    yzNoises = yzNoises.ravel()
    axs[0].hist(xzNoises, density=True, alpha=0.3, bins='auto', label=f"R{realizationNum}")
    axs[1].hist(yzNoises, density=True, alpha=0.3, bins='auto', label=f"R{realizationNum}")

    std_xz = np.std(xzNoises)
    std_yz = np.std(yzNoises)
    print(f"std_xz = {std_xz} Pa2")
    print(f"std_yz = {std_yz} Pa2")

    # Annotate in the top right, offsetting each annotation
    #axs[0].text(
    #    0.98, 0.95 - 0.05*i,
    #    f'R{realizationNum} = STD {std_xz:.3g} Pa2',
    #    transform=axs[0].transAxes,
    #    ha='right',
    #    va='top'
    #)
    #axs[1].text(
    #    0.98, 0.95 - 0.05*i,
    #    f'R{realizationNum} = STD {std_yz:.3g} Pa2',
    #    transform=axs[1].transAxes,
    #    ha='right',
    #    va='top'
    #)

    axs[0].grid(True)
    axs[0].set_title(f"noise xz")
    axs[0].set_xlabel(f"")
    axs[0].set_ylabel(f"")
    axs[0].legend()
    axs[1].grid(True)
    axs[1].set_title(f"noise yz")
    axs[1].set_xlabel(f"")
    axs[1].set_ylabel(f"")
    axs[1].legend()
    fig.savefig(f'sampledNoiseDistribution_R{realizationNum}.png')
    plt.close()

    return 0


if __name__ == "__main__":
    main()
