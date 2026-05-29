import numpy as np

class AUX:
    # Mesh node data
    meshNodes = np.empty([0, 0])

    # Gauss point data
    gaussPoints = np.empty([0, 0])
    sourceID = np.empty(0)
    sinkID = np.empty(0)
    qID = np.empty(0)
    r = np.empty([0, 0])
    j = np.empty(0)
    rl = np.empty([0, 0])
    stress = np.empty([0, 0, 0])
    pkForce = np.empty([0, 0])
    stackingFaultForce = np.empty([0, 0])
    lineTensionForce = np.empty([0, 0])
    velocity = np.empty([0, 0])
    elasticEnergyPerLength = np.empty(0)
    coreEnergyPerLength = np.empty(0)
    cCD = np.empty([0, 0])
    cDD = np.empty([0, 0])

    # Periodic patches (if any)
    periodicPatches = np.empty([0, 0])

def readAUXtxt(filename):
    try:
        with open(filename + '.txt', "r") as auxFile:
            numNodes = int(auxFile.readline().strip())
            numGPs = int(auxFile.readline().strip())
            numPGPP = int(auxFile.readline().strip())

            aux = AUX()

            aux.meshNodes = np.empty((numNodes, 4))
            aux.gaussPoints = np.empty((numGPs, 39))
            aux.sourceID = np.empty(numGPs)
            aux.sinkID = np.empty(numGPs)
            aux.qID = np.empty(numGPs)
            aux.r = np.empty((numGPs, 3))
            aux.j = np.empty(numGPs)
            aux.rl = np.empty((numGPs, 3))
            aux.stress = np.empty((numGPs, 3, 3))
            aux.pkForce = np.empty((numGPs, 3))
            aux.stackingFaultForce = np.empty((numGPs, 3))
            aux.lineTensionForce = np.empty((numGPs, 3))
            aux.velocity = np.empty((numGPs, 3))
            aux.elasticEnergyPerLength = np.empty(numGPs)
            aux.coreEnergyPerLength = np.empty(numGPs)
            aux.cCD = np.empty((numGPs, 3))
            aux.cDD = np.empty((numGPs, 3))
            aux.periodicPatches = np.empty((numPGPP, 4))

            # mesh nodes
            for k in range(numNodes):
                aux.meshNodes[k, :] = np.fromstring(auxFile.readline().strip(), sep=' ')

            # gauss points
            for k in range(numGPs):
                data = np.fromstring(auxFile.readline().strip(), sep=' ')
                aux.gaussPoints[k, :] = data

                aux.sourceID[k] = data[0]
                aux.sinkID[k] = data[1]
                aux.qID[k] = data[2]
                aux.r[k, :] = data[3:6]
                aux.j[k] = data[6]
                aux.rl[k, :] = data[7:10]
                aux.stress[k, :, :] = data[10:19].reshape(3, 3)
                aux.pkForce[k, :] = data[19:22]
                aux.stackingFaultForce[k, :] = data[22:25]
                aux.lineTensionForce[k, :] = data[25:28]
                aux.velocity[k, :] = data[28:31]
                aux.elasticEnergyPerLength[k] = data[31]
                aux.coreEnergyPerLength[k] = data[32]
                aux.cCD[k, :] = data[33:36]
                aux.cDD[k, :] = data[36:39]

            # periodic patches
            for k in range(numPGPP):
                aux.periodicPatches[k, :] = np.fromstring(auxFile.readline().strip(), sep=' ')

        return aux

    except Exception as e:
        print(f"Error reading file {filename}: {e}")
        return None
