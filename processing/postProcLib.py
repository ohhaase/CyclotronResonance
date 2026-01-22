import numpy as np
import matplotlib.pyplot as plt

import os

# Function for importing raw data from csv file
def importData(filePath):
    return np.genfromtxt(filePath, delimiter=',')

# Functions to import the different sim outputs
def importAverages(runFilePath, simParams):
    rawdata = importData(f"{runFilePath}/avgdata.csv")

    Nbins = simParams["Nbins"]

    omegas = np.reshape(rawdata[:, 0], [Nbins, Nbins])
    thetas = np.reshape(rawdata[:, 1], [Nbins, Nbins])
    avgOmegas = np.reshape(rawdata[:, 2], [Nbins, Nbins])
    avgThetas = np.reshape(rawdata[:, 3], [Nbins, Nbins])
    avgNum = np.reshape(rawdata[:, 4], [Nbins, Nbins])
    avgPol = np.reshape(rawdata[:, 5], [Nbins, Nbins])

    if np.shape(rawdata)[1] > 5:
        binTimings = np.reshape(rawdata[:, 6], [Nbins, Nbins])
    else:
        binTimings = np.zeros_like(omegas)

    avgData = {
        "w": omegas,
        "th": thetas,
        "w_avg": avgOmegas,
        "th_avg": avgThetas,
        "N_avg": avgNum,
        "pol_avg": avgPol,
        "timings": binTimings
    }

    return avgData

def importHistogram(runFilePath, name):
    rawData = importData(f"{runFilePath}/hist_{name}.csv")

    wallsData = rawData[0, :-1]
    parCountsData = rawData[1, :-2]
    perpCountsData = rawData[2, :-2]

    # Trick to deal with overflows lol. Didn't think I'd ever need this
    # Python doesn't care fortunately, because ints are implemented very smartly
    parNegMask = parCountsData < 0
    perpNegMask = perpCountsData < 0

    parCountsData[parNegMask] = 2147483647 + (parCountsData[parNegMask] + 2147483647 + 2)
    perpCountsData[perpNegMask] = 2147483647 + (perpCountsData[perpNegMask] + 2147483647 + 2)

    hist = {
        "walls": wallsData,
        "parCounts": parCountsData,
        "perpCounts": perpCountsData,
        "totalCounts": parCountsData + perpCountsData
    }

    return hist

def importHistogram2D(runFilePath, name):
    rawData = importData(f"{runFilePath}/hist2D_{name}.csv")

    xWallsData = rawData[0, 1:-1:2]
    yWallsData = rawData[1:, 0]
    parCountsData = rawData[1:-1, 1:-3:2]
    perpCountsData = rawData[1:-1, 2:-3:2]

    hist = {
        "xWalls": xWallsData,
        "yWalls": yWallsData,
        "parCounts": parCountsData,
        "perpCounts": perpCountsData,
        "totalCounts": parCountsData + perpCountsData
    }

    return hist

# Functions to import all different outputs
def importAllHists(runFilePath):
    betaHist = importHistogram(runFilePath, "beta")
    nrgHist = importHistogram(runFilePath, "nrg")
    thetaHist = importHistogram(runFilePath, "theta")
    numHist = importHistogram(runFilePath, "num",)
    escapeNRGHist = importHistogram(runFilePath, "esc_nrg")
    escapeThetaHist = importHistogram(runFilePath, "esc_theta")

    histograms = {
        "beta": betaHist,
        "nrg": nrgHist,
        "num": numHist,
        "theta": thetaHist,
        "esc_nrg": escapeNRGHist,
        "esc_theta": escapeThetaHist
    }

    return histograms

def importAll2DHists(runFilePath):
    nrgXnrgHist = importHistogram2D(runFilePath, "nrgXnrg")
    nrgXthetaHist = importHistogram2D(runFilePath, "nrgXtheta")
    thetaXnrgHist = importHistogram2D(runFilePath, "thetaXnrg")
    thetaXthetaHist = importHistogram2D(runFilePath, "thetaXtheta")
    finalValsHist = importHistogram2D(runFilePath, "finalVals")

    histograms = {
        "nrgXnrg": nrgXnrgHist,
        "nrgXtheta": nrgXthetaHist,
        "thetaXnrg": thetaXnrgHist,
        "thetaXtheta": thetaXthetaHist,
        "finalVals": finalValsHist
    }

    return histograms

# Function to import a single run (within the 4-run runs)
def importSubRunData(runFilePath, simParams):
    runData = {
        "avgData": importAverages(runFilePath, simParams),
        "hists": importAllHists(runFilePath),
        "hists2D": importAll2DHists(runFilePath)
    }

    return runData

# Function to get extra tags based on simulation type
def tagsFromSubFolder(folder, quadRunName):
    Theta = 0
    recoil = ""
    distb = ""

    if folder[0] == "0":
        Theta = 0.05
        if "lowertemps" in quadRunName: Theta = 0.01
    else:
        Theta = 0.025
        if "lowertemps" in quadRunName: Theta = 0.005

    if folder[1] == "0":
        recoil = "No Recoil"
    else:
        recoil = "Recoil"

    if ("MB" in quadRunName) or (quadRunName[:3].isdigit()):
        # I think this covers all the runs. Might be wrong tho with some of the unlabeled ones
        distb = "MB"
    else:
        distb = "MJ"

    simTags = {
        "Theta": Theta,
        "recoil": recoil,
        "distb": distb
    }

    return simTags

# Function to get all 4 runs from a single 4-run
def importQuadRunData(quadRunName):
    rawParams = importData(f"../data/{quadRunName}/simParams.csv")

    Nparticles = rawParams[0]
    Nbins = rawParams[1].astype(int)
    simType = rawParams[2]

    simParams = {
        "Nparticles": Nparticles,
        "Nbins": Nbins,
        "simType": simType
    }

    folderNames = ["00", "01", "10", "11"]
    runs = []

    for folder in folderNames:
        runData = importSubRunData(f"../data/{quadRunName}/{folder}", simParams)

        runTags = tagsFromSubFolder(folder, quadRunName)

        thisRun = {
            "params": simParams,
            "tags": runTags,
            "data": runData
        }

        runs.append(thisRun)
    
    return runs

# Function to get all 4 runs from a 4x scatternum run
def importScatterNumData(scatterLim):
    rawParams = importData(f"../data/scatterNum_{scatterLim}/simParams.csv")

    Nparticles = rawParams[0]
    Nbins = rawParams[1].astype(int)
    simType = rawParams[2]

    simParams = {
        "Nparticles": Nparticles,
        "Nbins": Nbins,
        "simType": simType
    }

    folderNames = ["0", "1", "2", "3"]
    folderTemps = [0.01, 0.007, 0.005, 0.001]
    runs = []

    for folder in folderNames:
        runData = importSubRunData(f"../data/scatterNum_{scatterLim}/{folder}", simParams)

        runTags = {
            "Theta": folderTemps[int(folder)],
            "recoil": "No Recoil",
            "distb": "MB",
            "scatterLim": scatterLim
        }

        thisRun = {
            "params": simParams,
            "tags": runTags,
            "data": runData
        }

        runs.append(thisRun)
    
    return runs

def plotAllKeys(keys, plotFunc):
    # Loop over every key (diff for hists vs avg data)
    for key in keys:
        
        # 1x4 subplot to show them all
        fig, axes = plt.subplots(1, 4, figsize=(22, 10))

        for i in range(4):

            plotFunc(key, i, axes[i])

        plt.show()

def test():
    print(os.getcwd())

myMap = 'inferno'