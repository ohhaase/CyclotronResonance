# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.0
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %%
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import postProcLib

# All runs are run with the following standard settings:
# Nparticles: 100
# Nbins: 20
# Simtype: 4 (to allow some comparison between temp and recoil)
# Ncores will change

# %%
mac_runs = [f"mac_{num}thread" for num in [1,4,8,16]]
win_runs = [f"win_{num}thread" for num in [1,4,8,16]]
mac_local_runs = [f"mac_{num}thread_local" for num in [1,4,8,16]]
win_local_runs = [f"win_{num}thread_local" for num in [1,4,8,16]]
mac_hists_runs = [f"mac_{num}thread_hists" for num in [1,4,8,16]]
mac_long_runs = [f"mac_{num}thread_longrun" for num in [1,2,8,16]] # Note: 16 thread run undoes the histogram changes

# testRuns = mac_runs + win_runs + mac_local_runs + win_local_runs
testRuns = mac_runs + mac_local_runs + mac_hists_runs + mac_long_runs

data = []

for testRun in testRuns:
    data.extend(postProcLib.importQuadRunData(f"testRuns/{testRun}"))
    

# %%
# Angular hists (This is what we're trying to make sure works!)
def plotHistData(key, i, ax):
    
    for j, testRun in enumerate(testRuns):
            
        # Get relevant average data
        thisHist = data[j*4 + i]["data"]["hists"][key]

        countType = "totalCounts"
        # countType = "parCounts"
        # countType = "perpCounts"

        xWalls = thisHist["walls"]
        xVals = np.linspace((xWalls[0]+xWalls[1])/2, (xWalls[-2] + xWalls[-1])/2, xWalls.size-1) # Cell centered xvals

        plotVals = thisHist[countType]/np.trapezoid(thisHist[countType], xVals)
        wallVals = xWalls

        if (key == "theta"):
            plotVals = plotVals / np.sin(xVals) - 0.5
            wallVals = np.cos(xWalls)

        ax.stairs(plotVals, wallVals, fill=False)

    ax.legend(testRuns)

    if (i == 0):
        match (key):
            case "beta":
                name = r"$\beta$"
            case "nrg" | "esc_nrg":
                name = r"$\omega$"
            case "num":
                name = "Num Scatterings"
            case "theta" | "esc_theta":
                name = r"$\theta$"
        ax.set_xlabel(name)
        ax.set_ylabel("Relative Counts")

    ax.set_box_aspect(1)

    match (key):
        case "nrg" | "esc_nrg":
            ax.set_xscale('log')
        case "num":
            ax.set_xscale('log')
            ax.set_yscale('log')
        case "theta":
            ax.hlines(0, -1, 1, "grey", "dashed")

    ax.set_title(postProcLib.titleFromFolder(data[j*4 + i]["name"]))

# postProcLib.plotAllKeys(["theta"], plotHistData)
postProcLib.plotAllKeys(["beta", "nrg", "theta", "num", "esc_nrg", "esc_theta"], plotHistData)


# %%
# Timing data
def plotTimingData(key, i, ax):
    # Key is just here so we can use the plotAllKeys func to get 4 graphs
    theseScatterPoints = []

    for j, testRun in enumerate(testRuns):
        timing = np.sum(data[j*4 + i]["data"]["avgData"]["timings"]/1e6) / (data[j*4 + i]["params"]["Nparticles"] * data[j*4 + i]["params"]["Nbins"]**2)

        theseScatterPoints.append(timing)

    ax.bar(testRuns, theseScatterPoints)

    if (i == 0):
        ax.set_ylabel("Average time per particle(sec)")
    
    plt.setp(ax.get_xticklabels(), rotation=30, ha='right')

postProcLib.plotAllKeys([""], plotTimingData)


# %%
# Plot avg data for one run

def plotAvgData(key, i, ax):
    # Get relevant average data
    thisAvgDict = data[12+i]["data"]["avgData"]

    omegas = thisAvgDict["w"]
    thetas = thisAvgDict["th"]
    if key == "timings":
        plotVal = thisAvgDict[key] / 1e6
        ax.text(0.085, 3.0, f"Total time elapsed: {np.sum(plotVal):.2f}", color="w")
    else:
        plotVal = thisAvgDict[f"{key}_avg"]

    if (key == "w"):
        image = ax.pcolormesh(omegas, thetas, plotVal, cmap=postProcLib.myMap, norm=colors.Normalize(7e-2, 1.5e-1))
    else:
        image = ax.pcolormesh(omegas, thetas, plotVal, cmap=postProcLib.myMap)

    if (i == 0):
        ax.set_xlabel(fr"$\varepsilon_i$")
        ax.set_ylabel(fr"$\theta_i$")

    cbar = plt.colorbar(image, ax=ax, fraction=0.046, pad=0.04)

    # if (ind != 3 and key == "w"):
    #     cbar.set_ticks([])

    ax.set_box_aspect(1)

    # figure out the title
    ax.set_title(postProcLib.titleFromFolder(data[i]["name"]))

postProcLib.plotAllKeys(["w", "th", "N", "pol", "timings"], plotAvgData)
