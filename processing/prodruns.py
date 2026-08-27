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

from scipy.optimize import curve_fit

import postProcLib

from importlib import reload
reload(postProcLib)

# %%
# data = postProcLib.importRun("prodrun2", 8)
data = postProcLib.importRun("prodrun_singlescatterbins_100000", 8)
# data = postProcLib.importRun("prodrun_temp_0025", 2)
# data = postProcLib.importRun("prodrun_cutoff2", 8)

# %%
fig, axs = plt.subplots(1, 2, figsize=(11, 4))

def fit(x, a, b, c, d):
    return x**a * np.exp(b * x) * (x - c)**d

for run in data:
    ax = axs[0]
    if run["info"]["Recoil"] == True:
        ax = axs[1]

    numHist = run["data"]["hists"]["num"]

    ax.stairs(numHist["totalCounts"], numHist["walls"], label=run["info"]["ElectronTemp"])

    xVals = numHist["centers"]
    yVals = numHist["totalCounts"]
    # yVals = np.abs(np.gradient(numHist["totalCounts"], numHist["centers"])

    # ax.scatter(xVals, yVals, marker=",", s=1, label=run["info"]["ElectronTemp"])

    # coeffs, _ = curve_fit(fit, xVals ,yVals, [-0.5, -1, 1, 1])

    # ax.scatter(xVals, fit(xVals, coeffs[0], coeffs[1], coeffs[2], coeffs[3]), label=f"{run["info"]["ElectronTemp"]} fit")

    # newCounts, newWalls = postProcLib.rebin(numHist["totalCounts"], numHist["walls"], 50, True)
    # ax.stairs(newCounts, newWalls, label=run["info"]["ElectronTemp"])
    ax.set_xscale("log")
    ax.set_yscale("log")
    # ax.set_xlim([0.1, 1000])

axs[0].set_title("No Recoil")
axs[1].set_title("Recoil")

axs[0].legend()
axs[1].legend()

# %%
reload(postProcLib)

def numPlotVals(theseParams, thisHist):
    # plotVals = thisHist["totalCounts"]
    plotVals = thisHist["totalNormalized"]

    # diffs = thisHist["perpNormalized"] - thisHist["parNormalized"]
    diffs = thisHist["perpCounts"]/thisHist["parCounts"] - 3

    xWalls = thisHist["walls"]

    return plotVals, xWalls, diffs

def muPlotVals(theseParams, thisHist):
    plotVals = thisHist["totalNormalized"]
    plotVals = plotVals / np.sin(thisHist["centers"]) - 0.5

    diffs = thisHist["perpNormalized"] - thisHist["parNormalized"]
    diffs = diffs / np.sin(thisHist["centers"])

    xWalls = np.cos(thisHist["walls"])

    return plotVals, xWalls, diffs

def escThetaPlotVals(theseParams, thisHist):
    # plotVals = thisHist["totalCounts"]
    plotVals = thisHist["totalNormalized"]

    diffs = thisHist["perpNormalized"] - thisHist["parNormalized"]

    xWalls = thisHist["walls"]

    return plotVals, xWalls, diffs

def nrgPlotVals(theseParams, thisHist):
    # plotVals = thisHist["totalCounts"]
    plotVals = thisHist["totalNormalized"]

    diffs = thisHist["perpNormalized"] - thisHist["parNormalized"]

    xWalls = thisHist["walls"] / theseParams["FieldStrength"]

    return plotVals, xWalls, diffs

def escNRGPlotVals(theseParams, thisHist):
    # plotVals = thisHist["totalCounts"]
    plotVals = thisHist["totalNormalized"]

    diffs = thisHist["perpNormalized"] - thisHist["parNormalized"]

    xWalls = thisHist["walls"] / theseParams["FieldStrength"]

    return plotVals, xWalls, diffs

keys = ["num", "theta", "nrg", "esc_theta", "esc_nrg"]
funcs = [numPlotVals, muPlotVals, nrgPlotVals, escThetaPlotVals, escNRGPlotVals]

for key, func in zip(keys, funcs):
    postProcLib.recoilComparisonDiffPlot(data, key, func)


# %%
def getAxesFromKey(key):
    if key[0:3] == "nrg":
        xlabel = fr"$\varepsilon_i$"
    else:
        xlabel = fr"$\theta_i$"
        
    if (key[-3:] == "nrg"):
        ylabel = fr"$\varepsilon_f$"
    else:
        ylabel = fr"$\theta_f$"

    if (key == "finalVals"):
        xlabel = fr"$\varepsilon_f$"
        ylabel = fr"$\theta_f$"

    return (xlabel, ylabel)

def plot2DHists(key, i, ax):
    tempVals = [0.05, 0.025, 0.01, 0.005]

    for run in data:
        recoil = True

        if (run["info"]["ElectronTemp"] == tempVals[i] and run["info"]["Recoil"] == recoil):
            # Get relevant 2D hist
            thisHist2D = run["data"]["hists2D"][key]

            xWalls = thisHist2D["xWalls"]
            yWalls = thisHist2D["yWalls"]

            # plotVals = thisHist2D["totalCounts"]
            # plotVals = thisHist2D["totalNormalized"]
            plotVals = thisHist2D["perpNormalized"] - thisHist2D["parNormalized"]

            image = ax.pcolormesh(xWalls, yWalls, plotVals, cmap='inferno')

            xlabel, ylabel = getAxesFromKey(key)

            if (i == 0):
                ax.set_xlabel(xlabel)
                ax.set_ylabel(ylabel)
            else:
                ax.get_yaxis().set_visible(False)
            
            if (ylabel==fr"$\varepsilon_f$"):
                ax.set_yscale('log')

            cbar = plt.colorbar(image, ax=ax, fraction=0.046, pad=0.04)

            if (i == 3):
                cbar.set_label("Relative Counts")

            ax.set_box_aspect(1)

            # figure out the title
            ax.set_title(run["info"]["ElectronTemp"])



postProcLib.plotAllKeys(["nrgXnrg", "nrgXtheta", "thetaXnrg", "thetaXtheta", "finalVals"], plot2DHists)

# %%
run = data[2]

thisAvgDict = run["data"]["avgData"]

myMap = 'inferno'

def plotAvgData(key, i, ax):
    # Get relevant average data
    omegas = thisAvgDict["w"]
    thetas = thisAvgDict["th"]
    plotVal = thisAvgDict[f"{key}_avg"]

    if (key == "w"):
        image = ax.pcolormesh(omegas, thetas, plotVal, cmap=myMap, norm=colors.Normalize(7e-2, 1.3e-1))
    else:
        image = ax.pcolormesh(omegas, thetas, plotVal, cmap=myMap)

    if (i == 0):
        ax.set_xlabel(fr"$\varepsilon_i$")
        ax.set_ylabel(fr"$\theta_i$")

        ax.set_title("Average escape energy")
    if (i == 1):
        ax.set_title("Average escape angle")

    if (i == 2):
        ax.set_title("Average scatter count")

    cbar = plt.colorbar(image, ax=ax, fraction=0.046, pad=0.04)

    # if (ind != 3 and key == "w"):
    #     cbar.set_ticks([])

    ax.set_box_aspect(1)

fig, axes = plt.subplots(1, 3, figsize=(18, 10))

keys = ["w", "th", "N"]

for i in range(3):

    plotAvgData(keys[i], i, axes[i])
