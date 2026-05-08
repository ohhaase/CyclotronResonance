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

# %%
data = postProcLib.importRun("prodrun1", 8)
# data = postProcLib.importRun("quicktests", 8)

# %%
def numPlotVals(thisHist):
    # plotVals = thisHist["totalCounts"]
    plotVals = thisHist["totalNormalized"]

    diffs = thisHist["perpNormalized"] - thisHist["parNormalized"]

    xWalls = thisHist["walls"]

    return plotVals, xWalls, diffs

def muPlotVals(thisHist):
    plotVals = thisHist["totalNormalized"]
    plotVals = plotVals / np.sin(thisHist["centers"]) - 0.5

    diffs = thisHist["perpNormalized"] - thisHist["parNormalized"]
    diffs = diffs / np.sin(thisHist["centers"])

    xWalls = np.cos(thisHist["walls"])

    return plotVals, xWalls, diffs

def escThetaPlotVals(thisHist):
    # plotVals = thisHist["totalCounts"]
    plotVals = thisHist["totalNormalized"]

    diffs = thisHist["perpNormalized"] - thisHist["parNormalized"]

    xWalls = thisHist["walls"]

    return plotVals, xWalls, diffs

def nrgPlotVals(thisHist):
    # plotVals = thisHist["totalCounts"]
    plotVals = thisHist["totalNormalized"]

    diffs = thisHist["perpNormalized"] - thisHist["parNormalized"]

    xWalls = thisHist["walls"]

    return plotVals, xWalls, diffs

def escNRGPlotVals(thisHist):
    # plotVals = thisHist["totalCounts"]
    plotVals = thisHist["totalNormalized"]

    diffs = thisHist["perpNormalized"] - thisHist["parNormalized"]

    xWalls = thisHist["walls"]

    return plotVals, xWalls, diffs

keys = ["num", "theta", "nrg", "esc_theta", "esc_nrg"]
funcs = [numPlotVals, muPlotVals, nrgPlotVals, escThetaPlotVals, escNRGPlotVals]

for key, func in zip(keys, funcs):
    # postProcLib.recoilComparisonPlot(data, key, func)
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
        recoil = False

        if (run["info"]["ElectronTemp"] == tempVals[i] and run["info"]["Recoil"] == recoil):
            # Get relevant 2D hist
            thisHist2D = run["data"]["hists2D"][key]

            xWalls = thisHist2D["xWalls"]
            yWalls = thisHist2D["yWalls"]

            plotVals = thisHist2D["totalCounts"]
            # plotVals = thisHist2D["totalNormalized"]
            # plotVals = thisHist2D["perpNormalized"] - thisHist2D["parNormalized"]

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
                cbar.set_label("Counts")

            ax.set_box_aspect(1)

            # figure out the title
            ax.set_title(run["info"]["ElectronTemp"])



postProcLib.plotAllKeys(["nrgXnrg", "nrgXtheta", "thetaXnrg", "thetaXtheta", "finalVals"], plot2DHists)
