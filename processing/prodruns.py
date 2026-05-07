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

import postProcLib as processing

# %%
data = processing.importRun("prodrun1", 8)

# %%
# Plot for the scattering numbers vs temps

fig, axes = plt.subplots(1, 2, figsize=(11, 5))

for i, ax in enumerate(axes):
    if (i == 0):
        recoil = "No Recoil"
        recoilBool = False
    else:
        recoil = "Recoil"
        recoilBool = True

    legendTitles = []
    for run in data:

        if (run["info"]["Recoil"] == recoilBool):
            # Get relevant numscatter hist data
            thisHist = run["data"]["hists"]["num"]
            # thisHist2 = data[j]["data"]["hists"]["esc_nrg"]

            xWalls = thisHist["walls"]
            # xWalls2 = thisHist2["walls"]
            xVals = np.linspace(xWalls[0], xWalls[-1], xWalls.size-1)
            # xVals2 = np.linspace(xWalls2[0], xWalls2[-1], xWalls2.size-1)
            
            # plotVals = thisHist["totalCounts"]
            plotVals = thisHist["totalCounts"]/np.trapezoid(thisHist["totalCounts"], xVals)
            # plotVals2 = thisHist2["totalCounts"]/np.trapezoid(thisHist2["totalCounts"], xVals2)
            
            # plotVals = plotVals / np.sin(xVals)
            # wallVals = np.cos(xWalls)

            
            ax.stairs(plotVals, xWalls, fill=False)
            # ax.stairs(plotVals2, xWalls2, fill=False)

            # ax.stairs(plotVals[:25], wallVals[:26], fill=False)
            # ax.stairs(thisHist["perpCounts"], thisHist["walls"], fill=False)
            # ax.stairs(thisHist["parCounts"], thisHist["walls"], fill=False)

            legendTitles.append(run["info"]["ElectronTemp"])    

    if (i == 0):
        # ax.set_ylabel("Counts")
        ax.set_ylabel("Relative Counts")
    
    
    ax.set_xlabel(r"$\theta$")

    ax.set_box_aspect(1)

    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.legend(legendTitles, title=r"$\mathcal{T}=kT/mc^2$")

    ax.set_title(recoil)


# %%
# Plot for the mu histograms

fig, axes = plt.subplots(1, 2, figsize=(11, 5))

for i, ax in enumerate(axes):
    if (i == 0):
        recoil = "No Recoil"
        recoilBool = False
    else:
        recoil = "Recoil"
        recoilBool = True

    # First 8 runs are the lognumhists
    legendTitles = []
    for run in data:

        if (run["info"]["Recoil"] == recoilBool):
            # Get relevant numscatter hist data
            thisHist = run["data"]["hists"]["theta"]

            xWalls = thisHist["walls"]
            xVals = np.linspace((xWalls[0]+xWalls[1])/2, (xWalls[-2] + xWalls[-1])/2, xWalls.size-1) # Cell centered xvals
            
            plotVals = thisHist["totalCounts"]/np.trapezoid(thisHist["totalCounts"], xVals)

            plotVals = plotVals / np.sin(xVals) - 0.5
            wallVals = np.cos(xWalls)
            # wallVals = xWalls 
            
            ax.stairs(plotVals, wallVals, fill=False)

            # ax.stairs(plotVals[:25], wallVals[:26], fill=False)
            # ax.stairs(thisHist["perpCounts"], thisHist["walls"], fill=False)
            # ax.stairs(thisHist["parCounts"], thisHist["walls"], fill=False)

            legendTitles.append(run["info"]["ElectronTemp"])    

    if (i == 0):
        ax.set_ylabel("Relative Counts")
    
    
    ax.set_xlabel(r"$\cos{\Theta}$")

    ax.set_box_aspect(1)

    # ax.set_xscale('log')
    # ax.set_yscale('log')

    ax.legend(legendTitles, title=r"$\mathcal{T}=kT/mc^2$")

    ax.set_title(recoil)
