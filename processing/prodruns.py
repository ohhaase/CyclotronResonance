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

            xWalls = thisHist["walls"]
            
            # plotVals = thisHist["totalCounts"]
            plotVals = thisHist["totalNormalized"]
            # plotVals = thisHist["perpNormalized"] - thisHist["parNormalized"]

            ax.stairs(plotVals, xWalls, fill=False)

            legendTitles.append(run["info"]["ElectronTemp"])    

    if (i == 0):
        # ax.set_ylabel("Counts")
        ax.set_ylabel("Relative Counts")
    
    
    ax.set_xlabel(r"Scatter Num")

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

    legendTitles = []
    for run in data:

        if (run["info"]["Recoil"] == recoilBool):
            # Get relevant numscatter hist data
            thisHist = run["data"]["hists"]["theta"]

            plotVals = thisHist["totalNormalized"]
            plotVals = plotVals / np.sin(thisHist["centers"]) - 0.5

            # plotVals = thisHist["perpNormalized"] - thisHist["parNormalized"]
            # plotVals = plotVals / np.sin(thisHist["centers"])

            wallVals = np.cos(thisHist["walls"])
            
            ax.stairs(plotVals, wallVals, fill=False)

            legendTitles.append(run["info"]["ElectronTemp"])    

    if (i == 0):
        ax.set_ylabel("Relative Counts")
    
    
    ax.set_xlabel(r"$\cos{\Theta}$")

    ax.set_box_aspect(1)

    # ax.set_xscale('log')
    # ax.set_yscale('log')

    ax.legend(legendTitles, title=r"$\mathcal{T}=kT/mc^2$")

    ax.set_title(recoil)

# %%
# Plot for the scatter nrg histograms

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
            thisHist = run["data"]["hists"]["nrg"]

            wallVals = thisHist["walls"]
            
            plotVals = thisHist["totalNormalized"]
            # plotVals = thisHist["perpNormalized"] - thisHist["parNormalized"]
            
            ax.stairs(plotVals, wallVals, fill=False)

            legendTitles.append(run["info"]["ElectronTemp"])    

    if (i == 0):
        ax.set_ylabel("Relative Counts")
    
    
    ax.set_xlabel(r"$\omega$")

    ax.set_box_aspect(1)

    # ax.set_xscale('log')
    # ax.set_yscale('log')

    ax.legend(legendTitles, title=r"$\mathcal{T}=kT/mc^2$")

    ax.set_title(recoil)

# %%
# Plot for the scatter escape theta histograms

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
            thisHist = run["data"]["hists"]["esc_theta"]

            wallVals = thisHist["walls"]
            
            plotVals = thisHist["totalNormalized"]
            # plotVals = thisHist["perpNormalized"] - thisHist["parNormalized"]
            
            ax.stairs(plotVals, wallVals, fill=False)

            legendTitles.append(run["info"]["ElectronTemp"])    

    if (i == 0):
        ax.set_ylabel("Relative Counts")
    
    
    ax.set_xlabel(r"$\Theta$")

    ax.set_box_aspect(1)

    # ax.set_xscale('log')
    # ax.set_yscale('log')

    ax.legend(legendTitles, title=r"$\mathcal{T}=kT/mc^2$")

    ax.set_title(recoil)

# %%
# Plot for the escape nrg histograms

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
            thisHist = run["data"]["hists"]["esc_nrg"]

            wallVals = thisHist["walls"]
            
            plotVals = thisHist["totalNormalized"]
            # plotVals = thisHist["perpNormalized"] - thisHist["parNormalized"]

            ax.stairs(plotVals, wallVals, fill=False)

            legendTitles.append(run["info"]["ElectronTemp"])    

    if (i == 0):
        ax.set_ylabel("Relative Counts")
    
    
    ax.set_xlabel(r"$\omega$")

    ax.set_box_aspect(1)

    ax.set_xscale('log')
    # ax.set_yscale('log')

    ax.legend(legendTitles, title=r"$\mathcal{T}=kT/mc^2$")

    ax.set_title(recoil)


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
