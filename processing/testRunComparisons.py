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
testRuns = ["mac_1thread"] #, "mac_4thread", "mac_8thread", "mac_16thread"]

data = []

for testRun in testRuns:
    data.extend(postProcLib.importQuadRunData(f"testRuns/{testRun}"))
    

# %%
def plotAvgData(key, i, ax):
    # Get relevant average data
    thisAvgDict = data[i]["data"]["avgData"]

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
