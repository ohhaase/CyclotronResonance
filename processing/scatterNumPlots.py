import numpy as np
import matplotlib.pyplot as plt
import postProcLib

# Import data
scatterLims = [50, 100, 500, 1000, 5000, 10000, 50000, 10001, 501, 502, 503, 504, 505]

data = []

for scatterLim in scatterLims:
    data.extend(postProcLib.importScatterNumData(scatterLim))


# Plot for the scattering numbers vs temps

keys = ["nrg", "esc_nrg", "theta", "esc_theta", "num"]

for key in keys:
    fig, axes = plt.subplots(1, 4, figsize=(16, 5))

    runTemps = [0.01, 0.007, 0.005, 0.001]

    for i, ax in enumerate(axes):
        temp = runTemps[i]

        legendTitles = []
        for runData in data[8:]: #only looking at relevant data

            if (runData["tags"]["Theta"] == temp):
                # Get relevant numscatter hist data
                thisHist = runData["data"]["hists"][key]

                xWalls = thisHist["walls"]
                xVals = np.linspace((xWalls[0]+xWalls[1])/2, (xWalls[-2] + xWalls[-1])/2, xWalls.size-1) # Cell centered xvals
                
                plotVals = thisHist["totalCounts"]
                # plotVals = thisHist["perpCounts"]
                # plotVals = thisHist["parCounts"]
                # plotVals = thisHist["perpCounts"] - thisHist["parCounts"]

                plotVals = plotVals/np.trapezoid(plotVals, xVals)

                if key == "theta":
                    plotVals = plotVals / np.sin(xVals) - 0.5
                    xWalls = np.cos(xWalls)
                    # print("")
                
                ax.stairs(plotVals, xWalls, fill=False)

                legendTitles.append(runData["tags"]["scatterLim"])

        if (i == 0):
            # ax.set_ylabel("Counts")
            ax.set_ylabel("Relative Counts")
            if key == "theta":
                ax.set_ylabel("Relative Counts - 0.5")
        
        if key == "nrg" or key == "esc_nrg":
            ax.set_xlabel(r"$\omega$")
            ax.set_xscale('log')
        
        if key == "theta":
            ax.set_xlabel(r"$\mu$")

        if key == "esc_theta":
            ax.set_xlabel(r"$\theta$")

        if key == "num":
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.set_xlabel("Num")

        ax.set_box_aspect(1)

        ax.legend(legendTitles, title=r"Limit")

        ax.set_title(r"$\Theta = $" + f"{temp}")


plt.show()