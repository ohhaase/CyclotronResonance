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

import postProcLib

from importlib import reload
reload(postProcLib)

# %%
data_without = postProcLib.importRun("outputhandler_comp", 8)

data_with = postProcLib.importRun("outputhandler", 8)

# %%
fig, axs = plt.subplots(1, 2, figsize=(11, 4))

colorMap = {0.05: "tab:blue",
            0.01: "tab:orange",
            0.025: "tab:green",
            0.005: "tab:red"}

for run in data_with:
    ax = axs[0]
    if run["info"]["Recoil"] == True:
        ax = axs[1]
    
    time = run["info"]["RunTime"]
    temp = run["info"]["ElectronTemp"]

    ax.scatter(temp, time, marker="x", label="New Outputs")

for run in data_without:
    ax = axs[0]
    if run["info"]["Recoil"] == True:
        ax = axs[1]
    
    time = run["info"]["RunTime"]
    temp = run["info"]["ElectronTemp"]

    ax.scatter(temp, time, marker="o", label="Old Outputs")


axs[0].set_title("No Recoil")
axs[1].set_title("Recoil")

axs[0].legend()
axs[1].legend()

axs[0].set_xlabel("Electron Temp")
axs[1].set_xlabel("Electron Temp")
