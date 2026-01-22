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
