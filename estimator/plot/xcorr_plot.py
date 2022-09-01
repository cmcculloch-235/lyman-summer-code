#!/usr/bin/env python3


import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np


# load the data
#assert(len(sys.argv) > 1)
data_file = "out/xcorr_delta_matter_delta_matter.dat"

# first argument specifies output format
EXTENSION = "pdf"
if len(sys.argv) > 1:
    EXTENSION = sys.argv[1]

with open(data_file) as data_fd:
    data_rows = np.array([[ float(j) for j in l.strip().split(" ")] for l in data_fd.readlines()])
    data = data_rows.transpose()
    ks = data[0]
    power = data[1]
    count = data[2]


# LaTeX mode
matplotlib.rcParams["text.usetex"] = True
matplotlib.rcParams["savefig.dpi"] = 300
matplotlib.rcParams["font.size"] = 18


fig_scale = 12
fig_aspect = 4/3
fig, axis = plt.subplots(1, 1, figsize=(fig_scale, fig_scale / fig_aspect))
stdev = power / np.sqrt(count)
axis.plot(ks, power, color="black", label="Extracted spectrum")
axis.plot(ks, power + stdev, color="red")
axis.plot(ks, power - stdev, color="red")
axis.plot(ks, power + 2 * stdev, "--", color="yellow")
axis.plot(ks, power - 2 * stdev, "--", color="yellow")




axis.legend()

axis.set_xlabel("$k/ (h/\mathrm{Mpc})$")
axis.set_title("Extracted matter power spectrum")

axis.set_xscale("log", basex=10)
axis.set_yscale("log", basey=10)

plt.tight_layout()
plt.savefig("out/pspec." + EXTENSION)


