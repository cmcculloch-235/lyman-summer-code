#!/usr/bin/env python3


import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

spec_type="test"
file_type="test"
# load the data
#assert(len(sys.argv) > 1)
data_file_z2 = "out/test/xcorr_delta_H_delta_H.dat"

# first argument specifies output format
EXTENSION = "pdf"
if len(sys.argv) > 1:
    EXTENSION = sys.argv[1]

with open(data_file_z2) as data_fd:
    data_rows = np.array([[ float(j) for j in l.strip().split(" ")] for l in data_fd.readlines()])
    data = data_rows.transpose()
    ks_z2 = data[0]
    power_z2 = data[1]
    count_z2 = data[2]


# LaTeX mode
matplotlib.rcParams["text.usetex"] = True
matplotlib.rcParams["savefig.dpi"] = 300
matplotlib.rcParams["font.size"] = 18


fig_scale = 12
fig_aspect = 4/3
fig, axis = plt.subplots(1, 1, figsize=(fig_scale, fig_scale / fig_aspect))
stdev = power_z2 / np.sqrt(count_z2)
axis.plot(ks_z2, power_z2 + stdev, linewidth=0.8, color="purple")
axis.plot(ks_z2, power_z2 - stdev, linewidth=0.8, color="purple")
axis.plot(ks_z2, power_z2 + 2 * stdev, linewidth=0.6, color="green")
axis.plot(ks_z2, power_z2 - 2 * stdev, linewidth=0.6, color="green")
axis.plot(ks_z2, power_z2, color="blue", label="$z=2$")


axis.set_xlim(left=2e-1)

axis.legend()

axis.set_xlabel("$k/ (h/\mathrm{Mpc})$")
axis.set_title("Extracted "+ spec_type +" power spectrum")

axis.set_xscale("log", basex=10)
axis.set_yscale("log", basey=10)

plt.tight_layout()
plt.savefig("out/xcorr_" + file_type + "." + EXTENSION)


