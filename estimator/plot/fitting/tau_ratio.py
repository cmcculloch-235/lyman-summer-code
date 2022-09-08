#!/usr/bin/env python3


import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np


# first argument specifies output format
EXTENSION = "pdf"
if len(sys.argv) > 1:
    EXTENSION = sys.argv[1]

# load the data
data_file_ics = "../out/80_ics/xcorr_delta_DM_delta_DM.dat"
data_file_tau = "../out/80_snap18_old/xcorr_delta_tau_delta_tau.dat"

def dat_parse(name):
    with open(name) as data_fd:
        data_rows = np.array([[ float(j) for j in l.strip().split(" ")] for l in data_fd.readlines()])
        data = data_rows.transpose()
        ks = list(data[0])
        power = list(data[1])
        count = list(datan2])
        return (ks, power, count)

ks_i, power_i, count_i = dat_parse(data_file_ics)
ks_t, power_t, count_t = dat_parse(data_file_tau)



# LaTeX mode
matplotlib.rcParams["text.usetex"] = True
matplotlib.rcParams["savefig.dpi"] = 300
matplotlib.rcParams["font.size"] = 18


fig_scale = 12
fig_aspect = 4/3
fig, axis = plt.subplots(1, 1, figsize=(fig_scale, fig_scale / fig_aspect))

axis.plot(ks_i, power_i)

axis.plot(ks_t, power_t)


#axis.set_xlim(left=2e-1)

axis.legend()

axis.set_xlabel("$k/ (h/\mathrm{cMpc})$")
axis.set_ylabel("Power$/ \\mathrm{cMpc}^3$")
axis.set_title("Extracted "+ spec_type +" power spectrum")

axis.set_xscale("log", basex=10)
axis.set_yscale("log", basey=10)

plt.tight_layout()
plt.savefig("out/xcorr_" + file_type + "." + EXTENSION)


