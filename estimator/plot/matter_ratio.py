#!/usr/bin/env python3


import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

# load the data
#assert(len(sys.argv) > 1)
data_file_z2 = "out/80_snap18/xcorr_delta_matter_delta_matter.dat"
data_file_dm_z2 = "out/80_snap18/xcorr_delta_DM_delta_DM.dat"
data_file_ic = "out/80_ics/xcorr_delta_matter_delta_matter.dat"

# first argument specifies output format
EXTENSION = "pdf"
if len(sys.argv) > 1:
    EXTENSION = sys.argv[1]



def dat_parse(name):
    with open(name) as data_fd:
        data_rows = np.array([[ float(j) for j in l.strip().split(" ")] for l in data_fd.readlines()])
        data = data_rows.transpose()
        #ks = list(data[0])
        #power = list(data[1])
        #count = list(data[2])
        ks = data[0]
        power = data[1]
        count = data[2]
        return (ks, power, count)

ks_z2, power_z2, count_z2 = dat_parse(data_file_z2)
ks_dm_z2, power_dm_z2, count_dm_z2 = dat_parse(data_file_dm_z2)
ks_ic, power_ic, count_ic = dat_parse(data_file_ic)


print(power_z2[0:20] / (power_ic[0:20] * ((1 + 99) / (1 + 2)) ** 2 ))

# LaTeX mode
matplotlib.rcParams["text.usetex"] = True
matplotlib.rcParams["savefig.dpi"] = 300
matplotlib.rcParams["font.size"] = 18


fig_scale = 12
fig_aspect = 4/3
fig, axis = plt.subplots(1, 1, figsize=(fig_scale, fig_scale / fig_aspect))

axis.plot(ks_ic, power_z2 / (power_ic * ((1 + 99) / (1 + 2))**2), label="Total matter")
axis.plot(ks_ic, power_dm_z2 / (power_ic * ((1 + 99) / (1 + 2))**2), label="DM")

axis.plot(ks_ic, [1.1 for k in ks_ic], "--", linewidth=0.4, color="black")
axis.plot(ks_ic, [0.9 for k in ks_ic], "--", linewidth=0.4, color="black")

axis.legend()

axis.set_xlim(left=ks_ic[0])
axis.set_ylim([0.8, 1.2])


axis.set_xlabel("$k/ (h/\mathrm{cMpc})$")
axis.set_ylabel("Power$/ (\\mathrm{cMpc}/h)^3$")
axis.set_title("Matter power ratio ($z=2$ / initial)")

axis.set_xscale("log", basex=10)

plt.tight_layout()
plt.savefig("out/ratio_matter." + EXTENSION)



