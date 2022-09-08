#!/usr/bin/env python3


import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np


# NOTE: reference data in Mpc^3/h; estimator data in Mpc^3
# data from estimator
data_file_es = "../80_ics/xcorr_delta_DM_delta_DM.dat"
with open(data_file_es) as data_fd:
    data_rows = np.array([[ float(j) for j in l.strip().split(" ")] for l in data_fd.readlines()])
    data = data_rows.transpose()
    ks_es = list(data[0])
    power_es = list(data[1])
    count_es = list(data[2])

while 0 in count_es:
    i = count_es.index(0)
    del ks_es[i]
    del power_es[i]
    del count_es[i]
ks_es = np.array(ks_es)
power_es = np.array(power_es)
count_es = np.array(count_es)


data_file_p = "../80_ics/posk_delta_DM_delta_DM.dat"
with open(data_file_p) as data_fd:
    data_rows = np.array([[ float(j) for j in l.strip().split(" ")] for l in data_fd.readlines()])
    data = data_rows.transpose()
    ks_p = list(data[0])
    power_p = list(data[1])
    count_p = list(data[2])

while 0 in count_p:
    i = count_p.index(0)
    del ks_p[i]
    del power_p[i]
    del count_p[i]
ks_p = np.array(ks_p)
power_p = np.array(power_p)
count_p = np.array(count_p)

# LaTeX mode
matplotlib.rcParams["text.usetex"] = True
matplotlib.rcParams["savefig.dpi"] = 300
matplotlib.rcParams["font.size"] = 18


fig_scale = 12
fig_aspect = 4/3
fig, axis = plt.subplots(1, 1, figsize=(fig_scale, fig_scale / fig_aspect))




#axis.plot(ks_p, count_p / ks_p**2, label="count / $k^2$")






axis.plot(ks_es, power_es, color="black", linewidth=3, label="Estimator")
axis.plot(ks_p, power_p, color="red", linewidth=1, label="Pk")


axis.set_xlabel("$k/ (h/\mathrm{cMpc})$")
axis.set_ylabel("Power$/ \\mathrm{cMpc}^3$")
axis.legend()

axis.set_xlim(left=ks_es[0], right=ks_es[-1])
axis.set_ylim(bottom=1e-4)

axis.set_xscale("log", basex=10)
axis.set_yscale("log", basey=10)

plt.tight_layout()
plt.savefig("comparison.pdf")


