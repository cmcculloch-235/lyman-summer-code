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

# input spectrum?
data_file_i = "planck1_matterpower_z=99_rightunits.dat"
with open(data_file_i) as data_fd:
    data_rows = np.array([[ float(j) for j in l.strip().split(" ")] for l in data_fd.readlines()])
    data = data_rows.transpose()
    # k in file: log_10 k/(h/Mpc) (?h?)
    ks_i = np.exp(np.log(10) * data[0])
    # power in file: log_10 (P(k) k^3 / (2pi^2)) the dimensionless pspec
    power_i = (2 * np.pi**2) * np.exp(np.log(10) * data[1]) / (ks_i**3)

interp_i = lambda x: np.interp(x, ks_i, power_i)


# spectrum from ics (Pylians)?
data_file_p = "Pk_CDM_z=99.000.dat"
with open(data_file_p) as data_fd:
    data_rows = np.array([[ float(j) for j in l.strip().split(" ")] for l in data_fd.readlines()])
    data = data_rows.transpose()
    # k in file: k/(h/Mpc) (¿h?)
    ks_p = data[0]
    # power in file: log_10 (P(k) k^3 / (2pi^2)) the dimensionless pspec
    power_p = data[1]

    count_p = data[4]


# LaTeX mode
matplotlib.rcParams["text.usetex"] = True
matplotlib.rcParams["savefig.dpi"] = 300
matplotlib.rcParams["font.size"] = 18


fig_scale = 12
fig_aspect = 4/3
fig, axis = plt.subplots(1, 1, figsize=(fig_scale, fig_scale / fig_aspect))



#axis.set_xlim(left=2e-1)
axis.plot(ks_i, power_i, color="blue", label="Input")

#axis.plot(ks_p, count_p / ks_p**2, label="count / $k^2$")



stdev = power_p / np.sqrt(count_p)
axis.plot(ks_p, power_p + stdev, linewidth=0.8, color="purple")
axis.plot(ks_p, power_p - stdev, linewidth=0.8, color="purple")
axis.plot(ks_p, power_p + 2 * stdev, linewidth=0.6, color="cyan")
axis.plot(ks_p, power_p - 2 * stdev, linewidth=0.6, color="cyan")
axis.plot(ks_p, power_p, color="green", label="Pylians")

# -ves
axis.plot(ks_p, -(power_p + stdev), "--", linewidth=0.8, color="purple")
axis.plot(ks_p, -(power_p - stdev), "--", linewidth=0.8, color="purple")
axis.plot(ks_p, -(power_p + 2 * stdev), "--", linewidth=0.6, color="cyan")
axis.plot(ks_p, -(power_p - 2 * stdev), "--", linewidth=0.6, color="cyan")
axis.plot(ks_p, -(power_p), "--", color="green")




axis.plot(ks_es, power_es, color="black", label="Estimator")

# also indicate where pylians samples are
# indicates that they are uniformly spaced
#axis.plot(ks_p, [1 for k in ks_p], ".")


axis.set_xlabel("$k/ (h/\mathrm{cMpc})$")
axis.set_ylabel("Power$/ \\mathrm{cMpc}^3$")
axis.legend()

axis.set_xlim(left=ks_es[0], right=ks_es[-1])
axis.set_ylim(bottom=1e-4)

axis.set_xscale("log", basex=10)
axis.set_yscale("log", basey=10)

plt.tight_layout()
plt.savefig("comparison.pdf")


# ratio

fig_scale = 12
fig_aspect = 4/3
fig, axis = plt.subplots(1, 1, figsize=(fig_scale, fig_scale / fig_aspect))

axis.plot(ks_es, power_es / interp_i(ks_es))
axis.set_xlabel("$k/ (h/\mathrm{cMpc})$")
axis.set_ylabel("estimated/input")

#axis.set_xlim(left=ks_es[0], right=ks_es[-1])
axis.set_xlim(left=ks_es[0], right=1)
axis.set_ylim(top=6)

axis.set_xscale("log", basex=10)
#axis.set_yscale("log", basey=10)

plt.tight_layout()
plt.savefig("comparison_ratio.pdf")


# fluctuation ratio

fig_scale = 12
fig_aspect = 4/3
fig, axis = plt.subplots(1, 1, figsize=(fig_scale, fig_scale / fig_aspect))

axis.plot(ks_es, power_es / interp_i(ks_es) - 1)
axis.set_xlabel("$k/ (h/\mathrm{cMpc})$")
axis.set_ylabel("estimated/input - 1")

axis.set_xlim(left=ks_es[0], right=ks_es[-1])

axis.set_xscale("log", basex=10)
#axis.set_yscale("log", basey=10)

plt.tight_layout()
plt.savefig("comparison_f_ratio.pdf")
