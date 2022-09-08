#!/usr/bin/env python3


import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

field_fnames = ["delta_H", "delta_H1", "delta_DM", "delta_matter", "delta_tau", "delta_flux", "delta_matter_et"]
field_lnames = ["$\\delta_{H}$", "$\\delta_{HI}$", "$\\delta_{DM}$", "$\\delta_{m}$",
        "$\\delta_{\\tau}$", "$\\delta_{F}$", "$\\delta_{m}$ (ET)"]
for i in range(0, len(field_fnames)):
    j = i - 1
    while j + 1 < len(field_fnames):
        j += 1

        file_type=field_fnames[i] + "_" + field_fnames[j]
        spec_type=field_lnames[i] + "--" + field_lnames[j]
        # load the data
        #assert(len(sys.argv) > 1)
        data_file_z2 = "out/80_snap18_new/xcorr_" + file_type + ".dat"
        data_file_z10 = "out/80_ics/xcorr_" + file_type + ".dat"

        # first argument specifies output format
        EXTENSION = "pdf"
        if len(sys.argv) > 1:
            EXTENSION = sys.argv[1]

        with open(data_file_z2) as data_fd:
            data_rows = np.array([[ float(j) for j in l.strip().split(" ")] for l in data_fd.readlines()])
            data = data_rows.transpose()
            ks_z2 = list(data[0])
            power_z2 = list(data[1])
            count_z2 = list(data[2])

        with open(data_file_z10) as data_fd:
            data_rows = np.array([[ float(j) for j in l.strip().split(" ")] for l in data_fd.readlines()])
            data = data_rows.transpose()
            ks_z10 = list(data[0])
            power_z10 = list(data[1])
            count_z10 = list(data[2])

        # assume the data fit together
        while 0 in count_z2:
            q = count_z2.index(0)
            del ks_z2[q]
            del power_z2[q]
            del count_z2[q]
            del ks_z10[q]
            del power_z10[q]
            del count_z10[q]
        ks_z2 = np.array(ks_z2)
        power_z2 = np.array(power_z2)
        count_z2 = np.array(count_z2)
        ks_z10 = np.array(ks_z10)
        power_z10 = np.array(power_z10)
        count_z10 = np.array(count_z10)


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

        # -ves
        axis.plot(ks_z2, -(power_z2 + stdev), "--", linewidth=0.8, color="purple")
        axis.plot(ks_z2, -(power_z2 - stdev), "--", linewidth=0.8, color="purple")
        axis.plot(ks_z2, -(power_z2 + 2 * stdev), "--", linewidth=0.6, color="green")
        axis.plot(ks_z2, -(power_z2 - 2 * stdev), "--", linewidth=0.6, color="green")
        axis.plot(ks_z2, -(power_z2), "--", color="blue")

        stdev = power_z10 / np.sqrt(count_z10)
        axis.plot(ks_z10, power_z10 + stdev, linewidth=0.8, color="red")
        axis.plot(ks_z10, power_z10 - stdev, linewidth=0.8, color="red")
        axis.plot(ks_z10, power_z10 + 2 * stdev, linewidth=0.6, color="yellow")
        axis.plot(ks_z10, power_z10 - 2 * stdev, linewidth=0.6, color="yellow")
        axis.plot(ks_z10, power_z10, color="black", label="Initial")

        # -ves
        stdev = power_z10 / np.sqrt(count_z10)
        axis.plot(ks_z10, -(power_z10 + stdev), "--", linewidth=0.8, color="red")
        axis.plot(ks_z10, -(power_z10 - stdev), "--", linewidth=0.8, color="red")
        axis.plot(ks_z10, -(power_z10 + 2 * stdev), "--", linewidth=0.6, color="yellow")
        axis.plot(ks_z10, -(power_z10 - 2 * stdev), "--", linewidth=0.6, color="yellow")
        axis.plot(ks_z10, -(power_z10), "--", color="black")

        #axis.set_xlim(left=2e-1)

        axis.legend()

        axis.set_xlabel("$k/ (h/\mathrm{cMpc})$")
        axis.set_ylabel("Power$/ \\mathrm{cMpc}^3$")
        axis.set_title("Extracted "+ spec_type +" power spectrum")

        axis.set_xscale("log", basex=10)
        axis.set_yscale("log", basey=10)

        plt.tight_layout()
        plt.savefig("out/xcorr_" + file_type + "." + EXTENSION)


# matter / matter_et ratio plot
file_type="delta_matter_et_delta_matter"
spec_type="matter (et) / matter"
# load the data
data_file_1 = "out/80_snap18_new/xcorr_delta_matter_et_delta_matter_et.dat"
data_file_2 = "out/80_snap18_new/xcorr_delta_matter_delta_matter.dat"


with open(data_file_1) as data_fd:
    data_rows = np.array([[ float(j) for j in l.strip().split(" ")] for l in data_fd.readlines()])
    data = data_rows.transpose()
    ks_1 = list(data[0])
    power_1 = list(data[1])
    count_1 = list(data[2])

with open(data_file_2) as data_fd:
    data_rows = np.array([[ float(j) for j in l.strip().split(" ")] for l in data_fd.readlines()])
    data = data_rows.transpose()
    ks_2 = list(data[0])
    power_2 = list(data[1])
    count_2 = list(data[2])

while 0 in count_1:
    i = count_1.index(0)
    del ks_1[i]
    del power_1[i]
    del count_1[i]
    del ks_2[i]
    del power_2[i]
    del count_2[i]
ks_1 = np.array(ks_1)
power_1 = np.array(power_1)
count_1 = np.array(count_1)
ks_2 = np.array(ks_2)
power_2 = np.array(power_2)
count_2 = np.array(count_2)

# LaTeX mode
matplotlib.rcParams["text.usetex"] = True
matplotlib.rcParams["savefig.dpi"] = 300
matplotlib.rcParams["font.size"] = 18


fig_scale = 12
fig_aspect = 4/3
fig, axis = plt.subplots(1, 1, figsize=(fig_scale, fig_scale / fig_aspect))


axis.plot(ks_z2, power_1 / power_2, color="green")
axis.plot(ks_z2, [1 for k in ks_z2], "--", color="black")

#axis.set_xlim(left=2e-1)


axis.set_xlabel("$k/ (h/\mathrm{cMpc})$")
axis.set_ylabel("Power$/ \\mathrm{cMpc}^3$")
axis.set_title("Extracted "+ spec_type + " ratio")

axis.set_xscale("log", basex=10)
axis.set_yscale("log", basey=10)

plt.tight_layout()
plt.savefig("out/ratio_" + file_type + "." + EXTENSION)


# ratio plot
file_type="delta_matter_delta_matter"
spec_type="matter spectrum"
# load the data
#assert(len(sys.argv) > 1)
data_file_z2 = "out/80_snap18_new/xcorr_" + file_type + ".dat"
data_file_z10 = "out/80_snap0/xcorr_" + file_type + ".dat"


with open(data_file_z2) as data_fd:
    data_rows = np.array([[ float(j) for j in l.strip().split(" ")] for l in data_fd.readlines()])
    data = data_rows.transpose()
    ks_z2 = list(data[0])
    power_z2 = list(data[1])
    count_z2 = list(data[2])

with open(data_file_z10) as data_fd:
    data_rows = np.array([[ float(j) for j in l.strip().split(" ")] for l in data_fd.readlines()])
    data = data_rows.transpose()
    ks_z10 = list(data[0])
    power_z10 = list(data[1])
    count_z10 = list(data[2])

while 0 in count_z2:
    i = count_z2.index(0)
    del ks_z2[i]
    del power_z2[i]
    del count_z2[i]
    del ks_z10[i]
    del power_z10[i]
    del count_z10[i]
ks_z2 = np.array(ks_z2)
power_z2 = np.array(power_z2)
count_z2 = np.array(count_z2)
ks_z10 = np.array(ks_z10)
power_z10 = np.array(power_z10)
count_z10 = np.array(count_z10)

# LaTeX mode
matplotlib.rcParams["text.usetex"] = True
matplotlib.rcParams["savefig.dpi"] = 300
matplotlib.rcParams["font.size"] = 18


fig_scale = 12
fig_aspect = 4/3
fig, axis = plt.subplots(1, 1, figsize=(fig_scale, fig_scale / fig_aspect))


axis.plot(ks_z2, power_z2/power_z10, color="green")
axis.plot(ks_z2, [((1 + 10) / (1 + 2))**2 for k in ks_z2], "--", color="black")

#axis.set_xlim(left=2e-1)


axis.set_xlabel("$k/ (h/\mathrm{cMpc})$")
axis.set_ylabel("Power$/ \\mathrm{cMpc}^3$")
axis.set_title("Extracted "+ spec_type + " ratio")

axis.set_xscale("log", basex=10)
axis.set_yscale("log", basey=10)

plt.tight_layout()
plt.savefig("out/ratio_" + file_type + "." + EXTENSION)
