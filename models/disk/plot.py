import numpy as np
import matplotlib.pyplot as plt
from glob import glob
from scipy.interpolate import interp1d
from plot_commons import species_names
from tqdm import tqdm
from sys import argv
import os, sys

istep = int(argv[1])

fname_png = "plot_%05d.png" % istep


if os.path.isfile(fname_png):
    sys.exit()

shading = "nearest"

n1 = 4  # rows
n2 = 7  # columns
fig, axs = plt.subplots(n1, n2, figsize=(16, 8), subplot_kw=dict(projection='polar'))

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

au2cm = 1.49598073e13

def tr(data):
    nx = len(np.unique(data[0]))
    data = [x[:len(x) - len(x) % nx] for x in data]
    return np.array([x.reshape((nx, -1), order="C") for x in data])

def loader(fname):
    wdata = np.loadtxt(fname).T
    wdata = tr(wdata)
    wdata[1] = np.log10(wdata[1] / au2cm + 1e0)
    # theta, r, ngas, tgas
    rmin = wdata[1].min()
    rmax = wdata[1].max()
    return rmin, rmax, wdata

_, _, data_cpu = loader("output_tcpu_%05d.dat" % istep)


# **********************************

_, _, data = loader("output_ngas_%05d.dat" % istep)


dax = axs[0, 0].pcolormesh(data[0], data[1], np.log10(data[2]), cmap="jet", shading=shading)
fig.colorbar(dax, ax = axs[0, 0])
axs[0, 0].title.set_text("ngas / cm-3")


# **********************************
_, _, data = loader("output_chem_%05d.dat" % istep)

vmax = max(np.log10(data[2].max()), 4.)
vmin = np.log10(5.)

dax = axs[1, 0].pcolormesh(data[0], data[1], np.log10(data[2]+1e-1), cmap="jet", shading=shading, vmin=vmin, vmax=vmax)
fig.colorbar(dax, ax = axs[1, 0])
axs[1, 0].title.set_text("Tgas / K")

dax = axs[2, 0].pcolormesh(data[0], data[1], np.log10(data[3]+1e-1), cmap="jet", shading=shading)
fig.colorbar(dax, ax = axs[2, 0])
axs[2, 0].title.set_text("Tdust / K")

ioff = 3

# **********************************
_, _, data = loader("output_cool_%05d.dat" % istep)
_, _, data_heat = loader("output_heat_%05d.dat" % istep)

diff_thermo = np.abs(np.sum(data_heat[3:], axis=0) + 1e-40) / np.abs(np.sum(data[3:], axis=0) + 1e-40) #) / (np.sum(data_cool[3:]) + 1e-60)

vmax = 0 #1e0 #max(data[3:].max(), data_cool[3:].max())
vmin = -4 # 1e-4

labs_cool = ["atomic", "chem", "dust", "H2", "CO"]
cool_tot = data[3:].sum(axis=0)

ff = np.ones_like(data_cpu[4])
ff[(data_cpu[4] / data_cpu[4].max()) > 1e-1] = 1e0
ff[(data_cpu[4] / data_cpu[4].max()) <= 1e-1] = -1e0


for i, ll in enumerate(labs_cool):
    ii = (ioff + i) % n1
    jj = (ioff + i) // n1
    dax = axs[ii, jj].pcolormesh(data[0], data[1], np.log10(data[i+3]*ff / cool_tot), vmax=vmax, vmin=vmin, shading=shading)
    fig.colorbar(dax, ax = axs[ii, jj])
    axs[ii, jj].title.set_text("COOL: %s / tot_cool" % ll)

ioff += len(labs_cool)

# **********************************
_, _, data = loader("output_heat_%05d.dat" % istep)

labs_heat = ["photo", "phe", "CR", "H2diss"]
heat_tot = data[3:].sum(axis=0)

for i, ll in enumerate(labs_heat):
    ii = (ioff + i) % n1
    jj = (ioff + i) // n1
    dax = axs[ii, jj].pcolormesh(data[0], data[1], np.log10(data[i+3]*ff / heat_tot), vmax=vmax, vmin=vmin, cmap="inferno", shading=shading)
    #cs = axs[ii, jj].contour(data[0], data[1], np.log10(diff_thermo), levels=[-4, -2, 0, 2, 4])
    #axs[ii, jj].clabel(cs, inline=True, fontsize=10)
    fig.colorbar(dax, ax = axs[ii, jj])
    axs[ii, jj].title.set_text("HEAT: %s / tot_heat" % ll)

ioff += len(labs_heat)

ii = (ioff + 0) % n1
jj = (ioff + 0) // n1
dax = axs[ii, jj].pcolormesh(data[0], data[1], np.log10(diff_thermo*ff), cmap="coolwarm", vmin=-5, vmax=5, shading=shading)
fig.colorbar(dax, ax = axs[ii, jj])
axs[ii, jj].title.set_text("|heat / cool|")


ioff += 1

# **********************************

only = ["H", "E", "H+", "H2", "CO", "C", "C+", "E", "O", "O+", "HCO+", "H2O", "CO_DUST", "H2O_DUST"]
_, _, data = loader("output_chem_%05d.dat" % istep)
for i, ll in enumerate(only):
    idx = species_names.index(ll)
    ii = (ioff + i) % n1
    jj = (ioff + i) // n1
    dax = axs[ii, jj].pcolormesh(data[0], data[1], np.log10(data[idx+4] / data[4:].sum(axis=0)),
                                 cmap="viridis", vmin=-8, vmax=0, shading=shading)
    #cs = axs[ii, jj].contour(data[0], data[1], np.log10(diff_thermo), levels=[-4, -2, 0, 2, 4])
    #axs[ii, jj].clabel(cs, inline=True, fontsize=10)
    fig.colorbar(dax, ax = axs[ii, jj])
    axs[ii, jj].title.set_text("frac %s" % ll)

ioff += len(only)


# **********************************
_, _, data = loader("output_tcpu_%05d.dat" % istep)
ii = (ioff + 0) % n1
jj = (ioff + 0) // n1
dax = axs[ii, jj].pcolormesh(data[0], data[1], np.log10(data[4]*ff / data[4].max()), vmin=-3, vmax=0, cmap="jet", shading=shading)
fig.colorbar(dax, ax = axs[ii, jj])
axs[ii, jj].title.set_text("cpu time / max")


# **********************************

for i in range(len(axs)):
    for j in range(len(axs[i])):
        ax = axs[i, j]
        #axs[i, j].grid(ls="--", alpha=0.1, color="k")
        ax.set_thetamin(0)
        ax.set_thetamax(45)
        #ax.set_rmin(rmin)
        #ax.set_rmax(0.5)

plt.tight_layout()
plt.savefig(fname_png, bbox_inches="tight", dpi=120)
#plt.show()


