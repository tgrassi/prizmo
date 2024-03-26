import numpy as np
import matplotlib.pyplot as plt
from glob import glob
from scipy.interpolate import interp1d
from plot_commons import species_names

fig, axs = plt.subplots(2, 3, figsize=(14, 10))

ref = "v1"
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
lss = ["-", "--", ":"]

# **********************
def plot_species(ax, ymin=1e-2, ymax=None):
    labs_chem = species_names

    labs = ["Av", "t", "Tgas", "Tdust"] + labs_chem
    data = np.loadtxt("fort.23").T

    labs_ref = ["r", "Av", "H", "H2", "C+", "C", "CO", "O", "O2", "CH", "OH", "E", "He+", "H+", "H3+"]

    xpos = labs.index("Av")
    vmax = -1e99
    icount = 0
    for i, l in enumerate(labs_chem):
        if l not in labs_ref:
            continue
        ypos = labs.index(l)
        ax.loglog(data[xpos], data[ypos], label=l, color=colors[icount % len(colors)], ls=lss[icount // len(colors)])
        vmax = max(vmax, data[ypos].max())
        icount += 1

    icount = 0
    for i, l in enumerate(labs_chem):
        if l not in labs_ref:
            continue
        sp_min = np.zeros_like(data[xpos]) + 1e99
        sp_max = np.zeros_like(data[xpos]) - 1e99
        for fname in glob("benchmark/density_*_%s.dat" % ref):
            data_ref = np.loadtxt(fname).T
            ypos_ref = labs_ref.index(l)
            f = interp1d(data_ref[1], data_ref[ypos_ref], fill_value="extrapolate")
            sp_min = np.minimum(sp_min, f(data[xpos]))
            sp_max = np.maximum(sp_max, f(data[xpos]))
        ax.fill_between(data[xpos], sp_min, sp_max, color=colors[icount % len(colors)], alpha=0.2, label=l)
        icount += 1


    ax.legend(loc="best", ncol=4, fontsize=4)
    if ymax is None:
        ymax = vmax * 10
    ax.set_ylim(ymin, ymax)

plot_species(axs[0, 0])
plot_species(axs[1, 1], ymin=1e-6, ymax=2e-1)

# **********************
labs = ["Av", "t", "Tgas", "Tdust"]

data = np.loadtxt("fort.23").T
xpos = labs.index("Av")

ypos = labs.index("Tgas")
axs[0, 1].loglog(data[xpos], data[ypos], label="Tgas")

ypos = labs.index("Tdust")
axs[0, 1].loglog(data[xpos], data[ypos], label="Tdust")

tgas_ref_min = np.zeros_like(data[xpos]) + 1e99
tgas_ref_max = np.zeros_like(data[xpos]) - 1e99
tdust_ref_min = np.zeros_like(data[xpos]) + 1e99
tdust_ref_max = np.zeros_like(data[xpos]) - 1e99
# z Av T Tdust
for fname in glob("benchmark/temperature_*_%s.dat" % ref):
    if "sternberg" in fname:
      continue
    data_ref = np.loadtxt(fname).T
    f = interp1d(data_ref[1], data_ref[2], fill_value="extrapolate")
    tgas_ref_min = np.minimum(tgas_ref_min, f(data[xpos]))
    tgas_ref_max = np.maximum(tgas_ref_max, f(data[xpos]))
    f = interp1d(data_ref[1], data_ref[3], fill_value="extrapolate")
    tdust_ref_min = np.minimum(tdust_ref_min, f(data[xpos]))
    tdust_ref_max = np.maximum(tdust_ref_max, f(data[xpos]))
axs[0, 1].fill_between(data[xpos], tgas_ref_min, tgas_ref_max, alpha=0.2)
axs[0, 1].fill_between(data[xpos], tdust_ref_min, tdust_ref_max, alpha=0.2)


axs[0, 1].legend(loc="best")
axs[0, 1].set_xlim(left=data[xpos].min()*0.3)

# **********************
labs_cool = ["Av", "t", "Tgas", "atomic", "chem", "dust", "H2", "CO"]
labs_heat = ["Av", "t", "Tgas", "photo", "phe", "CR", "H2diss"]
data_cool = np.loadtxt("fort.24").T
data_heat = np.loadtxt("fort.25").T

for l in labs_cool[3:]:
    ypos = labs_cool.index(l)
    axs[1, 0].loglog(data_cool[xpos], data_cool[ypos], label=l)
    if data_cool[ypos].min() < 0e0:
        axs[1, 0].loglog(data_cool[xpos], -data_cool[ypos], label=l, ls="--")

#axs[1, 0].loglog(data_cool[xpos], np.sum(data_cool[3:], axis=0), label="tot", color="k")


for l in labs_heat[3:]:
    ypos = labs_heat.index(l)
    axs[1, 0].loglog(data_heat[xpos], data_heat[ypos], label=l, ls="--")

#axs[1, 0].loglog(data_cool[xpos], np.sum(data_heat[3:], axis=0), label="tot", color="k", ls="--")


for ii, lref in zip([7, 10, 11, -1], ["phe", "H2vib", "gg_cool", "emission"]):
    #   z(cm)      A_V      CII158     OI63    OI146     CI610     CI370    dust_pe  totalheat totalcool   h2_vib   gg_cool   gg_heat
    cmin = np.zeros_like(data_cool[xpos]) + 1e99
    cmax = np.zeros_like(data_cool[xpos]) - 1e99
    for fname in glob("benchmark/heatcool_*_%s.dat" % ref):
        data_ref = np.loadtxt(fname).T
        if ii >= len(data_ref):
            continue
        if ii == -1:
            ydata_ref = np.sum(data_ref[2:7], axis=0)
        else:
            ydata_ref = data_ref[ii]

        f = interp1d(data_ref[1], ydata_ref, fill_value="extrapolate")
        fdata = np.abs(f(data_cool[xpos]))
        # axs[1, 0].loglog(data_cool[xpos], fdata, ls=":")
        if fdata.max() > 1e-30:
            cmin = np.minimum(cmin, fdata)
            cmax = np.maximum(cmax, fdata)

    axs[1, 0].fill_between(data_cool[xpos], cmin, cmax, color=colors[ii % len(colors)], alpha=0.2, label=lref)


axs[1, 0].legend(loc="best", ncol=3, fontsize=6)
#axs[1, 0].set_xlim(left=1e-6)
axs[1, 0].set_ylim(bottom=1e-28)

# **********************
def plot_cool_function(ax, fname_cool, fname_heat):
    labs_cool = ["Tgas", "atomic", "chem", "dust", "H2", "CO"]
    labs_heat = ["Tgas", "photo", "phe", "CR", "H2diss"]
    data_cool = np.loadtxt(fname_cool).T
    data_heats = np.loadtxt(fname_heat).T

    for l in labs_cool[1:]:
        ypos = labs_cool.index(l)
        ax.loglog(data_cool[0], data_cool[ypos], label=l)
        if data_cool[ypos].min() < 0e0:
            ax.loglog(data_cool[0], -data_cool[ypos], label=l, ls="--")

    ax.loglog(data_cool[0], np.sum(data_cool[1:], axis=0), label="tot", color="k")

    for l in labs_heat[1:]:
        ypos = labs_heat.index(l)
        ax.loglog(data_heats[0], data_heats[ypos], label=l, ls="--")

    ax.loglog(data_heats[0], np.sum(data_heats[1:], axis=0), label="tot", color="k", ls="--")


    ax.legend(loc="best", ncol=2, fontsize=8)
    ax.set_ylim(bottom=1e-28)

plot_cool_function(axs[1, 2], "function_cooling_av1m4.dat", "function_heating_av1m4.dat")
#plot_cool_function(axs[1, 2], "function_cooling_av1m0.dat", "function_heating_av1m0.dat")

# **********************************

data = np.loadtxt("fort.27").T
labs_photo = ["H2", "CO", "C"]

for i, ll in enumerate(labs_photo):
    axs[0, 2].loglog(data[0], data[3+i], label=ll)


#    z(cm)   A_V_eff      H2        CO       C 
for i, ll in enumerate(["H2", "CO", "C"]):
    cmin = np.zeros_like(data[0]) + 1e99
    cmax = np.zeros_like(data[0]) - 1e99
    #axs[0, 2].loglog(data[0], data[2+i], label=ll)
    for fname in glob("benchmark/photo_*_%s.dat" % ref):
        data_ref = np.loadtxt(fname).T
        f = interp1d(data_ref[1], data_ref[2+i], fill_value="extrapolate")
        fdata = f(data[0])
        cmin = np.minimum(cmin, fdata)
        cmax = np.maximum(cmax, fdata)
    print(ll, 0.5 * (cmax[5] + cmin[5]) / data[3+labs_photo.index(ll)][5])
    axs[0, 2].fill_between(data[0], cmin, cmax, color=colors[i % len(colors)], alpha=0.2, label=ll)


axs[0, 2].legend(loc="best", ncol=2, fontsize=8)
axs[0, 2].set_ylim(1e-14, 1e-8)


# ********************************
def aaaa():
    data = np.loadtxt("fort.28").T
    labs_acool = ["H", "C", "O", "C+", "O+"]

    axs[0, 2].cla()

    for i, ll in enumerate(labs_acool):
        axs[0, 2].loglog(data[0], data[3+i], label=ll)

    axs[0, 2].legend(loc="best", ncol=2, fontsize=8)
    axs[0, 2].set_ylim(1e-28, data[3:].max()*2)


# **********************************

for i in range(len(axs)):
    for j in range(len(axs[i])):
        axs[i, j].grid(ls="--", alpha=0.1, color="k")

plt.savefig("plot.png", bbox_inches="tight", dpi=120)
# plt.show()


