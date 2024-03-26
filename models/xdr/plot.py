import numpy as np
import matplotlib.pyplot as plt
from glob import glob
from scipy.interpolate import interp1d
from plot_commons import species_names

fig, axs = plt.subplots(2, 3, figsize=(14, 10))

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

def add_ncol_lines(ax):
  for ncol in [5e21, 1e22, 2e22]:
    ax.axvline(ncol, ls=":", alpha=0.3, color="k")

jxi_plot = 1e-4

# **********************************
ax = axs[0, 0]

# Jscale, Av, rad_Ncol_tot, Tgas, Tdust, x
data = np.loadtxt("fort.23").T

Ncols = np.array([5e21, 1e22, 2e22])
for js in np.unique(data[0]):
  idx = data[0] == js
  f = interp1d(data[2][idx], data[3][idx])
  ax.scatter(np.zeros_like(Ncols)+js, f(Ncols))

data_ref = np.loadtxt("ncol.dat").T
ax.plot(data_ref[0], data_ref[1])
ax.plot(data_ref[0], data_ref[2])
ax.plot(data_ref[0], data_ref[3])

ax.set_xscale("log")
ax.set_yscale("log")

# **********************************
ax = axs[1, 0]

data = np.loadtxt("fort.23").T

only = ["O", "O+", "O++", "O+++", "O++++", "E", "C", "C+", "C++", "C+++", "C++++", "H", "H+", "H2", "CO"]

jxi = np.unique(data[0])
this_jxi = jxi[np.argmin(np.abs(jxi - jxi_plot))]

lss = ["-", "--", ":"]
icount = 0
for i, sp in enumerate(species_names):
  if sp not in only:
    continue
  idx = data[0] == this_jxi
  ax.loglog(data[2][idx], data[5+i][idx], label=sp, ls=lss[icount // len(colors)])
  icount += 1

add_ncol_lines(ax)

ax.title.set_text(this_jxi)
ax.legend(loc="best", fontsize=6, ncol=2)
ax.set_ylim(bottom=1e-10)


# **********************************

ax = axs[0, 1]

# Jscale, Av, rad_Ncol_tot, Tgas, Tdust, x
idx = data[0] == this_jxi
ax.loglog(data[2][idx], data[3][idx], label="Tgas")
ax.loglog(data[2][idx], data[4][idx], label="Tdust")
ax.legend(loc="best", fontsize=6, ncol=2)

ax.title.set_text(this_jxi)
add_ncol_lines(ax)


# **********************
#  Jscale, Av, rad_Ncol_tot, Tgas, cools
labs_cool = ["atomic", "chem", "dust", "H2", "CO"]
labs_heat = ["photo", "phe", "CR", "H2diss"]
data_cool = np.loadtxt("fort.24").T
data_heat = np.loadtxt("fort.25").T

ax = axs[1, 1]

jxi = np.unique(data_cool[0])
this_jxi = jxi[np.argmin(np.abs(jxi - jxi_plot))]

idx = data_cool[0] == this_jxi
for i, sp in enumerate(labs_cool):
  ax.loglog(data_cool[2][idx], data_cool[4+i][idx], label=sp)

idx = data_heat[0] == this_jxi
for i, sp in enumerate(labs_heat):
  ax.loglog(data_heat[2][idx], data_heat[4+i][idx], label=sp, ls="--")

ax.title.set_text(this_jxi)

ax.legend(loc="best", fontsize=6, ncol=2)
ax.set_ylim(bottom=1e-30)
add_ncol_lines(ax)


# **********************************
ax = axs[0, 2]

#  Jscale, Av, rad_Ncol_tot, Tgas, cools
data = np.loadtxt("fort.28").T

only = ["O", "O+", "O++", "O+++", "O++++", "C", "C+", "C++", "C+++", "C++++", "H"]

jxi = np.unique(data[0])
this_jxi = jxi[np.argmin(np.abs(jxi - jxi_plot))]
idx = data[0] == this_jxi

lss = ["-", "--", ":"]
icount = 0
for i, sp in enumerate(species_names):
  if sp not in only:
    continue
  ax.loglog(data[2][idx], data[4+i][idx], label=sp, ls=lss[icount // len(colors)])
  icount += 1

add_ncol_lines(ax)

ax.title.set_text(this_jxi)
ax.legend(loc="best", fontsize=6, ncol=2)
ax.set_ylim(bottom=1e-28)

# **********************************
ax = axs[1, 2]

#  Jscale, Av, rad_Ncol_tot, Tgas, Tdust, x
data = np.loadtxt("fort.23").T

jxi = np.unique(data[0])
this_jxi = jxi[np.argmin(np.abs(jxi - jxi_plot))]
idx = data[0] == this_jxi


data_ions = np.vstack([data[5+i] * x.count("+") for i, x in enumerate(species_names)])
idx_e = species_names.index("E")
ax.semilogx(data[2][idx], np.sum(data_ions, axis=0)[idx] - data[5+idx_e][idx])
ax.axhline(0)
ax.set_yscale("symlog")

add_ncol_lines(ax)

ax.title.set_text(this_jxi)

# **********************************

for i in range(len(axs)):
    for j in range(len(axs[i])):
        axs[i, j].grid(ls="--", alpha=0.1, color="k")

plt.savefig("plot.png", bbox_inches="tight", dpi=120)
plt.show()


