import numpy as np
import matplotlib.pyplot as plt
from plot_commons import species_names

data = np.loadtxt('output.dat').T

spy = 365. * 24. * 3600.

# lss = ['-', '--', '-.', ':']
# # get colors in the color cycle
# colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

# polar plot
fig, axs = plt.subplots(3, 3, subplot_kw={'projection': 'polar'})

au2cm = 1.496e13

r_unique = np.unique(data[0]) / au2cm
theta_unique = np.unique(data[1])

nr = len(r_unique)
ntheta = len(theta_unique)

data = [x.reshape((ntheta, nr), order='C').T for x in data]

# density
ax = axs[0, 0]
dd = np.log10(data[3])
pc = ax.pcolor(np.pi / 2. - theta_unique, r_unique, dd, shading='nearest', cmap='jet')
fig.colorbar(pc, ax=ax)
ax.set_title('ngas')

# temperature
ax = axs[0, 1]
dd = np.log10(data[4])
pc = ax.pcolor(np.pi / 2. - theta_unique, r_unique, dd, shading='nearest', cmap='jet')
fig.colorbar(pc, ax=ax)
ax.set_title('Tgas')

for i, sp in enumerate(["H2", "CO", "C+", "O", "H", "H2O", "E"]):
    idx = species_names.index(sp) + 5
    ax = axs.flatten()[i + 2]
    dd = np.log10(data[idx] / data[3])
    pc = ax.pcolor(np.pi / 2. - theta_unique, r_unique, dd, shading='nearest', cmap='jet')
    fig.colorbar(pc, ax=ax)
    ax.set_title(sp)

for ax in axs.flatten():
    ax.set_thetamin(0)
    ax.set_thetamax(45)

#ax.set_rscale('log')


#plt.rscale('log')
#plt.yscale('log')
#plt.legend(ncols=4, fontsize=8, loc='lower left')
plt.tight_layout()
#plt.ylim(1e-6, 2e4)
plt.show()