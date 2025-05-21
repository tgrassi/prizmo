import numpy as np
import matplotlib.pyplot as plt
from plot_commons import species_names

data = np.loadtxt('output.dat').T

spy = 365. * 24. * 3600.

# polar plot
fig, axs = plt.subplots(3, 4, subplot_kw={'projection': 'polar'}, figsize=(12, 10))

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
fig.colorbar(pc, ax=ax, fraction=0.03)
ax.set_title('ngas / cm-3')

# temperature
ax = axs[0, 1]
dd = np.log10(data[4])
pc = ax.pcolor(np.pi / 2. - theta_unique, r_unique, dd, shading='nearest', cmap='jet')
fig.colorbar(pc, ax=ax, fraction=0.03)
ax.set_title('Tgas / K')

# dust temperature
ax = axs[0, 2]
dd = np.log10(data[5])
pc = ax.pcolor(np.pi / 2. - theta_unique, r_unique, dd, shading='nearest', cmap='jet')
fig.colorbar(pc, ax=ax, fraction=0.03)
ax.set_title('Tdust / K')


for i, sp in enumerate(["H2", "CO", "C+", "O", "H", "H2O", "E", "CO_DUST", "H2O_DUST"]):
    idx = species_names.index(sp) + 6
    ax = axs.flatten()[i + 3]
    dd = np.log10(data[idx] / data[3])
    pc = ax.pcolor(np.pi / 2. - theta_unique, r_unique, dd, shading='nearest', cmap='jet', vmin=-10, vmax=0)
    fig.colorbar(pc, ax=ax, fraction=0.03)
    ax.set_title('n_%s / ngas' % sp)

for ax in axs.flatten():
    ax.set_thetamin(0)
    ax.set_thetamax(45)
    ax.set_xlabel('r [AU]')

plt.subplots_adjust(hspace=0.05, wspace=0.08, left=0.05, right=0.95, top=0.95, bottom=0.05)
#plt.tight_layout()
plt.show()