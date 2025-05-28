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

# cpu time
ax = axs[0, 3]
dd = np.log10(data[6])
pc = ax.pcolor(np.pi / 2. - theta_unique, r_unique, dd, shading='nearest', cmap='jet')
fig.colorbar(pc, ax=ax, fraction=0.03)
ax.set_title('CPU / s')
#ax.set_rlim(0, 20)


for i, sp in enumerate(["H2", "CO", "C+", "H", "H2O", "E", "CO_DUST", "H2O_DUST"]):
    idx = species_names.index(sp) + 7
    ax = axs.flatten()[i + 4]
    dd = np.log10((data[idx] + 1e-20) / data[3])
    pc = ax.pcolor(np.pi / 2. - theta_unique, r_unique, dd, shading='nearest', cmap='jet', vmin=-10, vmax=1)
    fig.colorbar(pc, ax=ax, fraction=0.03)
    ax.set_title('n_%s / ngas' % sp)

for ax in axs.flatten():
    ax.set_thetamin(0)
    ax.set_thetamax(45)
    ax.set_xlabel('r [AU]')

plt.subplots_adjust(hspace=0.05, wspace=0.08, left=0.05, right=0.95, top=0.95, bottom=0.05)
#plt.tight_layout()
plt.savefig("disk_maps.png", dpi=300)
print("Saved disk_maps.png")

plt.close('all')

# ******************************
fig, axs = plt.subplots(3, 4, figsize=(12, 10))

step = 5

# density
ax = axs[0, 0]
for th in range(ntheta-1, 0, -step):
    ax.plot(r_unique, data[3][:, th], label='%.1f' % (90.-np.degrees(theta_unique[th])))
ax.legend(loc='upper right')
ax.set_title('ngas / cm-3')


# temperature
ax = axs[0, 1]
for th in range(ntheta-1, 0, -step):
    ax.plot(r_unique, data[4][:, th])
ax.set_title('Tgas / K')

# dust temperature
ax = axs[0, 2]
for th in range(ntheta-1, 0, -step):
    ax.plot(r_unique, data[5][:, th])
ax.set_title('Tdust / K')

# cpu time
ax = axs[0, 3]
for th in range(ntheta-1, 0, -step):
    ax.plot(r_unique, data[6][:, th])
ax.set_title('CPU / s')

for i, sp in enumerate(["H2", "CO", "C+", "H", "H2O", "E", "CO_DUST", "H2O_DUST"]):
    idx = species_names.index(sp) + 7
    ax = axs.flatten()[i + 4]
    dd = data[idx] / data[3]
    for th in range(ntheta-1, 0, -step):
        ax.plot(r_unique, dd[:, th])
    ax.set_ylim(1e-10, 1e1)
    ax.set_title('n_%s / ngas' % sp)


for ax in axs.flatten():
    ax.set_yscale('log')

#plt.subplots_adjust(hspace=0.1, wspace=0.15, left=0.05, right=0.95, top=0.95, bottom=0.05)
plt.tight_layout()

plt.savefig("disk_radial.png", dpi=300)
print("Saved disk_radial.png")