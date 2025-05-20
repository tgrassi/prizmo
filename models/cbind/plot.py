import numpy as np
import matplotlib.pyplot as plt
from plot_commons import species_names

data = np.loadtxt('output.dat').T

spy = 365. * 24. * 3600.

lss = ['-', '--', '-.', ':']
# get colors in the color cycle
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

plt.plot(data[0] / spy, data[1], label="Tgas", lw=2, color='k')

for i, sp in enumerate(species_names):
    plt.plot(data[0] / spy, data[i + 2], label=sp, linestyle=lss[i // len(colors)])

plt.xscale('log')
plt.yscale('log')
plt.legend(ncols=4, fontsize=8, loc='lower left')
plt.tight_layout()
plt.ylim(1e-6, 2e4)
plt.show()