import numpy as np
import matplotlib.pyplot as plt
from plot_commons import species_names


data = np.loadtxt("fort.22").T

xdata = data[0]

for i, s in enumerate(species_names):
    plt.loglog(xdata, data[2+i], label=s)

plt.legend(loc="best")
plt.show()


