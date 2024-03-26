import shutil
import numpy as np


def prepare():
    shutil.copyfile("../data/CO_cooling/cooling.dat", "../runtime_data/CO_cooling.dat")

    # data = np.loadtxt("../runtime_data/CO_cooling.dat").T
