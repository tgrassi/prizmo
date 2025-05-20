import os.path
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from prizmo_commons import amin, amax, pslope, clight, hplanck, rho_bulk, erg2ev, print_title, fuv_energy1, fuv_energy2, plotOn, refInd_file
from prizmo_preprocess import preprocess
from bhmie import bhmie_qabs
from scipy.interpolate import interp1d
try:
    from scipy.integrate import trapz
except ImportError:
    from scipy.integrate import trapezoid as trapz

def prepare(user_energy):

    print_title("dust opacity")

    fname = "../runtime_data/kappa_dust.dat"
    if os.path.isfile(fname):
        fk = np.loadtxt(fname).T
        kappa = compute_kabs_integral(user_energy, fk)
        print("kabs integral 912-1100 AA:", kappa)
        return

    # w(micron)  Re(eps-1)  Im(eps)    Re(m-1)    Im(m)
    #data_eps = np.loadtxt("../data/dust_refractive_index/eps_Sil.txt").T
    data_eps = np.loadtxt(refInd_file).T
    wl = data_eps[0] * 1e-4  # micron -> cm (NOTE: are already reversed)

    # e_real = data_eps[1] + 1e0
    # e_img = data_eps[2]
    nref = data_eps[3] + 1e0
    kref = data_eps[4]

    # wl = wl[wl>2.5e-5]
    # grain size range
    arange = np.logspace(np.log10(amin), np.log10(amax), 100)

    f_nref = interp1d(wl[::-1], nref[::-1])
    f_kref = interp1d(wl[::-1], kref[::-1])

    p2 = pslope + 2e0
    p4 = pslope + 4e0

    anorm = 4. / 3. * np.pi * rho_bulk * (amax**p4 - amin**p4) / p4

    # loop on wavelengths
    fk = np.zeros_like(user_energy)
    for i, energy in enumerate(tqdm(user_energy)):
        wlen = clight * hplanck / energy
        # convolve with dust distribution
        fa = np.zeros_like(arange)
        for j, a in enumerate(arange):
            fa[j] = np.pi * a**p2 * bhmie_qabs(2e0 * np.pi / wlen * a, complex(f_nref(wlen), f_kref(wlen)))
        fk[i] = trapz(fa, arange) / anorm

    kappa = compute_kabs_integral(user_energy, fk)
    print("kabs integral 912-1100 AA:", kappa)

    plt.loglog(user_energy * erg2ev, fk)
    if plotOn:
        plt.show()
    else:
        plt.close()

    np.savetxt(fname, fk.T)


def compute_kabs_integral(energy, fk):
    cond = (energy >= fuv_energy1) & (energy <= fuv_energy2)
    kappaint = trapz(fk[cond], energy[cond]) / (fuv_energy2 - fuv_energy1)

    kabs_integral = "! kabs integral in FUV range\n"
    kabs_integral += "real*8,parameter::kabs_integral=%s" % ("%.18e" % kappaint).replace("e", "d")

    preprocess("../prizmo_commons.f90", {"KABS_INTEGRAL": kabs_integral})

    return kappaint
