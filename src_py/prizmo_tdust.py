from prizmo_commons import clight, hplanck, amin, amax, kboltzmann, pmass, rho_bulk, pslope, print_title, plotOn, refInd_file
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy.interpolate import interp1d
from scipy.optimize import bisect, brentq
from bhmie import bhmie_qabs
from tqdm import tqdm
import os
import sys

# compute table of Tdust as function of Ea (absorption integral against Qabs), Tgas, and ngas
# from this compute also table cooling (should be multiplied by mu and by dust-to-gas mass ratio).
# it uses bisection of: emission(Td) = absorption - collisional(Tgas - Tdust)


# nt: number of tgas points
# ng: number of ngas points
# ne: number of Ea points
# total grid points = nt*ngne
def prepare(user_energy_eV, ne=100, nt=100, ng=50):

    print_title("dust temperature and cooling")

    fname_grid = "../runtime_data/tdust_grid.dat"
    if os.path.isfile(fname_grid):
        print("skipping, tdust grid file found", fname_grid)
        return

    # wls, qabs = load_qabs()
    loge, Cabs = load_eps()  # wls, qabs)
    Bint = make_Bint(1e1**loge, Cabs)

    out = "# energy/eV, Ea pre-integral\n"
    for e in user_energy_eV:
        out += "%.18e %.18e\n" % (e, 1e1**Cabs(np.log10(e)))
    fh = open("../runtime_data/tdust_Ea_preint.dat", "w")
    fh.write(out)
    fh.close()

    # temperature range for start bisection, K
    td_min = 1e-1
    td_max = 1e6

    # store tdust per Ea slice for plots
    td = np.zeros((nt, ng))
    cool = np.zeros((nt, ng))

    # grid points to loop
    ea_range = np.logspace(-5, 10, ne)
    tgas_range = np.logspace(0, 6, nt)
    ngas_range = np.logspace(-4, 25, ng)

    tdiff = np.zeros((ne, ng)) + 1e99

    out_tdust = ""  # string that will be saved to file with the final table
    out_cool = ""  # same but for cooling
    out_heat = ""  # same but for heating, to avoid negative numbers in the log
    for k, ea in enumerate(tqdm(ea_range)):
        for i, tgas in enumerate(tgas_range):
            for j, ngas in enumerate(ngas_range):
                fek = get_fek(tgas, ngas)
                if f(td_min, fek, Bint, ea, tgas) * f(td_max, fek, Bint, ea, tgas) > 0e0:
                    plot(td_min, td_max, fek, Bint, ea, tgas, ngas)
                    sys.exit("ERROR: bisection failed!")
                td[i, j] = brentq(f, td_min, td_max, args=(fek, Bint, ea, tgas), full_output=False)
                # print(rstat.iterations)
                cool[i, j] = get_cooling(fek, td[i, j], tgas)
                out_tdust += "%.18e %.18e %.18e %.18e\n" % (ea, tgas, ngas, td[i, j])
                out_cool += "%.18e %.18e %.18e %.18e\n" % (ea, tgas, ngas, max(cool[i, j], 1e-40))
                out_heat += "%.18e %.18e %.18e %.18e\n" % (ea, tgas, ngas, max(-cool[i, j], 1e-40))
                if i == nt // 2:
                    tdiff[k, j] = min(td[i, j] - tgas, tdiff[k, j])


        if plotOn:
            ngrid, tgrid = np.meshgrid(ngas_range, tgas_range)
            plt.pcolor(tgas_range, ngas_range, (td / tgrid).T, norm=matplotlib.colors.LogNorm(vmin=0.1, vmax=1e1), cmap="seismic")
            plt.show()
        else:
            plt.close()

        # plot Tdust at each Ea for debug
        if plotOn:
            plt.pcolor(tgas_range, ngas_range, cool.T, cmap="jet", norm=matplotlib.colors.LogNorm())
            # plt.pcolor(tgas_range, ngas_range, td.T, cmap="jet", norm=matplotlib.colors.LogNorm())
            plt.xscale("log")
            plt.yscale("log")
            plt.xlabel("T_gas")
            plt.ylabel("n_gas")
            plt.title("ea %e" % ea)
            plt.colorbar()
            plt.show()
        else:
            plt.close()

    plt.pcolor(np.log10(ea_range), np.log10(ngas_range), np.abs(tdiff.T), cmap="jet", norm=matplotlib.colors.LogNorm())
    plt.colorbar()
    if plotOn:
        plt.show()
    else:
        plt.close()

    # save Tdust results to file, (Ea, Tgas, ngas, Tdust)
    fh = open(fname_grid, "w")
    fh.write("# Ea, Tgas, ngas, Tdust\n")
    fh.write(out_tdust)
    fh.close()

    # save cooling results to file, (Ea, Tgas, ngas, cooling)
    # note: it should be multiplied by rho_dust and d2g in the main code
    fh = open("../runtime_data/dust_cooling_grid.dat", "w")
    fh.write("# Ea, Tgas, ngas, cooling/d2g/rho_d\n")
    fh.write(out_cool)
    fh.close()

    # save heating results to file, (Ea, Tgas, ngas, heating)
    # note: it should be multiplied by rho_dust and d2g in the main code
    fh = open("../runtime_data/dust_heating_grid.dat", "w")
    fh.write("# Ea, Tgas, ngas, heating/d2g/rho_d\n")
    fh.write(out_heat)
    fh.close()


# plot quantities as a function of Tdust for bisection debug
def plot(td_min, td_max, fek, Bint, ea, tgas, ngas):
    trange = np.logspace(np.log10(td_min), np.log10(td_max), 100)

    plt.semilogx(trange, [f(x, fek, Bint, ea, tgas) for x in trange], label="f")
    plt.yscale("symlog")
    plt.semilogx(trange, [em(x, Bint) for x in trange], label="em")
    plt.semilogx(trange, ea * np.ones_like(trange), label="ea")
    plt.semilogx(trange, [ek(x, tgas, fek) for x in trange], label="ek")
    plt.title("T=%.1e, ngas=%.1e" % (tgas, ngas))
    plt.legend(loc="best")
    if plotOn:
        plt.show()
    else:
        plt.close()


# bisection function per grain (hence it doesn't depend on the amount of dust, but only on their distribution)
def f(x, fek, Bint, ea, tgas):
    return em(x, Bint) - ea - ek(x, tgas, fek)


# collisional term, without (Tgas-Td)
def get_fek(tgas, ngas):
    p3 = pslope + 3.
    vgas = np.sqrt(8e0 * kboltzmann * tgas / 2 / pmass / np.pi)
    return 2e0 * vgas * kboltzmann * (amax**p3 - amin**p3) / p3 * ngas


# load Qabs from file for debug
def load_qabs():
    wls = []
    qabs = []
    for row in open("../data/dust_refractive_index/draine_Rv3.1.dat"):
        srow = row.strip()
        if srow == "" or srow.startswith("#"):
            continue
        arow = srow.split()
        wls.append(float(arow[0]))
        qabs.append(float(arow[4]))
    return clight * hplanck / np.array(wls) / 1e-4, np.array(qabs)


# load optical epsilon and compute Cabs, or load Cabs from file if already there
def load_eps():
    if os.path.isfile("Cabs.npy"):
        loge, Cabs = np.load("Cabs.npy", allow_pickle=True)
        return loge, Cabs

    # w(micron)  Re(eps-1)  Im(eps)    Re(m-1)    Im(m)
    #data_eps = np.loadtxt("../data/dust_refractive_index/eps_Sil.txt").T
    data_eps = np.loadtxt(refInd_file).T
    wl = data_eps[0] * 1e-4  # micron -> cm (NOTE: are already reversed)

    e_real = data_eps[1] + 1e0
    e_img = data_eps[2]

    nref = np.sqrt((np.sqrt(e_real**2 + e_img**2) + e_real) / 2e0)
    kref = np.sqrt((np.sqrt(e_real**2 + e_img**2) - e_real) / 2e0)

    # wl = wl[wl>2.5e-5]
    # grain size range
    arange = np.logspace(np.log10(amin), np.log10(amax), 100)

    fe = np.zeros_like(wl)
    # loop on wavelengths
    for i in tqdm(range(len(wl))):

        # convolve with dust distribution
        fa = np.zeros_like(arange)
        for j, a in enumerate(arange):
            fa[j] = a**2 * a**pslope * bhmie_qabs(2e0 * np.pi / wl[i] * a, complex(nref[i], kref[i]))
        # integrate along dust size
        fe[i] = np.trapz(fa, arange)

    # plot for debug
    plt.loglog(wl, fe)
    if plotOn:
        plt.show()
    else:
        plt.close()

    # check
    if fe.min() < 0e0:
        print("WARNING: negative qabs from bhmie!")

    # compute energy, erg
    energy = clight * hplanck / wl

    # obtain an interpolation function for Cabs
    loge = np.log10(energy)
    Cabs = interp1d(loge, np.log10(fe))
    np.save("Cabs.npy", np.array([loge, Cabs], dtype=object))

    return loge, Cabs


# emission factor
def em(x, Bint):
    return 4e0 * 1e1**Bint(np.log10(x)) / hplanck


# collisional term
def ek(x, tgas, fek):
    return (tgas - x) * fek


# compute integral of Qabs*Planck_function
def make_Bint(energy, Cabs):
    nt = 100
    loge = np.log10(energy)
    trange = np.logspace(-2, 6, nt)
    ft = np.zeros_like(trange)
    for i, t in enumerate(trange):
        ft[i] = np.trapz(B(energy, t) * 1e1 ** Cabs(loge), energy) + 1e-80
    print("4*Bint/h, (min/max)", 4 * ft.min() / hplanck, 4*ft.max() / hplanck)
    return interp1d(np.log10(trange), np.log10(ft))


# compute cooling after getting Tdust from bisection
def get_cooling(fek, td, tgas):
    p4 = pslope + 4.
    amass = 4. / 3. * np.pi * rho_bulk * (amax**p4 - amin**p4) / p4
    return np.pi * ek(td, tgas, fek) / amass


# planck function, erg/cm2
def B(energy, tbb):
    nu = energy / hplanck
    xexp = np.minimum(hplanck * nu / kboltzmann / tbb, 2e2)
    return 2e0 * hplanck * nu ** 3 / clight ** 2 / (np.exp(xexp) - 1e0)
