import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d, interp2d
from prizmo_commons import hplanck, echarge2, wpot, clight, amin, amax, pslope, rho_bulk, \
    ev2erg, Zmin, Zmax, kboltzmann, pmass, emass, print_title, plotOn, refInd_file
from bhmie import bhmie_qabs
from tqdm import tqdm


# this routine produces three files
# 1) jptot.dat: contains the photodisorbtion and photododetachment rates as a function of the energy range
#  provided by the user code, columns are: Z, Jptot
#  Jptot needs to be multiplied by rho_gas * dust_gas_ratio * integral(Jptot * Jflux, energy)
# 2) jion.dat: contains the H+ + grain recombination rates as a function of Z grain and energy,
#  columns are: Z, tgas, Jion
#  Jion needs to be multiplied as: n_H+ * rho_gas * dust_gas_ratio * Jion(Z, tgas)
# 3) jelectron.dat: analogous to (2) but for electrons
def prepare(user_energy, ndust=100, nt=30):
    print_title("photoelectric")

    # dust grains range for integration
    arange = np.logspace(np.log10(amin), np.log10(amax), ndust)
    trange = np.logspace(0e0, 8e0, nt)

    #plot_test()

    prepare_jpetot(user_energy, arange)
    prepare_jion(trange, arange)


# prepare a file with Jptot that takes into account the proper photoelectric effect and the photodetachment.
# note that to have the final value these arrays should be integrated against Jflux
def prepare_jpetot(user_energy, arange):
    # get dust optical properties, as the Qabs(a, E) fit in log-log coordinates
    # and the imaginary part of the dielectric fit
    get_qabs, get_eps2 = load_optical(arange, user_energy)

    get_y0_WD06 = load_y0_WD06()

    # loop on charge
    out = out_heating = ""
    for Z in range(Zmin, Zmax+1):
        jpe = np.zeros_like(user_energy)
        jpd = np.zeros_like(user_energy)
        jpe_heating = np.zeros_like(user_energy)
        jpd_heating = np.zeros_like(user_energy)
        for i, e in enumerate(tqdm(user_energy)):
            jpe[i] = get_jpe(e, Z, get_qabs, get_eps2, get_y0_WD06, arange) / hplanck
            jpd[i] = get_jpd(e, Z, arange) / hplanck
            jpe_heating[i] = get_jpe_heating(e, Z, get_qabs, get_eps2, get_y0_WD06, arange) / hplanck
            jpd_heating[i] = get_jpd_heating(e, Z, arange) / hplanck

        # jpe is the array that needs to be integrated against by Jflux * rho_gas * dust_to_gas in the main code, i.e.
        # Jpe = rho_gas * dust_to_gas * trapz(Jflux * jpe, energy)
        # Jpd = rho_gas * dust_to_gas * trapz(Jflux * jpd, energy)
        # these can be summed to a single value and integrated as
        # Jptot = rho_gas * dust_to_gas * trapz(Jflux * jptot, energy)
        jptot = jpe + jpd
        jptot_heating = jpe_heating + jpd_heating
        out += "\n".join(["%d %.18e" % (Z, x) for x in jptot]) + "\n"
        out_heating += "\n".join(["%d %.18e" % (Z, x) for x in jptot_heating]) + "\n"

    fh = open("../runtime_data/jptot.dat", "w")
    fh.write(out)
    fh.close()

    fh = open("../runtime_data/jptot_heating.dat", "w")
    fh.write(out_heating)
    fh.close()


# prepare two files, H+ and electron, with the temperature and Z dependent rates for ion-grain recombination
# and electron-grain attachment
def prepare_jion(trange, arange):

    # H+
    Zpartner = 1
    mpartner = pmass
    out = out_cool = ""
    for Z in range(Zmin, Zmax+1):
        for tgas in trange:
            jion = get_jion(Z, Zpartner, mpartner, tgas, arange)
            out += "%d %.18e %.18e\n" % (Z, tgas, jion)
            jion_cool = get_jion_cooling(Z, Zpartner, mpartner, tgas, arange)
            out_cool += "%d %.18e %.18e\n" % (Z, tgas, jion_cool)

    fh = open("../runtime_data/jion.dat", "w")
    fh.write(out)
    fh.close()

    fh = open("../runtime_data/jion_cooling.dat", "w")
    fh.write(out_cool)
    fh.close()

    # electron
    Zpartner = -1
    mpartner = emass
    out = out_cool = ""
    for Z in range(Zmin, Zmax+1):
        for tgas in trange:
            jion = get_jion(Z, Zpartner, mpartner, tgas, arange)
            out += "%d %.18e %.18e\n" % (Z, tgas, jion)
            jion_cool = get_jion_cooling(Z, Zpartner, mpartner, tgas, arange)
            out_cool += "%d %.18e %.18e\n" % (Z, tgas, jion_cool)

    fh = open("../runtime_data/jelectron.dat", "w")
    fh.write(out)
    fh.close()

    fh = open("../runtime_data/jelectron_cooling.dat", "w")
    fh.write(out_cool)
    fh.close()


def load_y0_WD06():
    energy_ev, yield0 = np.loadtxt("../data/photoelectric_heating/WD2006_y0.dat").T
    return interp1d(np.log10(energy_ev * ev2erg), np.log10(yield0), bounds_error=False, fill_value=-99.)


def get_jion(Z, Zpartner, mpartner, tgas, arange):
    p2 = pslope + 2
    p4 = pslope + 4

    fint = np.array([np.pi * a**p2 * np.sqrt(8e0 * kboltzmann * tgas / mpartner / np.pi)
                     * get_jtilde(a, Z, Zpartner, tgas) for a in arange])

    anorm = 4. / 3. * np.pi * rho_bulk * (amax**p4 - amin**p4) / p4
    return np.pi * np.trapz(fint, arange) / anorm


def load_optical(arange, user_energy):
    if os.path.isfile("Cabs_phe.npy"):
        qabs, kref = np.load("Cabs_phe.npy", allow_pickle=True)
        return qabs, kref

    # w(micron)  Re(eps-1)  Im(eps)    Re(m-1)    Im(m)
    #data_eps = np.loadtxt("../data/dust_refractive_index/eps_Sil.txt").T
    data_eps = np.loadtxt(refInd_file).T
    wl = data_eps[0] * 1e-4  # micron -> cm
    energy = hplanck * clight / wl  # cm -> erg (NOTE: are already reversed)

    if user_energy.min() < energy.min():
        print(user_energy.min(), "<", energy.min())
        sys.exit("ERROR: user energy min is below the optical tables energy lower limit!")
    if user_energy.max() > energy.max():
        print(user_energy.max(), ">", energy.max())
        sys.exit("ERROR: user energy max is above the optical tables energy upper limit!")

    # e_real = data_eps[1]
    e_img = data_eps[2]
    nref = data_eps[3] + 1e0
    kref = data_eps[4]

    qabs = np.zeros((len(arange), len(wl)))
    for i, a in enumerate(tqdm(arange)):
        for j, w in enumerate(wl):
            qabs[i, j] = bhmie_qabs(2 * np.pi * a / w, complex(nref[j], kref[j])) + 1e-40

    np.save("Cabs_phe.npy",
            [interp2d(np.log10(arange), np.log10(energy), np.log10(qabs).T),
             interp1d(np.log10(energy), np.log10(kref))])
    return interp2d(np.log10(arange), np.log10(energy), np.log10(qabs).T), interp1d(np.log10(energy), np.log10(kref))


def get_jpe(E, Z, get_qabs, get_eps2, get_y0_WD06, arange):

    p2 = pslope + 2
    p4 = pslope + 4

    fp = np.array([get_epet(a, Z) <= E for a in arange]).astype(int)
    fint = np.array([a**p2 * get_y(E, a, Z, get_eps2, get_y0_WD06)
                     * 1e1**get_qabs(np.log10(a), np.log10(E))[0] / E for a in arange]) * fp

    if fint.min() < 0e0:
        # plt.semilogx(arange, [get_y0(get_theta(E, a, Z)) for a in arange])
        plt.semilogx(arange, [get_y(E, a, Z, get_eps2, get_y0_WD06) for a in arange])
        if plotOn:
            plt.show()
        else:
            plt.close()
        plt.semilogx(arange, [get_y2(E, a, Z) for a in arange])
        if plotOn:
            plt.show()
        else:
            plt.close()
        sys.exit("min(Jpe) < 0, E: %e, Z: %d" % (E, Z))

    anorm = 4. / 3. * np.pi * rho_bulk * (amax**p4 - amin**p4) / p4
    return np.pi * np.trapz(fint, arange) / anorm

    # idx = energy >= get_epet(a, Z)
    # erange = energy[idx]
    # = np.array([get_y(x, a, Z, get_eps2) for x in erange])
    # fint = [1e1**get_qabs(np.log10(a), np.log10(x))[0] * y[i] * jflux[idx][i] / x for i, x in enumerate(erange)]
    # return np.pi * a**2 / hplanck * np.trapz(fint, erange)  # + int_pdt() / hplanck


def get_y(E, a, Z, get_eps2, get_y0_WD06):
    theta = get_theta(E, a, Z)

    return get_y2(E, a, Z) * np.minimum(1e1**get_y0_WD06(np.log10(E)) * get_y1(E, a, get_eps2), 1e0)
    # return get_y2(E, a, Z) * np.minimum(get_y0(theta) * get_y1(E, a, get_eps2), 1e0)


def get_y0(theta):
    return 0.5 * theta / wpot / (1e0 + 5e0 * theta / wpot)


def get_y1(E, a, get_eps2):
    bb = a / get_la(E, get_eps2)
    aa = bb + a / get_le(E)
    return (bb / aa)**2 * (aa**2 - 2 * aa + 2 - 2 * np.exp(-aa)) / (bb**2 - 2 * bb + 2 - 2 * np.exp(-bb))


def get_y2(E, a, Z, el=None, eh=None):
    if Z < 0:
        return 1e0
    else:
        if eh is None:
            eh = get_eh(E, a, Z)
        if el is None:
            el = get_el(a, Z)
        return eh**2 * (eh - 3 * el) / (eh - el)**3


def get_theta(E, a, Z):
    if Z < 0:
        return max(E - get_epet(a, Z), 0e0)
    else:
        return max(E - get_epet(a, Z) + (Z + 1) * echarge2 / a, 0e0)


def get_la(E, get_eps2):
    return hplanck * clight / E / (4e0 * np.pi * 1e1**get_eps2(np.log10(E)))


def get_le(E=None):
    if E is None:
        return 1e-7

    if E >= 211 * ev2erg:
        return 3.27e-3 * 1e-8 * (E / ev2erg)**1.5
    else:
        return 1e-7


def get_el(a, Z):
    if Z < 0:
        return get_emin(a, Z)
    else:
        return - (Z + 1) * echarge2 / a


def get_eh(E, a, Z):
    if Z < 0:
        return get_emin(a, Z) + E - get_epet(a, Z)
    else:
        return E - get_epet(a, Z)


def get_epet(a, Z):
    if Z >= -1:
        return get_ip(a, Z)
    else:
        return get_ip(a, Z) + get_emin(a, Z)


def get_ip(a, Z):
    return wpot + (Z + 0.5) * echarge2 / a + (Z + 2) * echarge2 / a * 3e-9 / a


def get_emin(a, Z):
    if Z < -1:
        nu = abs(Z + 1)
        theta_nu = nu / (1e0 + 1e0 / np.sqrt(nu))
        return theta_nu * (1 - 0.3 * (a / 1e-7)**(-0.45) * abs(Z + 1)**(-0.26))
        # return - (Z + 1) * echarge2 / a / (1e0 + (2.7e-7 / a)**0.75)
    else:
        return 0e0


def get_jpd(E, Z, arange):

    p4 = pslope + 4

    fp = np.array([get_epdt(E, a, Z) <= E for a in arange]).astype(int)
    fint = np.array([get_sigma_pdt(E, a, Z) * a**pslope for a in arange]) * fp

    if fint.min() < 0e0:
        plt.loglog(arange, fint)
        plt.loglog(arange, -fint)
        if plotOn:
            plt.show()
        else:
            plt.close()
        sys.exit("min(Jpd) < 0, E: %e, Z: %d" % (E, Z))

    anorm = 4. / 3. * np.pi * rho_bulk * (amax**p4 - amin**p4) / p4
    return np.pi * np.trapz(fint, arange) / anorm


def get_epdt(E, a, Z):
    if Z >= 0:
        return 0e0
    return get_ea(a, Z + 1) + get_emin(a, Z)


def get_sigma_pdt(E, a, Z):
    xx = (E - get_epdt(E, a, Z)) / (3 * ev2erg)
    return 1.2e-17 * np.abs(Z) * xx / (1e0 + xx**2 / 3.)**2


def get_ea(a, Z):
    # return wpot + (Z - 0.5) * echarge2 / a - echarge2 * 4e-7 / a / (a + 7e-7)
    return 3 * ev2erg + (Z - 0.5) * echarge2 / a


def get_Zmin(a):
    uait = 2.5 + 0.07 * (a / 1e-7) + 2 * (1e-7 / a)
    return int(-uait / 14.4 * a / 1e-7) + 1


def get_jtilde(a, Z, Zpartner, tgas):
    tau = a * kboltzmann * tgas / echarge2
    nu = Z / Zpartner

    if nu == 0e0:
        jt = 1e0 + np.sqrt(np.pi / 2e0 / tau)
    elif nu < 0e0:
        jt = (1e0 - nu / tau) * (1e0 + np.sqrt(2e0 / (tau - 2e0 * nu)))
    else:
        theta_nu = nu / (1e0 + 1e0 / np.sqrt(nu))
        jt = (1e0 + 1e0 / np.sqrt(4*tau + 3*nu))**2 * np.exp(-theta_nu / tau)

    if Zpartner != -1:
        stick = 1e0
    else:
        if Z == 0:
            stick = 0.5 * (1e0 - np.exp(-a / get_le()))
        elif Z < 0:
            if Z < get_Zmin(a):
                stick = 0e0
            else:
                stick = 0.5 * (1e0 - np.exp(-a / get_le()))
        else:
            stick = 0.5 * (1 - np.exp(-a / get_le()))

    return jt * stick


def get_fe0(E, a, Z, el=None, eh=None):
    if el is None:
        el = get_el(a, Z)
    if eh is None:
        eh = get_eh(E, a, Z)
    return 6 * (E - el) * (eh - E) / (eh - el)**3


def get_jpe_heating(E, Z, get_qabs, get_eps2, get_y0_WD06, arange):

    p2 = pslope + 2
    p4 = pslope + 4

    fp = np.array([get_epet(a, Z) >= 0 for a in arange]).astype(int)
    fint = np.array([a**p2 * get_y(E, a, Z, get_eps2, get_y0_WD06) * 1e1**get_qabs(np.log10(a), np.log10(E))[0]
                     * get_pet_heating_int(E, a, Z) / E for a in arange])

    anorm = 4. / 3. * np.pi * rho_bulk * (amax**p4 - amin**p4) / p4
    return np.pi * np.trapz(fint * fp, arange) / anorm


def get_pet_heating_int(E, a, Z):
    if E <= get_epet(a, Z):
        return 0e0
    if Z >= 0:
        emin = 0e0
    else:
        emin = get_emin(a, Z)
    emax = E - get_epet(a, Z) + emin
    emin = max(emin, 1e-6*ev2erg)
    if emin > emax:
        sys.exit("ERROR: emax < emin!")
    erange = np.logspace(np.log10(emin), np.log10(emax), 300)
    peint = np.trapz([get_fe0(x, a, Z, el=emin, eh=emax) / get_y2(x, a, Z, el=emin, eh=emax) * x for x in erange],
                     erange)
    if peint < 0e0:
        sys.exit("ERROR: negative get_pet_heating_int! %.18e %.18e %d" % (E, a, Z))
    return peint


def get_jpd_heating(E, Z, arange):

    p4 = pslope + 4

    fp = np.array([get_epdt(E, a, Z) >= 0 for a in arange]).astype(int)
    fint = np.array([a**pslope * get_sigma_pdt(E, a, Z) * (E - get_epdt(E, a, Z) + get_emin(a, Z)) for a in arange])

    anorm = 4. / 3. * np.pi * rho_bulk * (amax**p4 - amin**p4) / p4
    return np.pi * np.trapz(fint * fp, arange) / anorm


def get_jion_cooling(Z, Zpartner, mpartner, tgas, arange):
    p2 = pslope + 2
    p4 = pslope + 4

    fint = np.array([a**p2 * np.sqrt(8e0 * kboltzmann * tgas / mpartner / np.pi)
                     * get_lambda_tilde(a, Z, Zpartner, tgas) for a in arange])

    anorm = 4. / 3. * np.pi * rho_bulk * (amax**p4 - amin**p4) / p4
    return np.pi * np.trapz(fint, arange) / anorm * kboltzmann * tgas


def get_lambda_tilde(a, Z, Zpartner, tgas):
    tau = a * kboltzmann * tgas / echarge2
    nu = Z / Zpartner

    if nu < 0e0:
        jl = (2e0 - nu / tau) * (1e0 + 1e0 / np.sqrt(tau - nu))
    elif nu > 0e0:
        theta_nu = nu / (1e0 + 1e0 / np.sqrt(nu))
        jl = (2e0 + nu / tau) * (1e0 + 1e0 / np.sqrt(3e0 / 2e0 / tau + 3e0 * nu)) * np.exp(-theta_nu / tau)
    else:
        jl = 2e0 + 1.5 * np.sqrt(np.pi / 2e0 / tau)

    if Zpartner != -1:
        stick = 1
    else:
        if Z == 0:
            stick = 0.5 * (1e0 - np.exp(-a / get_le()))
        elif Z < 0:
            if Z < get_Zmin(a):
                stick = 0e0
            else:
                stick = 0.5 * (1e0 - np.exp(-a / get_le()))
        else:
            stick = 0.5 * (1 - np.exp(-a / get_le()))

    return jl * stick


# *********************

def plot_test():
    from matplotlib import pyplot as plt

    erange = np.logspace(0, 3, 3000) * ev2erg
    arange = np.array([4, 10, 30, 100, 300]) * 1e-8
    get_qabs, get_eps2 = load_optical(arange, erange)

    p2 = pslope + 2
    p4 = pslope + 4
    anorm = 4. / 3. * np.pi * rho_bulk * (amax ** p4 - amin ** p4) / p4

    get_y0_WD06 = load_y0_WD06()

    for a in arange:
        aa = a / 1e-8
        Z = 0
        y = [get_y(x, a, Z, get_eps2, get_y0_WD06) * (x >= get_epet(a, Z)).astype(int) for x in erange]
        # y = [get_pet_heating_int(x, a, Z) * (x >= get_epet(a, Z)).astype(int) for x in erange]
        #y = [a**p2 * a * get_y(E, a, Z, get_eps2, get_y0_WD06) * 1e1**get_qabs(np.log10(a), np.log10(E))[0]
        #             * get_pet_heating_int(E, a, Z) / E / anorm / hplanck for E in erange]
        plt.loglog(erange / ev2erg, y, label=aa)
    plt.legend(loc="best")
    if plotOn:
        plt.show()
    else:
        plt.close()
