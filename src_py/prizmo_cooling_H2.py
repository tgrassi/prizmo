import numpy as np
import matplotlib.pyplot as plt
from prizmo_commons import print_title, plotOn


def prepare():
    # number of temperature grid points
    npoints = 1000

    print_title("H2 cooling")

    trange = np.logspace(0, 8, npoints)
    coolH2H = np.zeros_like(trange)
    coolH2Hj = np.zeros_like(trange)
    coolH2H2 = np.zeros_like(trange)
    coolH2e = np.zeros_like(trange)
    coolHDL = np.zeros_like(trange)

    for i, temp in enumerate(trange):
        coolH2H[i] = np.clip(fH2H(temp), 1e-80, 1e99)
        coolH2Hj[i] = np.clip(fH2Hj(temp), 1e-80, 1e99)
        coolH2H2[i] = np.clip(fH2H2(temp), 1e-80, 1e99)
        coolH2e[i] = np.clip(fH2e(temp), 1e-80, 1e99)
        coolHDL[i] = np.clip(fHDL(temp), 1e-80, 1e99)

    np.savetxt("../runtime_data/cool_H2_H.dat", np.array([trange, coolH2H]).T)
    np.savetxt("../runtime_data/cool_H2_Hj.dat", np.array([trange, coolH2Hj]).T)
    np.savetxt("../runtime_data/cool_H2_H2.dat", np.array([trange, coolH2H2]).T)
    np.savetxt("../runtime_data/cool_H2_e.dat", np.array([trange, coolH2e]).T)
    np.savetxt("../runtime_data/cool_H2_HDL.dat", np.array([trange, coolHDL]).T)

    if plotOn:
        plt.loglog(trange, coolH2Hj, label="H+")
        plt.loglog(trange, coolH2H, label="H")
        plt.loglog(trange, coolH2H2, label="H2")
        plt.loglog(trange, coolH2e, label="e-")
        plt.loglog(trange, coolHDL, label="HDL")
        plt.legend()
        plt.show()
    else:
        plt.close()


def fHDL(temp):
    t3 = 1e-3 * temp
    logt3 = np.log10(temp * 1e-3)
    logt = np.log10(temp)
    logt32 = logt3 * logt3
    logt33 = logt32 * logt3
    logt34 = logt33 * logt3
    logt35 = logt34 * logt3
    logt36 = logt35 * logt3
    logt37 = logt36 * logt3
    logt38 = logt37 * logt3

    if temp < 2e3:
        HDLR = ((9.5e-22 * t3 ** 3.76) / (1. + 0.12 * t3 ** 2.1) * np.exp(-(0.13 / t3) ** 3)
                + 3.e-24 * np.exp(-0.51 / t3))
        HDLV = (6.7e-19 * np.exp(-5.86 / t3) + 1.6e-18 * np.exp(-11.7 / t3))
        HDL = HDLR + HDLV
    elif 2e3 <= temp <= 1e4:
        HDL = 1e1 ** (-2.0584225e1 + 5.0194035 * logt3
                      - 1.5738805 * logt32 - 4.7155769 * logt33
                      + 2.4714161 * logt34 + 5.4710750 * logt35
                      - 3.9467356 * logt36 - 2.2148338 * logt37
                      + 1.8161874 * logt38)
    else:
        # dump14 = 1e0 / (1e0 + np.exp(min((temp - 3e4) * 2e-4, 3e2)))
        w14 = wcool(logt, 1e0, 4e0)
        HDL = 5.531333679406485e-19 * w14

    return HDL


def fH2H(temp):
    logt3 = np.log10(temp * 1e-3)
    logt = np.log10(temp)

    logt32 = logt3 * logt3
    logt33 = logt32 * logt3
    logt34 = logt33 * logt3
    logt35 = logt34 * logt3

    if temp <= 1e2:
        f = 1e1 ** (-16.818342e0 + 3.7383713e1 * logt3
                    + 5.8145166e1 * logt32 + 4.8656103e1 * logt33
                    + 2.0159831e1 * logt34 + 3.8479610e0 * logt35)
    elif 1e2 < temp <= 1e3:
        f = 1e1 ** (-2.4311209e1 + 3.5692468e0 * logt3
                    - 1.1332860e1 * logt32 - 2.7850082e1 * logt33
                    - 2.1328264e1 * logt34 - 4.2519023e0 * logt35)
    elif 1e3 < temp <= 6e3:
        f = 1e1 ** (-2.4311209e1 + 4.6450521e0 * logt3
                    - 3.7209846e0 * logt32 + 5.9369081e0 * logt33
                    - 5.5108049e0 * logt34 + 1.5538288e0 * logt35)
    else:
        f = 1.862314467912518e-22 * wcool(logt, 1e0, np.log10(6e3))

    return f


def fH2Hj(temp):
    logt3 = np.log10(temp * 1e-3)
    logt = np.log10(temp)
    logt32 = logt3 * logt3
    logt33 = logt32 * logt3
    logt34 = logt33 * logt3
    logt35 = logt34 * logt3

    w14 = wcool(logt, 1e0, 4e0)
    if temp <= 1e4:
        f = 1e1 ** (-2.2089523e1 + 1.5714711e0 * logt3
                    + 0.015391166e0 * logt32 - 0.23619985e0 * logt33
                    - 0.51002221e0 * logt34 + 0.32168730e0 * logt35)
    else:
        f = 1.182509139382060E-021 * w14

    return f


def fH2H2(temp):
    logt3 = np.log10(temp * 1e-3)
    logt = np.log10(temp)
    logt32 = logt3 * logt3
    logt33 = logt32 * logt3
    logt34 = logt33 * logt3
    logt35 = logt34 * logt3
    w24 = wcool(logt, 2e0, 4e0)

    return w24 * 1e1 ** (-2.3962112e1 + 2.09433740e0 * logt3
                         - .77151436e0 * logt32 + .43693353e0 * logt33
                         - .14913216e0 * logt34 - .033638326e0 * logt35)


def fH2e(temp):
    logt3 = np.log10(temp * 1e-3)
    logt = np.log10(temp)
    logt32 = logt3 * logt3
    logt33 = logt32 * logt3
    logt34 = logt33 * logt3
    logt35 = logt34 * logt3
    logt36 = logt35 * logt3
    logt37 = logt36 * logt3
    logt38 = logt37 * logt3

    w24 = wcool(logt, 2e0, 4e0)

    if temp <= 5e2:
        f = 1e1 ** (min(-2.1928796e1 + 1.6815730e1 * logt3
                        + 9.6743155e1 * logt32 + 3.4319180e2 * logt33
                        + 7.3471651e2 * logt34 + 9.8367576e2 * logt35
                        + 8.0181247e2 * logt36 + 3.6414446e2 * logt37
                        + 7.0609154e1 * logt38, 3e1))
    else:
        f = 1e1 ** (-2.2921189e1 + 1.6802758e0 * logt3
                    + .93310622e0 * logt32 + 4.0406627e0 * logt33
                    - 4.7274036e0 * logt34 - 8.8077017e0 * logt35
                    + 8.9167183 * logt36 + 6.4380698 * logt37
                    - 6.3701156 * logt38)

    return f * w24


def sigmoid(x, x0, s):
    return 1e1 / (1e1 + np.exp(-s * (x - x0)))


def wcool(logTgas, logTmin, logTmax):
    x = (logTgas - logTmin) / (logTmax - logTmin)
    wCool = 1e1 ** (2e2 * (sigmoid(x, -2e-1, 5e1) * sigmoid(-x, -1.2e0, 5e1) - 1e0))
    if wCool < 1e-199:
        return 0e0
    else:
        return wCool
