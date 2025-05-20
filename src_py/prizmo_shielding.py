import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator as rgi
from prizmo_commons import print_title, plotOn


def prepare(H2_inc, CO_inc):
    if H2_inc:
        shielding_H2()
    if CO_inc:
        shielding_CO()

def shielding_H2(nt=100, nn=30):

    print_title("shielding H2")

    trange = np.logspace(0, 8, nt)
    nrange = np.logspace(12, 22, nn)
    b5 = 1e0

    out = ""
    for n in nrange:
        f = np.zeros_like(trange)
        for i, t in enumerate(trange):
            if t < 3e3:
                a = 1.4
                ncrit = 1.3e14 * (1. + (t / 6e2)**0.8)
            elif t >= 4e3:
                a = 1.1
                ncrit = 2e14
            else:
                a = (t / 4.5e3)**(-0.8)
                ncrit = 1e14 * (t / 4.76e3)**(-3.8)

            x = n / ncrit
            w = 0.013 * (1e0 + (t / 2.7e3)**1.3)**(1. / 1.3) * np.exp(-(t / 3.9e3)**14.6)
            f[i] = (1e0 - w) / (1e0 + x / b5)**a * np.exp(-5e-7 * (1e0 + x)) \
                + w / np.sqrt(1e0 + x) * np.exp(-8.5e-4 * np.sqrt(1e0 + x))

            out += "%.18e %.18e %.18e\n" % (n, t, f[i])
        plt.loglog(trange, f, label="%.2e" % n)

    plt.legend(loc="best")
    if plotOn:
        plt.show()
    else:
        plt.close()

    fout = open("../runtime_data/shielding_H2.dat", "w")
    fout.write(out)
    fout.close()


def shielding_CO(nco=50, nh2=50):

    print_title("shielding CO")

    in_NCO = in_NH2 = in_shield = False
    NCO = []
    NH2 = []
    shield = []
    for row in open("../data/CO_shielding/shield.03.5.69-557-36.dat"):
        srow = row.strip()
        if srow == "N(12CO)":
            in_NCO, in_NH2, in_shield = True, False, False
            continue
        if srow == "N(H2)":
            in_NCO, in_NH2, in_shield = False, True, False
            continue
        if srow == "12C16O":
            in_NCO, in_NH2, in_shield = False, False, True
            continue
        if srow == "12C17O":
            break

        if in_NCO:
            NCO.append(float(srow))
        elif in_NH2:
            NH2.append(float(srow))
        elif in_shield:
            shield += [float(x) for x in srow.split()]
        else:
            continue

    NCO = np.array(NCO)
    NH2 = np.array(NH2)
    shield = np.array(shield).reshape((len(NH2), len(NCO)))
    f_shield = rgi((np.log10(NH2), np.log10(NCO)), np.log10(shield))

    out = ""
    for xNCO in np.logspace(np.log10(NCO.min()), np.log10(NCO).max(), nco):
        for xNH2 in np.logspace(np.log10(NH2.min()), np.log10(NH2).max(), nh2):
            out += "%.18e %.18e %.18e\n" % (xNCO, xNH2, 1e1**f_shield((np.log10(xNH2), np.log10(xNCO) + 1e-40)))

    fh = open("../runtime_data/shielding_CO.dat", "w")
    fh.write(out)
    fh.close()
