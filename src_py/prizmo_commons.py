### Constants and Conversions
erg2ev = 6.24150647996e11  # erg -> eV
ev2erg = 1e0 / erg2ev  # eV -> erg
v2statv = 0.0033356405  # V -> statV
clight = 2.99792458e10  # cm/s
clight_nm = clight * 1e7  # nm/s
hplanck_eV = 4.135667696e-15  # eV*s
hplanck = hplanck_eV / erg2ev  # erg*s
kboltzmann = 1.380658e-16  # erg / K
kboltzmann_eV = kboltzmann * erg2ev  # eV / K
pmass = 1.6726219e-24  # g
emass = 9.10938e-28  # g
echarge = 4.80320425e-10  # statC
echarge2 = echarge**2  # stataC^2


### User Input Properties
import argparse
plotOn=False

# dust properties
amin = 5e-7  # cm
amax = 2.5e-5  # cm
pslope = -3.5  # MNR slope
rho_bulk = 3.  # g/cm3
wpot = 8e0 / erg2ev  # erg
refInd_file = "../data/dust_refractive_index/silD03.txt"

# chemical network file
chemNet = "../networks/network_approx.dat"

# atomic data file
atomData = "../data/atomic_cooling/krome_data.dat"

# number of photobins
nphoto = 1000

# photo energy range, erg
energy_min = 1e-2 * ev2erg
energy_max = 1e4 * ev2erg

# radiation grid spacing (excluded thresholds)
photo_logspacing = True

# radiation types are, "draine", "xdr", "BB@Tbb", "file@LX@X_lo-X_hi",
# xdr loads radiation from data/spectra/xdr_spectrum.dat
# file columns are: eV, eV/cm2/Hz/s
radiation_type = "E09_spec.dat@2.04e30@3e2-1e4"

# Get command line overwrites
parser = argparse.ArgumentParser()
parser.add_argument("--input", "-i", type=str, default=None, help='Specify input file (default: None)')
parser.add_argument("--chemNet", "-c", type=str, default=chemNet, help='Specify the chemical network (default: {})'.format(chemNet))
parser.add_argument("--atomData", "-a", type=str, default=atomData, help='Specify the atomic data file (default: {})'.format(atomData))
parser.add_argument("--radiation_type", "-r", type=str, default=radiation_type, help='Specify the radiation type (default: {})'.format(radiation_type))
parser.add_argument("--nphoto", "-n", type=int, default=nphoto, help='Specify the number of photobins (default: {})'.format(nphoto))
parser.add_argument("--energy_minmax", "-E", type=float, nargs='+', default=[energy_min/ev2erg, energy_max/ev2erg], help='Specify the maximum and minimum energies (default: {})'.format([energy_min/ev2erg, energy_max/ev2erg]))
parser.add_argument("--dust_minmax", "-d", type=float, nargs='+', default=[amin, amax], help='Specify the maximum and minimum dust grain sizes (default: {})'.format([amin,amax]))
parser.add_argument("--refInd_file", "-e", type=str, default=refInd_file, help='Specify the refactory index file (default: {})'.format(redInd_file))
parser.add_argument("--plot", "-p", action='store_true', help='Show plots produced by each stage')
args = parser.parse_args()
input_file = args.input
chemNet = args.chemNet
atomData = args.atomData
radiation_type = args.radiation_type
nphoto = args.nphoto
energy_min = min(args.energy_minmax) * ev2erg
energy_max = max(args.energy_minmax) * ev2erg
amin = min(args.dust_minmax)
amax = max(args.dust_minmax)
refInd_file = args.redInd_file
plotOn = args.plot

def parse_input_file(fname):
    opts = dict()
    for row in open(fname):
        srow = row.strip()
        if not srow or srow.startswith("#"):
            continue
        key, val = srow.split("=")
        key = key.strip()
        val = val.strip()
        if key == "plot":
            val = val.lower() in ["true", "1", "yes", "on", "t", "y"]
        elif key == "energy_minmax":
            energy_min, energy_max = [float(x) for x in val.replace(",", " ").split()]
            opts[energy_min] = energy_min * ev2erg
            opts[energy_max] = energy_max * ev2erg
        elif key == "dust_minmax":
            amin, amax = [float(x) for x in val.replace(",", " ").split()]
            opts[amin] = amin
            opts[amax] = amax
        elif key == "nphoto":
            val = int(val)

        opts[key] = val

    return opts

if input_file:
    print("Reading input file: {}".format(input_file))
    opts = parse_input_file(input_file)
    for k, v in opts.items():
        if k in globals():
            print("Overwriting {} with {} from file {}".format(k, v, input_file))
            globals()[k] = v
        else:
            print("Unknown option: {}".format(k))

### Limits

# charge limits (for photoelectric effect)
Zmin, Zmax = -2, 2

# FUV limits, erg
fuv_energy1 = 9.69e-12  # 2050 AA in erg (FUV lower limit)
fuv_energy2 = 2.178e-11  # 912 AA in erg (FUV upper limit)

### Atomic data

# use prototypes for vectorization
use_reaction_prototypes = True

natom2name = {1: "H", 2: "He", 3: "Li", 4: "Be", 5: "B", 6: "C", 7: "N", 8: "O", 9: "F", 10: "Ne", 11: "Na", 12: "Mg",
              13: "Al", 14: "Si", 15: "P", 16: "S", 17: "Cl", 18: "Ar", 19: "K", 20: "Ca", 21: "Sc", 22: "Ti", 23: "V",
              24: "Cr", 25: "Mn", 26: "Fe", 27: "Co", 28: "Ni", 29: "Cu", 30: "Zn", 31: "Ga", 32: "Ge", 33: "As",
              34: "Se", 35: "Br", 36: "Kr", 37: "Rb", 38: "Sr", 39: "Y", 40: "Zr", 41: "Nb", 42: "Mo", 43: "Tc",
              44: "Ru", 45: "Rh", 46: "Pd", 47: "Ag", 48: "Cd", 49: "In", 50: "Sn", 51: "Sb", 52: "Te", 53: "I",
              54: "Xe", 55: "Cs", 56: "Ba", 57: "La", 58: "Ce", 59: "Pr", 60: "Nd", 61: "Pm", 62: "Sm", 63: "Eu",
              64: "Gd", 65: "Tb", 66: "Dy", 67: "Ho", 68: "Er", 69: "Tm", 70: "Yb", 71: "Lu", 72: "Hf", 73: "Ta",
              74: "W", 75: "Re", 76: "Os", 77: "Ir", 78: "Pt", 79: "Au", 80: "Hg", 81: "Tl", 82: "Pb", 83: "Bi",
              84: "Po", 85: "At", 86: "Rn", 87: "Fr", 88: "Ra", 89: "Ac", 90: "Th", 91: "Pa", 92: "U", 93: "Np"}
name2natom = {el: Z for Z, el in natom2name.items()}

### Useful Functions

def init():
    import os
    if not os.path.exists("../runtime_data/"):
        os.makedirs("../runtime_data/")


def sp2spj(sp):
    return sp.replace("+", "j").replace("-", "q")


def idx2spj(idx):
    sp = idx2sp(idx)
    return sp.replace("+", "j").replace("-", "q")


def sp2idx(sp):
    return "idx_" + sp.replace("+", "j").replace("-", "q")


def idx2sp(idx):
    return idx.replace("idx_", "").replace("j", "+").replace("q", "-")


def idx2mass(idx):
    sp = idx2sp(idx)
    return sp2mass(sp)


def sp2mass(sp):
    if sp == "E":
        return emass

    mdict = {"H": 1*pmass,
             "He": 4*pmass,
             "C": 12*pmass,
             "O": 16*pmass}
    mass = 0e0
    for k, v in mdict.items():
        mass += count_X(sp, k) * v
    mass -= emass * sp.count("+")

    return mass


def print_title(message):
    print("******************")
    print(message.upper())


def count_X(sarg, Xin="H"):
    arg = sarg.replace("_DUST", "")
    arg = arg.replace("He", "W")
    if Xin == "He":
        X = "W"
    else:
        X = Xin
    if X not in arg:
        return 0

    hcount = arg.count(X)
    sold = ""
    for s in arg:
        if sold == X and is_int(s):
            hcount += int(s) - 1
        sold = s

    return hcount


def is_int(x):
    try:
        int(x)
        return True
    except ValueError:
        return False


def py2f90(arg):
    return("%e" % float(arg)).replace("e", "d")
