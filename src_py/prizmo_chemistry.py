import sys

import numpy as np
from prizmo_preprocess import preprocess
from prizmo_commons import nphoto, sp2idx, print_title, ev2erg, idx2sp, count_X, idx2spj, idx2mass, \
    amin, amax, pslope, rho_bulk, kboltzmann, sp2mass, py2f90, use_reaction_prototypes


def prepare(fname="../networks/test.dat", main=False, speciesList=None):
    print_title("chemistry")
    rows = ""
    for row in open(fname):
        srow = row.strip()
        if srow == "" or srow.startswith("#"):
            continue
        rows += srow + "\n"

    data_chemical_cooling = load_chemical_cooling()

    rows = rows.replace("&\n", "&").strip()
    fluxes = krates = krates_photo = load_xsecs = photo_thresholds = photoheating_rate = photoheating = ""
    recombination_cooling = attenuate = H2diss = chemical_cooling = ""

    reactants = []
    products = []
    species = []
    photo_limits = []
    verbatim_reactions = []
    prototype = prototype_vars = prototype_idxs = None
    prototype_pragma = prototype_define = ""
    i = 0
    for line in rows.split("\n"):

        if line.startswith("PROTOTYPE["):
            if use_reaction_prototypes:
                prototype = parse_prototype(line)
            continue
        if line == "}":
            prototype_pragma += get_prototype_pragma(prototype_vars, prototype_idxs, prototype)
            prototype_define += get_prototype_define(prototype_vars)
            prototype = prototype_vars = prototype_idxs = None
            continue

        rr_org, pp_org, krate = parse_line(line)
        rr = [sp2idx(x) for x in rr_org]
        pp = [sp2idx(x) for x in pp_org]

        verbatim_rate = get_verbatim_rate(rr_org, pp_org)
        check_charge_conservation(rr_org, pp_org, verbatim_rate)

        fluxes += "flux(%d) = %s\n" % (i + 1, get_rhs(i + 1, rr))
        rate = parse_krate(i + 1, krate, verbatim_rate, prototype)
        krates += rate
        prototype_vars, prototype_idxs = parse_prototype_rate(rate, prototype, prototype_vars, prototype_idxs, i+1)
        krates_photo += parse_krate_photo(i + 1, krate, verbatim_rate)
        photoheating_rate += parse_photoheating_rate(i + 1, krate, verbatim_rate)
        photo_limits.append(parse_photolimit(krate))
        photoheating += parse_photoheating(i + 1, krate, rr, verbatim_rate)
        H2diss += parse_H2diss(i + 1, krate, verbatim_rate)
        attenuate += parse_attenuation(i + 1, krate, rr, verbatim_rate)
        photo_thresholds += parse_threshold(krate)
        load_xsecs += parse_xsecs(i + 1, krate, rr_org, pp_org)
        recombination_cooling += parse_recombination_cooling(i + 1, rr_org, pp_org, verbatim_rate)
        chemical_cooling += parse_chemical_cooling(i + 1, verbatim_rate, data_chemical_cooling)

        verbatim_reactions.append(verbatim_rate)

        reactants.append(rr)
        products.append(pp)
        species += rr + pp
        i += 1

    # check_reverse(reactants, products)

    species = np.sort(np.unique(species))

    ode = ""
    for sp in species:
        ode += "dy(%s) = " % sp
        for i, rr in enumerate(reactants):
            for _ in range(rr.count(sp)):
                ode += " - flux(%d)" % (i + 1)
        for i, pp in enumerate(products):
            for _ in range(pp.count(sp)):
                ode += " + flux(%d)" % (i + 1)
        ode += "\n"

    # write python useful stuff
    species_py = "species_names = [%s]" % ", ".join(["\"%s\"" % idx2sp(x) for x in species])
    plot_commons = species_py
    if main:
        open("plot_commons.py", "w").write(plot_commons)
    else:
        open("../plot_commons.py", "w").write(plot_commons)

    # prepare get_electrons
    cations = ["%d * x(%s)" % (x.count("j"), x) for x in species if "j" in x]
    electron_sum = "ne = " + " + ".join([x.replace("1 * ", "") for x in cations])

    # prepare get_nuclei for different atoms
    xnuclei_pragma = dict()
    for xa in ["H", "C", "O", "He"]:
        xnuclei = ["%d * x(%s)" % (count_X(x, xa), x) for x in species if count_X(x, xa) > 0]
        nuclei = [x.replace("1 * ", "") for x in xnuclei]
        if len(nuclei) == 0:
            xnuclei = "n%s = 0d0" % xa
        else:
            xnuclei = ("n%s =" % xa) + " + ".join(nuclei)
        xnuclei_pragma[xa.upper() + "NUCLEI"] = xnuclei

    # save verbatim reactions to file
    verbatim_reactions = [x.ljust(100) for x in verbatim_reactions]
    if not main:
        open("../runtime_data/reactions.dat", "w").write("\n".join(verbatim_reactions))

    # useful PLUTO variables saved to file
    species_pluto = ["#define IDX_CHEM_%s (NFLX + NIONS + %d)" % (idx2spj(x), i) for i, x in enumerate(species)]
    pluto_txt = "#define NTRACER %d\n" % len(species)
    pluto_txt += "#define NPHOTO %d\n\n" % nphoto
    pluto_txt += "\n".join(species_pluto)
    if main:
        open("pluto_definitions.txt", "w").write(pluto_txt)
    else:
        open("../pluto_definitions.txt", "w").write(pluto_txt)

    # mass array in commons
    masses = "real*8,parameter::masses(nspecies) = (/" + \
             ", &\n\t".join([("%.18e" % idx2mass(x)).replace("e", "d") for x in species]) \
             + "/)\n"

    # prepare get_rho
    rhos = ["x(%s) * masses(%s)" % (x, x) for x in species if "DUST" not in x]
    get_rho = "rho = " + "&\n + ".join(rhos)

    # finalize pre-processor strings
    common_vars = "integer,parameter::nspecies=%d\n" % len(species)
    common_vars += "integer,parameter::nphoto=%d\n" % nphoto
    common_vars += "integer,parameter::nreactions=%d\n\n" % len(reactants)

    indexes = ""
    for i, sp in enumerate(species):
        common_vars += "integer,parameter::%s=%d\n" % (sp, i + 1)
        indexes += "integer,parameter::prizmo_%s=%s\n" % (sp, sp)

    if not main:
        preprocess("../prizmo_ode.f90", {"ODE": ode})
        preprocess("../prizmo_rates.f90", {"RATES": krates.replace("e", "d").replace("Hd", "He").replace("dxp", "exp").replace("usdr", "user").replace("ddbyd", "debye").replace("*invTgas","d0*invTgas"),
                                           "PROTOTYPES": prototype_pragma,
                                           "PROTOTYPES_DEFINE": prototype_define})
        preprocess("../prizmo_rates_photo.f90", {"PHOTORATES": krates_photo})
        preprocess("../prizmo_rates_heating.f90", {"PHOTOHEATING_RATE": photoheating_rate})
        preprocess("../prizmo_heating_photo.f90", {"PHOTOHEATING": photoheating})
        preprocess("../prizmo_heating_H2diss.f90", {"H2DISS": H2diss})
        preprocess("../prizmo_cooling_chemical.f90", {"RECOMBINATION": recombination_cooling,
                                                      "CHEMICAL": chemical_cooling})
        preprocess("../prizmo_loaders.f90", {"LOAD_XSECS": load_xsecs})
        preprocess("../prizmo_flux.f90", {"FLUXES": fluxes})
        preprocess("../prizmo_attenuate.f90", {"ATTENUATE": attenuate})
        preprocess("../prizmo_commons.f90", {"COMMON_VARS": common_vars,
                                             "MASSES": masses})
        preprocess("../prizmo_utils.f90", {"ELECTRONS": electron_sum,
                                           "GET_RHO": get_rho})
        preprocess("../prizmo_utils.f90", xnuclei_pragma)
        preprocess("../prizmo.f90", {"INDEXES": indexes})

        open("../runtime_data/energy_thresholds.dat", "w").write(photo_thresholds)

    return species, np.array(photo_limits)


def check_reverse(reactants, products):
    rrs = [[idx2sp(y) for y in x] for x in reactants]
    pps = [[idx2sp(y) for y in x] for x in products]
    verbs = []
    rev_verbs = []
    for rr, pp in zip(rrs, pps):
        verbs.append("_".join(sorted(rr)) + "__" + "_".join(sorted(pp)))
        rev_verbs.append("_".join(sorted(pp)) + "__" + "_".join(sorted(rr)))

    for v, rv in zip(verbs, rev_verbs):
        if rv not in verbs:
            print(v)


def get_prototype_define(prototype_vars):
    if prototype_vars is None:
        return ""

    defs = ["%s(%d)" % (v, len(k)) for v, k in prototype_vars.items()]
    return "real*8::" + ", ".join(defs) + "\n"


def get_prototype_pragma(prototype_vars, prototype_idxs, prototype):
    if prototype_vars is None:
        return ""

    s = ""
    for v, k in prototype_vars.items():
        nums = [py2f90(x) for x in k]
        s += "%s = (/%s/)" % (v, ", ".join(nums)) + "\n"

    imin = min(prototype_idxs)
    imax = max(prototype_idxs)
    s += "kall(%d:%d) = %s\n" % (imin, imax, prototype)

    s += "\n"

    return s


def parse_prototype_rate(rate, prototype, prototype_vars, prototype_idxs, idx):
    if prototype is None:
        return None, None
    crate = rate.lower().replace(" ", "").split("=")[1].strip()
    cproto = prototype.lower()

    vs = []
    for j in range(10):
        vs += ["var%d_%02d" % (i, j) for i in range(10) if "var%d_%02d" % (i, j) in prototype]

    for v in vs:
        cproto = cproto.replace(v, "|")

    parts_p = cproto.split("|")
    parts_p = [x for x in parts_p if x.strip() != ""]
    for p in parts_p:
        crate = crate.replace(p, "|")
    parts_val = crate.split("|")

    if prototype_vars is None:
        prototype_vars = {v: [] for v in vs}

    for v, val in zip(vs, parts_val):
        prototype_vars[v].append(val)

    if prototype_idxs is None:
        prototype_idxs = []
    prototype_idxs.append(idx)

    return prototype_vars, prototype_idxs


def parse_prototype(line):
    prototype = line.replace("]", "[").split("[")[1]
    prototype = prototype.lower().replace(" ", "")
    return prototype


def check_charge_conservation(rr_names, pp_names, verbatim):
    rr_charge = sum([x.count("+") for x in rr_names]) - sum([x == "E" for x in rr_names])
    pp_charge = sum([x.count("+") for x in pp_names]) - sum([x == "E" for x in pp_names])
    if rr_charge != pp_charge:
        print("problems with charge conservation in", verbatim)
        sys.exit()


def load_chemical_cooling():
    data = dict()
    for row in open("../data/chemical_cool_heat/cooling_heating.dat"):
        srow = row.strip()
        if srow == "" or srow.startswith("#"):
            continue
        verbatim, energy_eV = [x.strip() for x in srow.split(";")]
        energy = float(energy_eV) * ev2erg
        data[verbatim] = energy

    return data


def parse_chemical_cooling(i, verbatim, data_chemical_cooling):
    if verbatim in data_chemical_cooling:
        rate = "! %s\n" % verbatim
        rate += "coola(%d) = (%s) * fluxes(%d)\n\n" % (i, ("%.18e" % -data_chemical_cooling[verbatim]).replace("e", "d"), i)
        return rate
    else:
        return ""


def get_fname_rate(rr_names, pp_names):
    return "_".join(rr_names) + "__" + "_".join(pp_names)


def get_verbatim_rate(rr_names, pp_names):
    return (" + ".join(rr_names) + " -> " + " + ".join(pp_names)).strip()


def parse_recombination_cooling(i, rr_name, pp_name, verbatim):
    if len(rr_name) != 2:
        return ""
    if "E" not in rr_name:
        return ""
    if "".join(rr_name).count("+") - "".join(pp_name).count("+") != 1:
        return ""
    rec = "! %s\n" % verbatim
    rec += "coola(%d) = fluxes(%d)\n\n" % (i, i)
    return rec


def parse_threshold(krate):
    if krate.split(",")[0].strip() != "PHOTO":
        return "0e0\n"
    return krate.split(",")[1].strip() + "\n"


def parse_photolimit(krate):
    if krate.split(",")[0].strip() != "PHOTO":
        return -1e99

    phlim = float(krate.split(",")[1].strip())
    if phlim > 1e90:
        phlim = -1e99

    return phlim * ev2erg


def parse_attenuation(i, krate, rr, verbatim):
    if krate.split(",")[0].strip() != "PHOTO":
        return ""
    if verbatim.strip() == "H2 -> H + H":
        return ""
    if verbatim.strip() == "CO -> C + O":
        return ""
    att = "tau = tau + x(%s) * photo_xsecs(:, %d)\n" % (rr[0], i)
    return att


def parse_xsecs(i, krate, rr_names, pp_names):
    if krate.split(",")[0].strip() != "PHOTO":
        return ""
    fname_rate = get_fname_rate(rr_names, pp_names)
    k = "photo_xsecs(:, %d) = load_photo_xsecs(\"photo_xsecs_%s.dat\")\n" % (i, fname_rate)
    return k


def parse_krate_photo(i, krate, verbatim):
    if krate.split(",")[0].strip() != "PHOTO":
        return ""

    if verbatim.strip() == "H2 -> H + H":
        shielding = " * shielding_H2(log_NH2, log_tgas)"
    elif verbatim.strip() == "CO -> C + O":
        shielding = " * shielding_CO(log_NH2, log_NCO)"
    else:
        shielding = ""

    k = "! %s\n" % verbatim
    k += "f(:) = photo_xsecs(:, %d) * kernel\n" % i
    k += "kall(%d) = sum((f(2:nphoto) + f(1:nphoto-1)) * delta_energy) / 2d0 %s\n\n" % (i, shielding)

    return k


def parse_photoheating_rate(i, krate, verbatim):
    if krate.split(",")[0].strip() != "PHOTO":
        return ""

    if verbatim.strip() == "H2 -> H + H":
        k = "! %s\n" % verbatim
        k += "! skipped\n"
        return k

    k = "! %s\n" % verbatim
    k += "f(:) = photo_xsecs(:, %d) * max(energy - energy_threshold(%d), 0d0) * kernel\n" % (i, i)
    k += "kall_heat(%d) = sum((f(2:nphoto) + f(1:nphoto-1)) * delta_energy) / 2d0\n\n" % i

    return k


def parse_photoheating(i, krate, rr, verbatim):
    if krate.split(",")[0].strip() != "PHOTO":
        return ""

    heat = "! %s\n" % verbatim

    if verbatim.strip() == "H2 -> H + H":
        return "\n"

    heat += "heat = heat + kall_heat(%d) * x(%s)\n" % (i, rr[0])
    return heat


def parse_H2diss(i, krate, verbatim):
    if krate.split(",")[0].strip() != "PHOTO":
        return ""
    if verbatim.strip() == "H2 -> H + H":
        return "Rdiss = kall(%d)\n" % i
    return ""


def parse_krate(i, krate, verbatim, prototype):
    k = "! %s\n" % verbatim
    if krate.split(",")[0].strip() == "PHOTO":
        k += "! "
    if prototype is not None:
        k += "! "

    if "DUST" in verbatim:
        p3 = pslope + 3.
        p4 = pslope + 4.
        pre_dust = np.sqrt(8e0 * kboltzmann / np.pi) * (amax ** p3 - amin ** p3) / (amax ** p4 - amin ** p4) \
                   * p4 / p3 / (4. / 3. * rho_bulk)
        if is_dust_evaporation(verbatim):
            k += "kall(%d) = nu_debye * exp(-%s / Tdust)\n\n" % (i, ("%.18e" % float(krate)).replace("e", "d"))
        elif is_dust_freezing(verbatim):
            sp = verbatim.split("->")[0]
            inv_sqrt_mass = 1e0 / np.sqrt(sp2mass(sp))
            k += "kall(%d) = %.18e * rho_dust * sticking * sqrTgas * %.18e \n\n" % (i, pre_dust, inv_sqrt_mass)
            k += "kall(%d) = %s * rho_dust * sticking * sqrTgas * %s \n\n" % (i, ("%.18e" % pre_dust).replace("e", "d"), ("%.18e" % inv_sqrt_mass).replace("e", "d"))
        else:
            sys.exit("ERRROR: dust reaction uknown " + verbatim)
    else:
        k += "kall(%d) = %s\n\n" % (i, krate.replace("&", "&\n"))

    return k


def is_dust_evaporation(verbatim):
    rv, pv = verbatim.split("->")
    return "DUST" in rv and "DUST" not in pv


def is_dust_freezing(verbatim):
    rv, pv = verbatim.split("->")
    return "DUST" not in rv and "DUST" in pv


def parse_line(line):
    if "[" in line:
        line = line.replace("[", ";").replace("]", ";").replace("->", ";")
        rr, pp, _, kk = line.split(";")
    else:
        line = line.replace("->", ";")
        rr, pp, kk = line.split(";")
    rr = parse_species(rr)
    pp = parse_species(pp)
    return rr, pp, kk.strip()


def parse_species(specs):
    return [x.strip() for x in specs.split(" + ")]


def get_rhs(idx, rr):
    reacts = " * ".join(["x(%s)" % x for x in rr])
    return "kall(%d) * %s" % (idx, reacts)

def main():
    # Can be run just to extract species from a network
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--network", "-n", type=str, help='Network to use')
    args = parser.parse_args()

    prepare(fname=args.network, main=True)

if __name__ == "__main__":
    main()

