from prizmo_commons import print_title
from prizmo_preprocess import preprocess

def prepare(H2_inc, CO_inc):
    print_title("heating header")

    if H2_inc:
        preprocess("../prizmo_heating.f90", {"USE_H2DISS": "use prizmo_heating_H2diss"})
        preprocess("../prizmo_heating.f90", {"H2DISS_HEATING": "heats(4) = heating_H2diss(x, Tgas, ntot)"})
    else:
        preprocess("../prizmo_heating.f90", {"USE_H2DISS": ""})
        preprocess("../prizmo_heating.f90", {"H2DISS_HEATING": ""})

    print_title("cooling header")

    use_code = ""
    cools_code = ""
    if H2_inc:
        use_code   += "use prizmo_cooling_H2\n"
        cools_code += "cools(4) = cooling_H2(x, log_Tgas)\n"
    if CO_inc:
        use_code += "use prizmo_cooling_CO\n"
        cools_code += "cools(5) = cooling_CO(x, log_Tgas, log_Hnuclei)\n"

    preprocess("../prizmo_cooling.f90", {"USE_MOLECULES": use_code})
    preprocess("../prizmo_cooling.f90", {"MOLECULAR_COOLING": cools_code})

