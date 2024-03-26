from prizmo_commons import print_title
from prizmo_preprocess import preprocess

def prepare(H2_inc):
    print_title("cr heating")
    if H2_inc:
        preprocess("../prizmo_heating_CR.f90", {"CR_HEATING": "heat = user_cr * (5.5d-12 * x(idx_H) + 2.5d-11 * x(idx_H2))"})
    else:
        preprocess("../prizmo_heating_CR.f90", {"CR_HEATING": "heat = user_cr * (5.5d-12 * x(idx_H) )"})


