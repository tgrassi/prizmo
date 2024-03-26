from prizmo_commons import print_title
from prizmo_preprocess import preprocess
import os

def prepare(H2_inc, CO_inc):
    print_title("main")

    if not os.path.isfile("../main.f90"):
        print("skipping, main.f90 file not found")
        return

    update_code = ""

    if H2_inc:
        preprocess("../main.f90", {"INITIAL_H2": "x(prizmo_idx_H2) = 0d0"})
        update_code += "rad_Ncol_H2 = rad_Ncol_H2 + x(prizmo_idx_H2) * dr\n"
    else:
        preprocess("../main.f90", {"INITIAL_H2": ""})

    if CO_inc:
        update_code += "rad_Ncol_CO = rad_Ncol_CO + x(prizmo_idx_CO) * dr\n"
        update_code += "vert_Ncol_CO(ix) = vert_Ncol_CO(ix) + x(prizmo_idx_CO) * dz\n"

    preprocess("../main.f90", {"UPDATE_COLUMN": update_code})
