import prizmo_photo
import prizmo_tdust
import prizmo_cooling_atomic
import prizmo_chemistry
import prizmo_photoelectric
import prizmo_dust_opacity
import prizmo_cooling_H2
import prizmo_shielding
import prizmo_cooling_CO
import prizmo_heating_CR
import prizmo_heating_cooling
import prizmo_main
import numpy as np
import warnings
import os
from prizmo_commons import init, idx2sp, chemNet, atomData, radiation_type, plotOn, erg2ev, args
from prizmo_preprocess import preprocess

np.seterr(divide="raise", over="raise", invalid="raise")
warnings.simplefilter(action='ignore', category=FutureWarning)

# get current working directory
cwd = os.getcwd()

if not cwd.endswith("src_py"):
    print("Please run this script from the src_py directory.")
    exit(1)

init()

species, photo_limits = prizmo_chemistry.prepare(fname=chemNet)
species_names = [idx2sp(x) for x in species]
H2_inc = 'H2' in species_names
CO_inc = 'CO' in species_names

prizmo_cooling_atomic.prepare_atomic_cooling(species, H2_inc, fname=atomData)

user_energy = prizmo_photo.prepare(photo_limits, species)

prizmo_dust_opacity.prepare(user_energy)

# prizmo_photoelectric.prepare(user_energy)

prizmo_tdust.prepare(user_energy)

prizmo_shielding.prepare(H2_inc, CO_inc)
if H2_inc:
    prizmo_cooling_H2.prepare()
    preprocess("../prizmo.f90", {"H2": "call load_H2_cooling_tabs()\ncall load_shielding_H2_table()"})
else:
    preprocess("../prizmo.f90", {"H2": ""})
if CO_inc:
    prizmo_cooling_CO.prepare()
    preprocess("../prizmo.f90", {"CO": "call load_CO_cooling()\ncall load_shielding_CO_table()"})
else:
    preprocess("../prizmo.f90", {"CO": ""})

prizmo_heating_CR.prepare(H2_inc)

prizmo_heating_cooling.prepare(H2_inc, CO_inc)

prizmo_main.prepare(H2_inc, CO_inc)

f = open("../runtime_data/README.txt","w")
for arg, val in args.__dict__.items():
    f.write("{}: {}\n".format(arg, val))
f.close()

print("************************")
print("All preprocessing done!")


