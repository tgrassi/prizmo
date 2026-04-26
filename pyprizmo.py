
import numpy as np
import os
import uuid
from ctypes import cdll, POINTER, c_double, c_int
from glob import glob


class Prizmo:
    def __init__(self, preprocess=True, arguments=None, network_text=None, debug=False):

        self.lib_name = "libprizmo.so"

        if preprocess:
            self.lib_name = f"libprizmo_{uuid.uuid4().hex}.so"
            self.preprocessor(arguments, network_text=network_text)
            self.compile(debug=debug)

        self.lib = cdll.LoadLibrary(f"./{self.lib_name}")

        if preprocess:
            os.rename(self.lib_name, "libprizmo.so")

        self.init_bindings()

        self.load_variables()

        self.lib.prizmo_init_c()

    def cleaner(self):
        if not os.getcwd().endswith("/src_py"):
            os.chdir("src_py")
        os.system("python clean.py")

        if os.getcwd().endswith("/src_py"):
            os.chdir("..")
        print(os.getcwd())
        os.system("make clean")

    def preprocessor(self, arguments=None, network_text=None):
        if not os.getcwd().endswith("/src_py"):
            os.chdir("src_py")

        # to avoid key check issues, e.g. if "chemNet" is in arguments
        if arguments is None:
            arguments = {}

        if "energy_minmax" in arguments:
            if type(arguments["energy_minmax"]) == list:
                arguments["energy_minmax"] = " ".join(map(str, arguments["energy_minmax"]))

        if "chemNet" in arguments and network_text is not None:
            raise ValueError("Cannot specify both 'chemNet' in arguments and 'network_text'. Please choose one.")

        args = ""
        if arguments:
            args = " ".join(["--{} {}".format(k, v) for k, v in arguments.items()])

        if network_text is not None:
            # generate temporary file with random hash name to avoid conflicts
            fname = "temp_network_{}.dat".format(np.random.randint(1e6))
            tmp_path = os.path.join("../networks", fname)
            with open(tmp_path, "w") as f:
                f.write(network_text)
            args += " --chemNet {}".format(tmp_path)

        ret = os.system("python prizmo.py " + args)
        if ret != 0:
            raise RuntimeError("Preprocessing failed with return code {}".format(ret))

        # clean up temporary file if it was created
        if network_text is not None:
            os.remove(tmp_path)

    def compile(self, debug):

        os.environ["BUILD_LIB"] = self.lib_name

        if os.getcwd().endswith("/src_py"):
            os.chdir("..")
        print(os.getcwd())
        os.system("make clean")

        if debug:
            ret = os.system("make lib_debug")
        else:
            ret = os.system("make lib")

        if ret != 0:
            raise RuntimeError("Preprocessing failed with return code {}".format(ret))


    def init_bindings(self):

        self.lib.prizmo_init_c.argtypes = None
        self.lib.prizmo_init_c.restype = None

        self.lib.prizmo_evolve_c.argtypes = [POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_int), POINTER(c_int)]
        self.lib.prizmo_evolve_c.restype = None

        self.lib.prizmo_set_crate_c.argtypes = [POINTER(c_double)]
        self.lib.prizmo_set_crate_c.restype = None

        self.lib.prizmo_set_d2g_c.argtypes = [POINTER(c_double)]
        self.lib.prizmo_set_d2g_c.restype = None

        self.lib.prizmo_get_cooling_array_c.argtypes = [POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double)]
        self.lib.prizmo_get_cooling_array_c.restype = None

        self.lib.prizmo_get_heating_array_c.argtypes = [POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double)]
        self.lib.prizmo_get_heating_array_c.restype = None

        self.lib.prizmo_get_tdust_c.argtypes = [POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double)] #(x, Tgas, Tdust, jflux)
        self.lib.prizmo_get_tdust_c.restype = None

        self.lib.prizmo_get_electrons_c.argtypes = [POINTER(c_double), POINTER(c_double)]
        self.lib.prizmo_get_electrons_c.restype = None

        self.lib.prizmo_set_solve_thermo_c.argtypes = [POINTER(c_int)]
        self.lib.prizmo_set_solve_thermo_c.restype = None

        self.lib.prizmo_set_solve_chemistry_c.argtypes = [POINTER(c_int)]
        self.lib.prizmo_set_solve_chemistry_c.restype = None

        self.lib.prizmo_rt_c.argtypes = [POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double)]
        self.lib.prizmo_rt_c.restype = None

    def load_variables(self):
        self.energy = np.loadtxt("runtime_data/energy.dat")
        self.nphoto = len(self.energy)

        erg2ev = 6.241509074e11
        self.energy_ev = self.energy * erg2ev

        clight = 2.99792458e10  # cm/s
        hplanck = 6.62607015e-27  # erg s
        self.wavelength_cm = clight * hplanck / self.energy
        self.wavelength_micron = self.wavelength_cm * 1e4

        self.dust_opacity = np.loadtxt("runtime_data/kappa_dust.dat")

        self.reactions = open("runtime_data/reactions.dat").readlines()
        self.reactions = [x.strip() for x in self.reactions if x.strip()]
        self.nreactions = len(self.reactions)

        self.species = open("runtime_data/species.dat").readlines()
        self.species = [x.strip() for x in self.species if x.strip()]
        self.nspecies = len(self.species)

        self.xsecs = {}
        for g in glob("runtime_data/photo_xsecs_*.dat"):
            rr = g.replace("runtime_data/photo_xsecs_", "").replace(".dat", "")
            rr = rr.replace("__", " -> ").replace("_", " + ")
            self.xsecs[rr] = np.loadtxt(g)


    def species2index(self, sp):
        if sp in self.species:
            return self.species.index(sp)
        else:
            raise ValueError("Species {} not found in species list.".format(sp))


    def evolve(self, x, Tgas, jflux, dt, solve_thermo=True, solve_chemistry=True):
        x = (c_double * len(x))(*x)
        jflux = (c_double * len(jflux))(*jflux)

        tgas = c_double(Tgas)

        ierr = 0
        verboseChem = 0

        self.lib.prizmo_set_solve_thermo_c(c_int(1 if solve_thermo else 0))
        self.lib.prizmo_set_solve_chemistry_c(c_int(1 if solve_chemistry else 0))

        self.lib.prizmo_evolve_c(
            x,
            tgas,
            jflux,
            c_double(dt),
            c_int(verboseChem),
            c_int(ierr)
        )

        return np.array(list(x)), tgas.value

    def set_crate(self, crate):
        self.lib.prizmo_set_crate_c(c_double(crate))

    def set_d2g(self, d2g):
        self.lib.prizmo_set_d2g_c(c_double(d2g))

    def get_tdust(self, x, Tgas, jflux):
        x = (c_double * len(x))(*x)  # Convert to a C array
        jflux = (c_double * len(jflux))(*jflux)  # Convert to a C array

        Tdust = c_double(0.0)  # Initialize a variable to hold the dust temperature

        self.lib.prizmo_get_tdust_c(
            x,
            c_double(Tgas),
            Tdust,
            jflux
        )

        return Tdust.value  # Return the dust temperature as a Python float

    def get_cooling_array(self, x, Tgas, jflux):

        Tdust = self.get_tdust(x, Tgas, jflux)

        x = (c_double * len(x))(*x)  # Convert to a C array
        jflux = (c_double * len(jflux))(*jflux)  # Convert to a C array

        cools = [0.0] * 5  # Initialize a list for cooling values
        cools = (c_double * len(cools))(*cools)  # Create an array to hold cooling values

        self.lib.prizmo_get_cooling_array_c(
            x,
            c_double(Tgas),
            c_double(Tdust),
            jflux,
            cools
        )

        return np.array(list(cools))  # convert to a Python list

    def get_heating_array(self, x, Tgas, jflux):

        Tdust = self.get_tdust(x, Tgas, jflux)

        x = (c_double * len(x))(*x)  # Convert to a C array
        jflux = (c_double * len(jflux))(*jflux)  # Convert to a C array

        heats = [0.0] * 4  # Initialize a list for heating values
        heats = (c_double * len(heats))(*heats)  # Create an array to hold heating values

        self.lib.prizmo_get_heating_array_c(
            x,
            c_double(Tgas),
            c_double(Tdust),
            jflux,
            heats
        )

        return np.array(list(heats))  # convert to a Python list

    def get_electrons(self, x):
        x = (c_double * len(x))(*x)  # Convert to a C array

        ne = c_double(0.0)  # Initialize a variable to hold the electron density

        self.lib.prizmo_get_electrons_c(
            x,
            ne
        )

        return ne.value  # Return the electron density as a Python float

    def RT(self, x, Tgas, jflux, ds):
        x = (c_double * len(x))(*x)  # Convert to a C array
        jflux = (c_double * len(jflux))(*jflux)  # Convert to a C array

        self.lib.prizmo_rt_c(x,
                             c_double(Tgas),
                             jflux,
                             c_double(ds))

        return np.array(list(jflux))