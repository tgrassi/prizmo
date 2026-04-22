
import os
from ctypes import cdll, POINTER, c_double, c_int
import numpy as np


class Prizmo:
    def __init__(self, preprocess=True, arguments=None):

        if preprocess:
            self.preprocessor(arguments)
            self.compile()

        self.lib = cdll.LoadLibrary("./libprizmo.so")

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


    def preprocessor(self, arguments=None):
        if not os.getcwd().endswith("/src_py"):
            os.chdir("src_py")

        args = ""
        if arguments:
            args = " ".join(["--{} {}".format(k, v) for k, v in arguments.items()])

        os.system("python prizmo.py " + args)

    def compile(self):

        if os.getcwd().endswith("/src_py"):
            os.chdir("..")
        print(os.getcwd())
        os.system("make clean")
        os.system("make lib")

    def init_bindings(self):

        self.lib.prizmo_init_c.argtypes = None
        self.lib.prizmo_init_c.restype = None

        self.lib.prizmo_evolve_c.argtypes = [POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_int), POINTER(c_int)]
        self.lib.prizmo_evolve_c.restype = None

        self.lib.prizmo_set_crate_c.argtypes = [POINTER(c_double)]
        self.lib.prizmo_set_crate_c.restype = None

        self.lib.prizmo_get_cooling_array_c.argtypes = [POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double)]
        self.lib.prizmo_get_cooling_array_c.restype = None

        self.lib.prizmo_get_tdust_c.argtypes = [POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double)] #(x, Tgas, Tdust, jflux)
        self.lib.prizmo_get_tdust_c.restype = None

        self.lib.prizmo_get_electrons_c.argtypes = [POINTER(c_double), POINTER(c_double)]
        self.lib.prizmo_get_electrons_c.restype = None

    def load_variables(self):
        self.energy = np.loadtxt("runtime_data/energy.dat")
        self.nphoto = len(self.energy)

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

    def species2index(self, sp):
        if sp in self.species:
            return self.species.index(sp)
        else:
            raise ValueError("Species {} not found in species list.".format(sp))


    def evolve(self, x, Tgas, jflux, dt):
        x = (c_double * len(x))(*x)
        jflux = (c_double * len(jflux))(*jflux)

        tgas = c_double(Tgas)

        ierr = 0
        verboseChem = 0

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

    def get_cooling_array(self, x, Tgas, jflux):
        x = (c_double * len(x))(*x)  # Convert to a C array
        jflux = (c_double * len(jflux))(*jflux)  # Convert to a C array

        Tdust = 0e0 # dummy value
        # get the Tdust
        self.lib.prizmo_get_tdust_c(x, c_double(Tgas), c_double(Tdust), jflux)

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


    def get_electrons(self, x):
        x = (c_double * len(x))(*x)  # Convert to a C array

        ne = c_double(0.0)  # Initialize a variable to hold the electron density

        self.lib.prizmo_get_electrons_c(
            x,
            ne
        )

        return ne.value  # Return the electron density as a Python float