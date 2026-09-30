"""i2spec: open, physics-based model of the molecular iodine B-X spectrum."""

from .constants import MHZ_PER_CM, reduced_mass
from .local_nir import LocalNIRModel, load_local_nir
from .model import RovibronicModel
from .potentials import XRepPotential, load_potentials

__all__ = ["MHZ_PER_CM", "LocalNIRModel", "RovibronicModel", "XRepPotential", "load_local_nir",
           "load_potentials", "reduced_mass"]
__version__ = "0.0.1"
