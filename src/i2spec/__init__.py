"""i2spec: open, physics-based model of the molecular iodine B-X spectrum."""

from .constants import MHZ_PER_CM, reduced_mass
from .model import RovibronicModel
from .potentials import XRepPotential, load_potentials
from .research import ResearchModel

__all__ = ["MHZ_PER_CM", "ResearchModel", "RovibronicModel", "XRepPotential", "load_potentials", "reduced_mass"]
__version__ = "0.0.1"
