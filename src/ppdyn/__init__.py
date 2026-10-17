"""PPDyn: Plasma Particle Dynamics.

Molecular dynamics of charged grains (dusty plasma, debris) with Yukawa
interactions, neutral drag and external fields. Numba-parallel.
"""
__version__ = "2.0.0"

from .config import ConfigError, load, loads  # noqa: F401
from .simulation import Simulation, run  # noqa: F401
