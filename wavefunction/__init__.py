"""
Python package for calculating wavefunctions for 1D and 2D potentials.
"""

# utilities
from wavefunction.utils import *
from wavefunction.version import version as __version__

# 1D problems
import wavefunction.wavefunction1d as wf1d

# 2D problems
import wavefunction.wavefunction2d as wf2d
