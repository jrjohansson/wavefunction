"""
Python package for calculating wavefunctions for 1D and 2D potentials.
"""

# utilities
from wavefunction.misc import *
from wavefunction.nummath import *
from wavefunction.version import version as __version__

# 1D problems
from wavefunction.potential import *
from wavefunction.wavefunction1d import *

# 2D problems
from wavefunction.wavefunction2d import *

from wavefunction.matrix_recast import *
from wavefunction.matrix_generate import *
