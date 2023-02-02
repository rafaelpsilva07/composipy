import os
import sys

from ._version import VERSION
from .material import OrthotropicMaterial, IsotropicMaterial
from .property import LaminateProperty
from .panel import PlateStructure

__version__ = VERSION