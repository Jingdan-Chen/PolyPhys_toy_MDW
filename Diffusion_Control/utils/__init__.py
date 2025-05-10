"""
PolyPhys utility functions for polymer physics simulations.

This package contains utilities for:
- IO operations (configuration management)
- Kinetics calculations for polymerization models
"""

# Import functions from io module
from .io import (
    read_and_validate_config
)

# Import functions from kinetics module
from .kinetics import (
    pad_array_to_length,
    I_kinetic,
    M_kinetic,
    R1_kinetic,
    Rn_kinetic,
    Pn_kinetic,
    evolution
)

# Define what's available when using "from utils import *"
__all__ = [
    # IO functions
    'read_and_validate_config',
    
    # Kinetics functions
    'pad_array_to_length',
    'I_kinetic',
    'M_kinetic',
    'R1_kinetic',
    'Rn_kinetic',
    'Pn_kinetic',
    'evolution'
]
