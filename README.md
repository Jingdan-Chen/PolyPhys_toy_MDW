# PolyPhys: Polymer Physics Simulation Framework

PolyPhys is a comprehensive Python-based framework for simulating polymerization kinetics and polymer physics. This software enables researchers to model and analyze various polymerization processes with different reaction mechanisms and conditions.

## Overview

The framework consists of two main simulation modules:

1. **Simple_Simulation**: A basic implementation of polymerization kinetics focusing on core reaction mechanisms.
2. **Diffusion_Control**: An advanced implementation incorporating diffusion-controlled effects in polymerization processes.

Both simulation modules share similar structures but offer different levels of complexity for modeling polymerization reactions.

## Key Features

- Modeling of various polymerization mechanisms:
  - Free radical polymerization
  - Self-initiation
  - Depolymerization
  - Chain transfer to monomer
- Support for different monomers:
  - Methyl Methacrylate (MMA)
  - Styrene (Sty)
- Configurability through JSON parameters
- Comprehensive kinetics modeling with temperature-dependent rate constants
- Detailed simulation outputs for analysis
- Visualization and analysis tools through Jupyter notebooks

## Installation

1. Clone the repository:
```bash
git clone https://github.com/yourusername/PolyPhys_code.git
cd PolyPhys_code
```

2. Ensure you have the required dependencies:
```bash
pip install numpy matplotlib jupyter
```

## Usage

### Running Simulations

To run a basic polymerization simulation:

```bash
cd Simple_Simulation
python MMA_1.py  # For MMA monomer with config 1
```

Or for Styrene:

```bash
python Sty_1.py  # For Styrene monomer with config 1
```

### Configuration

Simulation parameters are controlled through JSON configuration files in the `exam_json` directory. These files allow you to specify:

- Initial concentrations
- Reaction conditions (temperature, time)
- Polymerization mechanisms to include (via flags)
- Model parameters

Example configuration:
```python
{
    "I0": 0.0258,      // Initial concentration of initiator
    "n_max": 300,      // Maximum polymerization degree considered
    "ending_thresh": 0.85, // Threshold for ending simulation
    "max_time": 200.0,  // Maximum simulation time (min)
    "dt": 0.000001,     // Time step for simulation (min)
    "flag1": true,      // Enable self-initiation
    "flag2": true,      // Enable depolymerization
    "flag3": true,      // Enable CT to monomer
    "flag4": false,     // Enable decay for kp
    "flag5": true,      // Enable ktc depreciation
    "x": 1,             // x=1 or 2 for self-initiation
    "y": 1,             // y=1 or 2 for self-initiation
    "T": 90.0           // Temperature (°C)
}
```

### Analysis

To analyze simulation results, use the included Jupyter notebook:

```bash
jupyter notebook analyze_result.ipynb
```

The notebook provides tools for:
- Visualizing conversion vs time
- Examining radical concentration profiles
- Analyzing polymer chain length distributions
- Calculating average polymerization degrees

## Code Structure

### Simple_Simulation

- `MMA_*.py`, `Sty_*.py`: Main simulation scripts for different monomers
- `*_AIBN_param.py`: Parameter files with temperature-dependent rate constants
- `utils/`: Core utilities
  - `io.py`: Configuration file handling
  - `kinetics.py`: Polymerization kinetics implementation
  - `__init__.py`: Package initialization

### Diffusion_Control

Similar to Simple_Simulation but with additional diffusion control mechanisms.

## Reaction Mechanisms

The framework models several key polymerization processes:

1. **Initiation**: Formation of radicals from initiator molecules
2. **Propagation**: Chain growth through monomer addition
3. **Termination**: Chain termination via combination or disproportionation
4. **Chain Transfer**: Transfer of the radical to monomer (when flag3 is true)
5. **Self-Initiation**: Spontaneous formation of radicals from monomers (when flag1 is true)
6. **Depolymerization**: Reverse of propagation (when flag2 is true)

## Mathematical Model

The system is modeled using a set of coupled differential equations:

- Initiator concentration kinetics: `d[I]/dt = -ki*[I]`
- Monomer concentration kinetics: `d[M]/dt = -kp*[M]*∑[R_n] - ...`
- Radical concentration kinetics: `d[R_n]/dt = ...`
- Polymer concentration kinetics: `d[P_n]/dt = ...`

These equations are solved numerically using a time-stepping approach.

## License

This project is licensed under the MIT License - see the LICENSE file for details.
