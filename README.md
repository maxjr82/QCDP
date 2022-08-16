# QCDP - Quantum Chemistry Data Processing

Python script to process and extract information from output files of quantum chemistry packages.
This script is currently working for Gaussian09.

## Requirements

To run this script, you need:

- python3 (tested with version 3.8.1)
- argparse
- numpy
- re

## How to use

You can simply run the script from the terminal as follows:

```
python3 get_properties.py -log output_from_g09 -prop "list of properties"
```

The list of properties should be separated by commas.

The list of available 'prop' arguments are:

- 'labels': atom labels
- 'coords': Cartesian coordinates
- 'forces': atomic forces
- 'Etot': Total energy
- 'e_hl': HOMO/LUMO energies
- 'alpha': polarizability
- 'dipole': dipole moment
  'quadrupole': quadrupole moment
- 'e_vib': Zero-point vibrational energy
- 'rot': Rotational constants
- 'r2': Electronic spatial extent
- 'e_thermal': Thermal Energies/Enthalpies/Free Energies
- 'freqs': Vibrational frequencies
- 'mulliken': atomic charges

