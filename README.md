# Fafoom - Flexible algorithm for optimization of molecules 

Fafoom is a tool for sampling the conformational space of organic molecules. Fafoom is intended to work with FHI-aims (Fritz Haber Institute ab initio molecular simulations package).

## Requirements

* functionality of the tool:
  * Python (used for testing: 2.7.6)
  * Numpy (used for testing: 1.8.2)
  * RDKit (used for testing: Release_2015_03_1)

* first-principles methods:
  * (recommended) FHI-aims (Fritz Haber Institute ab initio molecular simulations package)
  * (alternative) NWChem (NWChem: Open Source High-Performance Computational Chemistry)
  * (alternative) ORCA (- An ab initio, DFT and semiempirical SCF-MO package -)

## How to use

##### 1) Clone the fafoom repository
	git clone https://github.com/adrianasupady/fafoom

##### 2) Export the fafoom directory to you PYTHONPATH

##### 3) In python:

    import fafoom

## Example of usage

An implementation of a genetic algorithm is provided in the examples folder.
Depending on the used parameter file following genetic algorithm based searches can be run:

* parameters_aims.txt for first-principles (FHI-aims required)
* paremeters_ff.txt for force fields (force fields accessed from RDKit)
* parameters_nwchem.txt for first-principles via NWChem (NWChem required)
* parameters_orca.txt for first-principles via ORCA (ORCA required)

Get familiar with the provided manual to learn more about the tool and the parameters. 

## Outlook

Development goals:

* flexible formulation of the scoring function, e.g. allowing for optimizing a user-defined property
* adding more kinds of degrees of freedom, e.g. orientation of a molecule
* adding wrappers for different molecular simulations packages

Comments, feedback or development contributions are always welcome! 



## License

Fafoom is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Fafoom is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with fafoom.  If not, see <http://www.gnu.org/licenses/>.


Copyright 2015 Adriana Supady 
