# Fafoom - Flexible algorithm for optimization of molecules 

Fafoom is a tool for sampling the conformational space of organic molecules. Fafoom is intended to work with FHI-aims (Fritz Haber Institute ab initio molecular simulations package).

## Requirements

* functionality of the tool:
  * Python 2.7
  * Numpy
  * RDKit (used version: Release_2014_09_2)

* first-principles methods:
  * (recommended) FHI-aims (Fritz Haber Institute ab initio molecular simulations package)
  * (alternative) NWChem (NWChem: Open Source High-Performance Computational Chemistry)

## How to use

##### 1) Clone the fafoom repository
	git clone https://github.com/adrianasupady/fafoom

##### 2) Export the fafoom directory to you PYTHONPATH

##### 3) In python:

    import fafoom

## Examples of usage

Three examples are provided with the source code in the examples folder:

* genetic algorithm based search utilizing first-principles (FHI-aims required)
* genetic algorithm based search utilizing force fields (force fields accessed from RDKit)
* genetic algorithm based search utilizing first-principles via NWChem (NWChem required)

Get familiar with the provided manual to learn more about the tool and the parameters. 

## License

Fafoom is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Fafoom is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with fafoom.  If not, see <http://www.gnu.org/licenses/>.


Copyright 2015 Adriana Supady 
adriana.supady@gmail.com
