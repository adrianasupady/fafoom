#    Copyright 2015 Adriana Supady
#    adriana.supady@gmail.com
#
#    This file is part of fafoom.
#
#   Fafoom is free software: you can redistribute it and/or modify
#   it under the terms of the GNU Lesser General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   Fafoom is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#   along with fafoom.  If not, see <http://www.gnu.org/licenses/>.

'''Wrapper for RDKit force-field routines'''

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import ChemicalForceFields


class FFObject():
    """Create and handle force-field objects."""
    def __init__(self, parameter_file):
        """Get the parameters from file.
        Args(required):
            parameter_file (str): name of the parameter file
        Raises:
            KeyError: if the force field in not defined
            ValueError: if the force field is not 'uff' or 'mmff94'
        Optional parameters for the minimization method used for optimization
        that can be defined in the file:
        steps (default=1000)
        force_tol (default=1e-4)
        energy_tol (default=1e-6)
        """
        with open(parameter_file) as fin:
            parameter_dict = dict(line.strip().partition(' ')[::2] for line in fin)
        if 'force_field' not in parameter_dict:
            raise KeyError("The force field is not defined.")
        self.force_field = parameter_dict['force_field']
        if self.force_field not in ['uff', 'mmff94']:
            raise ValueError("Unknown force field.")
        s = [('steps', 1000), ('force_tol', 1.0e-4), ('energy_tol', 1.0e-6)]
        for k, v in s:
            if k not in parameter_dict:
                parameter_dict[k] = v
        self.steps = parameter_dict['steps']
        self.force_tol = parameter_dict['force_tol']
        self.energy_tol = parameter_dict['energy_tol']

    def run_ff(self, sdf_string):
        """Perform the force field minimization.

        Args:
            sdf_string (str)
        """
        mol = Chem.MolFromMolBlock(sdf_string, removeHs=False)
        if self.force_field == 'mmff94':
            molprop = ChemicalForceFields.MMFFGetMoleculeProperties(mol)
            ff = ChemicalForceFields.MMFFGetMoleculeForceField(mol, molprop)
        elif self.force_field == 'uff':
            ff = AllChem.UFFGetMoleculeForceField(mol)
        ff.Minimize(int(self.steps),
                    float(self.force_tol), float(self.energy_tol))
        self.sdf_string_opt = Chem.MolToMolBlock(mol)
        self.energy = float('{0:.4f}'.format(ff.CalcEnergy()))

    def get_energy(self):
        """Get the energy of the molecule.

        Returns:
            energy (float)
        Raises:
            AttributeError: if energy hasn't been calculated yet
        """
        if not hasattr(self, 'energy'):
            raise AttributeError("The calculation wasn't performed yet.")
        else:
            return self.energy

    def get_sdf_string_opt(self):
        """Get the optimized sdf string.

        Returns:
            optimized sdf string (str)
        Raises:
            AttributeError: if the optimization hasn't been performed yet
        """
        if not hasattr(self, 'sdf_string_opt'):
            raise AttributeError("The calculation wasn't performed yet.")
        else:
            return self.sdf_string_opt
